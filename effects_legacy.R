#### calculate various quantities needed for evaluating RMSTs ####

#### Libraries etc. ####
library(tidyverse)
library(statmod)
library(splines2)
library(survival)
nsk <- function(...)splines2::nsk(...)

# read surv_dat here

# create model matrix with possibly dropping given columns
my.D <- function(arg_list,data){
  if(is.null(arg_list[["drop"]])){
    model.matrix(arg_list[["formula"]],data)
  }else{
    model.matrix(arg_list[["formula"]],data)[,-arg_list[["drop"]]]
  }
}

# define function that calculates spline basis integrals with given time window history
bh_breaks <- c(0,quantile(surv_dat[["stop_time"]][surv_dat[["status"]]==1], seq(.1,.9,.1)),Inf)
kn_x <- c(1,3,5)
kn_z <- c(1,5)
Bk <- c(0,10)
# read posterior sample matrix
theta <- read_csv2("posterior_leg.csv") %>%
  select(starts_with(c("beta","gamma","zeta","xi","lambda","r[","sd_","sigma","L[")))
#set up data with each w_strata
wdat <- expand_grid(sex=0:1,age0f=factor(0:2),smoke=factor(0:2), score0=factor(0:2)) %>%
  mutate(w_strata = factor(1:nrow(.)))
wdat_u <- expand_grid(wdat,score=factor(0:2))

#### Calculate posterior sample of the probability distribution of U ####
# read posterior sample matrix for treatment-dependent confounder
theta_mult <- read_csv2("posterior_u.csv") %>%
  select(starts_with("psi"))
params_u <- theta_mult
dat_u <- expand_grid(gr=0:1,sex=0:1,age0f=factor(0:2), smoke=factor(0:2),score0=factor(0:2)) %>%
  left_join(wdat) #dataframe representing each stratum
X_u <- model.matrix(~gr+sex+age0f+smoke+score0,dat_u%>%select(gr,sex,age0f,smoke,score0)%>% distinct) # design matrix

## calculate probabilities P(U=u|W) for each stratum
pr_U_list <- lapply(1:nrow(params_u),function(i){
  partmp <- matrix(unlist(params_u[i,]), ncol=2, byrow=F)
  linpreds <- X_u %*% partmp
  pr_mat <- apply(linpreds,1,function(vec)c("pr_0"=1/(1+exp(vec[1])+exp(vec[2])),
                                            "pr_1"=exp(vec[1])/(1+exp(vec[1])+exp(vec[2])),
                                            "pr_2" = exp(vec[2])/(1+exp(vec[1])+exp(vec[2])))) %>% t
  cbind(dat_u%>%select(gr,sex,age0f,smoke,score0,w_strata)%>% distinct,pr_mat) %>%
    pivot_longer(cols=starts_with("pr_"), names_to=c(".value","score"), names_pattern="(pr)_(\\d)")
})
pr_U_tmp <- cbind(pr_U_list[[1]] %>% select(-pr),lapply(pr_U_list,function(ll)ll%>%select(pr))) %>%
  `colnames<-`(c("gr","sex","age0f","smoke","score0","w_strata","score", paste0("iter_",1:8000)))

## Check which iterations are not compatible with monotonicity assumption
system.time(
  bad_iters <- lapply(unique(wdat$w_strata),function(ws){
    dtmp <- pr_U_tmp %>% filter(w_strata==ws)
    iters <- dtmp %>% select(starts_with("iter"))
    
    outer <- with(dtmp, 1- unlist(iters[score==0&gr==0,])-unlist(iters[score==2&gr==1,]))
    inners <- rbind(with(dtmp, unlist(iters[score==1&gr==1,])), with(dtmp, unlist(iters[score==1&gr==0,]))) %>% apply(2,min)
    
    as.numeric(outer > inners)
  }) %>%
    do.call(rbind,.) %>%
    apply(2,function(x)any(x==1)) %>% which
)

## do the same for step monotonicity assumption
bad_iters_step <- lapply(unique(wdat$w_strata),function(ws){
  dtmp <- pr_U_tmp %>% filter(w_strata==ws)
  iters <- dtmp %>% select(starts_with("iter"))
  
  outer <- with(dtmp, 1- unlist(iters[score==0&gr==0,])-unlist(iters[score==2&gr==1,]))
  inners <- rbind(with(dtmp, unlist(iters[score==1&gr==1,])), with(dtmp, unlist(iters[score==1&gr==0,]))) %>% apply(2,min)

  as.numeric(outer > inners | outer < 0)
}) %>%
  do.call(rbind,.) %>%
  apply(2,function(x)any(x==1)) %>% which

## remove iterations not compatible with monotnoicity assumption
pr_U <- pr_U_tmp %>% select(-c(paste0("iter_",bad_iters)))

#P(U=u|A=1,W)
pr_u1 <- wdat_u %>% select(w_strata,sex,age0f,smoke,score0,score)%>%mutate(gr=1) %>% left_join(pr_U, by=c("sex","age0f","smoke","score","gr","score0","w_strata")) %>%
  select(w_strata,score,starts_with("iter_"))
#P(U=u|A=0,W)
pr_u0 <- wdat_u %>% select(w_strata,sex,age0f,smoke,score0,score)%>%mutate(gr=0) %>% left_join(pr_U, by=c("sex","age0f","smoke","score","gr","score0","w_strata")) %>%
  select(w_strata,score,starts_with("iter_"))

# function to calculate the joint distribution of the counterfactuals of U given sensitivity parameter rho
get_Ujoint <- function(rho){
  if(is.numeric(rho)){
  out <- pr_U %>% split(.$w_strata) %>%
    lapply(function(dd){
      dtmp <- data.frame(w_strata=dd$w_strata[1],u0 = c(0:2,0:2,0:2), u1=rep(0:2,each=3))
      iters <- dd %>% select(starts_with("iter_"))
      p11min <- cbind(with(dd, 1 - unlist(iters[score==0&gr==0,]) - unlist(iters[score==2&gr==1,])), 0) %>%
        apply(1,max)
      p11max <- cbind(with(dd, unlist(iters[score==1&gr==1,])), with(dd,unlist(iters[score==1&gr==0,]))) %>% apply(1,min)
      
      p00 <- unlist(with(dd, iters[score==0&gr==1,]))
      p22 <- with(dd, unlist(iters[score==2&gr==0,]))
      p11 <- (1-rho)*p11min + rho*p11max
      p01 <- with(dd, unlist(iters[score==1&gr==1,])-p11)
      p12 <- with(dd, unlist(iters[score==1&gr==0,]) - p11)
      p02 <- 1-p00-p22-p11-p01-p12
      p10 <- 0;p20<-0;p21<-0
      
      cbind.data.frame(dtmp %>% mutate(u0=factor(u0),u1=factor(u1)), rbind(p00,p10,p20,p01,p11,p21,p02,p12,p22))
      
    }) %>%
    do.call(rbind.data.frame,.)
  }
  if(rho=="step"){
    out <- pr_U %>% split(.$w_strata) %>%
      lapply(function(dd){
        dtmp <- data.frame(w_strata=dd$w_strata[1],u0 = c(0:2,0:2,0:2), u1=rep(0:2,each=3))
        iters <- dd %>% select(starts_with("iter_"))
        p11 <- with(dd, 1 - unlist(iters[score==0&gr==0,]) - unlist(iters[score==2&gr==1,]))
        p11 <- ifelse(p11>0,p11,NA)
        p00 <- unlist(with(dd, iters[score==0&gr==1,]))
        p22 <- with(dd, unlist(iters[score==2&gr==0,]))
        p01 <- with(dd, unlist(iters[score==1&gr==1,])-p11)
        p12 <- with(dd, unlist(iters[score==1&gr==0,]) - p11)
        p02 <- 0 # 1-p00-p22-p11-p01-p12
        p10 <- 0;p20<-0;p21<-0
        
        cbind.data.frame(dtmp %>% mutate(u0=factor(u0),u1=factor(u1)), rbind(p00,p10,p20,p01,p11,p21,p02,p12,p22))
        
      }) %>%
      do.call(rbind.data.frame,.)
  }
  out
}

# weights for each strata
pr_w_tmp <- surv_dat %>% mutate(sex=factor(sex),age0f=factor(age0f), smoke=factor(smoke),score0=factor(score0)) %>%
  left_join(wdat%>%mutate(sex=factor(sex))) %>%
  count(w_strata, .drop=F) %>%
  mutate(pr=n/sum(n))
iternames <- colnames(pr_U %>% select(starts_with("iter")))
pr_w <- cbind(pr_w_tmp[,"w_strata"],pr_w_tmp[,c(rep("pr",8000))]%>%`colnames<-`(c(paste0("iter_",1:8000)))) %>%
  select(w_strata,all_of(iternames))

#### Compute RMST ####

# simulate one realization from random effects vector with each estimated covariance matrix
set.seed(42)
simR <- sapply(1:nrow(theta),function(i){
  L_mat <- matrix(as.vector(unlist(theta[i,]%>%select(starts_with("L[")))),nrow=4)
  sds <- diag(as.vector(theta[i,] %>% select(starts_with("sd_"))))
  out <- MASS::mvrnorm(1,c(0,0,0,0),Sigma = sds%*%L_mat%*%t(L_mat)%*%sds) %>% `names<-`(c("r0","r1","r2","r3"))
  out
}) %>% t

# new parameter matrix, containing the above random effects
params <- (theta %>% select(starts_with(c("beta","gamma","zeta","xi","lambda"))) %>%
  cbind(.,simR))[-bad_iters,]

# set maximum time for RMST
max_time <- 15
# baseline interval at 'max_time'
max_interval <- max(which(bh_breaks<max_time))
nGL <- 3 # number of Gauss-Laguerre quadrature points
ww <- gauss.quad(nGL)$weights #GL weights
pp <- gauss.quad(nGL)$nodes #GL nodes

# to be passed to compute_rmst_const
argL <- list("M_t" = list("formula"=~nsk(time,knots=kn_x,Boundary.knots=Bk,integral=TRUE)+nsk(time,knots=kn_x,Boundary.knots=Bk,integral=TRUE):gr+
                            nsk(time,knots=kn_x,Boundary.knots=Bk,integral=TRUE):factor(score,levels=0:2)+
                            nsk(time,knots=kn_x,Boundary.knots=Bk,integral=TRUE):gr:factor(score,levels=0:2),"drop"=c(1)),
             "M_z" = list("formula"=~nsk(time,knots=kn_z,Boundary.knots=Bk,integral=TRUE),drop=c(1)),
             "G" = list("formula"=~factor(score,levels=0:2)+gr),
             "S_const" = list("formula"=~sex+factor(age0f,levels=0:2)+factor(smoke,levels=0:2)+factor(score0,levels=0:2)+factor(score,levels=0:2)+I(score==1&gr==1) + I(score==2&gr==1),
                              "drop"=1))
mtime <- 3

dat0 <- expand_grid(interval = 1:max_interval,
                    ind = 1:nGL) %>%
  mutate(intstart = bh_breaks[interval],
         intstop = pmin(bh_breaks[interval+1],max_time),
         intlen = intstop-intstart,
         time = (intlen/2)*pp[ind]+(intstop+intstart)/2,
         coeff = (intlen/2)*ww[ind],
         mtime=mtime)

## Compute the terms in the additive decomposition
source("compute_rmst_const.R")
t0 <- Sys.time()
beta_ac1 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=T,score=0),score=factor(0,levels=0:2),gr=1,
                               params=params, dat0=dat0, arg_list=argL,prW=pr_w,change=T)#8sek
beta_ac0 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=T,score=0),score=factor(0,levels=0:2),gr=0,
                               params=params, dat0=dat0, arg_list=argL,prW=pr_w,change=T)

tau_u1 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=T,score=c(0,0,0)),score=factor(0:2),gr=1,
                             params=params, dat0=dat0, arg_list=argL,prU=pr_u1,prW=pr_w,change=T)
tau_u0 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=T,score=c(0,0,0)),score=factor(0:2),gr=0,
                             params=params, dat0=dat0, arg_list=argL,prU=pr_u0,prW=pr_w,change=T)

tau_m11 <- compute_rmst_const(wdat,M_args=list(gr=1,ref=F,score=factor(0:2)),score=factor(rep(0,3),levels=0:2),gr=1,
                              params=params, dat0=dat0, arg_list=argL,
                         prU=pr_u1%>%rename(score_m:=score),prW=pr_w,change=T)#45sek
tau_m10 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(0:2)),score=factor(rep(0,3),levels=0:2),gr=1,
                              params=params, dat0=dat0, arg_list=argL,
                         prU=pr_u0%>%rename(score_m:=score),prW=pr_w,change=T)
tau_m00 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(0:2)),score=factor(rep(0,3),levels=0:2),gr=0,
                              params=params, dat0=dat0, arg_list=argL,
                         prU=pr_u0%>%rename(score_m:=score),prW=pr_w,change=T)
Sys.time()-t0

## Compute the delta-terms from additive interaction effects
t0 <- Sys.time()
tau_a1m1u1 <- compute_rmst_const(wdat,M_args=list(gr=1,ref=F, score=factor(0:2)),
                           score=factor(0:2),gr=1,
                           params=params,dat0=dat0,arg_list=argL,prU=pr_u1,prW=pr_w,change=T)
tau_a0m0u0 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F, score=factor(0:2)),
                           score=factor(0:2),gr=0,
                           params=params,dat0=dat0,arg_list=argL,prU=pr_u0,prW=pr_w,change=T)
Sys.time()-t0 

tau_a1m0ur <-  compute_rmst_const(wdat,M_args=list(gr=0,ref=F, score=factor(0:2)),
                            score=factor(c(0,0,0),levels=0:2),gr=1,
                            params=params,dat0=dat0,arg_list=argL,prU=pr_u0%>%rename(score_m:=score),prW=pr_w,change=T)

dat_taus <- cbind.data.frame(beta_ac0=unlist(beta_ac0),beta_ac1=unlist(beta_ac1),
                             tau_u0=unlist(tau_u0),tau_u1=unlist(tau_u1),
                             tau_m11=unlist(tau_m11),tau_m10=unlist(tau_m10),tau_m00=unlist(tau_m00),
                             tau_a0m0u0=unlist(tau_a0m0u0), tau_a1m1u1=unlist(tau_a1m1u1),
                             tau_a1m0ur=unlist(tau_a1m0ur))

## joint probability of counterfactuals of U with different choices of rho
pr_uj00 <- get_Ujoint(0)
pr_uj25 <- get_Ujoint(.25)
pr_uj50 <- get_Ujoint(.5)
pr_uj75 <- get_Ujoint(.75)
pr_uj100 <- get_Ujoint(1)
tau_mu_rho25 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(rep(0:2,each=3))),gr=1,score=factor(rep(0:2,3)),change=T,
                             params=params,dat0=dat0,arg_list=argL,
                             prU=pr_uj25%>%rename("score_m":=u0,"score":=u1),prW=pr_w)
tau_mu_rho50 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(rep(0:2,each=3))),gr=1,score=factor(rep(0:2,3)),change=T,
                             params=params,dat0=dat0,arg_list=argL,
                             prU=pr_uj50%>%rename("score_m":=u0,"score":=u1),prW=pr_w)
tau_mu_rho75 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(rep(0:2,each=3))),gr=1,score=factor(rep(0:2,3)),change=T,
                             params=params,dat0=dat0,arg_list=argL,
                             prU=pr_uj75%>%rename("score_m":=u0,"score":=u1),prW=pr_w)
tau_mu_rho0 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(rep(0:2,each=3))),gr=1,score=factor(rep(0:2,3)),change=T,
                            params=params,dat0=dat0,arg_list=argL,
                            prU=pr_uj00%>%rename("score_m":=u0,"score":=u1),prW=pr_w)
tau_mu_rho100 <- compute_rmst_const(wdat,M_args=list(gr=0,ref=F,score=factor(rep(0:2,each=3))),gr=1,score=factor(rep(0:2,3)),change=T,
                              params=params,dat0=dat0,arg_list=argL,
                              prU=pr_uj100%>%rename("score_m":=u0,"score":=u1),prW=pr_w)

dat_taus <- dat_taus %>% mutate(tau_mu_rho0=unlist(tau_mu_rho0),
                                tau_mu_rho25 = unlist(tau_mu_rho25),
                                tau_mu_rho50 = unlist(tau_mu_rho50),
                                tau_mu_rho75 = unlist(tau_mu_rho75),
                                tau_mu_rho100 = unlist(tau_mu_rho100))

#### Compute TE, DE and IE

delta_rho0 <- with(dat_taus,tau_mu_rho0 - tau_u1 - tau_a1m0ur + beta_ac1)
delta_rho25 <- with(dat_taus,tau_mu_rho25 - tau_u1 - tau_a1m0ur + beta_ac1)
delta_rho50 <- with(dat_taus,tau_mu_rho50 - tau_u1 - tau_a1m0ur + beta_ac1)
delta_rho75 <- with(dat_taus,tau_mu_rho75 - tau_u1 - tau_a1m0ur + beta_ac1)
delta_rho100 <- with(dat_taus,tau_mu_rho100 - tau_u1 - tau_a1m0ur + beta_ac1)

Delta_DE <- with(dat_taus, tau_a0m0u0 - tau_u0 - tau_m00 + beta_ac0)
Delta_IE <- with(dat_taus, tau_a1m1u1 - tau_u1 - tau_m11 + beta_ac1)

DE <- with(dat_taus, beta_ac0 - beta_ac1 + tau_u1 - tau_u0 +tau_m10-tau_m00)
IE <- with(dat_taus, tau_m11 - tau_m10)

eff_de_ie <- rbind(
  data.frame(eff="IE",
             rho = c(0,25,50,75,100),
             mean = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)mean(IE+Delta_IE-xx)),
             sd = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                         function(xx)sd(IE+Delta_IE-xx)),
             q025 = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)as.numeric(quantile(IE+Delta_IE-xx,.025))),
             q975 = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)as.numeric(quantile(IE+Delta_IE-xx,.975)))),
  data.frame(eff="DE",
             rho = c(0,25,50,75,100),
             mean = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)mean(DE-Delta_DE+xx)),
             sd = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                         function(xx)sd(DE-Delta_DE+xx)),
             q025 = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)as.numeric(quantile(DE-Delta_DE+xx,.025))),
             q975 = sapply(list(delta_rho0,delta_rho25,delta_rho50,delta_rho75,delta_rho100),
                           function(xx)as.numeric(quantile(DE-Delta_DE+xx,.975)))))


## Total effect
te <- with(dat_taus,tau_a1m1u1-tau_a0m0u0)
te_summ <- c("Est"=mean(te),"sd"=sd(te),"median"=median(te),quantile(te,c(.025,.975)))

# simple nonparametric
te_nonp <- survRM2::rmst2(surv_dat$stop_time,surv_dat$status,surv_dat$gr, tau=15)


## Throw away posterior samples giving negative effects
dd_c <- DE - Delta_DE + delta_rho50
ii_c <- IE + Delta_IE - delta_rho50
incl_c <- which(dd_c>0 & ii_c > 0)

mysumm <- function(x)c("mean"=mean(x), sd=sd(x), q025=as.numeric(quantile(x,.025)),q975=as.numeric(quantile(x,.975)))
eff_positive <- sapply(list(dd_c[incl_c],ii_c[incl_c],dd_c[incl_c]+ii_c[incl_c],ii_c[incl_c]/(dd_c[incl_c]+ii_c[incl_c])),mysumm) %>%
  `colnames<-`(c("DE","IE","TE","prc_mediated"))
