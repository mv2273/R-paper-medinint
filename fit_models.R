library(tidyverse)
library(statmod)
library(splines2)
library(readr)
library(rstan)

nsk <- function(...)splines2::nsk(...)

build_Z <- function(formula,data,block_id_var){
  allvars <- all.vars(formula)
  charvars <- allvars[which(data %>% ungroup %>% select(all_of(allvars)) %>% 
                              sapply(is.character))]
  dat_tmp <- data
  if(length(charvars)>0){
    for(vr in charvars){
      dat_tmp[[vr]] <- factor(dat_tmp[[vr]], levels=sort(unique(dat_tmp[[vr]])))
    }
  }
  sparse_mat <- split(dat_tmp,f=data[[block_id_var]]) %>%
    lapply(.,function(dd)model.matrix(formula,dd)) %>%
    do.call(Matrix::bdiag,.) 
  sparse_mat
}
#### Read datasets ####
#assign survival dataset to 'surv_dat'
#assign lngitudinal dataset to 'long_dat'

#### fit multinomial model for the time-dependent confounder ####
standat_mult <- list(
  y = surv_dat$score + 1,
  X = model.matrix(~gr+sex+factor(age0f)+factor(smoke)+factor(score0),surv_dat),
  N=nrow(surv_dat),
  n_cat = length(unique(surv_dat$score)),
  n_psi=ncol(model.matrix(~gr+sex+factor(age0f)+factor(smoke)+factor(score0),surv_dat))
)
mult_fit <- stan(file="multinomial.stan",
                 data=standat_mult,
                 iter=4000,chains=4,cores=4)
# save results
write_csv2(as.data.frame(mult_fit),file="posterior_u.csv")

#### fit model with current change parameterisation ####
source("build_standata.R")
stan_data_cc <- build_standata(data_long=long_dat, data_surv=surv_dat, long_var = "bmi_c",
                            formula_list = list("long_c" = ~gr*factor(score)+sex+factor(age0f)+factor(smoke)+factor(score0),
                                                "long_t" = ~nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10))+
                                                  nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10)):(gr+factor(score))+
                                                  nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10)):gr:factor(score),
                                                "long_zc" = ~1,
                                                "long_zt" = ~nsk(time,knots=c(1,5),Boundary.knots=c(0,10)),
                                                "surv_c" = ~sex+factor(age0f)+factor(smoke)+factor(score0)+factor(score)+I(score==1&gr==1)+I(score==2&gr==1),
                                                "G" = ~factor(score)+gr,
                                                "pred_Xt_gl" = ~my.nsk(time,kn=c(1,3,5),bk=c(0,10),tf="current")+
                                                  my.nsk(time,kn=c(1,3,5),bk=c(0,10),tf="current"):(gr+factor(score))+
                                                  my.nsk(time,kn=c(1,3,5),bk=c(0,10),tf="current"):gr:factor(score),
                                                "pred_Zt_gl" = ~my.nsk(time,kn=c(1,5),bk=c(0,10),tf="current"),
                                                "pred_Xt_exit" = ~my.nsk(stop_time,kn=c(1,3,5),bk=c(0,10),tf="current")+
                                                  my.nsk(stop_time,kn=c(1,3,5),bk=c(0,10),tf="current"):(gr+factor(score))+
                                                  my.nsk(stop_time,kn=c(1,3,5),bk=c(0,10),tf="current"):gr:factor(score),
                                                "pred_Zt_exit" = ~my.nsk(stop_time,kn=c(1,5),bk=c(0,10),tf="current")),
                            gfun="change",
                            lagtime=0,
                            kn_x = c(1,3,5),
                            kn_z = c(1,5),
                            Bk = c(0,10))

joint_fit_cc <- stan(file="jm_tv.stan",
                  data=stan_data_cc$standata,
                  iter=4000, chains=4, cores=4)
# save results
write_csv2(as.data.frame(joint_fit_cc) %>%select(starts_with(c("beta","gamma","zeta","xi","lambda","sigma","sd_r","r[","L["))), 
           "posterior_cc.csv")

#### fit joint model with three-year legacy parameterisation ####
source("build_standata_const.R")
stan_data_const <- build_standata_const(data_long=long_dat, data_surv=surv_dat%>%mutate(time_tmp = 3), long_var = "bmi_c",
                            formula_list = list("long_c" = ~gr*factor(score)+sex+factor(age0f)+factor(smoke)+factor(score0),
                                                "long_t" = ~nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10))+
                                                  nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10)):(gr+factor(score))+
                                                  nsk(time, knots=c(1,3,5), Boundary.knots=c(0,10)):gr:factor(score),
                                                "long_zc" = ~1,
                                                "long_zt" = ~nsk(time,knots=c(1,5),Boundary.knots=c(0,10)),
                                                "surv_c" = ~sex+factor(age0f)+factor(smoke)+factor(score0)+factor(score)+I(score==1&gr==1)+I(score==2&gr==1),
                                                "G" = ~factor(score)+gr,
                                                "pred_Xt_exit" = ~nsk(time_tmp,kn=c(1,3,5),Boundary.knots=c(0,10),integral=TRUE)+
                                                  nsk(time_tmp,knots=c(1,3,5),Boundary.knots=c(0,10),integral=TRUE):(gr+factor(score))+
                                                  nsk(time_tmp,knots=c(1,3,5),Boundary.knots=c(0,10),integral=TRUE):gr:factor(score),
                                                "pred_Zt_exit" = ~nsk(time_tmp,kn=c(1,5),Boundary.knots=c(0,10), integral=TRUE)),
                            kn_x = c(1,3,5),
                            kn_z = c(1,5),
                            Bk = c(0,10),
                            gfun="integral")

joint_fit_leg <- stan(file="jm_const.stan",
                  data=stan_data_const$standata,
                  iter=4000, chains=4, cores=4)
# save results
write_csv2(as.data.frame(joint_fit_leg), "posterior_leg.csv")
