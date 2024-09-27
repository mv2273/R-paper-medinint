## Build stan data for jm using time-constant functional of longitudinal outcome
# Arguments:
# data_long:    data containing the longitudinal follow up
# data_surv:    data containing the survival outcome
# long_var:     name of the longitudinal outcome in data_long
# formula_list: list of formulas;
#               long_c: formula for constant covariates in longitudinal submodel
#               surv_c: formula for constant covariates in survival submodel
#               long_zc & long_zt: time-constant and time-varying random effects for longitudinal model
#               G:      formula providing main effect + interactions containing g(M(t))
#               pred_Xt_exit: formula to evaluate g(M(t)) at exit times
#               pred_Zt_exit: random effects part
# bh_breaks:    break points for baseline hazard; if NULL, placed at .1, .2, ..., .9 quantiles of observed event times
# nGL:          number of Gauss-Legendre quadrature points per baseline hazard piece
# id_var:       name of the variable containing individual identifiers
# stop_var:     name of the variable containing stopping times in data_surv
# status_var:   name of the variable containing status at exit
# kn_x:         inner knots for population level spline
# kn_z:         inner knots for individual level spline
# Bk:           boundary knots for nsk
# shared_r:     intger vector indicating which ranefs are shared
# gfun:         how to transform the mediator trajectory into the survival submodel
build_standata_const <- function(data_long, data_surv, long_var, formula_list,
                           kn_x,kn_z,Bk,gfun,
                           bh_breaks=NULL,
                           id_var="id", stop_var="stop_time", status_var="status",
                           shared_r=1,time_var="time"){
  require(tidyverse)
  require(statmod)
  require(splines2)
  ## Check varnames and that the datasets match and all...
  stopifnot(all(data_surv[[id_var]]%in%data_long[[id_var]]))
  new_ids <- cbind.data.frame(sort(data_surv[[id_var]]), 1:length(data_surv[[id_var]])) %>% `colnames<-`(c(id_var,"id_new"))
  if(is.null(Bk)){
    Bk <- c(0,10)
  }
  
  ## Set up datasetes
  ds_tmp <- data_surv%>% mutate(status=data_surv[[status_var]], stop_time = data_surv[[stop_var]])
  dl_tmp <- data_long %>% mutate(time=data_long[[time_var]])
  
  if(is.null(bh_breaks)){
    evtimes <- ds_tmp %>% filter(status==1) %>% select(stop_time)%>%unlist
    bh_breaks <- c(0,quantile(evtimes, seq(.1,.9,.1)))
  }
  
  if("id_new"%in%colnames(data_surv)){
    ds_tmp <- ds_tmp %>% select(-id_new)
    message("original 'id_new' column replaced in survival data")
  }
  
  ds <- ds_tmp %>%
    left_join(new_ids,by=id_var) %>%
    arrange(id_new) %>%
    rowwise %>%
    mutate(end_interval = max(which(bh_breaks<stop_time))) %>%
    ungroup
  
  if("id_new"%in%colnames(data_long)){
    dl_tmp <- dl_tmp %>% select(-id_new)
    message("original 'id_new' column replaced in longitudinal data")
  }
  
  dl <- dl_tmp %>%
    left_join(new_ids,by=id_var) %>%
    arrange(id_new,time)
  
  source("my.nsk.R")
  gfun_num <- case_when(gfun=="integral" ~ 2,
                        gfun=="current" ~ 1,
                        gfun=="avg" ~ 3,
                        gfun=="shared_only" ~ 4,
                        gfun=="cut_integral" ~ 5,
                        gfun=="change" ~ 7)
  
  ## Set up design matrices
  X_s <- as.matrix(model.matrix(formula_list[["surv_c"]], ds)[,-1])
  X_c <- model.matrix(formula_list[["long_c"]], ds)
  X_t <- model.matrix(formula_list[["long_t"]], dl)[,-1]
  Z_c <- model.matrix(formula_list[["long_zc"]], ds)
  Z_t <- model.matrix(formula_list[["long_zt"]], dl)[,-1]
  
  G <- model.matrix(formula_list[["G"]], ds)
  
  X_t_exit <- model.matrix(formula_list[["pred_Xt_exit"]], ds)[,-1]
  Z_t_exit <- model.matrix(formula_list[["pred_Zt_exit"]], ds)[,-1]
  

  stan_data <- list(
    N = nrow(dl), # n.o. data rows
    N_i = nrow(ds), #n.o. individuals
    n_c = ncol(X_c), # n.o. covariates
    n_t = ncol(X_t), # n.o. time-varying covariates
    n_zc = ncol(Z_c), # n.o. time-constant random effects
    n_zt = ncol(Z_t), # n.o. time-varying random effects
    n_bh_breaks = length(bh_breaks), # n.o. break points for piecewise constant baseline hazard
    n_s = ncol(X_s), # n.o. covariates for survival submodel
    n_zeta = ncol(G), # how many parameters will be associated with g(M(.))
    nsr = length(shared_r), #n.o. shared random effects
    
    # longitudinal data
    y = dl[[long_var]], # longitudinal outcome
    X_c = X_c, # design matrix holding the time-independent population level covariates for longitudinal submodel
    X_t=X_t, # time-varying variables 
    id = dl[["id_new"]], # indexing individuals corresponding to the outcome values y
    Z_c = Z_c, # dense design matrix for time-independent random effects; i'th row corresponds to individual id[i]
    Z_t = Z_t, # time-varying variables
    
    # survival data
    id_surv = ds[["id_new"]], # id vector for survival data
    X_s = X_s, # design matrix for survival submodel
    stop_time = ds[["stop_time"]], # stopping times
    status = ds[["status"]], # 1 = event; 0 = censoring
    interval = ds[["end_interval"]], # in which piecewise h_0 interval exit occurred
    treatment = ds[["gr"]],
    
    # data matrices needed to evaluate M(t) at the exit time and
    # the evaluation times necessary for numerical integration
    G = G, # rows should correspond to 'design vectors'; M_i(t)*G_i*zeta (1 x 1)*(1 x n_zeta)*(n_zeta x 1)
    
    X_t_exit = X_t_exit, # design matrix at exit times (time-varying)
    Z_t_exit = Z_t_exit, # random effects part
    
    bh_breaks = bh_breaks, # breakpoints for piecewise constant baseline hazard
    shared_ranefs = array(shared_r), # which random effects are shared
    gfun_type=gfun_num
  )
  
  out <- list("standata"=stan_data,
              "parnames" = list(
                "beta_c"=colnames(X_c),
                "beta_t"=colnames(X_t),
                "gamma"=colnames(X_s),
                "zeta" = paste0(colnames(G),":gM")
              ),
              "id_key"=new_ids,
              nskfun = function(...)my.nsk(tf=gfun,...))
  return(out)
}
