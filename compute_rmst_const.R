## A function to compute RMST for a given strata
# w_strata: which strata of backgraound variables, as specified in wdat_all (defined outside this function), is used
# M_args:   gr: to which level is treatment set when computing M; 
#           ref: whether the reference trajectory only is needed; if TRUE, the 'gr' makes no difference
# score:    to which level is the treatment-dependent confounder set at
# gr:       to which level the treatment is set at

compute_rmst_const <- function(w_dat,M_args,score,gr,params,dat0,arg_list,prU=NULL,prW=NULL,change=F){
  if(!change){
    arg_M_c <- arg_list[["M_c"]]
  }
  arg_M_t <- arg_list[["M_t"]]
  arg_M_z <- arg_list[["M_z"]]
  arg_G <- arg_list[["G"]]
  arg_S_const <- arg_list[["S_const"]]
  
  dat_for_M <- w_dat %>% mutate(gr=M_args[["gr"]]) %>%
    expand_grid(score=M_args[["score"]])
  score_dat <- data.frame(score=score,score_m=M_args[["score"]])
  dat_w <- w_dat %>% expand_grid(score_dat)
  dat_w$gr <- gr
  
  u_int_arg <- case_when(
    is.data.frame(prU)&&"score"%in%colnames(prU)&&(!"score_m"%in%colnames(prU)) ~ "score",
    is.data.frame(prU)&&"score_m"%in%colnames(prU)&&(!"score"%in%colnames(prU)) ~ "score_m",
    is.data.frame(prU)&&"score"%in%colnames(prU)&&"score_m"%in%colnames(prU) ~ "both",
    TRUE ~ "none"
  )
  w_int_arg <- case_when(
    is.data.frame(prW)&&"w_strata"%in%colnames(prW) ~ "W",
    TRUE ~ "none"
  )
  
  # repeats rows or columns of a matrix
  pile <- function(mat,n,op,expr=NULL){
    switch(op,
           "expand_rows"=mat[rep(1:nrow(mat),each=n),],
           "expand_cols"=mat[,rep(1:ncol(mat),each=n)],
           "stack_long"=mat[rep(1:nrow(mat),n),],
           "stack_wide"=mat[,rep(1:ncol(mat),n)])
  }
  
  ## Calculates g(M)
  gMt <- function(t,data,change,ref=F,custom=NULL){
    #browser()
    stopifnot(length(unique(t))==1)
    nd <- data %>% mutate(time=unique(t))
    if(ref){ # if calculating the reference quantities, m will be fixed a constant zero
      out <- matrix(0,nrow=nrow(nd),ncol=nrow(params))
    }else if(is.function(custom)){
      out <- lapply(1:nrow(params), function(i)custom(nd$time)) %>% do.call(cbind,.)
    }else{
      if(length(t)==1 & nrow(data==1)){
        tmpfun <- function(...)c(...)
      }else{
        tmpfun <- function(...)cbind(...)
      }
      if(change){ # if interested in change, don't use the constant terms
        out <- my.D(arg_M_t,nd)%*%t(params%>%select(starts_with("beta_t")))+
          my.D(arg_M_z,nd) %*% t(model.matrix(~-1+r1+r2+r3,params))
      }else{
        out <- tmpfun(my.D(arg_M_c,nd),my.D(arg_M_t,nd))%*%t(params%>%select(starts_with("beta_c"),starts_with("beta_t")))+
          my.D(arg_M_z,nd) %*% t(model.matrix(~-1+r0+r1+r2+r3,params))
      }
    }
    out
  }
  # evaluate g(M(t)) at each 'inner time' given in dat0
  # gMt() evaluates mediator functional for each row in dat_for_M at each inner time in dat0,
  # in the order that first is all times from dat0 for the first row of dat_for_M, and so on

  zeta_params <- params %>% select(starts_with("zeta")) %>% as.matrix
  gM_G_zeta <- gMt(dat0$mtime, dat_for_M, ref=M_args[["ref"]],custom=M_args[["custom_fun"]],
                   change=change)*my.D(arg_G,pile(dat_w,1,"expand_rows"))%*%t(zeta_params)
  #this should hold, for each row of dat_w, 8000 posterior samples of g(M(t)) G \zeta evaluated at times specified in dat0
  # ordered such that the w_stratas change slowest
  
  # likewise the X\gamma + r_0\xi parts of the linear predictor and the baseline-hazard values
  # in X_gamma, one row per dat_w row since these are time-constant. Later, broadcast this in order to perform matrix operations with gM_G_zeta
  X_gamma <- my.D(arg_S_const,dat_w) %*% t(params%>%select(starts_with("gamma"))) + (params[["xi[1]"]] * params[["r0"]])
  
  br2 <- c(bh_breaks[bh_breaks<max_time],max_time)
  intlens <- br2[2:length(br2)]-br2[1:(length(br2)-1)]
  chaz_tmp <- sapply(dat0$time, function(tt){
    int_tmp <- max(which(bh_breaks<tt))
    if(int_tmp == 1){
      intl_tmp <- tt
    }else{
      intl_tmp <- c(intlens[1:(int_tmp-1)], tt-bh_breaks[int_tmp])
    }
    bh_names_tmp <- paste0("lambda",gr,"[",1:int_tmp,"]")
    apply(pile(matrix(as.numeric(intl_tmp),ncol=1),nrow(params),"expand_cols") * t(params[,bh_names_tmp]),2,sum)
  }) %>% t %>%
    `colnames<-`(paste0("iter_",1:ncol(.))) %>%
    `rownames<-`(NULL)
  const <- pile(exp(X_gamma + gM_G_zeta), nrow(dat0), "expand_rows")
  S_t <- exp(- const*pile(chaz_tmp,nrow(dat_w),"stack_long")) %>%
    `colnames<-`(paste0("iter_",1:ncol(.)))
  
  # attach these values to the GL-skeleton, and sum over the GL-rule
  result <- cbind(expand_grid(dat_w,dat0),S_t) %>%
    group_by(w_strata,score,score_m) %>%
    summarise(across(starts_with("iter_"), ~sum(coeff * .x)),.groups="drop")%>%
    ungroup
  # at this point, we should have the 8k posterior samples of RMST's from each given w_strata--score-combination
  # we may further summarise over P(u|w) and p(w) if we like
  
  # checks that the intended summin over -probabilities sum up to 1 (tl allowing for some numerical instability)
  check_prob_summ <- function(pr_mat, cols, iter_prefix,tol=1e-6){
    mm <- as.data.frame(pr_mat) %>% group_by(pick(all_of(cols))) %>%
      summarise(across(starts_with(iter_prefix), ~sum(.x)),.groups="drop") %>%
      select(starts_with(iter_prefix)) %>%
      unlist
    all(abs(mm-1)<tol)
  }
  # performs summing over
  sum_over <- function(res,pr, id_col, summ_col, iter_prefix, ind){
    #if(summ_col=="w_strata"){
    #  browser()
    #}
    if(!check_prob_summ(pr_mat=pr, cols=id_col, iter_prefix=iter_prefix)){
      warning(paste0("In summing over ",ind,": probabilities do not seem to add to 1; summing not performed"))
      out <- res
    }else{
      pr_new <- res%>%select(all_of(c(id_col,summ_col))) %>% 
        left_join(pr%>%select(all_of(c(id_col,summ_col)),starts_with(iter_prefix)), by=c(id_col,summ_col)) %>%
        select(starts_with(iter_prefix))
      out <- cbind(pr %>% select(all_of(id_col)),
                   (res %>% select(starts_with("iter_")))*pr_new) %>%
        group_by(pick(all_of(id_col))) %>%
        summarise(across(starts_with("iter_"), ~sum(.x)))
    }
    out
  }
  
  # sum over U
  result <- switch(
    u_int_arg,
    "score"=sum_over(res=result,pr=prU,id_col="w_strata",summ_col="score",iter_prefix="iter",ind="U"),
    "score_m"=sum_over(res=result,pr=prU,id_col="w_strata",summ_col="score_m",iter_prefix = "iter_",ind="U(M)"),
    "both" = sum_over(res=result,pr=prU,id_col="w_strata",summ_col=c("score","score_m"),iter_prefix="iter_",ind="U(joint)"),
    "none"=result
  )

  #sum over W
  result <- switch(
    w_int_arg,
    "W" = sum_over(res=result,pr=prW,id_col=NULL, summ_col="w_strata", iter_prefix="iter_", ind="W"),
    "none"=result
  )
  
  return(result) # return data vector holding the posterior sample of the RMST
}
