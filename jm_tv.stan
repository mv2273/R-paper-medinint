// A joint model using the longitudinal outcome as time-varying covariate and shared random effects

// random effects parametrization for longitudinal model copied from brms
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
// function to evaluate g(M(t))
  /*
  * which:    the function to be applied on the mediator trajectory
  * constant: additive term produced by all time-constant variables
  * X_t:      population level time-varying covariates
  * beta_t:   associated parameter vector
  * Z_t:      time-varying covariates for random effects
  * r:        random effects
  * time:     time at evaluation; not used in current value parameterisation
  */
  real my_gfun(int which, real constant, row_vector X_t, vector beta, row_vector Z_t, row_vector r, real time){
    real out;
    if(which==1){// current value
      out = constant + X_t*beta + Z_t*r';
    }
    if(which==2){ // integral
      out = constant*time + X_t*beta + Z_t*r';
    }
    if(which==3){//average; X_t and Z_t should evaluate to integrals
      out = constant + (X_t*beta + Z_t*r')/time;
    }
    if(which==4){//shared random effects only
      out = 0;
    }
    if(which==5){// integral over the latest x years
      out = constant*time + X_t*beta + Z_t*r';
    }
    if(which==6){ // value at year x; the expression is same as one but here, X_t and Z_t are required to be constructed with a constant time-argument!
      out = constant + X_t*beta + Z_t*r';
    }
    if(which==7){ // change compared to baseline
      out = X_t*beta + Z_t*r';
    }
    if(which==8){ // relative change
      out = (X_t*beta + Z_t*r')/constant;
    }
    return out;
  }
}
data {
  // required integers
  int<lower=0> N; // n.o. data rows
  int<lower=0> N_i; // n.o. individuals
  int<lower=0> n_c; // n.o. covariates
  int<lower=0> n_t; // n.o. time-varying covariates
  int<lower=0> n_zc; // n.o. time-constant random effects
  int<lower=0> n_zt; // n.o. time-varying random effects
  int<lower=1> n_bh_breaks; // n.o. break points for piecewise constant baseline hazard
  int<lower=0> n_s; // n.o. covariates for survival submodel
  int<lower=0> n_gl_evals; // number of times the mediator value must be evaluated during GL; should equal nGL*sum(interval)
  int<lower=0> n_zeta; // how many parameters will be associated with g(M(.))
  int nGL; // n.o. Gauss-Legendre quadrature points per baseline-hazard piece
  int<lower=0> nsr; //n.o. shared random effects
  
  // longitudinal data
  vector[N] y; // longitudinal outcome
  matrix[N_i,n_c] X_c; // design matrix holding the time-independent population level covariates for longitudinal submodel
  matrix[N,n_t] X_t; // time-varying variables 
  array[N] int id; // indexing individuals corresponding to the outcome values y
  matrix[N_i, n_zc] Z_c; // dense design matrix for time-independent random effects; i'th row corresponds to individual id[i]
  matrix[N, n_zt] Z_t; // time-varying variables
  
  // survival data
  array[N_i] int id_surv; // id vector for survival data
  matrix[N_i, n_s] X_s; // design matrix for survival submodel
  vector[N_i] stop_time; // stopping times
  vector[N_i] status; // 1 = event; 0 = censoring
  array[N_i] int interval; // in which piecewise h_0 interval exit occurred
  vector[N_i] treatment;
  
  // data matrices needed to evaluate M(t) at the exit time and
  // the evaluation times necessary for numerical integration
  matrix[N_i,n_zeta] G; // rows should correspond to 'design vectors'; M_i(t)*G_i*zeta (1 x 1)*(1 x n_zeta)*(n_zeta x 1)
  matrix[n_gl_evals, n_t] X_t_gl; // design matrix for evaluating M(.) at points as per gauss-legendre rules (time-varying part)
  matrix[n_gl_evals, n_zt] Z_t_gl; // time-varying random effects
  
  vector[n_gl_evals] pGL; // GL evaluation points
  vector[nGL] wGL; // weights
  
  matrix[N_i, n_t] X_t_exit; // design matrix at exit times (time-varying)
  matrix[N_i, n_zt] Z_t_exit; // random effects part
  
  vector[n_bh_breaks] bh_breaks; // breakpoints for piecewise constant baseline hazard
  
  real<lower=0> lag; // needed if using restricted integral
  array[nsr] int shared_ranefs; // which random effects are shared
  int<lower=1> gfun_type; // which function is applied to mediator trajectory
  
}
transformed data{
  int n_zz = n_zc + n_zt; // n.o. total random effects params
  //vector[N_i] zeros_rw = rep_vector(0,N_i);
}

parameters {
  // parameters for longitudinal submodel
  vector[n_c] beta_c; // population level parameters for longitudinal submodel; corresponding to X_c
  vector[n_t] beta_t; // corresponding to time-varying covariates 
  vector[n_s] gamma; // parameters for survival submodel
  vector[n_zeta] zeta; // association parameters
  vector[nsr] xi; // parameters for shared random effects for survival submodel
  vector<lower=0>[n_bh_breaks] lambda0;
  vector<lower=0>[n_bh_breaks] lambda1;
  real<lower=0> sigma; // standard deviation for the measured mediators
  vector<lower=0>[n_zz] sd_r;  // group-level standard deviations
  matrix[n_zz, N_i] U;  // standardized group-level effects
  cholesky_factor_corr[n_zz] L;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_i, n_zz] r;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r = scale_r_cor(U, sd_r, L);
  lprior += lkj_corr_cholesky_lpdf(L | 1);
  for(i in 1:n_zz){
    lprior += cauchy_lpdf(sd_r[i] | 0,10);
  }
  lprior += cauchy_lpdf(sigma | 0, 10); //prior for the residual variance
}

model {
  // longitudinal submodel
  vector[N] mu; // population level contributions to the longitudinal linear predictors
  vector[N_i] mu_tmp = X_c*beta_c; // add time-constant parts here first
  for(i in 1:N_i){
    mu_tmp[i] += Z_c[i,:]*r[i,1:n_zc]';
  }
  // add time-varying terms
  for(i in 1:N){
    mu[i] = mu_tmp[id[i]] + X_t[i,:]*beta_t + Z_t[i,:]*r[id[i],(n_zc+1):n_zz]';
  }
  y ~ normal(mu,sigma); // assume conditional normal distribution for outcome
  target += std_normal_lpdf(to_vector(U));
  target += lprior;
  
  // survival submodel
  vector[N_i] Xs_i = X_s*gamma;
  if(nsr>0){ // add random effects from the longitudinal submodel?
    Xs_i += r[:,shared_ranefs]*xi;
  }

  vector[N_i] M_const;
  for(i in 1:N_i){
    M_const[i] = X_c[i,:]*beta_c + Z_c[i,:]*r[i,1:n_zc]';
  }
    // compute log-hazards at event times
  vector[N_i] log_hi;
  for(i in 1:N_i){
    log_hi[i] = log(lambda0[interval[i]]*(1-treatment[i]) + lambda1[interval[i]]*treatment[i]) + Xs_i[i] + 
    my_gfun(gfun_type,M_const[i], X_t_exit[i,:],beta_t, Z_t_exit[i,:], r[i,(n_zc+1):n_zz], stop_time[i])*G[i,:]*zeta;
  }
// numerical integration of the hazard function
  int count = 1; // counter used to index correct row from ...gl-matrices
  for(i in 1:N_i){
    int int_tmp = interval[i];
    real H_i = 0; // cumulative hazard for individual i; summing over baseline pieces
    for(k in 1:int_tmp){
      real intlen; //length of the time interval in consideration
      if(k==interval[i]){
        intlen = stop_time[i] - bh_breaks[k];
      }else{
        intlen = bh_breaks[k+1] - bh_breaks[k];
      }
      real j_sum = 0; 
      for(j in 1:nGL){ 
        j_sum += wGL[j] * exp(my_gfun(gfun_type,M_const[i], X_t_gl[count,:], beta_t, Z_t_gl[count,:], r[i,(n_zc+1):n_zz], pGL[count]) * G[i,:] * zeta);
        count += 1;
      }
      H_i += (lambda0[k]*(1-treatment[i]) + lambda1[k]*treatment[i]) * (intlen/2) * j_sum;
    }
    target += status[i]*log_hi[i] - exp(Xs_i[i])*H_i; // add log-likelihood contribution from individual i
  }
  // priors for model params
  beta_c ~ normal(0,5);
  beta_t ~ normal(0,5);
  gamma ~ normal(0,5);
  zeta ~ normal(0,5);
  xi ~ normal(0,5);
  lambda0 ~ gamma(.5,.5);
  lambda1 ~ gamma(.5,.5);
}
