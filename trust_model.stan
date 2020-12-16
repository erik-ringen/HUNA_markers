functions{
  /* compute monotonic effects, function lifted from brms package
  * Args:
    *   scale: a simplex parameter
  *   i: index to sum over the simplex
  * Returns:
    *   a scalar between 0 and 1
  */
    real mo(vector scale, int i) {
      if (i > 0) return sum(scale[1:i]);
      else return 0;
    }
}

data{
  int N_obs;
  int N_id;
  int N_alter;
  int N_eth;
  int N_village;
  int trust[N_obs]; // trust, mistrust, or neutral?
  int id[N_obs];
  int alter_id[N_obs];
  int eth[N_obs];
  int village[N_obs];
  int edu[N_obs];
  int sex[N_obs];
  int sex_alter[N_obs];
  int age[N_obs]; // age category
  int magic[N_obs];
  int christian[N_obs];
  //// id-level covariates 
  int edu_id[N_id]; // education for each id
  int age_id[N_id]; // age for each id
}

parameters{
  vector[2] b0; // intercept
  matrix[7,2] b; // fixed effects
  simplex[max(edu)] scale_edu[2]; //  monotonic effect of education 
  
  matrix[2,N_id] id_z; // random effects for individual (ego), unscaled, uncentered, and uncorrelated
  matrix[2,N_eth] eth_z; // random effects for ethnicity
  matrix[2,N_village] village_z; // random effects for village
  matrix[2,N_alter] alter_z; // random effects for alter
  vector<lower=0>[2] sigma_eth; // variance components
  vector<lower=0>[2] sigma_village;
  vector<lower=0>[2] sigma_id;
  vector<lower=0>[2] sigma_alter;
  cholesky_factor_corr[2] L_id; // correlations between random effects
  cholesky_factor_corr[2] L_eth; // 
  cholesky_factor_corr[2] L_village; // 
  cholesky_factor_corr[2] L_alter; //
  
  ordered[max(edu)] c_edu; // cutpoints for education
  real<lower=0,upper=1> age_p; // prob of being adult
}

transformed parameters{
  // Now, we need to scale and correlate the random effects
  matrix[N_id,2] id_v;
  matrix[N_eth,2] eth_v;
  matrix[N_village,2] village_v;
  matrix[N_alter,2] alter_v;

  id_v = (diag_pre_multiply(sigma_id, L_id) * id_z)';
  eth_v = (diag_pre_multiply(sigma_eth, L_eth) * eth_z)';
  village_v = (diag_pre_multiply(sigma_village, L_village) * village_z)';
  alter_v = (diag_pre_multiply(sigma_alter, L_alter) * alter_z)';
}

model{
  // priors
  b0 ~ std_normal();
  to_vector(b) ~ std_normal();
  to_vector(id_z) ~ std_normal();
  to_vector(eth_z) ~ std_normal();
  to_vector(village_z) ~ std_normal();
  to_vector(alter_z) ~ std_normal();
  
  sigma_id ~ exponential(1);
  sigma_eth ~ exponential(1);
  sigma_village ~ exponential(1);
  sigma_alter ~ exponential(1);
  
  L_id ~ lkj_corr_cholesky(4);
  L_eth ~ lkj_corr_cholesky(4);
  L_village ~ lkj_corr_cholesky(4);
  L_alter ~ lkj_corr_cholesky(4);
  
  for (c in 1:2) {
  scale_edu[c] ~ dirichlet(rep_vector(2, max(edu)));
  }
  
  
  //// Imputation models ///////
    for (i in 1:N_id) {
      if (age_id[i] != -99) age_id[i] ~ bernoulli( age_p );
      if (edu_id[i] != -99) (edu_id[i] + 1) ~ ordered_logistic( 0, c_edu ); // adding one to satisfy ordinal constraint; this arbitrary change doesn't affect model
}
  
  for (i in 1:N_obs) {
  
    //// When no predictors missing
    if (age[i] != -99 && edu[i] != -99) {
      
      vector[3] pk; // probability of each response
      
      for (p in 1:2) pk[p] = b0[p] + id_v[id[i],p] + alter_v[alter_id[i],p] + eth_v[eth[i],p] + village_v[village[i],p] + b[1,p]*mo(scale_edu[p],edu[i]) + b[2,p]*sex[i] + b[3,p]*sex_alter[i] + b[4,p]*sex[i]*sex_alter[i] + b[5,p]*age[i] + b[6,p]*magic[i] + b[7,p]*christian[i];

    pk[3] = 0; // we get reference category (neutral, neither trust nor mistrust) for free
    
    trust[i] ~ categorical_logit( pk );
  }
  
    ///// When just age missing
    if (age[i] == -99 && edu[i] != -99) {
      
      matrix[3,2] pk; // probability of each response
      vector[2] lp; // log likelihood for each age possibiltiy 
      
      // loop over both possible age states (k)
      for (k in 1:2) {
      for (p in 1:2) { 
        pk[p,k] = b0[p] + id_v[id[i],p] + alter_v[alter_id[i],p] + eth_v[eth[i],p] + village_v[village[i],p] + b[1,p]*mo(scale_edu[p],edu[i]) + b[2,p]*sex[i] + b[3,p]*sex_alter[i] + b[4,p]*sex[i]*sex_alter[i] + b[5,p]*(k-1) + b[6,p]*magic[i] + b[7,p]*christian[i];
      }
    pk[3,k] = 0; // we get reference category (neutral, neither trust nor mistrust) for free
      }
    lp[1] = log( 1 - inv_logit( age_p ) ) + categorical_logit_lpmf(trust[i] | pk[,1] );
    lp[2] = log( inv_logit( age_p ) ) + categorical_logit_lpmf(trust[i] | pk[,2] );
    target += log_sum_exp( lp );
  }
  
      ///// When just edu missing
    if (age[i] != -99 && edu[i] == -99) {
      
      matrix[3,(max(edu)+1)] pk; // probability of each response
      vector[max(edu)+1] lp; // log likelihood for each age possibiltiy 
      vector[max(edu)+1] pr_edu; // probabiltiy of each education level
      vector[max(edu)] cum_pr_edu = inv_logit( c_edu ); // cumulative probability of each education level
      
      // loop over both possible edu states (k)
      for (k in 1:(max(edu)+1)) {
        
      if (k == 1) pr_edu[k] = cum_pr_edu[k];
      else if (k <= max(edu)) pr_edu[k] = cum_pr_edu[k] - cum_pr_edu[k-1];
      else pr_edu[k] = 1 - cum_pr_edu[k-1];
        
      for (p in 1:2) { 
        pk[p,k] = b0[p] + id_v[id[i],p] + alter_v[alter_id[i],p] + eth_v[eth[i],p] + village_v[village[i],p] + b[1,p]*mo(scale_edu[p],k-1) + b[2,p]*sex[i] + b[3,p]*sex_alter[i] + b[4,p]*sex[i]*sex_alter[i] + b[5,p]*age[i] + b[6,p]*magic[i] + b[7,p]*christian[i];
      }
    pk[3,k] = 0; // we get reference category (neutral, neither trust nor mistrust) for free
    
    lp[k] = log(pr_edu[k]) + categorical_logit_lpmf(trust[i] | pk[,k] );
      }
    target += log_sum_exp( lp );
  }
  
    ///// When both missing
    if (age[i] == -99 && edu[i] == -99) {
      
      matrix[3,(max(edu)+1)*2] pk; // probability of each response
      vector[(max(edu)+1)*2] lp; // log likelihood for each age possibiltiy 
      vector[max(edu)+1] pr_edu; // probabiltiy of each education level
      vector[max(edu)] cum_pr_edu = inv_logit( c_edu ); // cumulative probability of each education level
      
      // loop over both possible edu states (k) and age (k2)
      for (k in 1:(max(edu)+1))
      for (k2 in 1:2) {
        
      if (k == 1) pr_edu[k] = cum_pr_edu[k];
      else if (k <= max(edu)) pr_edu[k] = cum_pr_edu[k] - cum_pr_edu[k-1];
      else pr_edu[k] = 1 - cum_pr_edu[k-1];
        
      for (p in 1:2) { 
        pk[p,k + (max(edu)+1)*(k2-1)] = b0[p] + id_v[id[i],p] + alter_v[alter_id[i],p] + eth_v[eth[i],p] + village_v[village[i],p] + b[1,p]*mo(scale_edu[p],k-1) + b[2,p]*sex[i] + b[3,p]*sex_alter[i] + b[4,p]*sex[i]*sex_alter[i] + b[5,p]*(k2 - 1) + b[6,p]*magic[i] + b[7,p]*christian[i];
      }
    pk[3,k + (max(edu)+1)*(k2-1)] = 0; // we get reference category (neutral, neither trust nor mistrust) for free
    
    // We need pr(edu) * pr(adult) * pr(data|theta), additive on log scale
    if (k2 == 1) lp[k + (max(edu)+1)*(k2-1)] = log( 1-(inv_logit(age_p)) ) + categorical_logit_lpmf(trust[i] | pk[,k + (max(edu)+1)*(k2-1)] );
    else lp[k + (max(edu)+1)*(k2-1)] = log( (inv_logit(age_p)) ) + categorical_logit_lpmf(trust[i] | pk[,k + (max(edu)+1)*(k2-1)] );
      }
    target += log_sum_exp( lp );
  }
  
  } // end  loop over observations
} // end model block

generated quantities{
  matrix[2,2] Rho_id;
  matrix[2,2] Rho_eth;
  matrix[2,2] Rho_village;
  matrix[2,2] Rho_alter;
  matrix[N_obs,3] mu_pred; // categorical predictions
  int age_hat[N_obs] = age; // combination of observed and imputed age cats
  int edu_hat[N_obs] = edu; // combination of observed and imputed education

  Rho_id = L_id * L_id';
  Rho_alter = L_alter * L_alter';
  Rho_eth = L_eth * L_eth';
  Rho_village = L_village * L_village';
  
  for (i in 1:N_obs) {
    if (age[i] == -99) age_hat[i] = bernoulli_rng( age_p );
    if (edu[i] == -99) edu_hat[i] = ordered_logistic_rng( 0, c_edu ) - 1;
    
    for (p in 1:2) {
      mu_pred[i,p] = b0[p] + id_v[id[i],p] + alter_v[alter_id[i],p] + eth_v[eth[i],p] + village_v[village[i],p] + b[1,p]*mo(scale_edu[p],edu_hat[i]) + b[2,p]*sex[i] + b[3,p]*sex_alter[i] + b[4,p]*sex[i]*sex_alter[i] + b[5,p]*age_hat[i] + b[6,p]*magic[i] + b[7,p]*christian[i];
    }
   mu_pred[i,3] = 0; 
  }
}
