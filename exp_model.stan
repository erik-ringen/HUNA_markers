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
  int correct[N_obs]; // did they guess alter correctly?
  int id[N_obs];
  int alter_id[N_obs];
  int eth[N_obs];
  int eth_alter[N_obs];
  int village[N_obs];
  int edu[N_obs];
  int sex[N_obs];
  int sex_alter[N_obs];
  int age[N_obs]; // age category
  int magic[N_obs];
  int christian[N_obs];
}

parameters{
  real b0; // intercept
  vector[7] b; // fixed effects
  simplex[max(edu)] scale_edu; //  monotonic effect of education 
  
  matrix[4,N_id] id_z; // random effects for individual (ego), unscaled, uncentered, and uncorrelated
  matrix[5,N_eth] eth_z; // random effects for ethnicity
  matrix[4,N_village] village_z; // random effects for village
  vector[N_alter] alter_z; // random effects for alter
  vector<lower=0>[5] sigma_eth; // variance components
  vector<lower=0>[4] sigma_village;
  vector<lower=0>[4] sigma_id;
  real<lower=0> sigma_alter;
  cholesky_factor_corr[4] L_id; // correlations between random effects
  cholesky_factor_corr[5] L_eth; // correlations between random effects
  cholesky_factor_corr[4] L_village; // correlations between random effects
}

transformed parameters{
  // Now, we need to scale and correlate the random effects
  matrix[N_id,4] id_v;
  matrix[N_eth,5] eth_v;
  matrix[N_village,4] village_v;
  vector[N_alter] alter_v;

  id_v = (diag_pre_multiply(sigma_id, L_id) * id_z)';
  eth_v = (diag_pre_multiply(sigma_eth, L_eth) * eth_z)';
  village_v = (diag_pre_multiply(sigma_village, L_village) * village_z)';
  
  alter_v = alter_z * sigma_alter;
}

model{
  // priors
  b0 ~ std_normal();
  b ~ std_normal();
  to_vector(id_z) ~ std_normal();
  to_vector(eth_z) ~ std_normal();
  to_vector(village_z) ~ std_normal();
  alter_z ~ std_normal();
  
  sigma_id ~ exponential(1);
  sigma_eth ~ exponential(1);
  sigma_village ~ exponential(1);
  sigma_alter ~ exponential(1);
  
  L_id ~ lkj_corr_cholesky(4);
  L_eth ~ lkj_corr_cholesky(4);
  L_village ~ lkj_corr_cholesky(4);

  scale_edu ~ dirichlet(rep_vector(2, max(edu)));
  
  // likelihood
  for (i in 1:N_obs) {
  correct[i] ~ bernoulli_logit( b0 + id_v[id[i],1] + id_v[id[i],(eth_alter[i]+1)] + eth_v[eth[i],1] + eth_v[eth_alter[i],2] + eth_v[eth[i],(eth_alter[i]+2)] + village_v[village[i],1] + village_v[village[i],(eth_alter[i]+1)] + alter_v[alter_id[i]] + b[1]*mo(scale_edu,edu[i]) + b[2]*sex[i] + b[3]*sex_alter[i] + b[4]*sex[i]*sex_alter[i] + b[5]*age[i] + b[6]*magic[i] + b[7]*christian[i] );
  }
}

generated quantities{
  matrix[4,4] Rho_id;
  matrix[5,5] Rho_eth;
  matrix[4,4] Rho_village;
  
  Rho_id = L_id * L_id';
  Rho_eth = L_eth * L_eth';
  Rho_village = L_village * L_village';
}
