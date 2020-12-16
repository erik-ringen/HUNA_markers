data{
  int N;
  real CFST[N];
  int norm[N];
  int cat[N];
  int level[N];
}

parameters{
  real a; // intercept, latent scale
  real<lower=0> phi; // dispersion, latent scale
  
  vector[max(norm)] norm_z; // random effects for norm/question, unscaled
  vector[max(level)] level_z; // random effects for level, unscaled
  matrix[(max(level)+1),max(cat)] level_cat_z; // interaction of level and category, unscaled and uncorrelated
  real<lower=0> sigma_norm; // sd of norm random effects
  real<lower=0> sigma_level;
  vector<lower=0>[max(level) + 1] sigma_level_cat; // vector of sds
  cholesky_factor_corr[max(level)+1] L_Rho; // lower tri cholesky of the correlation between level random effects
}

transformed parameters{
  vector[max(norm)] norm_v;
  vector[max(level)] level_v;
  matrix[max(cat),(max(level) + 1)] level_cat_v; // note the transpose
  
  norm_v = norm_z * sigma_norm;
  level_v = level_z * sigma_level;
  level_cat_v = (diag_pre_multiply(sigma_level_cat, L_Rho) * level_cat_z)';
}

model{
  // Priors
  a ~ std_normal();
  phi ~ std_normal();
  norm_z ~ std_normal();
  level_z ~ std_normal();
  to_vector(level_cat_z) ~ std_normal();
  sigma_norm ~ exponential(1);
  sigma_level ~ exponential(1);
  sigma_level_cat ~ exponential(1);
  L_Rho ~ lkj_corr_cholesky(4);
  
  for (i in 1:N) {
    real mu; // latent-scale expectation  
    real alpha; // Beta dist parameters
    real beta; // same
  
    mu = a + norm_v[norm[i]] + level_v[level[i]] + level_cat_v[cat[i],1] + level_cat_v[cat[i],(level[i]+1)];
    alpha = inv_logit(mu)*(1/phi);
    beta = (1 - inv_logit(mu))*(1/phi);
    
    // Likelihood
    CFST[i] ~ beta( alpha, beta );
  }
}

generated quantities{
  matrix[max(level)+1,max(level)+1] Rho;
  Rho = L_Rho * L_Rho';
}
