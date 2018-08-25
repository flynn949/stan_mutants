data {
  int<lower=1> J; //number of mutants
  int<lower=1> N; //Number of observations
  int<lower=1,upper=J> mutant[N]; //mutant for observation n
  vector[N] x; //concentration for observation n
  int<lower=0> n_trials[N]; //number of trials for observation n
  int<lower=0> y[N]; //Number of survivors for observation n
} parameters {
  real mu_centre;
  real<lower=0> sigma_centre;
  vector[J] centre_raw;

  real mu_heightlogodds;
  real<lower=0> sigma_heightlogodds;
  vector[J] heightlogodds_raw;

  real<lower=0> mu_width_squared;
  real<lower=0> sigma_width_squared;
  vector[J] width_raw_squared;
  
} transformed parameters {
  vector[J] height;
  vector[J] centre;
  vector[J] heightlogodds;
  vector<lower=0>[J] width_squared;
  

  centre = mu_centre + sigma_centre*centre_raw;

  heightlogodds = mu_heightlogodds + sigma_heightlogodds * heightlogodds_raw;

  height = inv_logit(heightlogodds);

  width_squared = mu_width_squared + sigma_width_squared * width_raw_squared;


} model {

  vector[N] psurvive;
  for (n in 1:N)
    psurvive[n] = height[mutant[n]] * exp( -0.5* ( ((x[n] - centre[mutant[n]])^2) / (width_squared[mutant[n]]) ) ); 


  mu_centre ~ normal(65,20); 
  sigma_centre ~ cauchy(0,2);
  centre_raw ~ normal(0,1);

  mu_heightlogodds ~ normal(0.5,2);
  sigma_heightlogodds ~ cauchy(0,2);
  heightlogodds_raw ~ normal(0,1);

  mu_width_squared ~ cauchy(0,5);
  sigma_width_squared ~ cauchy(0,3);
  width_raw_squared ~ normal(0,1);

  y ~ binomial(n_trials, psurvive);

} generated quantities {

  vector<lower=0>[J] width;
  width = sqrt(width_squared);

}