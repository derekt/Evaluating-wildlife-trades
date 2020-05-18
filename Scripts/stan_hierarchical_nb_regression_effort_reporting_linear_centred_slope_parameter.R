hierarchical_nb_regression_effort_reporting_linear_centred_slope_parameter <- 
"data {
  int<lower=0> J; // Number of taxon-products
  int<lower=0> N; // Total number of observations
  vector[N] beta; // Legal trade observations
  int y[N]; // Seizures
  vector[N] gamma; // Yearly effort
  vector[N] reporting; // Yearly reporting level
  int timeseries_length[J]; // Length of time-series for each taxon-product
}

parameters {
  real overall_slope; // Overall meta-analysis mean
  real<lower = 0.00001> overall_slope_sd; // SD hyperparameter

  real effort_overall_slope; // Effort slope mean hyperparameter 
  real<lower = 0.00001> effort_overall_slope_sd; // Effort SD hyperparameter

  real alpha[J]; // Individual taxon-product intercepts
  real<lower = 0.00001, upper = 100> phi[J];
  real slope_raw[J]; // non-centred parameterization slope parameters
  real effort_slope_raw[J]; // non-centred effort slope parameters
}  

transformed parameters {
  real slope[J]; // Individual centred taxon-product slopes
  real effort_slope[J]; // Individual taxon-product effort slopes

  // implies: slope ~ normal(overall_slope, overall_slope_sd)
  for (taxa in 1:J){
    slope[taxa] = overall_slope + overall_slope_sd * slope_raw[taxa];

    // implies: effort_slope ~ normal(effort_overall_slope, effort_overall_slope_sd)
    effort_slope[taxa] = effort_overall_slope + effort_overall_slope_sd * effort_slope_raw[taxa];
  }  
}

model {
  int pos;
  pos = 1;

  // Prior distributions (i.e. distributions conditional on fixed hyperparameters)
  overall_slope ~ normal(0,1000);    //assumes a N(0,1000) prior on the overall slope
  overall_slope_sd ~ cauchy(0,5); //assumes a Cauchy(0,5) prior on the overall SD
  effort_overall_slope ~ normal(0,1000);    //assumes a N(0,1000) prior on the effort overall slope
  effort_overall_slope_sd ~ cauchy(0,5); //assumes a Cauchy(0,5) prior on the effort overall SD

  for (taxa in 1:J)
  {
    //alpha[taxa] ~ normal(0,1000); //N(0,1000) prior for the taxa-level intercept
    
    // Hierarchical distributions (i.e. distributions conditional on unknown parameters)
   slope_raw[taxa] ~ normal(0,1);
    effort_slope_raw[taxa] ~ normal(0,1);

   segment(y, pos, timeseries_length[taxa]) ~ neg_binomial_2(exp(alpha[taxa] + slope[taxa] * segment(beta, pos, timeseries_length[taxa]) + effort_slope[taxa] * segment(gamma, pos, timeseries_length[taxa])) .* segment(reporting, pos, timeseries_length[taxa]), phi[taxa]);  //likelihood of the data distributed Poisson

   pos = pos + timeseries_length[taxa];
  }

}

generated quantities {
vector[N] lambda;
vector[N] log_lik;
//vector[N] pred;
int location;

location = 0;

// Predicted values used for calculating the linear regression credible intervals
for (taxa in 1:J)
{
  for (year in 1:timeseries_length[taxa])
  {
    location = location + 1;

    // Predicted values used for calculating the linear regression credible intervals
    lambda[location] = exp(alpha[taxa] + beta[location] * slope[taxa] + effort_slope[taxa] * gamma[location]).* reporting[location];

    // Log-likelihood
    log_lik[location] = neg_binomial_2_lpmf(y[location] | lambda[location],phi[taxa]);

    // New draw for predictive confidence intervals
    //pred[location] = neg_binomial_2_rng(lambda[location], phi[taxa]);

  }
}
}"
