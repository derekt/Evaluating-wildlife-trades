#-----------------------------------------------------------
# Script to test models against simulated data
#
# prior to putting them into this script
# Derek Tittensor & Greg Britten
# with some elements derived from 
# http://tinyurl.com/j8mk8pl and http://tinyurl.com/hgua44a
# Load necessary libraries
#-----------------------------------------------------------

rm(list=ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(reshape2)

# Allows saving a Stan program so that it does't need to be recompiled
rstan_options(auto_write = TRUE)

# Lets us execute multiple Markov chains in parallel fashion
options(mc.cores = parallel::detectCores())

# Directory containing all scripts
script_directory = "c:/users/derekt/work/research/scripts/legal-illegal-analysis/"

# Load STAN models and prepare for modelling
source(paste(script_directory, "stan_hierarchical_nb_regression_linear_centred.R", sep="")) #static model 
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_linear_centred_slope_parameter.R", sep="")) #static model 
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_reporting_linear_centred_slope_parameter.R", sep="")) #static model 


# Pre-compile models, do this to avoid opening a new Rcpp DLL each model fit
mod_nb_centred <- stan_model(model_code=hierarchical_nb_regression_linear_centred) 
mod_nb_effort_centred <- stan_model(model_code=hierarchical_nb_regression_effort_linear_centred_slope_parameter) 
mod_nb_effort_reporting_centred <- stan_model(model_code=hierarchical_nb_regression_effort_reporting_linear_centred_slope_parameter) 

N_iter <- 2000 #number of MCMC samples

# Convert the data from matrices to vectors for fitting in STAN
linearize_data <- function(data_matrix, number_of_years)
{
  data_vector <- vector(length = sum(number_of_years))
  loc = 1
  for (ii in 1:length(number_of_years))
  {
    final_loc = sum(number_of_years[1:ii])
    data_vector[loc:final_loc] = data_matrix[ii,1:number_of_years[ii]]
    loc = final_loc + 1
  }
  data_vector
}

#------------------------------------------------------------------------------------------------------------------------
# Generate  test data
#------------------------------------------------------------------------------------------------------------------------
# Generate time-series for a number of taxa with an overall slope and standard deviaion (os, os_sd)
# and a variable number of years for each time-series (num_years)
num_taxa = 10
os = 2.0 # Overall slope
os_sd = 3 # SD of overall slope
alpha = seq(0,5,0.5) # Intercepts
slope = rnorm(10,os,os_sd) # Slopes linking legal to illegal trade
#num_years = c(20,20,20,20,20,20,20,20,20,20)
#num_years = c(50,50,50,50,50,50,50,50,50,50)
J = num_taxa
num_years = c(25,24,23,22,25,26,25,20,22,20)
effort_slope = 0.6

# Generate data matrices to hold the time-series, with a predictor variable (beta) and a 
# response variable (y) drawn from a neg binom distribution. 
legal <- matrix(nrow = J, ncol = max(num_years))
illegal <- matrix(nrow = J, ncol = max(num_years))
illegal_effort <- matrix(nrow = J, ncol = max(num_years))
illegal_effort_reporting <- matrix(nrow = J, ncol = max(num_years))
                         
effort <- rnorm(max(num_years), 5)
effort_vec = c()
reporting = rnorm(max(num_years),0.9, sd = 0.05)
reporting_vec = c()

for (i in 1:num_taxa)
{
  legal[i,1:num_years[i]] <- rnorm(num_years[i])
  illegal[i,1:num_years[i]] <- rnbinom(num_years[i], size = 0.5, mu = exp(alpha[i] + (slope[i] * legal[i,1:num_years[i]])))
  illegal_effort[i,1:num_years[i]] <- rnbinom(num_years[i], size = 0.5, mu = exp(alpha[i] + (slope[i] * legal[i,1:num_years[i]]) + effort_slope * effort[1:num_years[i]]))
  illegal_effort_reporting[i,1:num_years[i]] <- rnbinom(num_years[i], size = 0.5, mu = exp(alpha[i] + (slope[i] * legal[i,1:num_years[i]]) + effort_slope * effort[1:num_years[i]]) * reporting[1:num_years[i]])
  effort_vec = c(effort_vec, effort[1:num_years[i]])
  reporting_vec = c(reporting_vec, reporting[1:num_years[i]])
}

#------------------------------------------------------------------------------------------------------------------------
# Fit hierarchical model (converting input matrices to single long concatenated vector)
#------------------------------------------------------------------------------------------------------------------------

# Convert data into vectors for use in STAN
legal_vec <- linearize_data(legal, num_years)
illegal_vec <- linearize_data(illegal, num_years)
illegal_effort_vec = linearize_data(illegal_effort, num_years)
illegal_effort_reporting_vec = linearize_data(illegal_effort_reporting, num_years)

# # Non effort-corrected model - negative binomial
N = sum(num_years)
data_list <- list(J = J, N = N, y = illegal_vec, beta = legal_vec, timeseries_length = num_years)
fit_nb_inear_centred <- sampling(mod_nb_centred, data=data_list, iter=N_iter, show_messages=FALSE, control = list(adapt_delta = 0.99))

# Extract the slope estimates for each model
ii = extract(fit_nb_linear_centred)
print(paste("True overall slope: ", os), sep="")
print(paste("Mean slope from observations: ", mean(slope), sep=""))
#print(paste("Posterior mean: ", format(mean(ii$overall_slope), digits = 2), sep=""))
print(paste("Posterior median: ", format(median(ii$overall_slope), digits = 2), sep=""))
print(paste("True slope standard deviation ", os_sd), sep="")
print(paste("Posterior median: ", format(median(ii$overall_slope_sd), digits = 2), sep=""))
print(paste("True individual slopes: "))
print(format(slope, digits = 2))
print("Posterior medians: ")
print(format(apply(ii$slope, 2, quantile, probs = 0.5), digits = 2))
print(paste("True intercepts: "))
print(alpha)
print("Posterior medians: ")
print(format(apply(ii$alpha, 2, quantile, probs =0.5), digits = 2))
#summary(fit_nb_linear_centred)$summary[1:10,]



# # Effort-corrected model - negative binomial
N = sum(num_years)
data_list <- list(J = J, N = N, y = illegal_effort_vec, beta = legal_vec, gamma = effort_vec, timeseries_length = num_years)
fit_nb_linear_effort_centred <- sampling(mod_nb_effort_centred, data=data_list, iter=N_iter, show_messages=FALSE, control = list(adapt_delta = 0.99))

# Extract the slope estimates for each model
ii = extract(fit_nb_linear_effort_centred)
print(paste("True overall slope: ", os), sep="")
print(paste("Mean slope from observations: ", mean(slope), sep=""))
#print(paste("Posterior mean: ", format(mean(ii$overall_slope), digits = 2), sep=""))
print(paste("Posterior median: ", format(median(ii$overall_slope), digits = 2), sep=""))
print(paste("True slope standard deviation ", os_sd), sep="")
print(paste("Posterior median: ", format(median(ii$overall_slope_sd), digits = 2), sep=""))
print(paste("True effort parameter ", effort_slope), sep="")
print(paste("Posterior median: ", format(median(ii$effort_slope), digits = 2), sep=""))
print(paste("True individual slopes: "))
print(format(slope, digits = 2))
print("Posterior medians: ")
print(format(apply(ii$slope, 2, quantile, probs = 0.5), digits = 2))
print(paste("True intercepts: "))
print(alpha)
print("Posterior medians: ")
print(format(apply(ii$alpha, 2, quantile, probs =0.5), digits = 2))
#summary(fit_nb_linear_centred)$summary[1:10,]



# # Effort- and reporting corrected model - negative binomial
N = sum(num_years)
data_list <- list(J = J, N = N, y = illegal_effort_vec, beta = legal_vec, gamma = effort_vec, reporting = reporting_vec, timeseries_length = num_years)
fit_nb_linear_effort_reporting_centred <- sampling(mod_nb_effort_reporting_centred, data=data_list, iter=N_iter, show_messages=FALSE, control = list(adapt_delta = 0.99))

# Extract the slope estimates for each model
ii = extract(fit_nb_linear_effort_reporting_centred)
print(paste("True overall slope: ", os), sep="")
print(paste("Mean slope from observations: ", mean(slope), sep=""))
#print(paste("Posterior mean: ", format(mean(ii$overall_slope), digits = 2), sep=""))
print(paste("Posterior median: ", format(median(ii$overall_slope), digits = 2), sep=""))
print(paste("True slope standard deviation ", os_sd), sep="")
print(paste("Posterior median: ", format(median(ii$overall_slope_sd), digits = 2), sep=""))
print(paste("True effort parameter ", effort_slope), sep="")
print(paste("Posterior median: ", format(median(ii$effort_slope), digits = 2), sep=""))
print(paste("True individual slopes: "))
print(format(slope, digits = 2))
print("Posterior medians: ")
print(format(apply(ii$slope, 2, quantile, probs = 0.5), digits = 2))
print(paste("True intercepts: "))
print(alpha)
print("Posterior medians: ")
print(format(apply(ii$alpha, 2, quantile, probs =0.5), digits = 2))?
#summary(fit_nb_linear_centred)$summary[1:10,]

