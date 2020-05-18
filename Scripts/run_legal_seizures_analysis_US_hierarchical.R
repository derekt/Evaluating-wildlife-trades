#-----------------------------------------------------------
# Script to fit models to legal and seizure data for the US
#
# Models have been tested on simulated data prior to 
# being fitted in this script
# Derek Tittensor & Greg Britten
# with some elements derived from 
# http://tinyurl.com/j8mk8pl and http://tinyurl.com/hgua44a
# Load necessary libraries
#-----------------------------------------------------------


rm(list=ls())

library(loo)

# Directory containing all scripts
script_directory = "c:/users/derekt/work/research/scripts/legal-illegal-analysis/"

# Directory containing all data
data_directory = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/"

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(reshape2)

# Allows saving a Stan program so that it does't need to be recompiled
rstan_options(auto_write = TRUE)

# Lets us execute multiple Markov chains in parallel fashion
options(mc.cores = parallel::detectCores())

# Load STAN models and prepare for modelling
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_linear_centred_slope_parameter.R", sep=""))
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_both_linear_centred_slope_parameter.R", sep=""))

# Pre-compile models, do this to avoid opening a new Rcpp DLL each model fit
mod_nb_effort <- stan_model(model_code=hierarchical_nb_regression_effort_linear_centred_slope_parameter)
mod_nb_effort_both <- stan_model(model_code = hierarchical_nb_regression_effort_both_linear_centred_slope_parameter)

N_iter <- 4000 #number of MCMC samples (excludes burnin)

# Return the CITES Appendix for an individual species
get_cites_appendix <- function(data)
{
  
  if (sum(data$CITES == "I") == length(data$CITES))
    return("I")
  else if (sum(data$CITES == "II") == length(data$CITES))
    return("II")
  else if (sum(data$CITES == "I/II") == length(data$CITES))
    return("I/II")
  else
    return("Mixed")
}

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
# Model fitting  - US data
#------------------------------------------------------------------------------------------------------------------------
# Load, extract, and format the data
source(paste(script_directory, "read_data_US_robust.r", sep=""))

# US DATA
effort_seizures = read.csv(paste(data_directory, "total-annual-siezures.csv", sep=""))
l1 = which(effort_seizures[,1] == 1982)
effort_seizures = effort_seizures[l1:dim(effort_seizures)[1],]
inspections = read.csv(paste(data_directory, "Inspections.csv", sep=""))

# log_transform and z-transform
effort_seizures[,2] = log(as.numeric(effort_seizures[,2]))
effort_seizures[,2] = as.numeric(scale(effort_seizures[,2], center = TRUE, scale = TRUE))
inspections[,2] = log(as.numeric(inspections[,2]))
inspections[,2] = as.numeric(scale(inspections[,2], center = TRUE, scale = TRUE))

# Now go through and remove species with a mean of <5 seizures per year
retain = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (mean(DAT[[ii]]$illegal) >= 5)
    retain[ii] = TRUE
  else
    retain[ii] = FALSE
}

DAT_OLD = DAT
DAT = DAT[retain == TRUE]

# Remove species with a mean legal trade of <=3 per year
remove(retain)
retain = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (sum(DAT[[ii]]$legal_original == 0) <= 3)
    retain[ii] = TRUE
  else
    retain[ii] = FALSE
}

DAT = DAT[retain == TRUE]
print(length(DAT))

# Remove species with a mean original illegal trade of <= 5 per year
remove(retain)
retain = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (sum(DAT[[ii]]$illegal_original == 0) <= 5)
    retain[ii] = TRUE
  else
    retain[ii] = FALSE
}

DAT = DAT[retain == TRUE]
print(length(DAT))



# Remove species which are not CITES Appendix I or II for the entire period
retain = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (get_cites_appendix(DAT[[ii]]) == "I/II")
  {
    retain[ii] = FALSE
  } else if (get_cites_appendix(DAT[[ii]]) == "Mixed")
  {
    retain[ii] = FALSE
  } else
  {
    retain[ii] = TRUE
  }
}


DAT = DAT[retain == TRUE]
print(length(DAT))

J = length(DAT)

# Get number of years in each time-series
num_years_per_timeseries = vector(length = length(DAT))
for (jj in 1:length(DAT))
{
  num_years_per_timeseries[jj] = length(DAT[[jj]]$illegal_original)
}

# Get CITES Appendix for each taxa remaining in the model to feed into the analysis
cites_appendix = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (get_cites_appendix(DAT[[ii]]) == "I")
  {
    cites_appendix[ii] = 1
  } else
  {
    
    cites_appendix[ii] = 2
  }
}

# Get taxonomic group for each taxa remaining in the model to feed into the analysis
taxonomic_group = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (DAT[[ii]]$class == "REPTILIA")
  {
    taxonomic_group[ii] = 1
  } else if (DAT[[ii]]$class == "MAMMALIA") {
    
    taxonomic_group[ii] = 2
  } else if (DAT[[ii]]$class == "AVES") {
    taxonomic_group[ii] = 3
  } else if ((DAT[[ii]]$class == "GASTROPODA") || (DAT[[ii]]$class == "BIVALVIA")) {
    taxonomic_group[ii] = 4 
  }
    else {
      taxonomic_group[ii] = 5
  }
}

# Get product information to feed into the analysis
product_type = vector(length = length(DAT))
for (ii in 1:length(DAT))
{
  if (DAT[[ii]]$terms == "leather products & items")
  {
    product_type[ii] = 1
  } else if (DAT[[ii]]$terms == "live") {
    
    product_type[ii] = 2
  } else if ((DAT[[ii]]$terms == "trophies") || (DAT[[ii]]$terms == "skulls & trophies")) {
    product_type[ii] = 3
  } else {
    product_type[ii] = 4 
  }
}

# Create matrices for the hierarchical model
legal_effort <- matrix(0,nrow = J, ncol = max(num_years_per_timeseries))
illegal_effort <- matrix(0,nrow = J, ncol = max(num_years_per_timeseries))
legal_inspections <- matrix(0,nrow = J, ncol = dim(inspections)[1])
illegal_inspections <- matrix(0,nrow = J, ncol = dim(inspections)[1])
legal_short <- matrix(0,nrow = J, ncol = dim(inspections)[1])
illegal_short <- matrix(0,nrow = J, ncol = dim(inspections)[1])
years <- matrix(0, nrow = J, ncol = max(num_years_per_timeseries))
years_effort <- matrix(0,nrow = J, ncol = max(num_years_per_timeseries))
years_inspections <- matrix(0, nrow = J, ncol = dim(inspections)[1])
years_effort_short <- matrix(0,nrow = J, ncol = dim(inspections)[1])
effort_seizures_matrix <- matrix(nrow = J, ncol = max(num_years_per_timeseries))
inspections_matrix <- matrix(nrow = J, ncol = dim(inspections)[1])
effort_seizures_short_matrix <- matrix(nrow = J, ncol = dim(inspections)[1])
num_years_per_timeseries_effort <- vector(length = J)
num_years_per_timeseries_inspections <- vector(length = J)
num_years_per_timeseries_effort_short <- vector(length = J)

for (ii in 1:length(DAT))
{
  effort_years <- intersect(effort_seizures[,1], DAT[[ii]]$years)
  num_years_per_timeseries_effort[ii] <- length(effort_years)
  data_years_to_pick = DAT[[ii]]$years %in%effort_years
  legal_effort[ii,1:num_years_per_timeseries_effort[ii]] <- DAT[[ii]]$legal[data_years_to_pick]
  illegal_effort[ii,1:num_years_per_timeseries_effort[ii]] <- DAT[[ii]]$illegal[data_years_to_pick] - min(DAT[[ii]]$illegal[data_years_to_pick])
  
  effort_seizures_matrix[ii,1:length(effort_years)] <- effort_seizures[effort_seizures[,1] %in% DAT[[ii]]$years,2]
  years_effort[ii,1:length(effort_years)] = effort_years
  
  inspection_years <- intersect(inspections[,1], DAT[[ii]]$years)
  num_years_per_timeseries_inspections[ii] <- length(inspection_years)
  data_years_to_pick <- DAT[[ii]]$years %in% inspection_years
  legal_inspections[ii,1:num_years_per_timeseries_inspections[ii]] <- DAT[[ii]]$legal[data_years_to_pick]
  illegal_inspections[ii,1:num_years_per_timeseries_inspections[ii]] <- DAT[[ii]]$illegal[data_years_to_pick] - min(DAT[[ii]]$illegal[data_years_to_pick])
  inspections_matrix[ii,1:length(inspection_years)] <- inspections[inspections[,1] %in% DAT[[ii]]$years,2]
  years_inspections[ii,1:length(inspection_years)] = inspection_years
  
  effort_short_years <- intersect(inspection_years, DAT[[ii]]$years)
  num_years_per_timeseries_effort_short[ii] <- length(effort_short_years)
  data_years_to_pick <- DAT[[ii]]$years %in% effort_short_years
  legal_short[ii,1:num_years_per_timeseries_effort_short[ii]] <- DAT[[ii]]$legal[data_years_to_pick]
  illegal_short[ii,1:num_years_per_timeseries_effort_short[ii]] <- DAT[[ii]]$illegal[data_years_to_pick] - min(DAT[[ii]]$illegal[data_years_to_pick])
  effort_seizures_short_matrix[ii,1:length(effort_short_years)] <- effort_seizures[effort_seizures[,1] %in% effort_short_years,2]
  years_effort_short[ii,1:length(effort_short_years)] = effort_short_years
}

# Convert data into vectors for use in STAN
legal_effort_vec <- linearize_data(legal_effort, num_years_per_timeseries_effort)
illegal_effort_vec <- linearize_data(illegal_effort, num_years_per_timeseries_effort)
effort_vec <- linearize_data(effort_seizures_matrix, num_years_per_timeseries_effort)
years_vec_effort <- linearize_data(years_effort, num_years_per_timeseries_effort)

legal_inspections_vec <- linearize_data(legal_inspections, num_years_per_timeseries_inspections)
illegal_inspections_vec <- linearize_data(illegal_inspections, num_years_per_timeseries_inspections)
inspections_vec <- linearize_data(inspections_matrix, num_years_per_timeseries_inspections)
years_vec_inspections <- linearize_data(years_inspections, num_years_per_timeseries_inspections)

legal_short_vec <- linearize_data(legal_short, num_years_per_timeseries_effort_short)
illegal_short_vec <- linearize_data(illegal_short, num_years_per_timeseries_effort_short)
effort_short_vec <- linearize_data(effort_seizures_short_matrix, num_years_per_timeseries_effort_short)
years_vec_effort_short <- linearize_data(years_effort_short, num_years_per_timeseries_effort_short)

#--------------------------------------------------------------- 
# Effort as total number of seizures - negative binomial
data_list <- list(J = J, N = length(legal_effort_vec), y = illegal_effort_vec, beta = legal_effort_vec,
                  gamma = effort_vec, timeseries_length = num_years_per_timeseries_effort)
fit_nb_effort <- sampling(mod_nb_effort, data=data_list, iter=N_iter, show_messages=FALSE)
FITS_nb_effort <- list(nb = fit_nb_effort, data=data_list, years = years_vec_effort)

#--------------------------------------------------------------- 
# Effort as container searches - negative binomial
data_list2 <- list(J = J, N = length(legal_inspections_vec), y = illegal_inspections_vec, beta = legal_inspections_vec,
                   gamma = inspections_vec, timeseries_length = num_years_per_timeseries_inspections)
fit_nb_searches <- sampling(mod_nb_effort, data=data_list2, iter=N_iter, show_messages=FALSE)
FITS_nb_searches <- list(nb = fit_nb_searches, data=data_list2, years = years_vec_inspections)

#--------------------------------------------------------------- 
# EFfort as total number of seizures (short time-series to match container searches)
data_list3 <- list(J = J, N = length(legal_short_vec), y = illegal_short_vec, beta = legal_short_vec,
                   gamma = effort_short_vec, timeseries_length = num_years_per_timeseries_effort_short)
fit_nb_effort_short <- sampling(mod_nb_effort, data=data_list3, iter=N_iter, show_messages=FALSE)
FITS_nb_effort_short <- list(nb = fit_nb_effort_short, data=data_list3, years = years_vec_effort_short)

#--------------------------------------------------------------- 
# EFfort as container searches AND total number of seizures (short time-series to match container searches)
data_list4 <- list(J = J, N = length(legal_short_vec), y = illegal_short_vec, beta = legal_short_vec,
                   gamma = effort_short_vec, gamma2 = inspections_vec, timeseries_length = num_years_per_timeseries_effort_short)
fit_nb_effort_both <- sampling(mod_nb_effort_both, data=data_list4, iter=N_iter, control = list(adapt_delta = 0.95), show_messages=FALSE)
FITS_nb_effort_both <- list(nb = fit_nb_effort_both, data=data_list4, years = years_vec_effort_short)

save(FITS_nb_effort, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort.Rdata")
save(FITS_nb_searches, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_searches.Rdata")
save(FITS_nb_effort_short, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort_short.Rdata")
save(FITS_nb_effort_both, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort_both.Rdata")
save(DAT, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_DAT.Rdata")

load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort_short.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_nb_effort_both.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_FITS_searches.Rdata")

# Model comparison to assess differences between CITES Appendix, taxa, and products.

# Extract relevant posteriors and classify, then find 95% CI quantiles
effort_posteriors = extract(FITS_nb_effort$nb)
cites_appendix1_posteriors = sort(as.vector(effort_posteriors$slope[,cites_appendix == 1]))
cites_appendix2_posteriors = sort(as.vector(effort_posteriors$slope[,cites_appendix == 2]))
print(paste("Appendix 1 2.5% CI: ", quantile(cites_appendix1_posteriors,0.025), "97.5% CI: ",quantile(cites_appendix1_posteriors,0.975),sep=" "))
print(paste("Appendix 2 2.5% CI: ", quantile(cites_appendix2_posteriors,0.025), "97.5% CI: ",quantile(cites_appendix2_posteriors,0.975),sep=" "))

# Extract relevant posteriors and classify, then find 95% CI quantiles
product1_posteriors = sort(as.vector(effort_posteriors$slope[,product_type == 1]))
product2_posteriors = sort(as.vector(effort_posteriors$slope[,product_type == 2]))
product3_posteriors = sort(as.vector(effort_posteriors$slope[,product_type == 3]))
product4_posteriors = sort(as.vector(effort_posteriors$slope[,product_type == 4]))
print(paste("Product 1 2.5% CI: ", quantile(product1_posteriors,0.025), "97.5% CI: ",quantile(product1_posteriors,0.975),sep=" "))
print(paste("Product 2 2.5% CI: ", quantile(product2_posteriors,0.025), "97.5% CI: ",quantile(product2_posteriors,0.975),sep=" "))
print(paste("Product 3 2.5% CI: ", quantile(product3_posteriors,0.025), "97.5% CI: ",quantile(product3_posteriors,0.975),sep=" "))
print(paste("Product 4 2.5% CI: ", quantile(product4_posteriors,0.025), "97.5% CI: ",quantile(product4_posteriors,0.975),sep=" "))

# Extract relevant posteriors and classify, then find 95% CI quantiles
taxa1_posteriors = sort(as.vector(effort_posteriors$slope[,taxonomic_group == 1]))
taxa2_posteriors = sort(as.vector(effort_posteriors$slope[,taxonomic_group == 2]))
taxa3_posteriors = sort(as.vector(effort_posteriors$slope[,taxonomic_group == 3]))
taxa4_posteriors = sort(as.vector(effort_posteriors$slope[,taxonomic_group == 4]))
taxa5_posteriors = sort(as.vector(effort_posteriors$slope[,taxonomic_group == 5]))
print(paste("Taxa 1 2.5% CI: ", quantile(taxa1_posteriors,0.025), "97.5% CI: ",quantile(taxa1_posteriors,0.975),sep=" "))
print(paste("Taxa 2 2.5% CI: ", quantile(taxa2_posteriors,0.025), "97.5% CI: ",quantile(taxa2_posteriors,0.975),sep=" "))
print(paste("Taxa 3 2.5% CI: ", quantile(taxa3_posteriors,0.025), "97.5% CI: ",quantile(taxa3_posteriors,0.975),sep=" "))
print(paste("Taxa 4 2.5% CI: ", quantile(taxa4_posteriors,0.025), "97.5% CI: ",quantile(taxa4_posteriors,0.975),sep=" "))
print(paste("Taxa 5 2.5% CI: ", quantile(taxa5_posteriors,0.025), "97.5% CI: ",quantile(taxa5_posteriors,0.975),sep=" "))

# Model comparison for effort proxies
loo_effort <- loo(extract_log_lik(FITS_nb_effort$nb))

loo_effort_container <- loo(extract_log_lik(FITS_nb_searches$nb))
loo_effort_short <- loo(extract_log_lik(FITS_nb_effort_short$nb))
loo_effort_both <- loo(extract_log_lik(FITS_nb_effort_both$nb))

print(compare(loo_effort_container, loo_effort_short))
print(compare(loo_effort_container, loo_effort_both))
print(compare(loo_effort_short, loo_effort_both))

# Plot the results
#source(paste(script_directory, "plot_results_US_static_only.r", sep=""))