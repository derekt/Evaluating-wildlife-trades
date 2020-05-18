#-----------------------------------------------------------
# Script to fit models to legal and seizure data for the EU
#
# Models have been tested on simulated data prior to 
# being fitted in this script
# Derek Tittensor & Greg Britten
# with some elements derived from 
# http://tinyurl.com/j8mk8pl and http://tinyurl.com/hgua44a
# Load necessary libraries
#-----------------------------------------------------------

rm(list=ls())

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
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_linear_centred_slope_parameter.R", sep="")) #static model
source(paste(script_directory, "stan_hierarchical_nb_regression_effort_reporting_linear_centred_slope_parameter.R", sep="")) #static model 

# Pre-compile models, do this to avoid opening a new Rcpp DLL each model fit
mod_nb_effort <- stan_model(model_code=hierarchical_nb_regression_effort_linear_centred_slope_parameter)
mod_nb_effort_reporting <- stan_model(model_code=hierarchical_nb_regression_effort_reporting_linear_centred_slope_parameter)

N_iter <- 4000#1E4 #number of MCMC samples (excludes burnin)

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
# Model fitting  - EU data
#------------------------------------------------------------------------------------------------------------------------
# Load, extract, and format the data
source(paste(script_directory, "read_data_EU_robust.r", sep=""))
 
# Load the seizures effort proxy
effort_seizures_eu = read.csv(paste(data_directory, "total_annual_seizures_EU.csv", sep=""))

# Reporting scaling
#reporting_eu = c(0.82, 0.86, 0.86, 0.86, 0.79, 0.75, 0.71, 0.75, 0.57, 0.74)
reporting_eu = c(0.998081, 0.997845, 0.989204, 0.990523, 0.952611, 0.968155, 0.968777, 0.985148, 0.950096, 0.834571)

# log then z-transform seizure data
effort_seizures_eu[,2] = log(as.numeric(effort_seizures_eu[,2]))
effort_seizures_eu[,2] = as.numeric(scale(effort_seizures_eu[,2], center = TRUE, scale = TRUE))

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
   #print(get_cites_appendix(DAT[[ii]]))
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
   #print(get_cites_appendix(DAT[[ii]]) == "I/II")
 }
 
 
 DAT = DAT[retain == TRUE]
 print(length(DAT))
 
 

 J = length(DAT)
 
 
 
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
   } else if (DAT[[ii]]$terms == "trophies") {
     product_type[ii] = 3
   } else {
     product_type[ii] = 4 
   }
 }

 
# Create matrix for the hierarchical model
legal <- matrix(nrow = J, ncol = length(DAT[[1]]$year))
illegal <- matrix(nrow = J, ncol = length(DAT[[1]]$year))
years <- matrix(nrow = J, ncol = length(DAT[[1]]$year))
effort_seizures_matrix <- matrix(nrow = J, ncol = length(DAT[[1]]$year))
reporting_matrix <- matrix(nrow = J, ncol = length(DAT[[1]]$year))
num_years = rep(length(DAT[[1]]$year),J)

for (i in 1:length(DAT))
{
  legal[i,] <- DAT[[i]]$legal
  illegal[i,] <- DAT[[i]]$illegal - min(DAT[[i]]$illegal)
  years[i,] <- DAT[[i]]$years
  effort_seizures_matrix[i,] <- effort_seizures_eu[effort_seizures_eu[,1] %in% DAT[[i]]$years,2]
  reporting_matrix[i,] = reporting_eu
}

# Convert data into vectors for use in STAN
legal_vec <- linearize_data(legal, num_years)
illegal_vec <- linearize_data(illegal, num_years)
effort_vec <- linearize_data(effort_seizures_matrix, num_years)
reporting_vec <- linearize_data(reporting_matrix, num_years)
years_vec <- linearize_data(years, num_years)
N = sum(num_years)


# #--------------------------------------------------------------- 
# # Effort as total number of seizures - negative binomial
data_list <- list(J = J, N = N, y = illegal_vec, beta = legal_vec, gamma = effort_vec, timeseries_length = num_years)
fit_nb_effort <- sampling(mod_nb_effort, data=data_list, iter=N_iter, show_messages=FALSE) 
FITS_nb_effort <- list(nb = fit_nb_effort, data=data_list, years = years_vec)

#--------------------------------------------------------------- 
# Effort as total number of seizures plus reporting - negative binomial
data_list <- list(J = J, N = N, y = illegal_vec, beta = legal_vec, gamma = effort_vec, reporting = reporting_vec, timeseries_length = num_years)
fit_nb_effort_reporting <- sampling(mod_nb_effort_reporting, data=data_list, iter=N_iter, show_messages=FALSE) 
FITS_nb_effort_reporting <- list(nb = fit_nb_effort_reporting, data=data_list, years = years_vec)

save(FITS_nb_effort, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_FITS_nb_effort.Rdata")
save(FITS_nb_effort_reporting, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_FITS_nb_effort_reporting.Rdata")
save(DAT, file = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_DAT.Rdata")


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

# # Plot the results
#source(paste(script_directory, "plot_results_EU_static_only.r", sep=""))

