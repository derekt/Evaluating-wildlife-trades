#--------------------------------------------------------------------------------
# Plot the effort-adjusted fits
# model_selected = list()

source("c:/users/derekt/work/research/scripts/legal-illegal-analysis/plotting_helper_scripts.R")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_fits_nb_effort.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_fits_nb_effort_diff_slopes.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_fits_nb_effort_short.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarhical_fits_searches.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_dat.Rdata")

data_directory = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/"
setwd(data_directory)

effort_seizures = read.csv(paste(data_directory, "total-annual-siezures.csv", sep=""))
l1 = which(effort_seizures[,1] == 1982)
effort_seizures = effort_seizures[l1:dim(effort_seizures)[1],]
inspections = read.csv(paste(data_directory, "Inspections.csv", sep=""))


# Directory containing all scripts
script_directory = "c:/users/derekt/work/research/scripts/legal-illegal-analysis/"

#--------------------------------------------------------------------------------

# Plot the effort-only adjusted fits
pdf(
  'c:/users/derekt/work/research/dropbox/legalillegalanalogues/static_fits_US_hierarchical_effort.pdf', paper = "a4r", width = 11, height = 8.5
)

par(mfrow = c(3, 3), mar = c(3, 3, 2, 3), oma = c(1,1,1,1))

for (i in 1:length(DAT)) {
  print(i)
  
  # Plot the data
  plot_data(DAT[[i]]$years, DAT[[i]]$legal_original, DAT[[i]]$illegal, DAT[[i]]$illegal_original, DAT[[i]]$taxa, DAT[[i]]$unit, DAT[[i]]$terms)
  
  if ((i %% 3) == 1)
    legend("topright", legend = c("Reported legal trade", "Seizure volume"), col = c("black", "red"), lty = c(1,1), pch  = c(19,15), cex = 0.75)
  
  # Plot the model fit
  plot_model_nb_effort(i, DAT, extract(FITS_nb_effort$nb), FITS_nb_effort$data, FITS_nb_effort$years)
}
dev.off()

# Extract and save results
source(paste(script_directory, "save_and_plot_summaries_static_only.r", sep=""))
to_plot = extract_save_fits(FITS_nb_effort, DAT, model_selected, "US_results_effort_reporting.csv")
first_sort_order = sort(to_plot$PoissonMed, index.return = T)
to_plot = to_plot[first_sort_order$ix,]

DAT = DAT[first_sort_order$ix]
plot_caterpillar(to_plot, 1:length(DAT), extract(FITS_nb_effort$nb)$overall_slope, "US_summary_effort.pdf", show_appendix = FALSE)

plot_caterpillar_ratio(DAT, "US_summary_ratios.pdf")


#-------------------------------------------------------------------------------------

# Plot the model where effort is the number of container inspections per year
#--------------------------------------------------------------------------------
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_dat.Rdata")
# pdf(
#   'c:/users/derekt/work/research/dropbox/legalillegalanalogues/static_fits_effort_inspections_US.pdf', family="ArialMT"
# )
# 
# par(mfrow = c(4, 2), mar = c(3, 3, 2, 3))
# 
# for (i in 1:length(DAT)) {
#   print(i)
#   N <- DAT[[i]]$N
#   min_year = 2014 - N
#   years = (min_year + 1):2014
#   legal <- DAT[[i]]$legal
#   illegal <- DAT[[i]]$illegal
#   illegal_original = DAT[[i]]$illegal_original
#   legal_original = DAT[[i]]$legal_original
#   data_minimum_illegal = min(illegal)
#   
#   years2 = which(DAT[[i]]$years %in% inspections[,1] == TRUE)
#     
#   plot_data(years[years2], DAT[[i]]$legal_original[years2], DAT[[i]]$illegal[years2], DAT[[i]]$illegal_original[years2], DAT[[i]]$taxa, DAT[[i]]$unit, DAT[[i]]$terms)
#   
#   #taxon_number, raw_data, model_fit, model_fit_data, years, year_restriction = NA
#   # Plot the model fit
#   plot_model_nb_effort(i, DAT, extract(FITS_nb_searches$nb), FITS_nb_searches$data, FITS_nb_searches$years, year_restriction = years2)
#   }
# 
# dev.off()


# Extract and save results.
source(paste( script_directory, "save_and_plot_summaries_static_only.r", sep = ""))

to_plot = extract_save_fits(FITS_nb_searches, DAT, model_selected, "US_results_effort_inspections.csv")
first_sort_order = sort(to_plot$PoissonMed, index.return = T)
to_plot = to_plot[first_sort_order$ix,]

DAT = DAT[first_sort_order$ix]

plot_caterpillar(to_plot, 1:length(DAT), extract(FITS_nb_searches$nb)$overall_slope, "US_summary_effort_inspections.pdf", show_appendix = FALSE)

graphics.off()

#-------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------
# Plot the effort-adjusted fits where effort is the shortened seizures fits
model_selected = list()


load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_dat.Rdata")
# pdf(
#   'c:/users/derekt/work/research/dropbox/legalillegalanalogues/static_fits_effort_short_US.pdf', family="ArialMT"
# )
# 
# par(mfrow = c(4, 2), mar = c(3, 3, 2, 3))
# 
# for (i in 1:length(DAT)) {
#   print(i)
#   N <- DAT[[i]]$N
#   min_year = 2014 - N
#   years = (min_year + 1):2014
#   legal <- DAT[[i]]$legal
#   illegal <- DAT[[i]]$illegal
#   illegal_original = DAT[[i]]$illegal_original
#   legal_original = DAT[[i]]$legal_original
#   data_minimum_illegal = min(illegal)
#   
#   years2 = which(DAT[[i]]$years %in% inspections[,1] == TRUE)
#   
#   plot_data(years[years2], DAT[[i]]$legal_original[years2], DAT[[i]]$illegal[years2], DAT[[i]]$illegal_original[years2], DAT[[i]]$taxa, DAT[[i]]$unit, DAT[[i]]$terms)
#   plot_model_nb_effort(i, DAT, extract(FITS_nb_effort_short$nb), FITS_nb_effort_short$data, FITS_nb_effort_short$years, year_restriction = years2)
# }

  
  
  # Extract and save results.
  load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/US_hierarchical_dat.Rdata")
  source(paste( script_directory, "save_and_plot_summaries_static_only.r", sep = ""))
  
  to_plot = extract_save_fits(FITS_nb_effort_short, DAT, model_selected, "US_results_effort_shortened.csv")
  first_sort_order = sort(to_plot$PoissonMed, index.return = T)
  to_plot = to_plot[first_sort_order$ix,]
  #FITS = FITS[first_sort_order$ix]
  DAT = DAT[first_sort_order$ix]
  
  plot_caterpillar(to_plot, 1:length(DAT), extract(FITS_nb_effort_short$nb)$overall_slope, "US_summary_effort_shortened.pdf", show_appendix = FALSE)
  

graphics.off()
#-------------------------------------------------------------------------------------