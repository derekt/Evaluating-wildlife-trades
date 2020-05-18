# Plot results for EU legal-seizures analyses

remove(list = ls())
script_directory = "c:/users/derekt/work/research/scripts/legal-illegal-analysis/"

source("c:/users/derekt/work/research/scripts/legal-illegal-analysis/plotting_helper_scripts.R")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_fits_nb_effort.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_fits_nb_effort_reporting.Rdata")
load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_dat.Rdata")

reporting_eu = c(0.998081, 0.997845, 0.989204, 0.990523, 0.952611, 0.968155, 0.968777, 0.985148, 0.950096, 0.834571)

#--------------------------------------------------------------------------------
# Plot the effort- and reporting-adjusted fits
pdf(
  'c:/users/derekt/work/research/dropbox/legalillegalanalogues/static_fits_EU_hierarchical_effort_reporting.pdf', width = 11, height = 8.5, paper = "USr"
)

par(mfrow = c(3, 3), mar = c(3, 3, 2, 3), oma = c(1,1,1,1))

for (i in 1:length(DAT)) {
  print(i)

  # Plot the data
  plot_data(DAT[[i]]$years, DAT[[i]]$legal_original, DAT[[i]]$illegal, DAT[[i]]$illegal_original, DAT[[i]]$taxa, DAT[[i]]$unit, DAT[[i]]$terms)
  
  if ((i %% 3) == 1)
    legend("topright", legend = c("Legal trade", "Seizure volume"), col = c("black", "red"), lty = c(1,1), pch  = c(19,15), cex = 0.75)
  
  # Plot the model fit
  plot_model_nb_effort_reporting(i, DAT, extract(FITS_nb_effort_reporting$nb), FITS_nb_effort_reporting$data, FITS_nb_effort_reporting$years, NA, reporting_eu)
}
dev.off()



# Plot the effort-only adjusted fits

load("c:/users/derekt/work/research/dropbox/legalillegalanalogues/EU_hierarchical_dat.Rdata")

pdf(
  'c:/users/derekt/work/research/dropbox/legalillegalanalogues/static_fits_EU_hierarchical_effort.pdf'
)

par(mfrow = c(4, 2), mar = c(3, 3, 2, 3))

for (i in 1:length(DAT)) {
  print(i)
  
  # Plot the data
  plot_data(DAT[[i]]$years, DAT[[i]]$legal_original, DAT[[i]]$illegal, DAT[[i]]$illegal_original, DAT[[i]]$taxa, DAT[[i]]$unit, DAT[[i]]$terms)
  
  if ((i %% 4) == 1)
    legend("topright", legend = c("Legal trade", "Seizure volume"), col = c("black", "red"), lty = c(1,1), pch  = c(1,1), cex = 0.75)
  
  # Plot the model fit
  plot_model_nb_effort(i, DAT, extract(FITS_nb_effort$nb), FITS_nb_effort$data, FITS_nb_effort$years)
  
  
}
dev.off()

source(paste(script_directory, "save_and_plot_summaries_static_only.r", sep=""))
to_plot = extract_save_fits(FITS_nb_effort_reporting, DAT, model_selected, "EU_results_effort_reporting.csv")
first_sort_order = sort(to_plot$PoissonMed, index.return = T)
to_plot = to_plot[first_sort_order$ix,]
DAT = DAT[first_sort_order$ix]


# Make summary plots
plot_caterpillar(to_plot, 1:length(DAT), extract(FITS_nb_effort_reporting$nb)$overall_slope, "EU_summary_effort_reporting.pdf", show_appendix = FALSE)

plot_caterpillar_ratio(DAT, "EU_summary_ratios.pdf")
