# Helper scripts for the plotting functions

library(rstan)

# Extract concatenated values
devector <- function(data_vector, taxon_number, number_of_years)
{
  if (taxon_number == 1)
    point1 = 1
  else
    point1 = sum(number_of_years[1:(taxon_number - 1)]) + 1
  point2 = sum(number_of_years[1:taxon_number])

  return(data_vector[point1:point2])
}

# Extract concatenated values from model ouput
devector_model_results <- function(data_matrix, taxon_number, number_of_years)
{
  if (taxon_number == 1)
    point1 = 1
  else
    point1 = sum(number_of_years[1:(taxon_number - 1)]) + 1
  point2 = sum(number_of_years[1:taxon_number])

  return(data_matrix[,point1:point2])
}

# Function to calculate r-squared values
calculate_rsq <- function(y_data, estimates)
{
  return(cor(y_data, estimates)^2)
}

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


# Plot the raw data
plot_data <- function(years, legal_original, illegal, illegal_original, taxon, unit, term)
  
{
 
  # Plot the legal and illegal data
  plot(years, legal_original, type = 'b', ylim = c(min(c(legal_original)), max(c(legal_original))), xlab = "Year", pch = 19)
  mtext("Reported legal trade volume", side = 2, line = 2, cex = 0.5)
  lines(years, legal_original)
  par(new=TRUE)
  plot(illegal,col='red',type='b',xaxt='n',yaxt='n',ylab='n',xlab='n', pch = 15)
  lines(years, illegal, col = 'red')
  axis(4, ylim = c(min(illegal), max(illegal)))
  mtext("Seizure volume", side = 4, line = 2, cex = 0.5)
  
  #if ((mean(illegal) < 5) || (sum(legal_original == 0) > 7) || (sum(illegal_original == 0) > 5))
  #{
  #  print("ERROR! Trying to plot model fit for data removed due to insufficient legal/illegal")
  #} else
  #{
    mtext(bquote(paste(italic(.(DAT[[i]]$taxa)))), line = 1.03, cex = 0.7)
    if (DAT[[i]]$unit == "KIL")
    {  
      mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), " [kg])", sep="")), cex = 0.7, line = -0.1)
    } else
    {  
      mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), .(DAT[[i]]$unit), ")", sep="")), cex = 0.7, line = -0.1)
    }
      #}
}

# Plot model results - hierarchical negative binomial model
plot_model_nb_effort_reporting <- function(taxon_number, raw_data, model_fit, model_fit_data, years, year_restriction = NA, reporting)
{
  # Extract the raw data
  taxon_data <- DAT[[taxon_number]]
  
  # Extract the effort data
  gamma = devector(model_fit_data$gamma, taxon_number, model_fit_data$timeseries_length)

  # Get the minimum illegal data value to rescale back to original
  data_minimum_illegal = min(taxon_data$illegal)
  
  # Extract model results
  # Fitted values
  lambda <- devector_model_results(model_fit$lambda, taxon_number, model_fit_data$timeseries_length)
  lambda_median = apply(lambda, 2, quantile, probs = 0.5) + data_minimum_illegal
  lambda_lcb = apply(lambda, 2, quantile, probs = 0.025) + data_minimum_illegal
  lambda_ucb = apply(lambda, 2, quantile, probs = 0.975) + data_minimum_illegal
  lambda_mean = apply(lambda,2, mean) + data_minimum_illegal

  plot_model_fit(taxon_data, model_fit_data, taxon_number, lambda, lambda_median, lambda_mean, lambda_ucb, lambda_lcb, data_minimum_illegal, years, year_restriction)
  
   # Intercept
   alpha <- model_fit$alpha[,taxon_number]
   
   # Slope
   slope <- model_fit$slope[,taxon_number]
   effort_slope <- model_fit$effort_slope
  
  # Extract fitted values plus 95% CIs
  vals = matrix(nrow = dim(lambda)[1], ncol = dim(lambda)[2])

  for (ii in 1:dim(lambda)[1])
  {
    vals[ii,] = alpha[ii] + slope[ii] * taxon_data$legal
  }

  vals_median = apply(vals, 2, quantile, probs = 0.5)
  vals_lcb = apply(vals, 2, quantile, probs = 0.025)
  vals_ucb = apply(vals, 2, quantile, probs = 0.975)

  # Sort values from smallest to largest
  aa = sort(taxon_data$legal, index.return = T)
  legal = aa$x
  illegal = taxon_data$illegal[aa$ix]
  alpha = alpha[aa$ix]
  gamma = gamma[aa$ix]
  vals_median = vals_median[aa$ix]
  vals_lcb = vals_lcb[aa$ix]
  vals_ucb = vals_ucb[aa$ix]
  legal_original = taxon_data$legal_original[aa$ix]


  # Need to rescale back to original legal scale
  rescale_sd = sd(taxon_data$legal_original)
  rescale_centre = mean(taxon_data$legal_original / sd(taxon_data$legal_original))

  ylim_vals = c(min(min(vals_lcb), min(log((illegal - data_minimum_illegal)[(illegal - data_minimum_illegal) != 0]
                                           / reporting_eu[(illegal - data_minimum_illegal) != 0]))),
              max(max(log((illegal - data_minimum_illegal)[(illegal - data_minimum_illegal) != 0]
                          / reporting_eu[(illegal - data_minimum_illegal)!= 0])), max(vals_ucb)))

  
  # Plot the model fit
  plot((legal + rescale_centre) * rescale_sd, (log(illegal - data_minimum_illegal) * reporting - median(effort_slope) * gamma) , pch = 18, cex = 1, xlab = "Legal trade volume", ylab = "", yaxt = "n", ylim = ylim_vals)

  axis(4)
  mtext("Partial residuals (log scale)", side = 4, line = 2, cex = 0.5)
  mtext("Reported legal trade volume", side = 1, line = 2, cex = 0.65)

  lines((legal + rescale_centre) * rescale_sd, vals_median)
  lines((legal + rescale_centre) * rescale_sd, vals_ucb, lty = 2)
  lines((legal + rescale_centre) * rescale_sd, vals_lcb, lty = 2)

}


plot_model_fit <- function(taxon_data, model_fit_data, taxon_number, lambda, lambda_median, lambda_mean, lambda_ucb, lambda_lcb, data_minimum_illegal, years, year_restriction)
{
  # Plot results
  timeseries_years = years[(sum(model_fit_data$timeseries_length[0:(taxon_number - 1)]) + 1):(sum(model_fit_data$timeseries_length[1:(taxon_number)]))]
  ylimits = c(min(taxon_data$illegal, lambda_lcb),max(taxon_data$illegal, lambda_ucb))
  plot(timeseries_years, taxon_data$illegal, pch = 19, ylim = ylimits, yaxt = "n")
  lines(timeseries_years, lambda_median)
  lines(timeseries_years, lambda_ucb, lty = 2)
  lines(timeseries_years, lambda_lcb, lty = 2)

  
  axis(4)
  mtext("Seizure volume", side = 4, line = 2, cex = 0.5)
  # 
  mtext("Year", side = 1, line = 2, cex = 0.65)
  
  
  if (is.na(year_restriction))
  {
    if (cor(taxon_data$illegal[taxon_data$years %in%
                               devector(years, taxon_number, model_fit_data$timeseries_length)] - data_minimum_illegal,
            apply(lambda, 2, quantile, probs = 0.5)) < 0)
    {
      mtext(bquote(paste('Model fit:', 'r'^'2'*'= 0.0', sep = "")), cex = 0.7)
      
    } else
    {     mtext(bquote(paste('Model fit:', 'r'^'2'*'= ', .(round(calculate_rsq(taxon_data$illegal[taxon_data$years %in%
                                                                                                    devector(years, taxon_number, model_fit_data$timeseries_length)],
                                                                               lambda_median), digits = 2)), sep = "")), cex = 0.7)
    }
  } else
  {
    if (cor(taxon_data$illegal[taxon_data$years %in%
                               devector(years, taxon_number, model_fit_data$timeseries_length)],lambda_median) < 0)
    {
      mtext(bquote(paste('Model fit:', 'r'^'2'*'= 0.0', sep = "")), cex = 0.7)
      
    } else
    {
      mtext(bquote(paste('Model fit:', 'r'^'2'*'= ', .(round(calculate_rsq(taxon_data$illegal[taxon_data$years %in%
                                                                                                devector(years, taxon_number, model_fit_data$timeseries_length)],
                                                                           lambda_median), digits = 2)), sep = "")), cex = 0.7)
    } 
  }
}


#Plot model results - hierarchical negative binomial model
plot_model_nb_effort <- function(taxon_number, raw_data, model_fit, model_fit_data, years, year_restriction = NA)
{
  # Extract the raw data
  taxon_data <- DAT[[taxon_number]]

  # Extract the effort data
  gamma = devector(model_fit_data$gamma, taxon_number, model_fit_data$timeseries_length)
  
  
  # Get the minimum illegal data value to rescale back to original
  data_minimum_illegal = min(taxon_data$illegal)
  
  # Extract model results
  # Fitted values
  lambda <- devector_model_results(model_fit$lambda, taxon_number, model_fit_data$timeseries_length)

  lambda_median = apply(lambda, 2, quantile, probs = 0.5) + data_minimum_illegal
  lambda_lcb = apply(lambda, 2, quantile, probs = 0.025) + data_minimum_illegal
  lambda_ucb = apply(lambda, 2, quantile, probs = 0.975) + data_minimum_illegal
  lambda_mean = apply(lambda,2, mean) + data_minimum_illegal
  
  plot_model_fit(taxon_data, model_fit_data, taxon_number, lambda, lambda_median, lambda_mean, lambda_ucb, lambda_lcb, data_minimum_illegal, years, year_restriction)

  # PLOT PARTIAL RESIDUALS
  # Intercept
  alpha <- model_fit$alpha[,taxon_number]
  
  # Slope
  slope <- model_fit$slope[,taxon_number]
  
  # Effort slope
  effort_slope <- model_fit$effort_slope[,taxon_number]
  
  vals = matrix(nrow = dim(lambda)[1], ncol = dim(lambda)[2])
  
  # Extract partial residuals plus 95% CIs. Show on the LP scale
  if (is.na(year_restriction))
  {
    for (ii in 1:dim(lambda)[1])
    {
      vals[ii,] = alpha[ii] + slope[ii] * taxon_data$legal
    }
  } else
  {
    for (ii in 1:dim(lambda)[1])
    {
      vals[ii,] = alpha[ii] + slope[ii] * taxon_data$legal[year_restriction]
    }
  }
  
  vals_median = apply(vals, 2, quantile, probs = 0.5)
  vals_lcb = apply(vals, 2, quantile, probs = 0.025)
  vals_ucb = apply(vals, 2, quantile, probs = 0.975)
  
  # Sort values from smallest to largest
  if (is.na(year_restriction))
  {
    aa = sort(taxon_data$legal, index.return = T)
    legal = aa$x
    illegal = taxon_data$illegal[aa$ix]
    legal_original = taxon_data$legal_original[aa$ix]
  } else
  {
    legal_temp = taxon_data$legal[year_restriction]
    illegal_temp = taxon_data$illegal[year_restriction]
    legal_original_temp = taxon_data$legal_original[year_restriction]
    aa = sort(legal_temp, index.return = T)
    legal = aa$x
    illegal = illegal_temp[aa$ix]
    legal_original = legal_original_temp[aa$ix]
  }
  gamma = gamma[aa$ix]
  vals_median = vals_median[aa$ix]
  vals_lcb = vals_lcb[aa$ix]
  vals_ucb = vals_ucb[aa$ix]
  
  # Need to rescale back to original legal scale
  rescale_sd = sd(taxon_data$legal_original)
  rescale_centre = mean(taxon_data$legal_original / sd(taxon_data$legal_original))
  
  ylim_vals = c(min(min(vals_lcb), min(log((illegal - data_minimum_illegal)[(illegal - data_minimum_illegal) != 0]) - median(effort_slope) * gamma)),
                max(max(log((illegal - data_minimum_illegal)[(illegal - data_minimum_illegal) != 0])), max(vals_ucb)))
  
  # Plot the model fit
  plot((legal + rescale_centre) * rescale_sd, log(illegal - data_minimum_illegal) - median(effort_slope) * gamma, pch = 18, cex = 1, yaxt = "n", ylab = " ", ylim = ylim_vals)
  
  axis(4)
  mtext("Partial residuals (log scale)", side = 4, line = 2, cex = 0.5)
  
  mtext("Reported legal trade volume", side = 1, line = 2, cex = 0.65)
  
  lines((legal + rescale_centre) * rescale_sd, vals_median)
  lines((legal + rescale_centre) * rescale_sd, vals_ucb, lty = 2)
  lines((legal + rescale_centre) * rescale_sd, vals_lcb, lty = 2)
  
  mtext("Estimated slope", cex = 0.7)
}