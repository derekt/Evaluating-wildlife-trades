#-----------------------------------------------------------
# Script to load US legal-illegal data and perform robust
# statistics on outliers.
#
# See run_legal_illegal for full details and
# acknowledgements.
#-----------------------------------------------------------


# Directory containing all data. Will need to change this as appropriate
data_directory = "c:/users/derekt/work/research/dropbox/legalillegalanalogues/"
setwd(data_directory)


# Function to find the correct illegal row to match the legal row
find_illegal_row <- function()
{
  if (length(illegal_term_to_model) > 1)
  {
    to_use = which(illegal_term_to_model == legal_term_to_model)
    if (length(to_use) > 1)
    {
      to_use = to_use[which(illegal$Unit[illegal_row_to_fit[to_use]] == unit_to_fit)]
      if (length(to_use) == 1)
      {
        print("Units to use")
        print(illegal$Unit[illegal_row_to_fit[to_use]])
      }      else
        print("Comparable units not found")
    }
    illegal_row_to_fit = illegal_row_to_fit[to_use]
    illegal_term_to_model = illegal_term_to_model[to_use]
  }
  return(
    list(
      illegal_row_to_fit = illegal_row_to_fit,
      illegal_term_to_model = illegal_term_to_model
    )
  )
}

legal <-
  read.csv(paste(data_directory, "us_legal_sp_aggregations_taxa_term_agg-CAV_all-Terms.csv", sep = ""))
illegal <-
  read.csv(paste(data_directory, "lemis_sp_aggregations_taxa_term_agg-CAV_all-Terms.csv", sep = ""))

# Load the CITES Appendix information
CITES_appendix_us = read.csv(paste(data_directory, "us_legal_sp_aggregations_taxa_term_agg-CAV_CITES_Appendix.csv", sep=""))

# Go through each of the legal rows and find a matching row in the illegal data
DAT <- list()
k <- 1
for (ii in 1:dim(legal)[1])
{
  print(ii)
  # Get the data for a single taxon-product
  legal[] <- lapply(legal, as.character)
  illegal[] <- lapply(illegal, as.character)
  
  CITES_appendix_us[] <- lapply(CITES_appendix_us, as.character)
  
  # Legal data extraction
  taxa_to_model = legal$Taxon[ii]
  class_to_model = legal$Class[ii]
  legal_term_to_model = legal$Terms[ii]
  full_legal_term = legal$Full.term.name[ii]
  unit_to_fit = legal$Unit[ii]
  print(taxa_to_model)
  row_to_fit = ii
  
  if (taxa_to_model == CITES_appendix_us$Taxon[ii])
  {
    print(ii)
  } else {
    print("ERROR LOCATING CITES LISTING")
  }
  
  # Illegal data extraction
  illegal_row_to_fit = which(illegal$Taxon == taxa_to_model)
  illegal_term_to_model = illegal$Term[illegal_row_to_fit]
  
  # Loop through as some species have multiple terms and / or units, so find the correct row
  temp_list <- find_illegal_row()
  illegal_row_to_fit = temp_list[1]
  illegal_term_to_model = temp_list[2]
  
  # Note that the US data may have variable number of years of non-zero data so determine the years to be used for the analysis
   if (length(illegal_row_to_fit) > 0)
  {
    first_year = min(which(as.numeric(illegal[unlist(illegal_row_to_fit), 5:dim(illegal)[2]]) > 0))
    year = as.numeric(gsub('X', '', names(illegal)[first_year + 4]))
    
    cols_to_fit = (first_year + 4):dim(illegal)[2]
    illegal_to_fit <-
      as.numeric(illegal[unlist(illegal_row_to_fit), cols_to_fit])
    years_to_model = year:as.numeric(gsub('X', '', names(illegal)[dim(illegal)[2]]))
    
    cols_to_fit = (which(colnames(legal) == paste("X", min(years_to_model), sep =
                                                    ""))):(which(colnames(legal) == paste("X", max(years_to_model), sep = "")))
    legal_to_fit <- as.numeric(legal[row_to_fit, cols_to_fit])
    
    
    CITES_appendix <- CITES_appendix_us[row_to_fit, cols_to_fit]
    
    # Check to make sure that there is non-zero legal trade
    if (sum(legal_to_fit) > 0 &&
        (sum(is.na(legal_to_fit)) < length(legal_to_fit)) &&
        (sum(is.na(illegal_to_fit)) != length(illegal_to_fit))  && (sum(legal_to_fit != 0) > 1))
    {
      # Number of time_points
      N <- length(legal_to_fit)

      remove(outside_range)
      
      # Retain original data for plotting
      legal_to_fit_original = legal_to_fit
      
      if ((sum(legal_to_fit == 0) <= (length(legal_to_fit)/2)))
      {
        # Now use robust statistics to replace any points more than 4x outside the MAD. Legal data.
        median_legal = median(legal_to_fit)
        mad = median(abs(legal_to_fit - median_legal))
        for (jj in 1:N)
        {
          if (abs(legal_to_fit[jj] - median_legal) > (mad * 4))
            if (!exists("outside_range"))
              outside_range = jj
            else
              outside_range = c(outside_range, jj)
        }
        if (exists("outside_range"))
        {
          print("Outside range")
          max_value = max(abs(legal_to_fit[setdiff(1:length(legal_to_fit), outside_range)]))
          num_points_outside_range = length(outside_range)
          min_to_max = outside_range[sort(legal_to_fit[outside_range],
                                          decreasing = FALSE,
                                          index.return = TRUE)$ix]

          for (jj in 1:num_points_outside_range)
          {
            if ((max_value + max_value / (length(outside_range)) * jj) < legal_to_fit[min_to_max[jj]])
              legal_to_fit[min_to_max[jj]] = max_value + max_value / (length(outside_range)) * jj
          }
        }
      } else
      {
        print("Less than half of legal years non-zero")
        print(sum(legal_to_fit == 0))
        print(length(legal_to_fit))
      }
      remove(outside_range)
      
      legal_to_fit_robust = legal_to_fit
      
      # Z-scale the legal data
      sd_to_keep = sd(legal_to_fit)
      legal_to_fit = scale(legal_to_fit, center = TRUE, scale = TRUE)
      
    
      # Now the illegal data
      #if (max(illegal_to_fit) > 10  && (sum(illegal_to_fit == 0) <= (length(illegal_to_fit)/2)))
      illegal_original = illegal_to_fit
      if ((sum(illegal_to_fit == 0)) <= (length(illegal_to_fit)/2))
      {
        # Now use robust statistics to replace any points more than 4x outside the MAD. Illegal data.
        median_illegal = median(illegal_to_fit)
        mad = median(abs(illegal_to_fit - median_illegal))
        for (jj in 1:N)
        {
          if (abs(illegal_to_fit[jj] - median_illegal) > (mad * 4))
            if (!exists("outside_range"))
              outside_range = jj
            else
              outside_range = c(outside_range, jj)
        }
        if (exists("outside_range"))
        {
          print("Outside range")
          max_value = max(illegal_to_fit[setdiff(1:length(illegal_to_fit), outside_range)])
          num_points_outside_range = length(outside_range)
          min_to_max = outside_range[sort(illegal_to_fit[outside_range],
                                          decreasing = FALSE,
                                          index.return = TRUE)$ix]
          for (jj in 1:num_points_outside_range)
          {
            if ((max_value + max_value / (length(outside_range)) * jj) < illegal_to_fit[min_to_max[jj]])
            illegal_to_fit[min_to_max[jj]] = max_value + max_value / (length(outside_range)) * jj
          }
        }
        
        remove(outside_range)
      } else
      {
        print("Less than half of illegal years non-zero")
        print(sum(legal_to_fit == 0))
        print(length(legal_to_fit))
      }
      
      
      DAT[[k]] <-
        list(
          N = N,
          illegal = round(illegal_to_fit),
          illegal_original = as.numeric(illegal_original),
          legal = as.numeric(legal_to_fit),
          legal_original = as.numeric(legal_to_fit_robust),
          legal_original_non_robust = as.numeric(legal_to_fit_original),
          sd_to_keep = sd_to_keep,
          taxa = taxa_to_model,
          class = class_to_model,
          terms = full_legal_term, unit =unit_to_fit, years = years_to_model, CITES = CITES_appendix
        )
      k <- k + 1
      
      
    } else
    {
      print("No legal trade or illegal trade for these years")
    }
    
  } else
  {
    years_to_model = -1
    print("No equivalent illegal data taxa/products/units")
  }
}