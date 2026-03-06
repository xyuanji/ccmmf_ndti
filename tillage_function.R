library(data.table)
library(dplyr)
library(tidyr)
library(zoo)
library(readr)
library(ggplot2)

ndti  <- fread("/projectnb/dietzelab/XinyuanJi/ndti_subsample_n=400_2018-2023.csv")
phenology <- fread("/projectnb/dietzelab/XinyuanJi/chosen_pairs_subsample_n=400_2018-2023.csv")

tillage_metrics <- function(ndti, phenology){

  # STEP 1: extract phenology dates.
  
  # helper function - convert Julian dates
  convert_julian <- function(julian_days, year) {
    # Round up the dates
    julian_integers = ceiling(julian_days)
    
    # Define the start of the target year
    origin <- as.Date(paste0(year, "-01-01"))
    
    # In R, as.Date(0, origin="2023-01-01") is Jan 1st.
    # Usually, Julian day 1 is Jan 1st, so we subtract 1 from the input
    # to align the indices correctly.
    standard_dates <- as.Date(julian_integers - 1, origin = origin)
    
    return(standard_dates)
  }
  
  pheno_date <- phenology %>%
    transmute(
      parcel_id,
      year,
      season,
      OGI_date  = convert_julian(mslsp_OGI, year),
      OGMn_date = convert_julian(mslsp_OGMn, year)
    )
  
  
  
  
  
  # STEP 2: NDTI
  # 1. join the two table
  # Ensure ndti dates are in Date format
  ndti[, date := as.Date(date)]
  
  # Join phenology windows to the NDTI observations
  ndti_pheno <- ndti %>%
    inner_join(pheno_date, by = c("parcel_id", "year")) %>%
    arrange(parcel_id, date)
  
  # 2. Fill and Smooth 
  ndti_smooth <- ndti_pheno %>%
    group_by(parcel_id, year, PFT, season) %>%
    # Create daily sequence to ensure the Moving Average has a consistent time step
    complete(date = seq.Date(min(date), max(date), by = "day")) %>%
    fill(OGMn_date, OGI_date, .direction = "downup") %>%
    
    # A. Linear Interpolation (Your na.approx logic)
    mutate(mean_ndti_filled = na.approx(ndti_mean, x = date, na.rm = FALSE)) %>%
    
    # B. Moving Average Smoothing (Your filter logic: w=4)
    mutate(smoothed = as.numeric(stats::filter(
      mean_ndti_filled, 
      filter = rep(1/4, 4), 
      sides = 2
    ))) %>%
    ungroup()
  
  # 3. Create fallow windows by looking ahead to the next OGI
  fallow_periods <- pheno_date %>%
    arrange(parcel_id, OGI_date) %>%
    group_by(parcel_id) %>%
    mutate(
      fallow_start = OGMn_date,
      fallow_end   = lead(OGI_date)
    ) %>%
    # Keep only cases where a subsequent season exists to define the end date
    filter(!is.na(fallow_end)) %>%
    ungroup()
  
  # 4. Extract final metrics including the n_valid on the minNDTI date
  final_output <- ndti_smooth %>%
    # Join smoothed daily data with the fallow window definitions
    inner_join(fallow_periods, by = "parcel_id") %>%
    
    # Keep only the days falling within the fallow window (OGMn to next OGI)
    filter(date >= fallow_start & date <= fallow_end) %>%
    
    group_by(parcel_id, fallow_start) %>%
    summarize(
      year               = first(year.y),
      PFT                = first(PFT),
      season             = first(season.y),
      OGMn_date          = first(fallow_start),
      
      # 1. Values on the OGMn date
      ndti_on_OGMn       = smoothed[date == first(fallow_start)],
      n_valid_on_OGMn    = n_valid[date == first(fallow_start)],
      
      # 2. Minimum values found during the fallow period
      min_idx            = which.min(smoothed),
      minNDTI_date       = date[min_idx],
      ndti_on_minNDTI    = smoothed[min_idx],
      
      # 3. n_valid on the specific date of minNDTI
      # This will be the original value if an observation existed, or NA if interpolated
      n_valid_on_minNDTI = n_valid[min_idx], 
      
      .groups = "drop"
    ) %>%
    
    # remove helper columns
    select(-fallow_start, -min_idx)
  
  return(final_output)
  
}

# Testing
# example <- tillage_metrics(ndti,phenology)
# head(example)



# Save the result
# write_csv(example, "/projectnb/dietzelab/XinyuanJi/tillage_test.csv", na = "NA")



# Check one parcel
# check_parcel <- ndti_smooth %>%
#    filter(parcel_id == "1193555" & date == as.Date("2020-12-01")) %>%
#    select(date, ndti_mean, mean_ndti_filled, n_valid, smoothed)
# 
# print(check_parcel)

# # Result: date with no ndti_mean can also have ndti_filled and ndti_smoothed




# OPTIONAL: plotting

# # 1. Prepare data: Filter for the parcel and year
# target_parcel <- 1149047
# target_year   <- 2022
# 
# annual_data <- ndti_smooth %>% 
#   filter(parcel_id == target_parcel, year == target_year)
# 
# # 2. Get PFT for the title
# plot_pft <- unique(annual_data$PFT)[1] 
# 
# # 3. Get the specific dates for the vertical lines
# # We filter the final_output for the relevant parcel and year
# line_dates <- final_output %>% 
#   filter(parcel_id == target_parcel, year == target_year)
# 
# # 4. Create the plot
# ggplot(annual_data, aes(x = date, y = smoothed)) +
#   geom_line(color = "darkgreen", size = 1) +
#   # Vertical line 1: OGMn date (Red)
#   geom_vline(data = line_dates, aes(xintercept = OGMn_date), 
#              color = "red", linetype = "dashed", size = 0.8) +
#   # Vertical line 2: Next OGI date (Green) 
#   # Note: You may need to join the fallow_end back to final_output 
#   # or pull from fallow_periods if 'next OGI' isn't explicitly in final_output
#   geom_vline(data = line_dates, aes(xintercept = minNDTI_date), 
#              color = "blue", linetype = "dotted", size = 0.8) +
#   # Vertical line 3: min_ndti date (Blue)
#   labs(
#     title = paste("Annual NDTI Profile", target_year, "PFT: ", plot_pft),
#     subtitle = paste("Parcel_id =", target_parcel),
#     x = "Date",
#     y = "Smoothed NDTI"
#   ) +
#   theme_minimal()




