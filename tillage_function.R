# Library
librarian::shelf(
  data.table, dplyr, tidyr, zoo, readr, ggplot2, arrow, progressr
)

# parallel & progress
terraOptions(threads = 16)
handlers("txtprogressbar")

ndti  <- fread("/projectnb/dietzelab/XinyuanJi/ndti_subsample_n=400_2018-2023.csv")
phenology <- fread("/projectnb/dietzelab/XinyuanJi/chosen_pairs_subsample_n=400_2018-2023.csv")

tillage_metrics <- function(ndti_table, phenology_table) {
  
  # 1. Use as_tibble or copy to avoid data.table reference issues
  ndti_work <- as_tibble(ndti_table) %>% mutate(date = as.Date(date))
  
  # Helper: Julian converter
  convert_julian <- function(julian_days, year) {
    as.Date(ceiling(julian_days) - 1, origin = paste0(year, "-01-01"))
  }
  
  pheno_date <- phenology_table %>%
    transmute(parcel_id, year, season,
              OGI_date  = convert_julian(mslsp_OGI, year),
              OGMn_date = convert_julian(mslsp_OGMn, year))
  
  # 2. Smooth NDTI
  ndti_smooth <- ndti_work %>%
    inner_join(pheno_date, by = c("parcel_id", "year")) %>%
    arrange(parcel_id, date) %>%
    group_by(parcel_id, year, PFT, season) %>%
    complete(date = seq.Date(min(date), max(date), by = "day")) %>%
    fill(OGMn_date, OGI_date, .direction = "downup") %>%
    mutate(
      mean_ndti_filled = zoo::na.approx(ndti_mean, x = date, na.rm = FALSE),
      smoothed = as.numeric(stats::filter(mean_ndti_filled, rep(1/4, 4), sides = 2))
    ) %>%
    ungroup()
  
  # 3. Fallow periods
  fallow_periods <- pheno_date %>%
    arrange(parcel_id, OGI_date) %>%
    group_by(parcel_id) %>%
    mutate(fallow_start = OGMn_date, fallow_end = lead(OGI_date)) %>%
    filter(!is.na(fallow_end)) %>%
    ungroup()
  
  final_output <- ndti_smooth %>%
    inner_join(fallow_periods, by = "parcel_id") %>%
    filter(date >= fallow_start & date <= fallow_end) %>%
    group_by(parcel_id, fallow_start) %>%
    summarize(
      # 1. Define the index first so it is available for all following lines
      min_idx            = which.min(smoothed),
      
      # 2. Now you can safely use min_idx
      year               = first(year.y),
      PFT                = first(PFT),
      season             = first(season.y),
      OGMn_date          = first(fallow_start),
      
      ndti_on_OGMn       = smoothed[date == first(fallow_start)],
      n_valid_on_OGMn    = n_valid[date == first(fallow_start)],
      ndti_sd_on_OGMn    = ndti_sd[date == first(fallow_start)],
      
      minNDTI_date       = date[min_idx],
      ndti_on_minNDTI    = smoothed[min_idx],
      n_valid_on_minNDTI = n_valid[min_idx],
      ndti_sd_on_minNDTI = ndti_sd[min_idx],
      
      .groups = "drop"
    )
  
  # 5. External helper to find neighbors (defined OUTSIDE the pipe)
  get_neighbors_vec <- function(p_id, m_date) {
    valid_data <- ndti_smooth %>% 
      filter(parcel_id == p_id, !is.na(n_valid), n_valid > 0)
    
    prev <- valid_data %>% filter(date < m_date) %>% slice_max(date, n = 1, with_ties = FALSE)
    foll <- valid_data %>% filter(date > m_date) %>% slice_min(date, n = 1, with_ties = FALSE)
    
    tibble(
      valid_date_before_min_ndti = if(nrow(prev) > 0) prev$date else as.Date(NA),
      prev_n_valid = if(nrow(prev) > 0) prev$n_valid else NA,
      valid_date_after_min_ndti = if(nrow(foll) > 0) foll$date else as.Date(NA),
      foll_n_valid = if(nrow(foll) > 0) foll$n_valid else NA
    )
  }
  
  # 6. Apply search safely
  neighbors <- purrr::map2_dfr(final_output$parcel_id, final_output$minNDTI_date, get_neighbors_vec)
  
  final_output <- bind_cols(final_output, neighbors) %>%
    select(-fallow_start, -min_idx)
  
  return(final_output)
}




# # Load all final statewide outputs
# ndti_files  <- Sys.glob("/projectnb/dietzelab/ccmmf/management/tillage/ndti_by_parcel_id_v3/year=*/ndti_year=*_month=*.parquet")
# mslsp_files <- Sys.glob("/projectnb/dietzelab/ccmmf/management/phenology/raw_mslsp_v4.1/year=*/mslsp_year=*.parquet")
# 
# ndti  <- rbindlist(lapply(ndti_files, read_parquet),  use.names = TRUE, fill = TRUE)
# mslsp <- rbindlist(lapply(mslsp_files, read_parquet), use.names = TRUE, fill = TRUE)
# 
# 
# 
# ndti[,  parcel_id := as.character(parcel_id)]
# mslsp[, parcel_id := as.character(parcel_id)]
# 
# # keep only parcels present in both
# common_ids <- intersect(unique(ndti$parcel_id), unique(mslsp$parcel_id))
# ndti  <- ndti[parcel_id %in% common_ids]
# mslsp <- mslsp[parcel_id %in% common_ids]
# 
# # NDTI as monthly time series
# ndti[, date := as.Date(date)]
# setorder(ndti, parcel_id, date)
# 
# # MSLSP cycle timing (DOY -> Date)
# mslsp <- mslsp[, .(
#   parcel_id,
#   PFT,
#   year,
#   cycle,                # MSLSP cycle (1/2)
#   OGI  = OGI_mean,
#   OGMn = OGMn_mean
# )]
# 
# mslsp[, OGI_date  := as.Date(pmax(1, pmin(366, round(OGI)))  - 1L, origin = paste0(year, "-01-01"))]
# mslsp[, OGMn_date := as.Date(pmax(1, pmin(366, round(OGMn))) - 1L, origin = paste0(year, "-01-01"))]
# 
# setorder(mslsp, parcel_id, year, cycle)
# 
# 
# 
# # 1. Get unique parcel-PFT pairs
# # We assume 'PFT' exists in 'ndti'. If not, verify where PFT is stored.
# parcel_pft_map <- unique(ndti[, .(parcel_id, PFT)])
# 
# # 2. Sample 10 from each target PFT
# target_pfts <- c("woody", "row", "rice", "hay")
# sampled_ids <- parcel_pft_map[PFT %in% target_pfts, .SD[sample(.N, min(.N, 10))], by = PFT]$parcel_id
# 
# # 3. Subset the tables
# ndti_subset <- ndti[parcel_id %in% sampled_ids]
# mslsp_subset <- mslsp[parcel_id %in% sampled_ids]
# 
# # Rename MSLSP columns to match function expectations
# mslsp_formatted <- copy(mslsp_subset)
# setnames(mslsp_formatted, 
#          old = c("cycle", "OGI_date", "OGMn_date"), 
#          new = c("season", "mslsp_OGI", "mslsp_OGMn"), 
#          skip_absent = TRUE)
# 
# # Note: Your current tillage_metrics function expects 'mslsp_OGI' and 'mslsp_OGMn' 
# # as inputs (which it then converts). If your data is already converted, 
# # you might need to adjust the function slightly. 
# # Based on your previous code, let's assume you pass the raw DOY columns.
# 
# # Run the analysis on the 40-parcel subset
# results_subset <- tillage_metrics(ndti_subset, mslsp_formatted)
# 
# # View the result
# head(results_subset)





# Testing
example <- tillage_metrics(ndti,phenology)
head(example)



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

