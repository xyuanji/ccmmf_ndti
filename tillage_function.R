# Library
librarian::shelf(
  data.table, dplyr, tidyr, zoo, readr, ggplot2, arrow, progressr, terra, ggplot2
)

# parallel & progress
terraOptions(threads = 16)
handlers("txtprogressbar")

# ndti  <- fread("/projectnb/dietzelab/XinyuanJi/ndti_subsample_n=400_2018-2023.csv")
# phenology <- fread("/projectnb/dietzelab/XinyuanJi/chosen_pairs_subsample_n=400_2018-2023.csv")

tillage_metrics <- function(ndti_table, phenology_table) {
  
  # 1. Data Prep
  ndti_work <- as_tibble(ndti_table) %>% 
    mutate(date = as.Date(date),
           ss = ifelse(!is.na(n_valid) & n_valid > 0, n_valid * (ndti_sd^2), 0))
  
  convert_julian <- function(julian_days, year) {
    as.Date(ceiling(julian_days) - 1, origin = paste0(year, "-01-01"))
  }
  
  pheno_date <- phenology_table %>%
    transmute(parcel_id, year, OGI_date, OGMn_date)
              # OGI_date  = convert_julian(OGI, year),
              # OGMn_date = convert_julian(OGMn, year))
  
  # 2. Smooth NDTI
  ndti_smooth <- ndti_work %>%
    inner_join(pheno_date, by = c("parcel_id", "year")) %>%
    arrange(parcel_id, date) %>%
    group_by(parcel_id, year, PFT) %>%
    complete(date = seq.Date(min(date), max(date), by = "day")) %>%
    fill(OGMn_date, OGI_date, .direction = "downup") %>% 
    mutate(
      mean_ndti_filled = zoo::na.approx(ndti_mean, x = date, na.rm = FALSE),
      smoothed = as.numeric(stats::filter(mean_ndti_filled, rep(1/4, 4), sides = 2))
    ) %>%
    ungroup()
  
  # 3. Base Metrics Identification
  fallow_periods <- pheno_date %>%
    arrange(parcel_id, OGI_date) %>%
    group_by(parcel_id) %>%
    mutate(fallow_start = OGMn_date, fallow_end = lead(OGI_date)) %>%
    filter(!is.na(fallow_end)) %>%
    ungroup()
  
  base_metrics <- ndti_smooth %>%
    inner_join(fallow_periods, by = "parcel_id") %>%
    filter(date >= fallow_start & date <= fallow_end) %>%
    group_by(parcel_id, fallow_start) %>%
    summarize(
      min_idx          = which.min(smoothed),
      minNDTI_date     = date[min_idx],
      ndti_on_minNDTI  = smoothed[min_idx],
      max_pre_idx      = which.max(ifelse(date <= minNDTI_date, smoothed, -Inf)),
      maxNDTI_pre_date = date[max_pre_idx],
      maxNDTI_pre_min  = smoothed[max_pre_idx],
      ndti_pct_change  = ((maxNDTI_pre_min - ndti_on_minNDTI) / maxNDTI_pre_min) * 100,
      year = first(year.y), PFT = first(PFT),
      OGMn_date = first(fallow_start),
      .groups = "drop"
    )
  
  # 4. Final Robust Helper
  get_metric_details <- function(p_id, target_date, limit_date = NULL) {
    # 1. Get ALL data for this parcel (including the 0/NA n_valid days)
    all_days <- ndti_smooth %>% 
      filter(parcel_id == p_id) %>%
      distinct(date, .keep_all = TRUE)
    
    # 2. Get ONLY valid observation days for neighbor searching
    relevant <- all_days %>% filter(n_valid > 0)
    
    # Identify the specific target row
    target_row <- all_days %>% filter(date == target_date)
    
    # Neighbors for pooling
    prev <- relevant %>% filter(date < target_date) %>% slice_max(date, n = 1, with_ties = FALSE)
    foll <- relevant %>% filter(date > target_date) %>% slice_min(date, n = 1, with_ties = FALSE)
    
    if (!is.null(limit_date)) {
      prev <- prev %>% filter(date >= limit_date)
      foll <- foll %>% filter(date >= limit_date)
    }
    
    # LOGIC CHANGE: 
    # n_on_day: The actual n_valid from the input data for that specific date
    n_on_day <- if(nrow(target_row) > 0) coalesce(target_row$n_valid[1], 0) else 0
    
    # Calculate SD (Direct if n > 0, Pooled if n == 0)
    if (n_on_day > 0) {
      sd_val <- target_row$ndti_sd[1]
    } else {
      comb <- bind_rows(prev, foll)
      sd_val <- if(nrow(comb) > 0 && sum(comb$n_valid, na.rm=TRUE) > 0) {
        sqrt(sum(comb$ss) / sum(comb$n_valid))
      } else {
        NA_real_
      }
    }
    
    list(
      sd       = sd_val, 
      n_total  = n_on_day, # This will now show 0 for interpolated days
      d_before = if(nrow(prev) > 0) prev$date else as.Date(NA),
      n_before = if(nrow(prev) > 0) prev$n_valid else NA_real_,
      d_after  = if(nrow(foll) > 0) foll$date else as.Date(NA),
      n_after  = if(nrow(foll) > 0) foll$n_valid else NA_real_
    )
  }
  
  # 5. Assembly
  final_output <- base_metrics %>%
    rowwise() %>%
    mutate(
      res_max = list(get_metric_details(parcel_id, maxNDTI_pre_date, OGMn_date)),
      res_min = list(get_metric_details(parcel_id, minNDTI_date))
    ) %>%
    ungroup() %>%
    # unnest_wider creates columns like res_max_sd, res_max_n_total, etc.
    unnest_wider(res_max, names_sep = "_") %>%
    unnest_wider(res_min, names_sep = "_") %>%
    select(
      parcel_id, year, PFT, OGMn_date,
      
      # Keep ONLY these for Max
      max_date = maxNDTI_pre_date, 
      max_ndti = maxNDTI_pre_min, 
      
      # Keep everything for Min (as per your previous setup)
      min_date = minNDTI_date, 
      min_ndti = ndti_on_minNDTI, 
      min_n_valid = res_min_n_total, 
      min_sd = res_min_sd,
      
      # Pct Change
      ndti_pct_change,
      
      # Keep validation for Min, but skip them for Max
      min_val_date_before = res_min_d_before, min_val_n_before = res_min_n_before,
      min_val_date_after  = res_min_d_after,  min_val_n_after  = res_min_n_after
    )
  
  return(final_output)
}




library(data.table)
library(arrow)

# Load all final statewide outputs
ndti_files  <- Sys.glob("/projectnb/dietzelab/ccmmf/management/tillage/ndti_v4.1/year=*/ndti_year=*_month=*.parquet")
mslsp_files <- Sys.glob("/projectnb/dietzelab/ccmmf/management/phenology/matched_landiq_mslsp_v4.1/assigned_year=2020.parquet")

ndti  <- rbindlist(lapply(ndti_files, read_parquet),  use.names = TRUE, fill = TRUE)
mslsp <- rbindlist(lapply(mslsp_files, read_parquet), use.names = TRUE, fill = TRUE)

ndti[,  parcel_id := as.character(parcel_id)]
mslsp[, parcel_id := as.character(parcel_id)]

ndti <- ndti[year == 2020]
mslsp <- mslsp[year == 2020]

# keep only parcels present in both
common_ids <- intersect(unique(ndti$parcel_id), unique(mslsp$parcel_id))
ndti  <- ndti[parcel_id %in% common_ids]
mslsp <- mslsp[parcel_id %in% common_ids]

# NDTI as monthly time series
ndti[, date := as.Date(date)]
setorder(ndti, parcel_id, date)

# MSLSP cycle timing (DOY -> Date)
mslsp <- mslsp[, .(
  parcel_id,
  year,
  cycle = mslsp_cycle,                # MSLSP cycle (1/2)
  OGI  = mslsp_OGI,
  OGMn = mslsp_OGMn,
  PFT = landiq_PFT
)]

mslsp[, OGI_date  := as.Date(pmax(1, pmin(366, round(OGI)))  - 1L, origin = paste0(year, "-01-01"))]
mslsp[, OGMn_date := as.Date(pmax(1, pmin(366, round(OGMn))) - 1L, origin = paste0(year, "-01-01"))]

setorder(mslsp, parcel_id, year, cycle)





# Testing
example <- tillage_metrics(ndti,mslsp)
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



# # Box Plot
# 
# # Create the plot object
# clean_output <- example %>%
#   filter(ndti_pct_change > -50, 
#          ndti_pct_change < 200)
# 
# p <- ggplot(clean_output %>% filter(!is.na(ndti_pct_change)), 
#             aes(x = ndti_pct_change, y = PFT, fill = PFT)) +
#   geom_boxplot(outlier.alpha = 0.5) +
# #   geom_jitter(width = 0.1, alpha = 0.2) +
#   theme_minimal() +
#   labs(title = "NDTI Pct Change vs PFT")
# 
# # Force it to show up
# print(p)

