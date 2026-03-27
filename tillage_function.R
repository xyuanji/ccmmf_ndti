# Library
librarian::shelf(
  data.table, dplyr, tidyr, zoo, readr, ggplot2, arrow, progressr, terra
)

# parallel & progress
terraOptions(threads = 16)
handlers("txtprogressbar")

ndti  <- fread("/projectnb/dietzelab/XinyuanJi/ndti_subsample_n=400_2018-2023.csv")
phenology <- fread("/projectnb/dietzelab/XinyuanJi/chosen_pairs_subsample_n=400_2018-2023.csv")

tillage_metrics <- function(ndti_table, phenology_table) {
  
  # 1. Data Prep
  ndti_work <- as_tibble(ndti_table) %>% 
    mutate(date = as.Date(date),
           ss = ifelse(!is.na(n_valid) & n_valid > 0, n_valid * (ndti_sd^2), 0))
  
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
      year = first(year.y), PFT = first(PFT), season = first(season.y),
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
      parcel_id, year, PFT, season,
      OGMn_date,
      # Max NDTI
      max_date = maxNDTI_pre_date, max_ndti = maxNDTI_pre_min, 
      max_n_valid = res_max_n_total, max_sd = res_max_sd,
      # Min NDTI
      min_date = minNDTI_date, min_ndti = ndti_on_minNDTI, 
      min_n_valid = res_min_n_total, min_sd = res_min_sd,
      # Pct Change
      ndti_pct_change,
      # Validation Max
      max_val_date_before = res_max_d_before, max_val_n_before = res_max_n_before,
      max_val_date_after  = res_max_d_after,  max_val_n_after  = res_max_n_after,
      # Validation Min
      min_val_date_before = res_min_d_before, min_val_n_before = res_min_n_before,
      min_val_date_after  = res_min_d_after,  min_val_n_after  = res_min_n_after
    )
  
  return(final_output)
}




# # Load all final statewide outputs
# ndti_files  <- Sys.glob("/projectnb/dietzelab/ccmmf/management/tillage/ndti_by_parcel_id_v3/year=*/ndti_year=*_month=*.parquet")
# mslsp_files <- Sys.glob("/projectnb/dietzelab/ccmmf/management/phenology/raw_mslsp_v4.1/year=*/mslsp_year=*.parquet")
# 
# ndti  <- rbindlist(lapply(ndti_files, read_parquet),  use.names = TRUE, fill = TRUE)
# mslsp <- rbindlist(lapply(mslsp_files, read_parquet), use.names = TRUE, fill = TRUE)
# 
# colnames(ndti)
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




