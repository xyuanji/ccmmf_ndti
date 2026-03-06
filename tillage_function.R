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
  
  # 4. Extract base metrics
  final_output <- ndti_smooth %>%
    inner_join(fallow_periods, by = "parcel_id") %>%
    filter(date >= fallow_start & date <= fallow_end) %>%
    group_by(parcel_id, fallow_start) %>%
    summarize(
      year = first(year.y), PFT = first(PFT), season = first(season.y),
      OGMn_date = first(fallow_start),
      ndti_on_OGMn = smoothed[date == first(fallow_start)],
      n_valid_on_OGMn = n_valid[date == first(fallow_start)],
      min_idx = which.min(smoothed),
      minNDTI_date = date[min_idx],
      ndti_on_minNDTI = smoothed[min_idx],
      n_valid_on_minNDTI = n_valid[min_idx],
      .groups = "drop"
    )
  
  # 5. External helper to find neighbors (defined OUTSIDE the pipe)
  get_neighbors_vec <- function(p_id, m_date) {
    valid_data <- ndti_smooth %>% 
      filter(parcel_id == p_id, !is.na(n_valid), n_valid > 0)
    
    prev <- valid_data %>% filter(date < m_date) %>% slice_max(date, n = 1, with_ties = FALSE)
    foll <- valid_data %>% filter(date > m_date) %>% slice_min(date, n = 1, with_ties = FALSE)
    
    tibble(
      prev_date = if(nrow(prev) > 0) prev$date else as.Date(NA),
      prev_n_valid = if(nrow(prev) > 0) prev$n_valid else NA,
      foll_date = if(nrow(foll) > 0) foll$date else as.Date(NA),
      foll_n_valid = if(nrow(foll) > 0) foll$n_valid else NA
    )
  }
  
  # 6. Apply search safely
  neighbors <- purrr::map2_dfr(final_output$parcel_id, final_output$minNDTI_date, get_neighbors_vec)
  
  final_output <- bind_cols(final_output, neighbors) %>%
    select(-fallow_start, -min_idx)
  
  return(final_output)
}

# Testing
example <- tillage_metrics(ndti,phenology)
head(example)

# write_csv(example, "/projectnb/dietzelab/XinyuanJi/tillage_test.csv", na = "NA")
