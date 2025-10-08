# This code reads a list of polygons, and it calculates the percentage change for
# NDTI between max and min across a range of years (2018-2019 for testing). In
# addition, it'll fit a logistic regression to the data to determine if the calculated
# NDTI change can predict the tillage type (till vs. no till)

# --- STEP 1. Setup ---
# Core libraries
library(zoo)
librarian::shelf(
readr, dplyr, stringr, tidyr, purrr, ggplot2,
leaflet, sf, terra, progressr, lubridate
)

# parallel & progress
terraOptions(threads = 32)
handlers("txtprogressbar")



# --- STEP 2. Load Anchor CSV & Pick a Site ---
# Read and clean the CSV
path <- "/projectnb/dietzelab/XinyuanJi/anchor_sites_locations.csv"
raw  <- read_lines(path)
cleaned <- str_replace_all(raw,
                         ",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
                         ",\"c(\\1,\\2)\","
)
anchors <- read_csv(paste(cleaned, collapse = "\n"), show_col_types = FALSE) %>%
select(id, site_name,
       upper_left_lon, upper_left_lat,
       lower_right_lon, lower_right_lat) %>%
mutate(across(starts_with(c("upper_left","lower_right")), as.numeric),
       center_lon = (upper_left_lon + lower_right_lon)/2,
       center_lat = (upper_left_lat + lower_right_lat)/2)  %>%
filter(!is.na(center_lon), !is.na(center_lat))



# ── CHUNK 1 : read tillage sheet, keep Treatment column ─────────────────────
path <- "/projectnb/dietzelab/XinyuanJi/design_points_polygons.csv"
# path <- "/projectnb/dietzelab/XinyuanJi/till_treatment_polygons.csv"
# path <- "/projectnb/dietzelab/XinyuanJi/anchor_sites_locations.csv"

raw  <- readr::read_lines(path)
cleaned <- stringr::str_replace_all(
raw,
",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
",\"c(\\1,\\2)\","               # quote the “c(x,y)” blobs (same fix as before)
)

anchors <- readr::read_csv(paste(cleaned, collapse = "\n"),
                         show_col_types = FALSE) %>%

# filter anchors by years (remember to check the year name in the .csv)
dplyr::filter(year == c(2016,2017,2018,2019,2020,2021,2022,2023,2024)) %>%



# ── 1 standardise column names so downstream code is unchanged
dplyr::rename(
  id        = UniqueID,          # ← tillage file’s ID
  site_name = site_id      # ← adjust if your column is ProjectName
) %>%

# ── 2 include Treatment_Control
dplyr::select(
  id, site_name, pft,
  upper_left_lon, upper_left_lat,
  lower_right_lon, lower_right_lat
) %>%

dplyr::mutate(
  dplyr::across(
    dplyr::starts_with(c("upper_left", "lower_right")),
    as.numeric
  ),
  center_lon = (upper_left_lon + lower_right_lon) / 2,
  center_lat = (upper_left_lat + lower_right_lat) / 2
) %>%
dplyr::filter(!is.na(center_lon), !is.na(center_lat))



# Make the whole code a function
process_ndti_for_site <- function(selected_anchor) {
  site_id    <- selected_anchor$id
  site_name  <- selected_anchor$site_name
  
  

  # --- STEP 3. Visualize All Anchor Sites ---
  # Please refer to Step 3 in:
  # //projectnb/dietzelab/XinyuanJi/ndti-evi-timeseries(Sarah)-L+S.qmd
  
  
  
  # --- STEP 4. Build Exact tile-polygon ROI ---  
  # A) Build a closed polygon in WGS84 from the CSV bbox
  corner_df <- data.frame(
    lon = c(
      selected_anchor$upper_left_lon,
      selected_anchor$upper_left_lon,
      selected_anchor$lower_right_lon,
      selected_anchor$lower_right_lon,
      selected_anchor$upper_left_lon  # close ring
    ),
    lat = c(
      selected_anchor$upper_left_lat,
      selected_anchor$lower_right_lat,
      selected_anchor$lower_right_lat,
      selected_anchor$upper_left_lat,
      selected_anchor$upper_left_lat  # close ring
    )
  )
  
  # B) Make an sf and then a terra SpatVector
  bbox_sf      <- st_sfc(st_polygon(list(as.matrix(corner_df))), crs = 4326)
  roi_vect_wgs <- terra::vect(bbox_sf)
  
  
  
  # --- STEP 5. Determine Which HLS Tile Covers This Site ---
  ## Determine which HLS tile covers this site  ────────────────────────────
  ## (robust: uses ROI polygon → largest-overlap → nearest fallback)
  
  # ── Load the MS-LSP grid (already in WGS-84 / EPSG:4326) ───────────────────
  grid <- st_read(
    "/projectnb/dietzelab/skanee/ccmmf-phenology/MSLSP_tileGrid.geojson",
    quiet = TRUE
  )
  
  # ── Build the anchor-site polygon in the grid CRS (4326) -------------------
  bbox_anchor <- st_transform(bbox_sf, st_crs(grid))   # bbox_sf from Step 5 later
  
  # ── 1) tiles whose polygons intersect the ROI -----------------------------
  hits <- st_intersects(bbox_anchor, grid) |> unlist()
  
  if (length(hits) == 1) {
    
    tileID <- grid$tileID[hits]                        # exactly one tile
    
  } else if (length(hits) > 1) {
    
    # 2) choose the one with the largest overlapping area
    overlap_areas <- st_area(st_intersection(grid[hits, ], bbox_anchor))
    tileID <- grid$tileID[hits[which.max(overlap_areas)]]
    
  } else {
    
    # 3) ROI outside every grid polygon – pick the nearest grid feature
    idx    <- st_nearest_feature(st_centroid(bbox_anchor), grid)
    tileID <- grid$tileID[idx]
  }
  
  site_id   <- selected_anchor$id
  site_name <- selected_anchor$site_name
  crop_type <- site_name     # keep your original variable
  
  cat("Selected site:", site_id, "-", site_name, "\n")
  cat("Using HLS tile:", tileID, "\n")
  
  
  
  # --- STEP 6. List & Filter HLS NetCDFs for That Tile ---
  hls_dir   <- "/projectnb/dietzelab/malmborg/CARB/HLS_data"
  hls_files <- list.files(hls_dir, pattern = "MSLSP.*\\.nc$", full.names = TRUE)
  
  # Filter by tileID in filename (fast, no raster I/O)
  hls_tile_files <- hls_files[str_detect(basename(hls_files), fixed(tileID))]
  
  cat("Found", length(hls_tile_files), "HLS files for tile", tileID, "\n")
  
  
  
  # --- STEP 7. Precompute Projected ROI Extent ---
  # Load one file to get its CRS
  metrics     <- c("50PCGI","50PCGD")#,"OGMn", "OGI_2", "Peak_2", "OGMn_2")
  template_r  <- rast(hls_tile_files[1], subds = metrics)
  roi_proj    <- project(roi_vect_wgs, crs(template_r))
  roi_ext     <- ext(roi_proj)
  
  
  
  # --- STEP 8. Extract phenology metrics ---
  extract_one <- function(nc) {
    yr <- as.integer(str_extract(basename(nc), "\\d{4}"))
    r  <- try(rast(nc, subds = metrics), silent = TRUE)
    if (inherits(r, "try-error")) {
      message("Skipping load: ", basename(nc))
      return(tibble(year = yr, OGI = NA, Peak = NA, OGMn = NA))
    }
    cr <- try(crop(r, roi_ext), silent = TRUE)
    if (inherits(cr, "try-error")) {
      message(" Crop failed for ", basename(nc))
      return(tibble(year = yr, OGI = NA, Peak = NA, OGMn = NA))
    }
    ms  <- mask(cr, roi_proj)
    vals <- global(ms, mean, na.rm = TRUE)
    out  <- as_tibble(t(vals)) %>% setNames(names(ms))
    out %>% mutate(year = yr, .before = 1)
  }
  
  # Sequential extraction with a text progress bar
  ts_list <- vector("list", length(hls_tile_files))
  pb      <- txtProgressBar(0, length(hls_tile_files), style = 3)
  for (i in seq_along(hls_tile_files)) {
    ts_list[[i]] <- extract_one(hls_tile_files[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ts_raw <- bind_rows(ts_list) %>% arrange(year)
  
  
  
  # --- STEP 9. Pivot long & plot phenology ---
  ts_long <- ts_raw %>%
    pivot_longer(
      cols      = all_of(metrics),
      names_to  = "metric",
      values_to = "doy"
    ) %>%
    # ensure doy is integer, then add to the first of each year
    mutate(
      doy = as.integer(doy),
      date = as.Date(paste0(year, "-01-01")) + (doy - 1)
    )
  
  
  
  ## --- 10. Extract Mean NDTI over same tile polygon ---
  
  # .tif directory
  # tif_dir <- c("/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2018",
  #              "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2019",
  #              "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2018",
  #              "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2019"
  # )
  tif_dir <- c("/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2016",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2017",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2018",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2019",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2020",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2021",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2022",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2023",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL/2024",
               
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2016",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2017",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2018",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2019",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2020",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2021",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2022",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2023",
               "/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS/2024")
  
  # .Fmask directory
  # fmask_folder <- c("~/projectnb/XinyuanJi/State_of_California_HLSL/2018_Fmask",
  #                   "~/projectnb/XinyuanJi/State_of_California_HLSL/2019_Fmask",
  #                   "~/projectnb/XinyuanJi/State_of_California_HLSS/2018_Fmask",
  #                   "~/projectnb/XinyuanJi/State_of_California_HLSS/2019_Fmask"
  # )
  fmask_folder <- c("~/projectnb/XinyuanJi/State_of_California_HLSL/2016_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2017_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2018_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2019_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2020_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2021_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2022_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2023_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSL/2024_Fmask",
                    
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2016_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2017_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2018_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2019_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2020_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2021_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2022_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2023_Fmask",
                    "~/projectnb/XinyuanJi/State_of_California_HLSS/2024_Fmask")
  
  # read both Landsat & Sentinel
  ref_file <- list.files(
    tif_dir,
    pattern = "B(06|07|11|12).*\\.tif$",
    full.names = TRUE
  )[1]
  
  # setup the ROI
  r0       <- rast(list.files(tif_dir, pattern="B11.*\\.tif$", full.names=TRUE)[1])
  roi_ndti <- project(roi_vect_wgs, crs(r0)) |> ext()
  
  # --------------------------------------------------------------------------
  # 1) KEEP ONLY THIS SITE’S TILE  (e.g. T10SGF), *not* every B06 in 2018
  # --------------------------------------------------------------------------
  # fine relevant files
  b11_all <- list.files(
    tif_dir,
    pattern    = paste0("T", tileID, ".*B11.*\\.tif$"),
    full.names = TRUE
  )
  # ...file cleaning, etc (copy your existing cleaning/overlap code here) ...
  
  # extracting dates from file name
  b11_dates <- as.Date(stringr::str_extract(b11_all, "\\d{7}"), "%Y%j")
  
  # Date filtering & NA-removal
  keep_idx <- !is.na(b11_dates) &
    b11_dates >= as.Date("2016-01-01") &
    b11_dates <= as.Date("2024-12-31")
  
  b11_sub   <- b11_all[keep_idx]
  b11_dates <- b11_dates[keep_idx]
  # ...date filtering etc as before...
  
  # Pre-allocating Results Vectors
  n_sen <- length(b11_sub)
  dates_sen  <- as.Date(rep(NA, n_sen))
  ndvals_sen <- numeric(n_sen)
  
  # for-loop processing Sentinel-2 data
  pb <- txtProgressBar(0, n_sen, style = 3)
  for (i in seq_along(b11_sub)) {
    b11 <- b11_sub[i]
    b12 <- sub("B11", "B12", b11)
    if (!file.exists(b12)) { setTxtProgressBar(pb, i); next }
    scene      <- basename(b11)
    scene_year <- stringr::str_extract(scene, "\\d{4}")
    fmask_dir  <- file.path("/projectnb/dietzelab/XinyuanJi/State_of_California_HLSS", paste0(scene_year, "_Fmask"))
    fmask_file <- file.path(fmask_dir, sub("B11.tif$", "Fmask.tif", scene))
    
    # Skip if missing
    if (!file.exists(b12) || !file.exists(fmask_file)) {
      setTxtProgressBar(pb, i)
      next
    }
    
    r11   <- try(crop(rast(b11), roi_ndti), silent = TRUE)
    r12   <- try(crop(rast(b12), roi_ndti), silent = TRUE)
    fmask <- try(crop(rast(fmask_file), roi_ndti), silent = TRUE)
    if (inherits(r11, "try-error") || inherits(r12, "try-error") || inherits(fmask, "try-error")) {
      setTxtProgressBar(pb, i)
      next
    }
    
    # calculate NDTI
    ndti <- (r11 - r12) / (r11 + r12)
    
    # cloud/shadow masking
    # 1 = cloud; 3 = cloud shadow; 4 = snow
    badmask <- (bitwAnd(values(fmask), bitwShiftL(1, 1)) != 0) |
      (bitwAnd(values(fmask), bitwShiftL(1, 3)) != 0) |
      (bitwAnd(values(fmask), bitwShiftL(1, 4)) != 0)
    ndti[badmask] <- NA
    
    # calculate the mean of all pixels in the raster
    nd_mean <- global(ndti, fun = mean, na.rm = TRUE)[1, 1]
    
    # storing the result
    if (is.na(nd_mean)) { setTxtProgressBar(pb, i); next }
    ndvals_sen[i] <- nd_mean
    dates_sen[i]  <- b11_dates[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ndti_sent <- tibble(
    date      = dates_sen,
    mean_ndti = ndvals_sen
  ) %>%
    filter(!is.na(date) & !is.na(mean_ndti)) %>% # remove invalid rows
    arrange(date) # sort the final result by date
  
  # --------------------------------------------------------------------------
  # LANDSAT: B06/B07; a copy of the code above but for Band 6 & 7
  # --------------------------------------------------------------------------
  b06_all <- list.files(
    tif_dir,
    pattern    = paste0("T", tileID, ".*B06.*\\.tif$"),
    full.names = TRUE
  )
  # ...copy file cleaning, overlap check code from above, but with B06 instead of B11 ...
  
  b06_dates <- as.Date(stringr::str_extract(b06_all, "\\d{7}"), "%Y%j")
  keep_idx <- !is.na(b06_dates) &
    b06_dates >= as.Date("2016-01-01") &
    b06_dates <= as.Date("2024-12-31")
  b06_sub   <- b06_all[keep_idx]
  b06_dates <- b06_dates[keep_idx]
  # ...date filtering...
  
  n_lan <- length(b06_sub)
  dates_lan  <- as.Date(rep(NA, n_lan))
  ndvals_lan <- numeric(n_lan)
  
  pb <- txtProgressBar(0, n_lan, style = 3)
  for (i in seq_along(b06_sub)) {
    b06 <- b06_sub[i]
    b07 <- sub("B06", "B07", b06)
    if (!file.exists(b07)) { setTxtProgressBar(pb, i); next }
    
    scene      <- basename(b06)
    scene_year <- stringr::str_extract(scene, "\\d{4}")
    fmask_dir  <- file.path("/projectnb/dietzelab/XinyuanJi/State_of_California_HLSL", paste0(scene_year, "_Fmask"))
    fmask_file <- file.path(fmask_dir, sub("B06.tif$", "Fmask.tif", scene))
    
    if (!file.exists(fmask_file)) {
      setTxtProgressBar(pb, i)
      next
    }
    
    r6    <- try(crop(rast(b06), roi_ndti), silent = TRUE)
    r7    <- try(crop(rast(b07), roi_ndti), silent = TRUE)
    fmask <- try(crop(rast(fmask_file), roi_ndti), silent = TRUE)
    if (inherits(r6, "try-error") || inherits(r7, "try-error") || inherits(fmask, "try-error")) {
      setTxtProgressBar(pb, i)
      next
    }
    
    # calculate NDTI
    ndti <- (r6 - r7) / (r6 + r7)
    
    badmask <- (bitwAnd(values(fmask), bitwShiftL(1, 1)) != 0) |
      (bitwAnd(values(fmask), bitwShiftL(1, 3)) != 0) |
      (bitwAnd(values(fmask), bitwShiftL(1, 4)) != 0)
    ndti[badmask] <- NA
    
    nd_mean <- global(ndti, fun = mean, na.rm = TRUE)[1, 1]
    
    
    if (is.na(nd_mean)) { setTxtProgressBar(pb, i); next }
    ndvals_lan[i] <- nd_mean
    dates_lan[i]  <- b06_dates[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ndti_land <- tibble(
    date      = dates_lan,
    mean_ndti = ndvals_lan
  ) %>%
    filter(!is.na(date) & !is.na(mean_ndti)) %>%
    arrange(date)
  
  # Check for NA values
  if (nrow(ndti_land) == 0 & nrow(ndti_sent) == 0) {
    message("No valid NDTI found for site: ", selected_anchor$id, ". Returning NA row.")
    # You may also want to fill pft if you use it downstream
    na_row <- tibble(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      Treatment_Control = selected_anchor$Treatment_Control,
      pft = selected_anchor$pft,             # if pft is used/selected
      year = NA_integer_,
      min_ndti = NA_real_,
      min_ndti_date = as.Date(NA),
      max_before_min = NA_real_,
      max_before_min_date = as.Date(NA),
      percentage_change = NA_real_
    )
    return(na_row)
  }
  # --------------------------------------------------------------------------
  # COMBINE AND (OPTIONALLY) INDICATE SENSOR
  # --------------------------------------------------------------------------
  ndti_df <- bind_rows(ndti_land, ndti_sent)
  
  # Optional: double check for NA values
  if (nrow(ndti_df) == 0 || all(is.na(ndti_df$mean_ndti))) {
    message("No NDTI data present after binding for site: ", selected_anchor$id, ". Returning NA row.")
    na_row <- tibble(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      Treatment_Control = selected_anchor$Treatment_Control,
      pft = selected_anchor$pft,
      year = NA_integer_,
      min_ndti = NA_real_,
      min_ndti_date = as.Date(NA),
      max_before_min = NA_real_,
      max_before_min_date = as.Date(NA),
      percentage_change = NA_real_
    )
    return(na_row)
  }
  
  
  # fill the gap
  ndti_df$mean_ndti_filled <- na.approx(ndti_df$mean_ndti, x = ndti_df$date, na.rm = FALSE)
  # smoothing with a Moving Average
  w <- 4
  k <- rep(1/w, w)
  ndti_df$smoothed <- as.numeric(stats::filter(ndti_df$mean_ndti_filled, k, sides = 2))
  
  # calculate "max before min" logic for each year
  ndti_smooth <- ndti_df[!is.na(ndti_df$smoothed), ]
  ndti_smooth <- ndti_smooth %>% mutate(year = year(date))
  
  annual_stats <- ndti_smooth %>%
    mutate(year = lubridate::year(date)) %>%
    filter(year > 2016) %>% # remove first year
    group_by(year) %>%
    group_modify(~{
      yrdata <- .x
      min_i <- which.min(yrdata$smoothed)
      min_val <- yrdata$smoothed[min_i]
      min_date <- yrdata$date[min_i]
      
      before_min <- yrdata %>% filter(date < min_date)
      if (nrow(before_min) > 0 && !all(is.na(before_min$smoothed))) {
        max_before_min_i <- which.max(before_min$smoothed)
        max_before_min   <- before_min$smoothed[max_before_min_i]
        max_before_min_date <- before_min$date[max_before_min_i]
      } else {
        max_before_min <- NA_real_
        max_before_min_date <- NA
      }
      
      tibble(
        min_ndti = min_val,
        min_ndti_date = min_date,
        max_before_min = max_before_min,
        max_before_min_date = max_before_min_date,
        percentage_change = ifelse(!is.na(max_before_min) && max_before_min != 0,
                                   (max_before_min - min_val) / max_before_min * 100,
                                   NA_real_)
      )
    }) %>%
    
    mutate(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft,
    ) %>%
    select(id, site_name, pft, year,
           min_ndti, min_ndti_date,
           max_before_min, max_before_min_date, percentage_change)
  
  return(annual_stats)
  
}
  
  

# Testing
# anchor_num <- 11
# selected_anchor <- anchors[anchor_num, ]
# process_ndti_for_site(selected_anchor)



# summary
summary_table_raw <- anchors %>%
  split(.$id) %>%
  map_dfr(process_ndti_for_site)

# drop the NA rows
summary_table <- summary_table_raw %>%
  filter(!is.na(min_ndti) & !is.na(percentage_change) & !is.na(year))

cat("Rows before:", nrow(summary_table_raw), "after:", nrow(summary_table),
    "dropped:", nrow(summary_table_raw) - nrow(summary_table), "\n")

summary(summary_table)

# save the result
# write.csv(summary_table, "/projectnb/dietzelab/XinyuanJi/ndti_design_points_summary.csv", row.names = FALSE)



# Plot Treatment Control vs. Percentage_change
# Create a new till type column
# summary_table <- summary_table %>%
#   mutate(
#     till_type = case_when(
#       grepl("(?i)no[ _-]*till", Treatment_Control) ~ "no till",
#       grepl("(?i)till", Treatment_Control) ~ "till",
#       TRUE ~ NA_character_
#     )
#   )

# Plot the data
ggplot(summary_table, aes(x = percentage_change, y = pft, color = pft)) +
  geom_jitter(width = 0, height = 0.1, size = 2, alpha = 0.7) +
  labs(x = "Percentage Change", y = "Crop Type", title = "Annual vs. Woody") +
  scale_color_manual(values = c("no till" = "blue", "till" = "red")) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 100))



# Logistic regression
# Filter out NA values and set the factor levels explicitly
summary_table_filtered <- summary_table %>%
  filter(!is.na(pft), !is.na(percentage_change)) %>%
  mutate(
    pft_factor = factor(pft, levels = c("annual", "woody")),
    pft_numeric = as.numeric(pft_factor) - 1
  )



# --- Prepare the plotting data ---
plot_data <- summary_table_filtered %>%
  filter(!is.na(percentage_change), !is.na(pft)) %>%
  mutate(
    pft_factor = case_when(
      grepl("(?i)annual", pft) ~ "annual",
      grepl("(?i)woody", pft)  ~ "woody",
      TRUE ~ NA_character_
    ),
    pft_factor = factor(pft_factor, levels = c("annual", "woody")),
    pft_numeric = as.numeric(pft_factor) - 1
  ) %>%
  filter(!is.na(pft_factor)) %>%       # drop anything that didn't match
  filter(percentage_change >= 20, percentage_change <= 70)

ggplot(plot_data, aes(x = percentage_change, y = pft_numeric)) +
  geom_jitter(aes(color = pft_factor), height = 0.1, alpha = 0.7, size = 2) +
  stat_smooth(
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE,
    color = "black",
    size = 1.2
  ) +
  scale_color_manual(values = c("annual" = "blue", "woody" = "red")) +
  scale_y_continuous(
    limits = c(-0.1, 1.1),
    breaks = c(0, 1),
    labels = c("Annual", "Woody")
  ) +
  labs(
    title = "Logistic Regression: PFT vs. Percentage Change",
    x = "Percentage Change",
    y = "Predicted Probability of 'Woody'",
    color = "PFT"
  ) +
  theme_minimal()



# Prepare the data for logistic regression
logit_data <- summary_table_filtered %>%
  filter(!is.na(percentage_change), !is.na(pft)) %>%
  mutate(
    pft_factor = case_when(
      grepl("(?i)annual", pft) ~ "annual",
      grepl("(?i)woody", pft)  ~ "woody",
      TRUE ~ NA_character_
    ),
    pft_factor = factor(pft_factor, levels = c("annual", "woody")),
    pft_numeric = as.numeric(pft_factor) - 1
  ) %>%
  filter(!is.na(pft_factor))  # remove unmatched rows

# Fit the logistic regression model
logistic_model <- glm(
  pft_factor ~ percentage_change,
  data = logit_data,
  family = "binomial"
)

# Print a summary of the model
cat("\n\n--- Logistic Regression Model Summary ---\n")
print(summary(logistic_model))

# Calculate and print the odds ratio for easier interpretation
coef_pc <- coef(logistic_model)["percentage_change"]
odds_ratio_pc <- exp(coef_pc)

cat("\n\n--- Interpretation ---\n")
cat("Odds Ratio for percentage_change:", round(odds_ratio_pc, 3), "\n")
cat("Interpretation: For every 1-unit increase in 'percentage_change',\n")
cat("the odds of a field being 'Woody perennial crop' are multiplied by", 
    round(odds_ratio_pc, 3), ".\n")


