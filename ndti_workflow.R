# Library
librarian::shelf(
  readr, dplyr, stringr, tidyr, purrr, ggplot2,
  leaflet, sf, terra, progressr, lubridate, zoo
)

# parallel & progress
terraOptions(threads = 32)
handlers("txtprogressbar")



# BLOCK 1: EXTRACT LOCATIONS

# Step 1: Load the shapefile (polygons)
polygons <- st_read("~/projectnb/ccmmf/LandIQ_data/LandIQ_shapefiles/Spatial_Joins/all_crops_2018-2023_same_uids_try5.shp")

# Step 2: Load the point data
points_df <- read_csv("~/projectnb/XinyuanJi/design_points.csv")

# Step 3: Convert points to sf object
# Make sure that the CSV has 'lon' and 'lat' columns
points_sf <- st_as_sf(points_df, coords = c("lon", "lat"), crs = 4326)

# Step 4: Match CRS between polygons and points
if (st_crs(polygons) != st_crs(points_sf)) {
  polygons <- st_transform(polygons, crs = st_crs(points_sf))
}

# Step 5: Spatial join – find which polygon each point lies in
joined <- st_join(points_sf, polygons, join = st_intersects)
# Must check if the ID column is "UniqueID" or "UniqueID.x"
names(joined)

# Step 6: Extract UniqueIDs of matched polygons
# If the point .csv has UniqueID column
matched_ids <- joined$UniqueID.y[!is.na(joined$UniqueID.y)]

# If the point .csv doesn't have UniqueID column
# matched_ids <- joined$UniqueID[!is.na(joined$UniqueID)]

# Step 7: Subset original polygons based on matched IDs
matched_polygons <- polygons %>% 
  filter(UniqueID %in% matched_ids)

# Step 8: Compute bounding boxes for each matched polygon
bbox_df <- matched_polygons %>%
  rowwise() %>%
  mutate(
    bbox = list(st_bbox(geometry))
  ) %>%
  mutate(
    upper_left_lon  = bbox[["xmin"]],
    upper_left_lat  = bbox[["ymax"]],
    lower_right_lon = bbox[["xmax"]],
    lower_right_lat = bbox[["ymin"]]
  ) %>%
  select(UniqueID, upper_left_lon, upper_left_lat, lower_right_lon, lower_right_lat) %>%
  ungroup()

# Step 9: Merge bounding boxes with joined point data
# If the .csv doesn't have a 'UniqueID' column
# final_result <- left_join(joined, st_drop_geometry(bbox_df), by = "UniqueID")

# If the .csv has a 'UniqueID' column
final_result <- left_join(joined, st_drop_geometry(bbox_df), by = c("UniqueID.x" = "UniqueID"))
final_result <- final_result %>%
  rename(UniqueID = UniqueID.x)

# Quick view
# print(head(final_result))

# Optional: save result
# write.csv(final_result, "~/projectnb/XinyuanJi/design_points_polygons.csv", row.names = FALSE)





# BLOCK 2: EXTRACT PHENOLOGY DATES

# Path to your CSV with site bounding boxes
# path <- "/projectnb/dietzelab/XinyuanJi/design_points_polygons.csv"
# 
# # --- Read and clean ---
# raw <- readr::read_lines(path)
# cleaned <- stringr::str_replace_all(
#   raw,
#   ",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
#   ",\"c(\\1,\\2)\","  # quote "c(x,y)"
# )
# 
# anchors <- readr::read_csv(paste(cleaned, collapse = "\n"), show_col_types = FALSE) %>%
#   dplyr::rename(
#     id        = UniqueID,   # adjust if necessary
#     site_name = site_id
#   ) %>%
#   dplyr::select(id, site_name,
#                 upper_left_lon, upper_left_lat,
#                 lower_right_lon, lower_right_lat) %>%
#   dplyr::mutate(
#     across(starts_with(c("upper_left", "lower_right")), as.numeric),
#     center_lon = (upper_left_lon + lower_right_lon)/2,
#     center_lat = (upper_left_lat + lower_right_lat)/2
#   ) %>%
#   dplyr::filter(!is.na(center_lon), !is.na(center_lat))
anchors <- final_result %>%
  dplyr::rename(
    id        = UniqueID, # adjust if necessary
    site_name = site_id
  ) %>%
  dplyr::select(id, site_name,
                upper_left_lon, upper_left_lat,
                lower_right_lon, lower_right_lat) %>%
  dplyr::mutate(
    across(starts_with(c("upper_left", "lower_right")), as.numeric),
    center_lon = (upper_left_lon + lower_right_lon)/2,
    center_lat = (upper_left_lat + lower_right_lat)/2
  ) %>%
  dplyr::filter(!is.na(center_lon), !is.na(center_lat))



# --- Convert each site bounding box to polygon ---
make_roi_vect <- function(site_row) {
  corner_df <- data.frame(
    lon = c(site_row$upper_left_lon, site_row$upper_left_lon,
            site_row$lower_right_lon, site_row$lower_right_lon,
            site_row$upper_left_lon),
    lat = c(site_row$upper_left_lat, site_row$lower_right_lat,
            site_row$lower_right_lat, site_row$upper_left_lat,
            site_row$upper_left_lat)
  )
  st_polygon(list(as.matrix(corner_df)))
}

# Create sf object directly
anchors_sf <- anchors %>%
  rowwise() %>%
  mutate(geometry = st_sfc(make_roi_vect(cur_data()), crs = 4326)) %>%
  ungroup() %>%
  st_as_sf()



# ── Load the MS-LSP grid (already in WGS-84 / EPSG:4326) ───────────────────
grid <- st_read(
  "/projectnb/dietzelab/skanee/ccmmf-phenology/MSLSP_tileGrid.geojson",
  quiet = TRUE
)

# Now anchors_tiles has a column like "Name" or "tile_id" (e.g., "T10SGF")
anchors_tiles <- sf::st_join(
  anchors_sf,
  grid,
  join = sf::st_intersects,
  left = TRUE
) %>%
  # if both anchors_sf and grid have an "id" column
  dplyr::select(id = id.x, site_name, tileID, geometry)
# otherwise
# dplyr::select(id, site_name, tileID, geometry)



hls_dir <- "/projectnb/dietzelab/malmborg/CARB/HLS_data"
hls_files <- list.files(hls_dir, pattern = "MSLSP.*\\.nc$", full.names = TRUE)

# Example: only the files for the tile we care about
tileID <- "T10SGF"  # for example
hls_tile_files <- hls_files[str_detect(basename(hls_files), fixed(tileID))]

# function that extract the phenology
extract_phenology <- function(hls_files, roi_sfc) {
  # Metrics to extract
  metrics <- c("50PCGI", "50PCGD")
  
  # Use first file to get CRS
  template_r <- rast(hls_files[1], subds = metrics)
  
  # Convert sfc_POLYGON directly to SpatVector
  roi_vect <- terra::vect(roi_sfc)
  crs(roi_vect) <- "EPSG:4326"  # explicitly assign WGS84 CRS
  
  # Project to raster CRS
  roi_proj <- terra::project(roi_vect, crs(template_r))
  roi_ext <- terra::ext(roi_proj)
  
  # Function to extract one file
  extract_one <- function(nc_file) {
    yr <- as.integer(str_extract(basename(nc_file), "\\d{4}"))
    r <- try(rast(nc_file, subds = metrics), silent = TRUE)
    if (inherits(r, "try-error")) return(tibble(year = yr, `50PCGI` = NA, `50PCGD` = NA))
    
    cr <- try(crop(r, roi_ext), silent = TRUE)
    if (inherits(cr, "try-error")) return(tibble(year = yr, `50PCGI` = NA, `50PCGD` = NA))
    
    ms <- mask(cr, roi_proj)
    vals <- global(ms, mean, na.rm = TRUE)
    out <- as_tibble(t(vals)) %>% setNames(names(ms))
    out %>% mutate(year = yr, .before = 1)
  }
  
  # Apply to all files
  ts_list <- lapply(hls_files, extract_one)
  ts_raw <- bind_rows(ts_list) %>% arrange(year)
  
  # Convert DOY to dates
  phenology_dates <- ts_raw %>%
    mutate(
      greenup_date = as.Date(paste0(year, "-01-01")) + (`50PCGI` - 1),
      greendown_date = as.Date(paste0(year, "-01-01")) + (`50PCGD` - 1)
    ) %>%
    select(year, greenup_date, greendown_date)
  
  return(phenology_dates)
}



all_phenology <- map_dfr(1:nrow(anchors_tiles), function(i) {
  roi_vect_wgs <- anchors_tiles$geometry[[i]]       # extract polygon
  tileID <- anchors_tiles$tileID[i]                # the HLS tile
  hls_tile_files <- hls_files[str_detect(basename(hls_files), fixed(tileID))]
  
  df <- extract_phenology(hls_tile_files, roi_vect_wgs)
  df <- df %>%
    mutate(
      id = anchors_tiles$id[i],          # create id column
      site_name = anchors_tiles$site_name[i]
    )
})

summary(all_phenology)

# save as csv
# write.csv(all_phenology, "/projectnb/dietzelab/XinyuanJi/green_dates(test).csv", row.names = FALSE)





# BLOCK 3: NDTI



# --- STEP 2. Load Anchor CSV & Pick a Site ---
# Read and clean the CSV
# path <- "/projectnb/dietzelab/XinyuanJi/anchor_sites_locations.csv"
# raw  <- read_lines(path)
# cleaned <- str_replace_all(raw,
#                            ",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
#                            ",\"c(\\1,\\2)\","
# )
# anchors <- read_csv(paste(cleaned, collapse = "\n"), show_col_types = FALSE) %>%
#   select(id, site_name,
#          upper_left_lon, upper_left_lat,
#          lower_right_lon, lower_right_lat) %>%
#   mutate(across(starts_with(c("upper_left","lower_right")), as.numeric),
#          center_lon = (upper_left_lon + lower_right_lon)/2,
#          center_lat = (upper_left_lat + lower_right_lat)/2)  %>%
#   filter(!is.na(center_lon), !is.na(center_lat))



# ── CHUNK 1 : read tillage sheet, keep Treatment column ─────────────────────
# path <- "/projectnb/dietzelab/XinyuanJi/design_points_polygons.csv"
# # path <- "/projectnb/dietzelab/XinyuanJi/till_treatment_polygons.csv"
# 
# raw  <- readr::read_lines(path)
# cleaned <- stringr::str_replace_all(
#   raw,
#   ",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
#   ",\"c(\\1,\\2)\","               # quote the “c(x,y)” blobs (same fix as before)
# )
# 
# anchors <- readr::read_csv(paste(cleaned, collapse = "\n"),
#                            show_col_types = FALSE) %>%
#   
#   # filter anchors by years (remember to check the year name in the .csv)
#   dplyr::filter(year == c(2016,2017,2018,2019,2020,2021,2022,2023,2024)) %>%
#   
#   
#   
#   # ── 1 standardise column names so downstream code is unchanged
#   dplyr::rename(
#     id        = UniqueID,          # ← tillage file’s ID
#     site_name = site_id      # ← adjust if your column is ProjectName
#   ) %>%
#   
#   # ── 2 include Treatment_Control
#   dplyr::select(
#     id, site_name, pft,
#     upper_left_lon, upper_left_lat,
#     lower_right_lon, lower_right_lat
#   ) %>%
#   
#   dplyr::mutate(
#     dplyr::across(
#       dplyr::starts_with(c("upper_left", "lower_right")),
#       as.numeric
#     ),
#     center_lon = (upper_left_lon + lower_right_lon) / 2,
#     center_lat = (upper_left_lat + lower_right_lat) / 2
#   ) %>%
#   dplyr::filter(!is.na(center_lon), !is.na(center_lat))
anchors <- final_result %>%
  # optional year filter if final_result includes a 'year' column
  filter(year %in% 2016:2024) %>%
  rename(
    id        = UniqueID,          # adjust to match your column names
    site_name = site_id
  ) %>%
  select(
    id, site_name, pft,
    upper_left_lon, upper_left_lat,
    lower_right_lon, lower_right_lat
  ) %>%
  mutate(
    across(starts_with(c("upper_left", "lower_right")), as.numeric),
    center_lon = (upper_left_lon + lower_right_lon) / 2,
    center_lat = (upper_left_lat + lower_right_lat) / 2
  ) %>%
  filter(!is.na(center_lon), !is.na(center_lat))

# read Phenology data for future masking-out
# phenology_path <- "/projectnb/dietzelab/XinyuanJi/green_dates.csv"

phenology_dates <- all_phenology %>%
  mutate(
    greenup_date   = as.Date(greenup_date),
    greendown_date = as.Date(greendown_date)
  )#  %>%
#   rename(id = site_id)




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
  
  # # Check for NA values
  # if (nrow(ndti_land) == 0 & nrow(ndti_sent) == 0) {
  #   message("No valid NDTI found for site: ", selected_anchor$id, ". Returning NA row.")
  #   # You may also want to fill pft if you use it downstream
  #   na_row <- tibble(
  #     id = selected_anchor$id,
  #     site_name = selected_anchor$site_name,
  #     Treatment_Control = selected_anchor$Treatment_Control,
  #     pft = selected_anchor$pft,             # if pft is used/selected
  #     year = NA_integer_,
  #     min_ndti = NA_real_,
  #     min_ndti_date = as.Date(NA),
  #     max_before_min = NA_real_,
  #     max_before_min_date = as.Date(NA),
  #     percentage_change = NA_real_
  #   )
  #   return(na_row)
  # }
  # --------------------------------------------------------------------------
  # COMBINE AND (OPTIONALLY) INDICATE SENSOR
  # --------------------------------------------------------------------------
  ndti_df <- bind_rows(ndti_land, ndti_sent)
  
  # Optional: double check for NA values
  # if (nrow(ndti_df) == 0 || all(is.na(ndti_df$mean_ndti))) {
  #   message("No NDTI data present after binding for site: ", selected_anchor$id, ". Returning NA row.")
  #   na_row <- tibble(
  #     id = selected_anchor$id,
  #     site_name = selected_anchor$site_name,
  #     Treatment_Control = selected_anchor$Treatment_Control,
  #     pft = selected_anchor$pft,
  #     year = NA_integer_,
  #     min_ndti = NA_real_,
  #     min_ndti_date = as.Date(NA),
  #     max_before_min = NA_real_,
  #     max_before_min_date = as.Date(NA),
  #     percentage_change = NA_real_
  #   )
  #   return(na_row)
  # }
  
  # --- EARLY EXIT: Check for empty or all-NA NDTI data ---
  if (nrow(ndti_df) == 0 || all(is.na(ndti_df$mean_ndti))) {
    message("No valid NDTI for site ", selected_anchor$id, " — skipping.")
    return(tibble(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft,
      year = NA_integer_,
      max_ndti_growing = NA_real_,
      max_ndti_growing_date = as.Date(NA),
      min_ndti_fallow = NA_real_,
      min_ndti_fallow_date = as.Date(NA),
      ndti_range = NA_real_,
      percentage_change = NA_real_
    ))
  }
  
  # fill the gap
  ndti_df$mean_ndti_filled <- na.approx(ndti_df$mean_ndti, x = ndti_df$date, na.rm = FALSE)
  # smoothing with a Moving Average
  w <- 4
  k <- rep(1/w, w)
  ndti_df$smoothed <- as.numeric(stats::filter(ndti_df$mean_ndti_filled, k, sides = 2))
  
  # --- STEP: Mask NDTI by phenology season ---
  
  # Filter and prepare the smoothed NDTI time series
  ## ---- REPLACE_WITH_THIS: join phenology per-site and compute annual stats ----
  
  # prepare ndti_smooth (this already exists right above)
  ndti_smooth <- ndti_df %>%
    filter(!is.na(smoothed)) %>%
    mutate(
      year = lubridate::year(date),
      id = as.character(selected_anchor$id)   # make sure id is present and comparable
    ) %>%
    filter(year > 2016)
  
  # If no valid smoothed NDTI remain, return NA row for all years
  if (nrow(ndti_smooth) == 0) {
    warning(paste0("No valid NDTI for site ", selected_anchor$id))
    return(tibble(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft,
      year = NA_integer_,
      max_ndti_growing = NA_real_,
      max_ndti_growing_date = as.Date(NA),
      min_ndti_fallow = NA_real_,
      min_ndti_fallow_date = as.Date(NA),
      ndti_range = NA_real_,
      percentage_change = NA_real_
    ))
  }
  
  # ensure phenology_dates has id column (character) and date types
  phenology_dates <- phenology_dates %>%
    mutate(id = as.character(id),
           greenup_date = as.Date(greenup_date),
           greendown_date = as.Date(greendown_date))
  
  # get only this site's phenology rows (faster, avoids accidental cross-site joins)
  pheno_site <- phenology_dates %>% filter(id == as.character(selected_anchor$id))
  
  # --- SANITY CHECK ---
  message("Site ID: ", selected_anchor$id)
  message("Phenology rows found: ", nrow(pheno_site))
  message("Years available in phenology data: ", paste(unique(pheno_site$year), collapse = ", "))
  
  if (nrow(pheno_site) == 0) {
    message("No phenology data for site ", selected_anchor$id)
    return(tibble(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft,
      year = NA_integer_,
      max_ndti_growing = NA_real_,
      max_ndti_growing_date = as.Date(NA),
      min_ndti_fallow = NA_real_,
      min_ndti_fallow_date = as.Date(NA),
      ndti_range = NA_real_,
      percentage_change = NA_real_
    ))
  }
  
  # left-join by year only because pheno_site already filtered to the site
  annual_stats <- ndti_smooth %>%
    left_join(pheno_site %>% select(year, greenup_date, greendown_date), by = "year") %>%
    group_by(year) %>%
    group_modify(~{
      yrdata <- .x
      # if phenology dates are missing for the year, warn and return NA row
      if (is.na(yrdata$greenup_date[1]) | is.na(yrdata$greendown_date[1])) {
        return(tibble(
          year = yrdata$year[1],
          max_ndti_growing = NA_real_,
          max_ndti_growing_date = as.Date(NA),
          min_ndti_fallow = NA_real_,
          min_ndti_fallow_date = as.Date(NA),
          ndti_range = NA_real_
        ))
      }
      
      gu <- yrdata$greenup_date[1]
      gd <- yrdata$greendown_date[1]
      
      # define seasons
      # define seasons
      growing <- yrdata %>% filter(date >= gu & date <= gd)
      fallow  <- yrdata %>% filter(date < gu | date > gd)
      
      message("Year: ", yrdata$year[1])
      message("NDTI points total: ", nrow(yrdata))
      message("Growing season points: ", nrow(growing))
      message("Fallow season points: ", nrow(fallow))
      
      # find the minimum NDTI in fallow season
      min_ndti_fallow <- if (nrow(fallow) > 0) min(fallow$smoothed, na.rm = TRUE) else NA_real_
      min_ndti_fallow_date <- if (!is.na(min_ndti_fallow)) {
        fallow$date[which.min(fallow$smoothed)]
      } else {
        as.Date(NA)
      }
      
      # find max BEFORE the min_ndti date (within fallow season)
      max_before_min <- NA_real_
      max_before_min_date <- as.Date(NA)
      
      if (!is.na(min_ndti_fallow_date)) {
        fallow_before <- fallow %>% filter(date < min_ndti_fallow_date)
        
        if (nrow(fallow_before) > 0) {
          max_before_min <- max(fallow_before$smoothed, na.rm = TRUE)
          max_before_min_date <- fallow_before$date[which.max(fallow_before$smoothed)]
        }
      }
      
      tibble(
        max_ndti_growing       = max_before_min,
        max_ndti_growing_date  = max_before_min_date,
        min_ndti_fallow        = min_ndti_fallow,
        min_ndti_fallow_date   = min_ndti_fallow_date,
        ndti_range             = if (!is.na(max_before_min) && !is.na(min_ndti_fallow))
          max_before_min - min_ndti_fallow
        else NA_real_
      )
      
    }) %>%
    ungroup() %>%
    mutate(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft
    ) %>%
    select(id, site_name, pft, year,
           max_ndti_growing, max_ndti_growing_date,
           min_ndti_fallow, min_ndti_fallow_date, ndti_range)
  
  ## ---- end replacement ----
  print(head(ndti_smooth))
  
  
  # Compute percentage change between growing-season max and fallow-season min
  annual_stats <- annual_stats %>%
    mutate(
      ndti_range = max_ndti_growing - min_ndti_fallow,
      percentage_change = ifelse(
        !is.na(max_ndti_growing) & max_ndti_growing != 0,
        100 * (max_ndti_growing - min_ndti_fallow) / abs(max_ndti_growing),
        NA_real_
      )
    )
  
  # Final return — include site info
  summary_out <- annual_stats %>%
    mutate(
      id = selected_anchor$id,
      site_name = selected_anchor$site_name,
      pft = selected_anchor$pft
    ) %>%
    select(
      id, site_name, pft, year,
      max_ndti_growing, max_ndti_growing_date,
      min_ndti_fallow, min_ndti_fallow_date,
      ndti_range, percentage_change
    )
  
  return(summary_out)
}



# Testing
# test_summary <- map_dfr(10:15, function(i) {
#   process_ndti_for_site(anchors[i, ])
# })
# 
# test_summary





# valid_anchor <- anchors %>%
#   mutate(id = as.character(id)) %>%
#   filter(id == "1024636")
# # Single site test
# result <- process_ndti_for_site(valid_anchor)
# View(result)



# summary
summary_table_raw <- map_dfr(1:nrow(anchors), function(i) {
  process_ndti_for_site(anchors[i, ])
})

# drop the NA rows
summary_table <- summary_table_raw %>%
  filter(!is.na(min_ndti_fallow) & !is.na(percentage_change) & !is.na(year))

cat("Rows before:", nrow(summary_table_raw), "after:", nrow(summary_table),
    "dropped:", nrow(summary_table_raw) - nrow(summary_table), "\n")

summary(summary_table)

# save the result
# write.csv(summary_table_raw, "/projectnb/dietzelab/XinyuanJi/ndti_design_points_summary(test).csv", row.names = FALSE)

