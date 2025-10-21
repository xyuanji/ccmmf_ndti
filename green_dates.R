library(sf)
library(dplyr)
library(stringr)
library(purrr)
library(terra)
library(progressr)

# parallel & progress
terraOptions(threads = 32)
handlers("txtprogressbar")



# Path to your CSV with site bounding boxes
path <- "/projectnb/dietzelab/XinyuanJi/design_points_polygons.csv"

# --- Read and clean ---
raw <- readr::read_lines(path)
cleaned <- stringr::str_replace_all(
  raw,
  ",c\\(([0-9\\.-]+), *([0-9\\.-]+)\\),",
  ",\"c(\\1,\\2)\","  # quote "c(x,y)"
)

anchors <- readr::read_csv(paste(cleaned, collapse = "\n"), show_col_types = FALSE) %>%
  dplyr::rename(
    id        = UniqueID,   # adjust if necessary
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

# Create sf object directly (no need for map inside mutate)
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
  df$site_id <- anchors_tiles$id[i]
  df$site_name <- anchors_tiles$site_name[i]
  df
})

summary(all_phenology)

# save as csv
write.csv(all_phenology, "/projectnb/dietzelab/XinyuanJi/green_dates.csv", row.names = FALSE)
