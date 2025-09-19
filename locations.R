# This code takes in a list of points on the map based on lon/lat. It finds the
# polygons that each point lies within and return a list of the polygons' location

library(sf)
library(dplyr)
library(readr)

# Step 1: Load the shapefile (polygons)
polygons <- st_read("~/projectnb/ccmmf/LandIQ_data/LandIQ_shapefiles/Spatial_Joins/all_crops_2018-2023_same_uids_try5.shp")

# Step 2: Load the point data
points_df <- read_csv("~/projectnb/XinyuanJi/new_sites.csv")

# Step 3: Convert points to sf object
# Ensure the CSV has 'lon' and 'lat' columns
points_sf <- st_as_sf(points_df, coords = c("lon", "lat"), crs = 4326)

# Step 4: Match CRS between polygons and points
if (st_crs(polygons) != st_crs(points_sf)) {
  polygons <- st_transform(polygons, crs = st_crs(points_sf))
}

# Step 5: Spatial join – find which polygon each point lies in
joined <- st_join(points_sf, polygons, join = st_within)

# Step 6: Extract UniqueIDs of matched polygons
matched_ids <- joined$UniqueID[!is.na(joined$UniqueID)]

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
    upper_left_lon  = bbox["xmin"],
    upper_left_lat  = bbox["ymax"],
    lower_right_lon = bbox["xmax"],
    lower_right_lat = bbox["ymin"]
  ) %>%
  select(UniqueID, upper_left_lon, upper_left_lat, lower_right_lon, lower_right_lat) %>%
  ungroup()

# Step 9: Merge bounding boxes with joined point data
final_result <- left_join(joined, st_drop_geometry(bbox_df), by = "UniqueID")

print(tail(final_result))

# Optional: save result
# write.csv(final_result, "~/projectnb/XinyuanJi/new_sites_polygons.csv", row.names = FALSE)

