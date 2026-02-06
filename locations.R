# This code finds the big bounding box that contains all the testing fields, and convert the bounding box to a sf geometry object.
library(sf)
library(dplyr)
library(terra)
library(progressr)

# parallel & progress
terraOptions(threads = 16)
handlers("txtprogressbar")

# Load the shapefile (polygons)
polygons <- st_read("~/projectnb/ccmmf/LandIQ_data/LandIQ_shapefiles/Spatial_Joins/all_crops_2018-2023_same_uids_try5.shp")

# Load the point data
points_df <- read.csv("~/projectnb/XinyuanJi/design_points.csv")

# Find the matching polygon
field_matched <- polygons %>%
  filter(UniqueID %in% points_df$UniqueID)

# The big bounding box that contains all the fields (sf object)
big_bbox <- st_bbox(field_matched)

# convert to a sf geometry column
big_bbox_sf <- st_as_sfc(big_bbox)

# keep CRS
st_crs(big_bbox_sf) <- st_crs(field_matched)


