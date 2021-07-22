library(magrittr)
library(raster)
library(sf)
library(terra)
library(rgeos)


# Load in crowns - needed for both the single loop example, and the loop further down
Sep_crowns=st_read("~/lustre_scratch/sepilok_data/thomas_sepilok_crowns/thomas_no_intersection.shp")

splitRaster(r, nx, ny, buffer, path, cl)