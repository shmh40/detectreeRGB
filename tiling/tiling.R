library(magrittr)
library(raster)
library(sf)
library(terra)
library(rgeos)

library(parallel)

ncores <- detectCores()
cl <- makeCluster(ncores)

# Load in crowns - needed for both the single loop example, and the loop further down
crowns <- st_read("~/lustre_scratch/sepilok_data/thomas_sepilok_crowns/thomas_no_intersection.shp")


raster <- 
origin <- origin(raster)
tiledRaster <- splitRaster(r, nx, ny, buffer, path, cl)


stopCluster(cl)