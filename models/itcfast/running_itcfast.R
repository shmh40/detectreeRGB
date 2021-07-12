# This script extracts CHM values for the tree crown polygons. 
#library(profvis)
# Packages  =========
library(dplyr); library(tidyr); library(magrittr); library(rgdal)
library(ggplot2); library(RColorBrewer); library(readxl)

# LiDAR packages
library(raster); library(sf)
library(lidR);  library(UAVforestR) # This is just for plotting, you probably don't need it
library(rasterVis) # for gplot simplicity with rasters
library(spatialEco) # For Sobal transform
library(velox) # This is the key package for speed, try the following installations.
#devtools::install_github("hunzikp/velox")
#devtools::install_version("velox", version = "0.2.0")

# Make Danum difference tiff
Dan_2020=raster(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/danum_lidar/Danum_2020_CHM_g10_sub0.2_1m.tif"))
Dan_2014=raster(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/danum_lidar/Danum_2014_CHM_g10_sub0.2_1m.tif"))

Dan_diff= Dan_2020 - Dan_2014

# Set up look up table
alm<-read_xlsx("C:/Users/sebhi/ai4er/mres_project/work/data/GlobalAllometricDatabase.xlsx", sheet ="Data")
alm %<>%  filter(Biogeographic_zone == "Indo-Malaya", Biome == "Tropical forests")
lut<-rq_lut(x=alm$H, y=alm$CD/2, log=TRUE)


#writeRaster(Dan_diff,"C:/Users/sebhi/ai4er/mres_project/work/data/danum_lidar/Danum_CHM_diff_1m", format="GTiff")

# Load in the LiDAR data as a stack
# ALS is just the path to CHMs
Sep = raster::stack(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2014_coarse_CHM_g10_sub0.01_0.5m.tif"), paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sepilok_CHM_diff_0.5m.tif"))
#Dan = raster::stack(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/danum_lidar/Danum_2014_CHM_g10_sub0.2_1m.tif"), paste0("C:/Users/sebhi/ai4er/mres_project/work/data/danum_lidar/Danum_CHM_diff_1m.tif"))



# Load in the crown shapefiles
# These detectrons do not work! 
#Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/full_sepilok_deploy/606800_column.shp")
Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/full_sepilok_deploy/603800_column.shp")
#Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/sepilok_603450_648000.shp")
plot(Sep_detectron_crowns)
Sep_detectron_crowns$perimeter=as.numeric(polyPerimeter(Sep_detectron_crowns))
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$pixelvalue>100,] 
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$perimeter<300,]
plot(Sep_detectron_crowns)


### CROP the LiDAR data
bound = st_bbox(Sep_detectron_crowns)
Sep_cropped = raster::crop(Sep, bound)
#plot(Sep_cropped)

### This is the cropped area we care about. We now want to run ITCfast on it.


cropped_chm = Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m
plot(cropped_chm)

# Toms UAVforestR ############


# Set up imagery
chm_blur<-blur(cropped_chm)
#plot(chm_blur)
chm_sobel<-sobel_edge(cropped_chm)
#plot(chm_sobel)
imagery=chm_blur; im_sobel=chm_sobel;  pypath = NULL

# Ok, this works thus far. Now let's try to include the search params.

#out<-vector(mode="list", length=nrow(params))




# Separate treetops and segmentation
#THRESHSeed=0.6 # keep this set

#THRESHCrown=0.9 # optimize this between 0.6-0.95
tau = 30 # Optimize between 40-99
#SOBELstr=60 # Optimize this between 20-60
lm.searchwin = NULL
specT=30

# Extracts the image data as a matrix:
img <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[,], byrow = FALSE)
Sobel <- matrix(dim(imagery)[2], dim(imagery)[1], data = im_sobel[,], byrow = FALSE)
img <- img[1:dim(imagery)[2], dim(imagery)[1]:1]
Sobel <- Sobel[1:dim(imagery)[2], dim(imagery)[1]:1]
img[is.na(img)] <- 0 # Sets any nas to 0s.
img[img < specT] <- 0 # Any values beneath the minimum height are set to 0s.

# Finds the local maxima:
coordSeeds <-  detect.maxima(img, res(imagery)[1], 
                             lm.searchwin = lm.searchwin, lut = lut, tau = tau)

# define table for parameter search

THRESHSeed_vec<-0.6
THRESHCrown_vec<-0.6
SOBELstr_vec<-42.5
tau_vec<-91.5
#searchwin_vec<-seq(from=3, to=9, by=6)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
                    THRESHCrown=THRESHCrown_vec,
                    SOBELstr=SOBELstr_vec,
                    tau=tau_vec)
# searchwin=searchwin_vec)

# Segment crowns from local maxima
coordCrown <- segment.crowns(
  x = img,
  x.Sobel = Sobel,
  coordSeeds,
  THRESHSeed = params$THRESHSeed[1],
  THRESHCrown = params$THRESHCrown[1],
  SOBELstr = params$SOBELstr[1],
  scale = res(imagery)[1],
  lut = lut,
  tau = params$tau[1] # This one is now set and not related to tau in coordSeeds
)

crowns_sp<-crowns_to_spatial(coordCrown, imagery, pypath, convex.hull = FALSE)
#maxima_sp<-maxima_to_spatial(coordCrown, imagery, pypath)
#crowns_sp<-crown_max_combine(crowns_sp, maxima_sp)
plot(crowns_sp)

file_name<-paste('C:/Users/sebhi/ai4er/mres_project/work/data/shape/sepilok_deploy/603800_column_seed_', params[1,1],
                 '_crown_', params[1,2],
                 '_sobel_', params[1,3],
                 '_tau_', params[1,4],
                 '_searchwin_NULL',
                 '_specT_30', sep='')
writeOGR(crowns_sp,
         dsn = paste(file_name, '.shp', sep=''),
         layer = basename(file_name),
         drive = 'ESRI Shapefile')










######################### Now doing evaluation!



### Load these itcfast crowns back in
Sep_itc_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/shape/sepilok_deploy/left_plot_itc.shp")
plot(Sep_itc_crowns)

### Load in our manual crowns!
Sep_manual_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/manual_crowns/eval_sepilok_plots/sepilok_603450_648000_manuals.shp")
plot(Sep_manual_crowns)
#Dan_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/manual_crowns/yujie_crowns_danum/Danum_2014_RGB_manual_crowns_Yujie.shp")


# Extract height info from CHMs (function is below)
dfs_itc=get_poly_CHM_info_velox(rpoly=Sep_itc_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)
dfs_manual=get_poly_CHM_info_velox(rpoly=Sep_manual_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)

# Slow method
dfs=get_poly_CHM_info_stars(rpoly=Sep_itc_crowns,CHM_old=Sep$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep$Sepilok_CHM_diff_0.5m)
dfd=get_poly_CHM_info_stars(rpoly=Dan_crowns,CHM_old=Dan$Danum_2014_CHM_g10_sub0.2_1m,CHM_diff=Dan$Danum_CHM_diff_1m)


# Combine data frames
dfs_itc$Site="Sepilok"
dfs_manual$Site="Sepilok"

#Get tall trees
dfs_itc_no_nans_H_max = dfs_itc[!is.na(dfs_itc$Org_H_max), ]
dfs_itc_tall = dfs_itc_no_nans_H_max[dfs_itc_no_nans_H_max$Org_H_max>60,] 
plot(dfs_itc_tall)

dfs_manual_no_nans_H_max = dfs_manual[!is.na(dfs_manual$Org_H_max), ]
dfs_manual_tall = dfs_manual_no_nans_H_max[dfs_manual_no_nans_H_max$Org_H_max>60,]
plot(dfs_manual_tall)

### Now let's do our naughty evaluation to yield the F1 score
### CROWN_OVERLAP function below

crs(dfs_itc_tall)<-crs(dfs_manual_tall)

matched<-crown_overlap(auto_trees=dfs_itc_tall, manual_trees=dfs_manual_tall, buffer_by=0)
#matched<-crown_overlap(auto_trees=ut_large, manual_trees=manual_large, buffer_by=0)

##### evaluation statistics

no_of_preds <- length(dfs_itc_tall)
no_of_manuals <- length(dfs_manual_tall)
matched<-matched[!is.na(matched$id_auto_tree),]
auto<-dfs_itc_tall[matched$id_auto_tree,]
tp = length(matched)
fp = no_of_preds - tp
fn = no_of_manuals - length(matched)
precision = tp/(tp+fp)
recall = tp/(tp + fn)
fscore = (2*precision*recall)/(precision+recall)

matched<-matched[!is.na(matched$id_auto_tree),]
auto<-ut[matched$id_auto_tree,]
auto@data$shpbnd<-rowSums(auto@data[,c('hbnd_mn', 'hbnd_mx', 'soblbnd')])
auto@data$spcbnd<-rowSums(auto@data[,c('allmbnd', 'crwnbnd')])
auto@data<-auto@data[,c('trHghts_mn', 'trHghts_mx', 'R', 'shpbnd', 'spcbnd', 'sobl_mn')]
matched@data<-cbind(matched@data, auto@data)

# par(mfrow=c(1,1), mar=c(0,0,0,0))
length(ut)
plot(dfs_manual_tall)
plot(dfs_itc_tall, add=TRUE, border='red')
plot(auto, add=TRUE, border='blue')




#################################
### Now looking at plots of dead trees, growth rate etc. for merged itc crowns

### Load up and crop the LiDAR to the whole left hand plot of Sepilok
Sep = raster::stack(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2014_coarse_CHM_g10_sub0.01_0.5m.tif"), paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sepilok_CHM_diff_0.5m.tif"))

Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/full_sepilok_deploy/left_plot.shp")
plot(Sep_detectron_crowns)
Sep_detectron_crowns$perimeter=as.numeric(polyPerimeter(Sep_detectron_crowns))
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$pixelvalue>100,]
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$perimeter<300,]
plot(Sep_detectron_crowns)

bound = st_bbox(Sep_detectron_crowns)
Sep_cropped = raster::crop(Sep, bound)
plot(Sep_cropped)

### Load these itcfast crowns back in
Sep_itc_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/shape/sepilok_deploy/left_plot_itc.shp")
plot(Sep_itc_crowns)

### Load in our manual crowns!
Sep_manual_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/manual_crowns/eval_sepilok_plots/sepilok_603450_648000_manuals.shp")
plot(Sep_manual_crowns)
#Dan_crowns=shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/manual_crowns/yujie_crowns_danum/Danum_2014_RGB_manual_crowns_Yujie.shp")

# Extract height info from CHMs (function is below)
dfs_detectron=get_poly_CHM_info_velox(rpoly=Sep_detectron_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)
dfs_itc=get_poly_CHM_info_velox(rpoly=Sep_itc_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)
dfs_manual=get_poly_CHM_info_velox(rpoly=Sep_manual_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)

# Slow method
#dfs=get_poly_CHM_info_stars(rpoly=Sep_itc_crowns,CHM_old=Sep$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep$Sepilok_CHM_diff_0.5m)
#dfd=get_poly_CHM_info_stars(rpoly=Dan_crowns,CHM_old=Dan$Danum_2014_CHM_g10_sub0.2_1m,CHM_diff=Dan$Danum_CHM_diff_1m)


# Combine data frames
dfs_detectron$Site="Sepilok"
dfs_itc$Site="Sepilok"
dfs_manual$Site="Sepilok"

# Plot the results
ggplot(dfs_detectron@data,aes(Org_H_max,Change_H_mean,size=area,color=Site))+geom_point()

# Summarize by 10 m bins
dfs_detectron$Org_H_max_round=round(dfs_detectron$Org_H_max,digits=-1)
dfs_detectron$Org_H_mean_round=round(dfs_detectron$Org_H_mean,digits=-1)

df_summary=dfs_detectron@data %>% group_by(Site,Org_H_mean_round) %>% 
  summarize(num_tot=n(),num_died=sum(Change_H_mean<=-5, na.rm=TRUE),pct_died=100*num_died/num_tot)

# Plot the results of summary
ggplot(df_summary,aes(Org_H_mean_round,num_died,color=Site))+geom_point()


### Now let's just select the tallest trees...

dfs_no_nans_H_max = dfs[!is.na(dfs$Org_H_max), ]
dfs_tall = dfs_no_nans_H_max[dfs_no_nans_H_max$Org_H_max>60,] 


### crown_overlap function
crown_overlap<-function(auto_trees, manual_trees, buffer_by, verbose='off'){
  # out<-matrix(0, nrow=nrow(manual_trees), ncol=4)
  out<-matrix(0, nrow=nrow(manual_trees), ncol=6)
  out[,1]<-1:nrow(manual_trees)
  out[,3]<-1
  
  auto_trees$id<-1:nrow(auto_trees)
  sum(sapply(auto_trees@polygons, function(x) x@area)<0.001)
  
  for(i in 1:nrow(manual_trees)){
    i_esc<<-i
    #    print(i)
    poly<-manual_trees[i,] # selects the current polygon
    poly_area<-poly@polygons[[1]]@Polygons[[1]]@area # the area of the polygon
    poly_buffer<-gBuffer(poly, width=buffer_by) # makes a buffer around it
    #cropped_trees<-raster::crop(auto_trees, poly_buffer) # crops the auto trees to the buffer
    cropped_trees<-raster::crop(auto_trees, poly) # crops the auto trees to any that intersect with poly
    cropped_trees<-auto_trees[auto_trees$id %in% cropped_trees$id,]
    cropped_trees <- gBuffer(cropped_trees, byid=TRUE, width=0) # deals with self-intersection
    if(!is.null(cropped_trees)){
      for (j in 1:nrow(cropped_trees)) {
        #      print(j)
        overlap <- gIntersection(poly, cropped_trees[j, ])# extracts the intersecting area.
        overlap_area <- overlap@polygons[[1]]@Polygons[[1]]@area # the area of the intersection
        union1 <- gUnion(poly, cropped_trees[j, ])
        union1_area <- union1@polygons[[1]]@Polygons[[1]]@area
        iou <- overlap_area/union1_area
        if (iou > 0.49) {
          auto_area <- cropped_trees@polygons[[j]]@Polygons[[1]]@area
          overlap_area <- overlap@polygons[[1]]@Polygons[[1]]@area # the area of the intersection
          
          # now need to try and work out the metrics - precision/recall/F1
          
          
          #overseg <- (auto_area - overlap_area) / poly_area
          overseg <- 1-(overlap_area / auto_area) # false positive rate
          underseg <- (poly_area - overlap_area) / poly_area # false negative rate
          # overlap_percent<-overlap_area/poly_area # the percentage area
          size_ratio<-auto_area/poly_area # the percentage area
          
          if (out[i, 3] == 1) {
            out[i, 2] <- cropped_trees$id[j]
            out[i, 3] <- 0
            out[i, 4] <- overseg
            out[i, 5] <- underseg
            out[i, 6] <- size_ratio
          }
          else if ((overseg + underseg) < sum(out[i, 4:5])) {
            # stores the result
            out[i, 2] <- cropped_trees$id[j]
            out[i, 4] <- overseg
            out[i, 5] <- underseg
            out[i, 6] <- size_ratio
          }
          # if(verbose=='on')
          #   cat('j: ', j, 'overlap: ', overlap_percent, '\n')
          # if(overlap_percent>out[i,3]){ # stores the result
          #   out[i,1]<-i
          #   out[i,2]<-cropped_trees$id[j]
          #   out[i,3]<-overlap_percent
          #   out[i,4]<-size_ratio
          # }
        }
      }
      if (verbose == 'on')
        cat('out: ', out[i, ], '\n')
    }
  }
  out[out[,2]==0,2]<-NA
  # Loads these as additional columns for the manual trees:
  manual_trees$id_auto_tree<-out[,2]
  # manual_trees$overlap_auto_tree <- out[, 3]
  # manual_trees$size_ratio <- out[, 4]
  manual_trees$tree_match <- out[, 3]
  manual_trees$overseg <- out[, 4]
  manual_trees$underseg <- out[, 5]
  manual_trees$cost <- rowSums(out[,4:5])/2 + out[,3]
  manual_trees$size_ratio<-out[,6]
  
  if (verbose == 'on')
    cat('Cost: ', manual_trees$cost, '\n')
  
  return(manual_trees)
}




# This is the fast way of extracting raster info for a list of polygons - requires velox
get_poly_CHM_info_velox=function(rpoly,CHM_org,CHM_diff){
  
  CHM_orgv=velox(CHM_org)
  CHM_diffv=velox(CHM_diff)
  
  # Get info from original raster
  # does this introduce NaNs?
  poly_original_raster_data_list= CHM_orgv$extract(rpoly, fun = NULL)
  rpoly$area=sapply(poly_original_raster_data_list,FUN=length) # poly area in raster units
  rpoly$Org_H_mean=sapply(poly_original_raster_data_list,FUN=mean,na.rm=TRUE)
  rpoly$Org_H_max =sapply(poly_original_raster_data_list,FUN=max)
  rpoly$Org_H_min =sapply(poly_original_raster_data_list,FUN=min)
  rpoly$Org_H_var =sapply(poly_original_raster_data_list,FUN=var)
  
  
  poly_change_raster_data_list=CHM_diffv$extract(rpoly, fun = NULL)
  rpoly$Change_H_mean=sapply(poly_change_raster_data_list,FUN=mean)
  rpoly$Change_H_max =sapply(poly_change_raster_data_list,FUN=max)
  rpoly$Change_H_min =sapply(poly_change_raster_data_list,FUN=min)
  rpoly$Change_H_var =sapply(poly_change_raster_data_list,FUN=var)
  
  rpoly$perimeter=as.numeric(polyPerimeter(rpoly)) # perimeter of each polygon
  rpoly$shape_complexity = as.numeric(rpoly$perimeter/(2*sqrt(rpoly$area*pi)))
  rpoly$shape_circleness=as.numeric(4*pi*(rpoly$area)/((rpoly$perimeter)^2))
  return(rpoly)
}


rpoly$perimeter=as.numeric(polyPerimeter(rpoly))

##### SLOW OLD WAY TO EXTRACT INFO WITHOUT VELOX
get_poly_CHM_info_stars=function(rpoly,CHM_old,CHM_diff){
  
  # Get info from original raster
  #https://gis.stackexchange.com/questions/130522/increasing-speed-of-crop-mask-extract-raster-by-many-polygons-in-r
  poly_original_raster_data_list=list()
  for (i in 1:nrow(rpoly)) { #this is the number of polygons to iterate through
    #foreach (i = 1:nrow(rpoly),.packages= c("raster","foreach")) %do% {
    single <- rpoly[i,] #selects a single polygon
    clip1 <- raster::crop(CHM_old, extent(single)) #crops the raster to the extent of the polygon, I do this first because it speeds the mask up
    clip2 <- raster::rasterize(single, clip1, mask=TRUE) #crops to polygon edge & converts to raster
    poly_original_raster_data_list[[i]] <- raster::getValues(clip2) #much faster than extract
    #return(poly_original_raster_data_list)
    #ext[[i]]<-extract(clip2,single) #extracts data from the raster based on the polygon bound
  }
  
  #poly_original_raster_data_list=raster::extract(CHM_old,rpoly,fun=NULL)
  #test=CHM2016v$extract(rpoly,fun=mean)
  #test2=raster::extract(CHM2016,rpoly,fun=mean)
  
  #test=raster::extract(CHM_old,rpoly,fun=max)
  rpoly$area=sapply(poly_original_raster_data_list,FUN=length) # poly area in raster units
  rpoly$Org_H_mean=sapply(poly_original_raster_data_list,FUN=mean)
  rpoly$Org_H_max =sapply(poly_original_raster_data_list,FUN=max)
  rpoly$Org_H_min =sapply(poly_original_raster_data_list,FUN=min)
  rpoly$Org_H_var =sapply(poly_original_raster_data_list,FUN=var)
  
  # Get info from original raster
  poly_change_raster_data_list=list()
  #mypb <- tkProgressBar(title = "R progress bar", label = "", min = 0, max = nrow(rpoly), initial = 0, width = 300) 
  for (i in 1:nrow(rpoly)) { #this is the number of polygons to iterate through
    #foreach (i = 1:nrow(rpoly),.packages= c("raster","foreach")) %do% { 
    single <- rpoly[i,] #selects a single polygon
    clip1 <- raster::crop(CHM_diff, extent(single)) #crops the raster to the extent of the polygon, I do this first because it speeds the mask up
    clip2 <- raster::rasterize(single, clip1, mask=TRUE) #crops to polygon edge & converts to raster
    poly_change_raster_data_list[[i]] <- raster::getValues(clip2) #much faster than extract
    #setTkProgressBar(mypb, i, title = "number complete", label = i)
    #ext[[i]]<-extract(clip2,single) #extracts data from the raster based on the polygon bound
  }
  
  #poly_change_raster_data_list=raster::extract(CHM_diff,rpoly,fun=NULL)
  rpoly$Change_H_mean=sapply(poly_change_raster_data_list,FUN=mean)
  rpoly$Change_H_max =sapply(poly_change_raster_data_list,FUN=max)
  rpoly$Change_H_min =sapply(poly_change_raster_data_list,FUN=min)
  rpoly$Change_H_var =sapply(poly_change_raster_data_list,FUN=var)
  
  rpoly$shape_complexity = as.numeric(rpoly$perimeter/(2*sqrt(rpoly$area*pi)))
  rpoly$shape_circleness=as.numeric(4*pi*(rpoly$area)/((rpoly$perimeter)^2))
  
  return(rpoly)
}
