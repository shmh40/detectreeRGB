# Load up some packages!

rm(list=ls())

install.packages('bbplot')

library(bbplot)
library(reshape2)
library(raster)
library(dplyr); library(tidyr); library(magrittr)
library(ggplot2); library(RColorBrewer); library(rgdal)

# LiDAR packages
library(raster); library(sf); library(spatialEco)
library(rasterVis) # This is just for plotting, you probably don't need it

library(velox)

Sep = raster::stack(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2014_coarse_CHM_g10_sub0.01_0.5m.tif"), paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2020_CHM_g10_sub0.2_0.5m.tif"), paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sepilok_CHM_diff_0.5m.tif"))
Sep_2014 = raster(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2014_coarse_CHM_g10_sub0.01_0.5m.tif"))
Sep_2020 = raster(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sep_2020_CHM_g10_sub0.2_0.5m.tif"))
Sep_diff = raster(paste0("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sepilok_CHM_diff_0.5m.tif"))


DTM = raster("C:/Users/sebhi/ai4er/mres_project/work/data/sepilok_lidar/Sepilok_2014_DTM_g10_sub0.2_1m.tif")
plot(DTM)

DTM_agg = raster::aggregate(DTM, 20)
plot(DTM_agg)

slope_agg   <- terrain(DTM_agg, opt=c('slope'), unit='degrees',na.rm=TRUE)
aspect_agg  <- terrain(DTM_agg, opt=c('aspect'), na.rm=TRUE)
TPI_agg     <- terrain(DTM_agg, opt=c('TPI'), na.rm=TRUE)
plot(TPI_agg)

## Load in Mask RCNN crowns and remove the grid boxes and ones off the edge of the RGB

Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/full_sepilok_deploy/left_plot.shp")
Sep_detectron_crowns$perimeter=as.numeric(polyPerimeter(Sep_detectron_crowns))
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$pixelvalue>100,] 
Sep_detectron_crowns = Sep_detectron_crowns[Sep_detectron_crowns$perimeter<220,]
#file_name<-paste("C:/Users/sebhi/ai4er/mres_project/work/data/detectron_predicted_crowns/full_sepilok_deploy/maskrcnn_preds_sepilok")
#writeOGR(Sep_detectron_crowns,dsn = paste(file_name, '.shp', sep=''),layer = basename(file_name),drive = 'ESRI Shapefile')

plot(Sep_detectron_crowns)

### Load in itcfast predictions

Sep_detectron_crowns = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/shape/sepilok_deploy/left_plot_itc.shp")

Sep_plot_outline = shapefile("C:/Users/sebhi/ai4er/mres_project/work/data/shape/sepilok_2014_preds_outline.shp")

### CROP the LiDAR data
#bound = st_bbox(Sep_detectron_crowns)

Sep_2014_cropped = raster::crop(Sep_2014, Sep_plot_outline)
Sep_2020_cropped = raster::crop(Sep_2020, Sep_plot_outline) 
Sep_diff_cropped = raster::crop(Sep_diff, Sep_plot_outline)

Sep_2014_masked = raster::mask(Sep_2014_cropped, Sep_plot_outline)
Sep_2020_masked = raster::mask(Sep_2020_cropped, Sep_plot_outline) 
Sep_diff_masked = raster::mask(Sep_diff_cropped, Sep_plot_outline)

plot(Sep_2014_masked)
plot(Sep_2020_masked)
plot(Sep_diff_masked)

# stack them up...currently failing
Sep = raster::stack(Sep_2014_masked, Sep_2020_masked, Sep_diff_masked)

## CROP and MASK

### CROP the LiDAR data
TPI_cropped = raster::crop(TPI_agg, Sep_plot_outline)
TPI_masked = raster::mask(TPI_cropped, Sep_plot_outline)
plot(TPI_masked)
#plot(dfs, add=TRUE)

slope_cropped = raster::crop(slope_agg, Sep_plot_outline)
slope_masked = raster::mask(slope_cropped, Sep_plot_outline)
plot(slope_masked)
#plot(dfs, add=TRUE)

aspect_cropped = raster::crop(aspect_agg, Sep_plot_outline)
aspect_masked = raster::mask(aspect_cropped, Sep_plot_outline)
plot(aspect_masked)

#dfs=get_poly_CHM_info_velox(rpoly=Sep_detectron_crowns,CHM_org=Sep_cropped$Sep_2014_coarse_CHM_g10_sub0.01_0.5m,CHM_diff=Sep_cropped$Sepilok_CHM_diff_0.5m)
dfs=get_poly_CHM_info_velox(rpoly=Sep_detectron_crowns,CHM_org=Sep_2014_masked,CHM_2020=Sep_2020_masked,CHM_diff=Sep_diff_masked)

# Give the site of the dfs
dfs$Site="Sepilok"

### COUNT THE NUMBER OF TREES ABOVE CERTAN HEIGHT THRESHOLDS
dfs_2020 = dfs@data[dfs@data$H_max_2020>70,]
dfs_2020_drop_na = dfs_2020[!is.na(dfs_2020$H_max_2020), ]

dfs_30_drop_na$H_max_2020 = dfs_30_drop_na$Org_H_max + dfs_30_drop_na$Change_H_max
dfs_30_drop_na = dfs_30_drop_na[!is.na(dfs_30_drop_na$H_max_2020), ]
dfs_H = dfs_30_drop_na[dfs_30_drop_na$H_max_2020>100,]

### now doing analysis of not dead trees
ggplot(dfs@data,aes(Org_H_max,Change_H_mean,size=area,color=Site))+geom_point()

ggplot(dfs@data,aes(Org_H_max, Change_H_mean))+
  geom_point(color='red',aes(Org_H_max,Change_H_mean, size=area))+
  ggtitle('Change in tree height predicted by Mask R-CNN')+
  theme_classic()+
  theme(plot.title=element_text(size=10,face='bold',margin=margin(10,10,10,0), hjust = 0.5))+
  labs(x='2014 tree height / m', y='Change in tree height / m')+
  theme(legend.title = element_text(colour="black", size=8, face="bold"))+
  labs(size='Crown Area')+
  guides(size = guide_legend(reverse=TRUE))+
  xlim(30,80)+
  ylim(-60, 20)

dfs_drop_na = dfs[!is.na(dfs$Change_H_mean), ]
dfs_drop_dead = dfs_drop_na[dfs_drop_na$Change_H_mean>=-5,]

ggplot(dfs_drop_dead@data,aes(Org_H_max,Change_H_mean,size=area,color=Site))+geom_point()

ggplot(dfs_drop_dead@data,aes(Org_H_max, Change_H_mean))+
  geom_point(color='red',aes(Org_H_max,Change_H_mean, size=area))+
  ggtitle('Mask R-CNN predicts decrease in tree growth with height')+
  theme_classic()+
  theme(plot.title=element_text(size=10,face='bold',margin=margin(10,10,10,0), hjust = 0.5))+
  labs(x='2014 tree height / m', y='Change in tree height / m')+
  theme(legend.title = element_text(colour="black", size=8, face="bold"))+
  labs(size='Crown Area')+
  guides(size = guide_legend(reverse=TRUE))+
  geom_abline(intercept = model_intercept, slope = model_slope, color = 'black')+
  xlim(30,80)+
  ylim(-5,12)



model_linear1 <- lm(Change_H_mean ~ Org_H_max, data = dfs_drop_dead@data)
summary(model_linear1)

model_intercept <- coef(model_linear1)[1]
model_slope <- coef(model_linear1)[2]

ggplot(data = dfs_drop_dead@data, aes(Org_H_max, Change_H_mean, size=area, color=Site)) +
  geom_point() +
  geom_abline(intercept = model_intercept, slope = model_slope, color = 'black')+xlim(30,80)



### analysis of percentage deaths in each height bin

dfs$Org_H_max_round=round(dfs$Org_H_max,digits=-1)
dfs$Org_H_mean_round=round(dfs$Org_H_mean,digits=-1)

df_summary=dfs@data %>% group_by(Site,Org_H_mean_round) %>% 
  summarize(num_tot=n(),num_died=sum(Change_H_mean<=-5, na.rm=TRUE),pct_died_per_year=(100*num_died/num_tot)^(1/6))

# Plot the results of summary
ggplot(df_summary,aes(Org_H_mean_round,pct_died_per_year,color=Site))+geom_point()+xlim(20,60)

ggplot(df_summary,aes(Org_H_mean_round, pct_died_per_year))+
  geom_point(color='red',aes(Org_H_mean_round,pct_died_per_year))+
  ggtitle('Annual mortality rate of Mask R-CNN trees')+
  theme_classic()+
  theme(plot.title=element_text(size=10,face='bold',margin=margin(10,10,10,0), hjust = 0.5))+
  labs(x='Binned tree heights / m', y='Annual tree mortality / %')+
  theme(legend.title = element_text(colour="black", size=12, face="bold"))+
  xlim(20,60)+
  ylim(1, 2)

### Now do a carbon calculation

dfs_50 = dfs@data[dfs@data$Org_H_max>0,]
dfs_50_drop_na = dfs_50[!is.na(dfs_50$Org_H_max), ]

dfs_50_drop_na$diameter = 2*sqrt((dfs_50_drop_na$area)/pi)

dfs_50_drop_na$agb = 0.136 * (dfs_50_drop_na$Org_H_max * dfs_50_drop_na$diameter)^1.52

dfs_50_drop_na$carbon = 0.5 * dfs_50_drop_na$agb

sum_carbon = sum(dfs_50_drop_na$carbon, na.rm=TRUE)

dfs_50_drop_na$agb_change = (0.136 * ((dfs_50_drop_na$Org_H_max + dfs_50_drop_na$Change_H_max) * dfs_50_drop_na$diameter)^1.52) - (0.136 * (dfs$Org_H_max * dfs$diameter)^1.52)

dfs_50_drop_na$agb_change = (0.136 * ((dfs_50$Org_H_mean + dfs_50$Change_H_mean) * dfs_50$diameter)^1.52) - (0.136 * (dfs_50$Org_H_mean * dfs_50$diameter)^1.52)

dfs_50$carbon_change = 0.5 * dfs_50$agb_change 

sum_carbon_change = sum(dfs_50$carbon_change, na.rm=TRUE)

ggplot(dfs_50_drop_na,aes(Org_H_max, carbon))+
  geom_point(color='red',aes(Org_H_max,carbon, size=area))+
  ylim(0, 40000)+ggtitle('Carbon stored in each tree predicted by Mask R-CNN')+
  theme_classic()+
  theme(plot.title=element_text(size=10,face='bold',margin=margin(10,10,10,0), hjust = 0.5))+
  labs(x='Tree height / m', y='Carbon / kg')+
  theme(legend.title = element_text(colour="black", size=8, face="bold"))+
  labs(size='Crown Area')+
  guides(size = guide_legend(reverse=TRUE))


### Locate the DEAD trees
#dfs = dfs@data[dfs@data$Org_H_max>0,]
dfs_drop_na = dfs[!is.na(dfs$Change_H_mean), ]
dfs_dead = dfs_drop_na[dfs_drop_na$Change_H_mean<=-5,]


plot(TPI_masked)
plot(dfs_dead, add=TRUE)

hist(TPI_masked,
     main = "Distribution of TPI in Sepilok",
     xlab = "TPI", ylab = "Frequency",
     col = "springgreen")

hist(slope_masked,
     main = "Distribution of Slope in Sepilok",
     xlab = "Slope", ylab = "Frequency",
     col = "springgreen")

hist(aspect_masked,
     main = "Distribution of Aspect in Sepilok",
     xlab = "Aspect", ylab = "Frequency",
     col = "springgreen")

dfs_dead_velox_info=get_poly_TPI_slope_aspect_info_velox(rpoly=dfs_dead,TPI=TPI_masked,slope=slope_masked,aspect=aspect_masked)
plot(dfs_dead_velox_info)

hist(dfs_dead_velox_info$TPI,
     main = "Distribution of TPI in Sepilok",
     xlab = "TPI", ylab = "Frequency",
     col = "black")

hist(dfs_dead_velox_info$slope,
     main = "Distribution of TPI in Sepilok",
     xlab = "Slope", ylab = "Frequency",
     col = "black")

hist(dfs_dead_velox_info$aspect,
     main = "Distribution of TPI in Sepilok",
     xlab = "Slope", ylab = "Frequency",
     col = "black")


### Now trying to do histograms with ggplot2

dfs_dead_velox_df = dfs_dead_velox_info@data
dfs_dead_velox_df_rm_na = dfs_dead_velox_df[!is.na(dfs_dead_velox_df$TPI), ]
#dfs_dead_velox_df[['TPI','slope','aspect']]
#dfs_dead_velox_df.m = melt(dfs_dead_velox_df)

TPI_df = as.data.frame(TPI_masked)
slope_df = as.data.frame(slope_masked)
aspect_df = as.data.frame(aspect_masked)

ggplot(TPI_df, aes(x=tpi))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('TPI Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))

ggplot(dfs_dead_velox_df, aes(x=TPI))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('TPI Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))

ggplot(slope_df, aes(x=slope))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('Slope Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))

ggplot(dfs_dead_velox_df, aes(x = slope))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('Slope Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))

ggplot(aspect_df, aes(x=aspect))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('Aspect Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))

ggplot(dfs_dead_velox_df, aes(x = aspect))+
  geom_vline(xintercept=0, color="black", size=1)+
  stat_density(alpha = 1,size = 1, geom="line", position = "identity")+
  xlab('Aspect Distribution')+ylab("")+xlim(-1.5,1.5)+theme_light()+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.title = element_blank(),legend.position = c(0.2,0.8))+
  guides(color=guide_legend(nrow=4))+
  scale_color_manual(values=c(brewer.pal(6,"Paired")[c(1,2,5,6)]))


###


get_poly_TPI_slope_aspect_info_velox=function(rpoly,TPI,slope,aspect){
  
  TPI_orgv=velox(TPI)
  slope_orgv=velox(slope)
  aspect_orgv=velox(aspect)
  
  # Get info from the rasters
  poly_TPI_raster_data_list= TPI_orgv$extract(rpoly, fun = NULL)
  rpoly$TPI=sapply(poly_TPI_raster_data_list,FUN=mean,na.rm=TRUE)
  
  poly_slope_raster_data_list= slope_orgv$extract(rpoly, fun = NULL)
  rpoly$slope=sapply(poly_slope_raster_data_list,FUN=mean,na.rm=TRUE) # poly area in raster units
  
  poly_aspect_raster_data_list= aspect_orgv$extract(rpoly, fun = NULL)
  rpoly$aspect=sapply(poly_aspect_raster_data_list,FUN=mean,na.rm=TRUE) # poly area in raster units
  
  return(rpoly)
}


get_poly_CHM_info_velox=function(rpoly,CHM_org,CHM_2020,CHM_diff){
  
  CHM_orgv=velox(CHM_org)
  CHM_2020v=velox(CHM_2020)
  CHM_diffv=velox(CHM_diff)
  
  # Get info from original raster
  # does this introduce NaNs?
  poly_original_raster_data_list= CHM_orgv$extract(rpoly, fun = NULL)
  rpoly$area=sapply(poly_original_raster_data_list,FUN=length) # poly area in raster units
  rpoly$Org_H_mean=sapply(poly_original_raster_data_list,FUN=mean,na.rm=TRUE)
  rpoly$Org_H_max =sapply(poly_original_raster_data_list,FUN=max)
  rpoly$Org_H_min =sapply(poly_original_raster_data_list,FUN=min)
  rpoly$Org_H_var =sapply(poly_original_raster_data_list,FUN=var)
  
  poly_2020_raster_data_list= CHM_2020v$extract(rpoly, fun = NULL)
  rpoly$area=sapply(poly_2020_raster_data_list,FUN=length) # poly area in raster units
  rpoly$H_mean_2020=sapply(poly_2020_raster_data_list,FUN=mean,na.rm=TRUE)
  rpoly$H_max_2020 =sapply(poly_2020_raster_data_list,FUN=max)
  rpoly$H_min_2020 =sapply(poly_2020_raster_data_list,FUN=min)
  rpoly$H_var_2020 =sapply(poly_2020_raster_data_list,FUN=var)
  
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


## old velox with the diff

get_poly_CHM_info_velox_old=function(rpoly,CHM_org,CHM_diff){
  
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

