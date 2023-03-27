#### Aerial coverage of phenology
library(tidyverse)
library(terra)
library(sf)
library(brms)


myEPSG = 3488
#### Import data ----
## NAIP image available on earth engine
## 3DEP resampled to 30m; download from earth engine
## Contact CDFW for ROI file
NAIP_img = rast("data/NAIP/NAIP_30m.tif")
DEM = rast("data/DEP_30m/DEM_30m.tif")
ASP = terra::terrain(DEM,v="aspect", neighbors=8, unit="degrees")
ROI = st_read("data/PRIVATE/HerdUnits_UpdatedMetadata.shp") %>%
  st_transform(st_crs(DEM))
ROIvect = vect(ROI)

HUbuff = ROI %>%
  st_union() %>%
  st_buffer(5000) %>%
  st_transform(st_crs(DEM)) 
HUbuffvect = vect(HUbuff)

DEM_cropped = DEM  %>%
  crop(HUbuffvect) %>%
  mask(HUbuffvect) %>%
  project(paste0("epsg:",myEPSG))
plot(DEM_cropped)

ASP_deg_cropped = ASP %>%
  crop(HUbuffvect) %>%
  mask(HUbuffvect) %>%
  project(DEM_cropped)
plot(ASP_deg_cropped)

ASP_cropped = app(ASP_deg_cropped,
                  function(x){
                    cos(((x*pi)/180))
                  })
plot(ASP_cropped)

ROI_alb = st_transform(ROI,
                       myEPSG) %>%
  mutate(ID = 1:nrow(.))
ROI_alb_vect = vect(ROI_alb)

HU_alb = st_transform(HUbuff,
                      myEPSG)
HU_alb_vect = vect(HU_alb)

## Stan model
terrainModel = readRDS("stan_models/terrainModel.RDS")
coefYrMod = coefficients(terrainModel)

coefInt = mean(coefYrMod$Year[,,"Intercept"][,1])
coefElv = coefYrMod$Year[,,"Elevation"][1,1]
coefAsp = coefYrMod$Year[,,"cosAsp"][1,1]


#### Data wrangling ----
elevation = extract(DEM_cropped, ROI_alb_vect)
elevation$HerdUnit = ROI_alb$NAME[match(elevation$ID,
                                        as.numeric(ROI_alb$ID))]
elevation
aspect = extract(ASP_cropped, ROI_alb_vect)
aspect$HerdUnit = ROI_alb$NAME[match(aspect$ID,
                                     as.numeric(ROI_alb$ID))]
aspect


## Generate phenology plot:
myDOYRast = DEM_cropped*coefElv +
  ASP_cropped*coefAsp +
  coefInt
names(myDOYRast) <- "DOY"

# Take a day of year, calculate the temporal distance from peak green-up timing
# Identify pixels within 5 days of peak green-up timing
plotDistFromDOY = function(myDOY = 150,
                           temporalDistance = 5,
                           doyRast = myDOYRast,
                           makePlot = TRUE){
  closePhenology = doyRast - myDOY
  closePhenology[closePhenology > temporalDistance] <- NA
  closePhenology[closePhenology < -1*temporalDistance] <- NA
  closePhenology[!is.na(closePhenology)] <- 1
  # plot(closePhenology)
  
  if(makePlot){
    doyRastDF = doyRast %>%
      crop(HU_alb_vect) %>%
      mask(HU_alb_vect) %>%
      as.data.frame(xy=TRUE) 
    
    closePlot = doyRastDF %>%
      filter(!is.na(DOY)) %>%
      mutate(diffDOY = DOY - myDOY) %>%
      ggplot() +
      geom_raster(aes(x = x,
                      y = y,
                      fill = diffDOY)) +
      geom_sf(data = HU_alb,
              fill = "transparent") +
      geom_sf(data = ROI_alb,
              fill = "transparent") +
      scale_fill_distiller("Time from\npeak green-up",
                           type = "div",
                           lim = c(-100,100),
                           direction = 1) +
      ggtitle(paste0("Day ", myDOY)) +
      theme_void()
    return(list(closePhenology,closePlot))
    
  }else{
    return(list(closePhenology,NULL))
  }
}
plotInsideDOY = function(myDOY = 170,
                         temporalDistance = 5,
                         doyRast = myDOYRast,
                         makePlot = TRUE){
  closePhenology = doyRast - myDOY
  closePhenology[closePhenology > temporalDistance] <- NA
  closePhenology[closePhenology < -1*temporalDistance] <- NA
  closePhenology[!is.na(closePhenology)] <- 1
  # plot(closePhenology)
  
  doyRastDF = closePhenology %>%
    # aggregate(16) %>%
    crop(HU_alb_vect) %>%
    mask(HU_alb_vect) %>%
    as.data.frame(xy=TRUE) 
  
  closePlot = doyRastDF %>%
    filter(!is.na(DOY)) %>%
    mutate(DOY = as.factor(DOY)) %>%
    ggplot() +
    geom_raster(aes(x = x,
                    y = y,
                    fill = DOY)) +
    geom_sf(data = HU_alb,
            fill = "transparent") +
    geom_sf(data = ROI_alb,
            fill = "transparent") +
    scale_fill_manual("Peak green-up",
                      values = c("orange4")) +
    ggtitle(paste0("Day ", myDOY)) +
    theme_void()
  return(list(closePhenology,closePlot))
}

# Take a temporal distance from peak green-up raster and return
# a data.frame with herd unit and 1 or NA for within window or not
summarizeByPheno = function(distRast = day115[[1]]){
  within5 = extract(distRast,ROI_alb_vect)
  within5$HerdUnit = ROI_alb$NAME[match(within5$ID,
                                        as.numeric(ROI_alb$ID))]
  return(within5)  
}
# Apply across time period of interest (takes ~12 sec per frame)
distFromDOY_extracts = lapply(X = c(50:250),
                              FUN = function(x){
                                print(x)
                                print(Sys.time())
                                myDistRast = plotDistFromDOY(x,
                                                             makePlot = FALSE)
                                summaryDF = summarizeByPheno(myDistRast[[1]])
                                summaryDF$targetDOY = x
                                summaryDF$elevation = elevation$CA_DEM_UTM
                                return(summaryDF)
                              })
# Investigate proximity to green-up by elevation 
# over the course of the growing season
head(distFromDOY_extracts[[1]])

## Summarize by herd unit
distFromDOY_summaries = lapply(
  distFromDOY_extracts,
  function(x){
    x %>%
      group_by(HerdUnit, targetDOY) %>%
      summarize(
        available = n(),
        count = sum(DOY,
                    na.rm=T),
        area = count*res(DEM)[1]*res(DEM)[2])
  })

# Visualize phenology be herd unit across time
HU_propGreenupProximityDF = do.call("rbind",
                                    distFromDOY_summaries) %>%
  mutate(HerdUnit = case_when(HerdUnit == "Mt.Baxter" ~
                                "Mt. Baxter",
                              TRUE ~ HerdUnit),
         HUfact = factor(HerdUnit,
                         levels = c("Mt. Warren",
                                    "Mt. Gibbs",
                                    "Cathedral Range",
                                    "Convict Creek",
                                    "Wheeler Ridge",
                                    "Taboose Creek",
                                    "Sawmill Canyon",
                                    "Mt. Baxter",
                                    "Bubbs Creek",
                                    "Mt. Williamson",
                                    "Big Arroyo",
                                    "Mt. Langley",
                                    "Laurel Creek",
                                    "Olancha Peak"))) %>%
  filter(!is.na(HUfact))
HU_propGreenupProximity = HU_propGreenupProximityDF %>%
  mutate("DateChr" = paste0("2003-",targetDOY),
         "Date" = as.Date(DateChr, format = "%Y-%j")) %>%
  ggplot() +
  geom_raster(aes(x = Date,
                  y = HUfact,
                  fill = count/available)) +
  scale_x_date(date_labels = "%b") +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis_c(
    "Proportional area\nwithin 5 days of\npeak green-up") +
  xlab("Date") +
  ylab("Herd unit") +
  CJsBasics::BasicTheme

ggsave("plots/greenupProximity_herdunit.jpg",
       HU_propGreenupProximity,
       width = 6, height = 6, units = "in", dpi = 300)

#### Plot ----
#### Perspective plotting
library(rayshader)
aggFactor = 2

NAIP_img = NAIP_img %>%
  resample(DEM_cropped) %>%
  mask(HU_alb_vect)
NAIP_str = stretch(NAIP_img) %>%
  aggregate(aggFactor)
# plotRGB(NAIP_str)
NAIP_mat = as.array(NAIP_str)
NAIP_mat = NAIP_mat/255
# hist(NAIP_mat[,,1])
# hist(NAIP_mat[,,2])
# hist(NAIP_mat[,,3])

DEM_clipped = mask(DEM_cropped,
                   HU_alb_vect)
DEM_agg = DEM_clipped %>%
  aggregate(aggFactor)
makeRayshade = function(x,
                        DOY,
                        NPm){
  insideRast = x[[1]] %>%
    aggregate(aggFactor)
  
  DEM_agg_mat = t(as.array(DEM_agg)[,,1])
  ins_agg_mat = as.array(insideRast)[,,1]
  # image(ins_agg_mat)
  ins_agg_phot = NPm
  # ins_agg_phot[,,][is.na(ins_agg_mat)] <- 1
  ins_agg_phot[,,1][ins_agg_mat==1] <- 238/255
  ins_agg_phot[,,2][ins_agg_mat==1] <- 75/255
  ins_agg_phot[,,3][ins_agg_mat==1] <- 43/255
  
  # str(DEM_agg_mat)
  # str(ins_agg_mat)
  # str(ins_agg_phot)
  zscale = 40
  extDEM = ext(DEM_clipped)
  r = raster::raster(xmn = extDEM[1],
                     xmx = extDEM[2],
                     ymn=extDEM[3],
                     ymx=extDEM[4])
  re = raster::extent(r)
  DEM_agg_mat %>% 
    sphere_shade(texture = "imhof4") %>% 
    add_overlay(ins_agg_phot) %>%
    plot_3d(DEM_agg_mat, zscale = zscale, 
            theta = 55, zoom = 0.25, phi = 15,
            solidcolor = "grey40",
            windowsize = c(3200,2400))
  render_snapshot(paste0("plots/RegionalPhenologyWindow_55_",
                         DOY,
                         ".png"),
                  software_render = TRUE,
                  width = 3200,
                  height = 2400)
}

# STOP: NEED TO RUN THESE 1 AT A TIME
day1 = plotInsideDOY(130)
day2 = plotInsideDOY(150)
day3 = plotInsideDOY(170)
day4 = plotInsideDOY(190)

makeRayshade(day1, 130, NAIP_mat)
makeRayshade(day2, 150, NAIP_mat)
makeRayshade(day3, 170, NAIP_mat)
makeRayshade(day4, 190, NAIP_mat)
