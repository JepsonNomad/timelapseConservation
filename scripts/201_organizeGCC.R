#### Import packages ----
library(tidyverse)
library(sf)
source("scripts/001_cameraFunctions.R")


set.seed(1917) # Grinnell CA Thrasher

#### Define parameters ----
minYr = 2019
maxYr = 2020
plotDir = "plots"
dataDir = "data/stanPrep"

#### Create output directories ----
if(!dir.exists(plotDir)){
  dir.create(plotDir)
}
if(!dir.exists(dataDir)){
  dir.create(dataDir)
}


#### Import data ----
## Greenness data
fileListExtract = list.files("data/extract","*.csv", full.names = TRUE)
length(fileListExtract)
# 132 cameras
GCC = lapply(fileListExtract,
             read_csv) %>%
  do.call("rbind",.)
length(unique(GCC$filename))
# 69748 pictures

## Exif data
exifListExtract = list.files("data/exif","*.csv", full.names = TRUE)
exifs = lapply(exifListExtract,
               read_csv) %>%
  do.call("rbind",.)

exifsDF = exifs %>%
  mutate(Position = paste0(site,"-",plot),
         Date = lubridate::ymd_hms(dt),
         Year = lubridate::year(Date),
         Hour = lubridate::hour(Date)) %>%
  filter(Hour == 12) %>%
  select(Position, Date, Year, exposureTime)

## A function to find the most common exposure settings within all years for a cam position
exposureFinder = function(pos, depth = 1){
  # Narrow down to position
  posDF = exifsDF %>%
    filter(Position == pos) %>%
    filter(lubridate::yday(Date) > 100 & lubridate::yday(Date) <200)
  val = sort(table(posDF$exposureTime),decreasing=TRUE)[1:depth]
  return(names(val))
}

## Snow data
filePathSnow = "data/snowCover.csv"
snowWide = read_csv(filePathSnow, skip = 1)

## Alignment metadata and camera-level covariates
metadata = read_csv("data/metadata.csv", skip = 1)
camTopo = read_csv("data/camTopo.csv") # Topography etc
camLat = read_csv("data/camLat.csv") # Latitude
camNLCD = read_csv("data/camNLCD.csv") # Land cover

#### Data wrangling ----
## Compile metadata
positionData = camTopo %>%
  left_join(camLat,
            by = c("Position")) %>%
  left_join(camNLCD,
            by = c("Position"))

## Extract info from filenames in cam data
GCC = GCC %>%
  separate(filename, sep = "_",
           into = c("Position","Year","Month","DayOfMonth","TimeEtc"),
           remove = FALSE) %>%
  mutate(date_Date = as.Date(paste0(Year,"-",Month,"-",DayOfMonth)))

## Convert raw snow data into something useful
snowLong = snowWide %>%
  mutate("SnowEvent" = c(1:nrow(.)),
         SnowStart = as.Date(paste0(SnowStartYear,"-",SnowStartMonth,"-",SnowStartDay)),
         SnowStop = as.Date(paste0(SnowStopYear,"-",SnowStopMonth,"-",SnowStopDay)),
         snowDuration = SnowStop - SnowStart)  %>%
  group_by(Position, SnowEvent) %>%
  do(data.frame(SnowEvent = .$SnowEvent,
                date_Date = seq.Date(from = .$SnowStart,
                                     to = .$SnowStop,
                                     by = "1 day"))) %>%
  mutate("Snow" = 1)

## Identify camera positions to retain:
keeperCams = metadata %>%
  filter(quality > 2) %>%
  pull(Position)
length(keeperCams)
# 118 cameras

## Compile greenness and snow data into single data.frame. 
## Retain only cameras in keeperCams (i.e. quality rating was 3 or better)
## BA-01 gets bad lighting effects between mid-October and early March
## BA-02 region 4 gets overgrown by other plant spp.
## BA-06 had camera door open starting Sept 2020.
## CN-05 had a catastrophic frameshift during winter 2020-2021. Remove all GCC from Dec 2020 onward
## LN-07 fell between 16 Apr 2020 and 23 June 2020
## LN-10 fell between 25 Oct 2019 and 25 Apr 2020
## LN-11 gets bad lighting effects during the month of October
## LN-14 is offset by a month but I'm unsure which direction and can't fix that now
## OL-05 fell during Aug 2019. Remove those dates.
## OL-06 has branches that seasonally get in the way of regions 1, 3, and 4.
## Note that in WA-01, year was offset by 1 prior to Aug 2020. Fix here. Note that the logical test of 
# date_Date in case_when must test the final column to modify so that R doesn't give up on changing 
# things after changing date_Date.
## WA-13 fell from 2019-08-20 through 2020-06-11
## WH-01 gets bad lighting effects between early November and early February
## WH-06 and WH-07 fallen data were removed pre-alignment
## WL-10 regions 2, 3, 4, and 5 get seasonally sparse and visible rock messes up illumination
camdata = GCC %>%
  left_join(snowLong,
            by = c("Position", "date_Date")) %>%
  left_join(exifs,
            by = "filename") %>%
  mutate(doy = as.numeric(format(date_Date,"%j"))) %>%
  filter(Position %in% unique(snowLong$Position)) %>%
  filter(Position %in% keeperCams) %>%
  filter(!(Position == "BA-01" & 
             (doy < 75 |
                doy > 274))) %>%
  filter(!(Position == "BA-02" & regionID == 4)) %>%
  filter(!(Position == "BA-06" & date_Date > as.Date("2020-09-01"))) %>%
  mutate(species = ifelse(Position == "CN-01" & regionID == 4,
                          "Ribes sp",
                          species)) %>%
  filter(!(Position == "CN-05" & date_Date > as.Date("2020-12-01"))) %>%
  filter(!(Position == "LN-07" & 
             !(date_Date > as.Date("2020-06-23") | date_Date < as.Date("2020-04-16")))) %>%
  filter(!(Position == "LN-10" &
             !(date_Date > as.Date("2020-04-25") | date_Date < as.Date("2019-10-25")))) %>%
  filter(!(Position == "LN-11" & 
             (doy >= 274 & doy <= 304))) %>%
  filter(Position != "LN-14") %>%
  filter(!(Position == "OL-05" & 
             !(date_Date > as.Date("2019-08-25") | date_Date < as.Date("2019-08-09")))) %>%
  filter(!(Position == "OL-06" & (regionID %in% c(1,3,4,6,15,16)))) %>%
  filter(!(Position == "OL-07" & 
             (doy < 63 |
                doy > 262))) %>%
  filter(!(Position == "SW-11" & 
             !(date_Date > as.Date("2020-04-26") | date_Date < as.Date("2019-12-10")))) %>%
  mutate(Year = case_when(Position == "WA-01" & date_Date < as.Date("2020-08-01") ~ 
                            as.character(as.integer(Year) + 1),
                          TRUE ~ Year),
         date = case_when(Position == "WA-01" & date_Date < as.Date("2020-08-01") ~ 
                            (date + lubridate::years(1)),
                          TRUE ~ date),
         date_Date = case_when(Position == "WA-01" & date_Date < as.Date("2020-08-01") ~ 
                                 (date_Date + lubridate::years(1)),
                               TRUE ~ date_Date)
  ) %>%
  filter(!(Position == "WA-13" & 
             !(date_Date > as.Date("2020-06-11") | 
                 date_Date < as.Date("2019-08-20")))) %>%
  filter(!(Position == "WH-01" & 
             (doy < 36 |
                doy > 312))) %>%
  filter(!(Position == "WL-10" & (regionID %in% c(2,3,4,5)))) 

## Filter to only include consistent exposure settings
camdata_ExposureControl = lapply(unique(camdata$Position),
                      function(p){
                        pos = camdata %>%
                          filter(Position == p)
                        acceptableExp = exposureFinder(p, depth = 20)
                        pos %>%
                          filter(exposureTime %in% acceptableExp)
                      }) %>%
  bind_rows()


#### Apply time series filters to camdata ----
## Create useful date/time columns
## Remove snowy data
## Find the winterGCC value and apply across position/region/year combinations.
gccWinterVals = camdata_ExposureControl %>%
  mutate(Date = date_Date,
         Year = format(Date, "%Y"),
         DOY = as.numeric(format(Date, "%j"))) %>%
  filter(is.na(Snow),
         Year >= minYr,
         Year <= maxYr) %>%
  # For each position, apply filtering functions
  group_by(Position,regionID,Year) %>%
  mutate(rawGCC = GCC,
         wGCC = winterGCC(GCC))

## Re-join the snow data and assign snowy days the wGCC value
snowRejoiner = function(p){
  snowDat = snowLong %>%
    filter(Position == p) %>%
    rename(Date = date_Date) %>%
    mutate(Year = as.character(lubridate::year(Date))) %>%
    filter(Year >= minYr,
           Year <= maxYr)
  gccDat = gccWinterVals %>%
    filter(Position == p) %>%
    select(-Snow)
  snowDatROI = do.call("rbind",
                       lapply(unique(gccDat$regionID),
                              function(x){
                                d = snowDat
                                d$regionID = x
                                d$species = gccDat$species[gccDat$regionID==x][1]
                                return(d)
                              }))
  fullDat = full_join(gccDat,
                      snowDatROI,
                      by = c("Position","Date","regionID", "Year","species"))
  
}
gccSnowRejoined = do.call("rbind",
                       lapply(unique(gccWinterVals$Position),
                              snowRejoiner)) 

## Convert days that are missing due to snow cover to the winter GCC value
gccPreSmoothing = gccSnowRejoined %>%
  mutate(DOY = as.numeric(lubridate::yday(Date))) %>%
  group_by(Position, regionID, Year) %>%
  dplyr::arrange(DOY) %>%
  mutate(Snow = ifelse(is.na(Snow),
                       0,
                       Snow)) %>%
  mutate(wGCC = ifelse(Snow == 1,
                       min(wGCC, na.rm = T),
                       wGCC)) 

## Preview
gccWinterVals %>%
  filter(Position == "GB-06") %>%
  ggplot(aes(x = Date, y = wGCC, col = as.factor(regionID))) +
  geom_point()
gccPreSmoothing %>%
  filter(Position == "GB-06") %>%
  ggplot(aes(x = Date, y = wGCC, col = as.factor(regionID))) +
  geom_point()

gccSmoothed = gccPreSmoothing %>%
  group_by(Position,regionID,Year) %>%
  mutate(mwGCC = quantilewindow(wGCC,windowsize = 3),
         scGCC = ts_scale(mwGCC)) %>%
  ungroup() %>%
  # Retain only relevant columns
  select(scGCC, mwGCC, wGCC, rawGCC, DOY, Snow,
         Position, species, regionID, 
         totalPix, non0Pix,
         Year, Date, exposureTime) %>%
  # Rename as needed
  rename(GCC = scGCC) %>%
  separate(Position, into = c("Site","Plot"), sep = "-", remove = FALSE)

## Test the phenology function below
x = gccSmoothed %>%
  filter(Position == "GB-06",
         Year == 2020,
         regionID == 1)
plot(x$GCC ~ x$DOY)
gccp = gccPheno(x)
abline(v = gccp$xmidS, col = "green3", lwd = 2)
abline(v = gccp$xmidA, col = "orange", lwd = 2)

gccSmoothed

#### Save outputs ----
write_csv(gccSmoothed, file = paste0(dataDir,"/gccSmoothed.csv"))
write_csv(positionData, file = paste0(dataDir,"/positionData.csv"))
