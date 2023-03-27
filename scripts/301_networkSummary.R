library(tidyverse)
library(CJsBasics)

#### Import data ----
# Topographic and geographic info for camera locations
positionData = read_csv("data/stanPrep/positionData.csv")
# Greenness data
gccData = read_csv("data/stanPrep/gccSmoothed.csv")
# Which cam/ROI's had phenology data generated
gccPhenoSp = read_csv("data/stanPrep/gccPhenology.csv") %>%
  left_join(positionData, by = "Position") %>% 
  filter(!is.na(xmidA))

#### Data checks ----
gccData %>%
  mutate(posDate = paste0(Position, "_", as.character(Date))) %>% 
  pull(posDate) %>%
  unique() %>%
  length()
# 53488

#### Data wrangling ----
## Retain time series that had phenology data
goodPhenoCodes = paste0(gccPhenoSp$Position,"-", gccPhenoSp$Year,"-",gccPhenoSp$regionID)
## Find camera-wise mean daily GCC (aggregate across regions)
## Create a new date column for overlapped plotting in Fig 1.
gccSmoothed = gccData %>%
  left_join(positionData, by = "Position") %>%
  mutate(posRegYr = paste0(Position,"-",Year,"-",regionID)) %>%
  filter(posRegYr %in% goodPhenoCodes) %>%
  filter(!is.na(species),
         species != "ROCK ROCK",
         species != "TALUS TALUS",
         species != "DG DG") %>%
  group_by(regionID, Year, Position, Elevation, Latitude) %>%
  filter(all(sum(diff(GCC), na.rm = T) > 0)) %>%
  mutate(date_Date = Date) %>%
  ungroup() %>%
  group_by(Position, date_Date, Year, Elevation, Latitude) %>%
  summarize(GCC = mean(GCC,na.rm=T)) %>%
  ungroup() %>%
  mutate(DOY = strftime(date_Date, format = "%j"),
         date_reform = as.Date(paste0("2013-",DOY),
                               format = "%Y-%j")) %>%
  select(-DOY)

#### Analysis ----
## > Compare phenology across years ----
## gccPhenoCam is used to compare green-up timing across 
gccPhenoCam = gccPhenoSp %>%
  filter(species != "ROCK ROCK",
         species != "TALUS TALUS",
         species != "DG DG") %>%
  group_by(Position, Year,Elevation,Slope,Aspect,TPI,Latitude,
           NLCD,LandCoverHigherOrder,LandCover) %>%
  dplyr::summarize(meanMidS = mean(xmidS, na.rm = TRUE),
                   meanMidA = mean(xmidA, na.rm = TRUE)) %>%
  ungroup()
greenupComp = gccPhenoCam %>%
  select(Year,Position,meanMidS) %>%
  mutate(Year = paste0("Yr_",Year)) %>%
  pivot_wider(names_from = "Year",
              values_from = "meanMidS")
senescComp = gccPhenoCam %>%
  select(Year,Position,meanMidA) %>%
  mutate(Year = paste0("Yr_",Year)) %>%
  pivot_wider(names_from = "Year",
              values_from = "meanMidA")
## > Green-up ----
t.test(greenupComp$Yr_2019,
       greenupComp$Yr_2020,
       paired=T)
# Paired t-test
# 
# data:  greenupComp$Yr_2019 and greenupComp$Yr_2020
# t = 6.5009, df = 14, p-value = 1.4e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   10.10082 20.04741
# sample estimates:
#   mean of the differences 
# 15.07412 

## > Senescence ----
t.test(senescComp$Yr_2019,
       senescComp$Yr_2020,
       paired=T)
# Paired t-test
# 
# data:  senescComp$Yr_2019 and senescComp$Yr_2020
# t = -0.98754, df = 14, p-value = 0.3401
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -12.227545   4.517501
# sample estimates:
#   mean of the differences 
# -3.855022 

#### Generate plots ----
## > Network summary figure ----
## > - Camera-level time series
## Elevation
elGreenness = gccSmoothed %>%
  ggplot(aes(x = date_reform, y = GCC)) +
  geom_line(aes(col = Elevation, group = Position)) +
  facet_wrap(~Year, nrow = 2) +
  scale_color_viridis_c("Elevation (m)") +
  scale_x_date(date_labels = "%b") +
  ylab("Greenness") +
  xlab("Date") +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.2,0.875),
        legend.title.align = 0,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.background = element_blank(),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in"))
## Latitude
latGreenness = gccSmoothed %>%
  ggplot(aes(x = date_reform, y = GCC)) +
  geom_line(aes(col = Latitude, group = Position)) +
  facet_wrap(~Year, nrow = 2) +
  scale_color_viridis_c("Latitude (°N)",
                        option = "A") +
  scale_x_date(date_labels = "%b") +
  ylab("Greenness") +
  xlab("Date") +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.2,0.875),
        legend.title.align = 0,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.background = element_blank(),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in"))

## > - Compile
seasonalGreenness = cowplot::plot_grid(elGreenness,
                                       latGreenness,
                                       nrow = 1,
                                       ncol = 2,
                                       labels = c("A","B"),
                                       align = "hv",
                                       axis = "trbl") +
  theme(plot.margin = margin(0,0.25,0,0,"in"))

#### Save ----
ggsave("plots/fig2.jpg",
       seasonalGreenness,
       width = 6, height = 4, units = "in",dpi=300)
