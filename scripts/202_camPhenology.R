#### Import packages ----
library(tidyverse)
library(sf)
source("scripts/001_cameraFunctions.R")



set.seed(1957) # Hutchinson concluding remarks

#### Define parameters ----
minYr = 2019
maxYr = 2020
plotDir = "plots"
dataDir = "data/stanPrep"


#### Import data ----
gccSmoothed = read_csv("data/stanPrep/gccSmoothed.csv")

# For each year at each position, apply the gccPheno function
# Takes about 13 minutes for 2 full years of timelapse network data
startTime = Sys.time()
gccPhenology = gccSmoothed %>%
  group_by(Position, regionID, Year) %>%
  do(gccPheno(x = .))
stopTime = Sys.time()
print(stopTime-startTime)
# BRRR::skrrrahh(14)

gccPhenology = gccPhenology %>%
  left_join(gccSmoothed %>%
              select(Position,regionID,species) %>%
              filter(!duplicated(paste0(Position,"-",regionID))),
            by=c("Position","regionID"))

# For each year at each position, apply the snowPheno function
snowPhenoData = gccSmoothed %>%
  separate(Position, 
           into = c("Site","Plot"), 
           sep = "-", remove = FALSE) %>%
  mutate(Year = format(Date, "%Y"),
         DOY = as.numeric(format(Date, "%j"))) %>%
  # For each position, apply filtering functions
  group_by(Position,regionID) %>%
  ungroup() %>%
  # Rename as needed
  mutate(DOY = as.numeric(DOY)) %>%
  select(Position,Date,Snow) %>%
  mutate(Snow = replace_na(Snow, 0))  %>%
  mutate(snowYear = ifelse(Date < "2018-08-15",
                           "2017-2018",
                           ifelse(Date < "2019-08-15",
                                  "2018-2019",
                                  ifelse(Date < "2020-08-15",
                                         "2019-2020",
                                         "2020-2021")))) %>%
  group_by(Position, Date, snowYear) %>%
  summarize(Snow = mean(Snow))

## Takes 1 min for 2 yrs
startTime = Sys.time()
snowPhenology = snowPhenoData  %>%
  mutate(Year = format(Date, "%Y"),
         DOY = as.numeric(format(Date, "%j"))) %>%
  group_by(Position, snowYear) %>%
  do(snowPheno2(x = .))
stopTime = Sys.time()
print(stopTime-startTime)
# BRRR::skrrrahh(14)

#### Inspect input data ----
## Timing of spring green-up by plot
camLevelChange = gccPhenology %>%
  mutate(Year = as.integer(Year)) %>%
  ggplot(aes(x = Year,
             y = xmidS)) +
  geom_point(aes(group = paste0(Position,regionID),
                 col = species)) +
  geom_line(aes(group = paste0(Position,regionID),
                col = species)) +
  stat_smooth(aes(col = species),
              method = "lm", col = "black", size = 2) +
  scale_x_continuous(breaks = c(2019,2020),
                     limits = c(2018.95,2020.05)) +
  CJsBasics::BasicTheme +
  theme(legend.position = "none") +
  ylab("Green-up timing (day of year)")
camLevelChange


## Look at a single plant species for a single year
elymusExample = ggplot() +
  geom_vline(data = gccPhenology %>%
               filter(species == "Elymus elymoides",
                      !is.na(xmidS),
                      Year == 2020),
             aes(xintercept = as.Date(paste0(Year, "-",
                                             xmidS),
                                      format = "%Y-%j"),
                 col = Position),
             size = 0.75,
             lty = 2,
             alpha = 0.5) +
  geom_line(data = gccSmoothed %>%
              filter(Position %in% 
                       gccPhenology$Position[gccPhenology$Year==2020 &
                                               !is.na(gccPhenology$xmidS)],
                     species == "Elymus elymoides",
                     Year == 2020),
            aes(x = Date,
                y = GCC,
                col = Position,
                group = paste0(Position,species,regionID)),
            size = 0.5) +
  ggtitle("Squirreltail") +
  CJsBasics::BasicTheme
elymusExample


## Track change across years
phenoChange = gccPhenology %>%
  filter(species != "ROCK ROCK",
         species != "TALUS TALUS",
         species != "DG DG",
         Year >= minYr,
         Year <= maxYr) %>%
  pivot_wider(id_cols = c(Position, regionID, species), 
              names_from = Year,
              values_from = c(xmidS, xmidA, scalS, scalA)) %>%
  ggplot() +
  geom_abline(lty = 2) +
  geom_point(aes(x = xmidS_2019,
                 y = xmidA_2019,
                 col = species),
             arrow = arrow(length = unit(0.1,"cm")),
             size = 0.25) +
  geom_point(aes(x = xmidS_2020,
                 y = xmidA_2020,
                 col = species),
             arrow = arrow(length = unit(0.1,"cm")),
             size = 0.25) +
  geom_segment(aes(x = xmidS_2019,
                   xend = xmidS_2020,
                   y = xmidA_2019,
                   yend = xmidA_2020,
                   col = species),
               arrow = arrow(length = unit(0.1,"cm")),
               size = 0.25) +
  stat_smooth(aes(x = xmidS_2019,
                  y = xmidA_2019),
              col = "orange2",
              method = "lm") +
  stat_smooth(aes(x = xmidS_2020,
                  y = xmidA_2020),
              col = "blue3",
              method = "lm") +
  xlim(0,366) +
  ylim(0,366) +
  xlab("Green-up timing (Day of year)") +
  ylab("Senescence timing (Day of year)") +
  CJsBasics::BasicTheme +
  theme(legend.position = "none")
phenoChange

#### Save outputs ----
write_csv(snowPhenology, file = paste0(dataDir,"/snowPhenology.csv"))
write_csv(gccPhenology, file = paste0(dataDir,"/gccPhenology.csv"))
