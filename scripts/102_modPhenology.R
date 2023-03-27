#### Load packages
library(tidyverse)
source("scripts/000_modisFunctions.R")



set.seed(1927) # Elton Animal Ecology

#### Import data ----
myNDVI = read_csv("data/modis/timelapseCams_NDVI_tidy.csv")
mySites = myNDVI %>%
  separate(Position, into = c("Site","Plot"), 
           sep = "-", remove = FALSE) %>%
  select(Site) %>%
  arrange(Site) %>%
  unique() %>%
  pull()

# Preview the data
myNDVI %>%
  filter(Position == "GB-02") %>%
  ggplot(aes(x = Date,
             y = NDVI,
             col = Position)) +
  geom_line() +
  geom_point() +
  theme(legend.position = "none")

#### Data wrangling ----
myNDVI = myNDVI %>%
  separate(Position, 
           into = c("Site","Plot"), 
           sep = "-", remove = FALSE) %>%
  mutate(Year = format(Date, "%Y"),
         DOY = as.numeric(format(Date, "%j")))


#### Calculate phenology data ----
phenoOperator = function(x){
  y = x  %>%
    # For each position, apply filtering functions
    group_by(Position) %>%
    mutate(wNDVI = winterNDVI(NDVI),
           mwNDVI = medianwindow(wNDVI),
           scNDVI = ts_scale(mwNDVI)) %>%
    ungroup() %>%
    # Retain only relevant columns
    select(scNDVI, mwNDVI, wNDVI, DOY, Position, Year, Date) %>%
    # Rename as needed
    mutate(DOY = as.numeric(DOY)) %>%
    rename(NDVI = scNDVI) %>%
    # For each year at each position, apply the yearphenology funciton
    group_by(Position, Year) %>%
    do(yearphenology(x = .)) 
}
# Break the data into 3
NDVIdat1 = myNDVI %>%
  filter(Site %in% mySites[1:3])
NDVIdat2 = myNDVI %>%
  filter(Site %in% mySites[4:6])
NDVIdat3 = myNDVI  %>%
  filter(Site %in% mySites[7:9])
# Apply phenoOperator to the three datasets
NDVIpheno1 = phenoOperator(NDVIdat1)
NDVIpheno2 = phenoOperator(NDVIdat2)
NDVIpheno3 = phenoOperator(NDVIdat3)
# bind the data
NDVIpheno = rbind(ungroup(NDVIpheno1),
                  ungroup(NDVIpheno2), 
                  ungroup(NDVIpheno3))


#### Save outputs ----
write.table(NDVIpheno,
            file = "data/modis/modisPhenology.csv",
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

