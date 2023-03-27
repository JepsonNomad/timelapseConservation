#### Cleaning MOD data from pointwise extraction in Earth Engine ----
library(tidyverse)


#### Import data ----
modNDVIwide = read_csv("data/modis/timelapseCams_MOD13Q1_NDVI.csv")
modCDOYwide = read_csv("data/modis/timelapseCams_MOD13Q1_cDOY.csv")


#### Data wrangling ----
# Remove unnecessary columns
drops = c("system:index", ".geo", "PlotNum", "Site")
modNDVIwide = modNDVIwide[, !(names(modNDVIwide) %in% drops)]
modCDOYwide = modCDOYwide[, !(names(modCDOYwide) %in% drops)]
names(modCDOYwide)

# Transpose data.frame
modNDVI = modNDVIwide %>%
  rename("Position" = PlotID) %>%
  pivot_longer(-Position,
               names_to = "Date",
               values_to = "NDVI_raw") %>%
  # Replace "No data" with NA
  mutate(NDVI = ifelse(NDVI_raw == -99999,
                       NA,
                       NDVI_raw)) %>%
  # Convert NDVI to numeric
  mutate(NDVI = as.numeric(NDVI)) %>%
  # Convert Date column to date
  mutate(posDT = as.POSIXct(Date,
                            format = "%Y-%m-%dT%H:%M:%S")) %>%
  # Convert datetime to date 
  mutate(TSymd = lubridate::as_date(posDT)) %>%
  select(-Date, -NDVI_raw, -posDT) %>%
  rename(Date = TSymd) 
modCDOY = modCDOYwide %>%
  rename("Position" = PlotID) %>%
  pivot_longer(-Position,
               names_to = "Date",
               values_to = "cDOY_raw") %>%
  # Replace "No data" with NA
  mutate(cDOY = ifelse(cDOY_raw == -99999,
                       NA,
                       cDOY_raw)) %>%
  # Convert NDVI to numeric
  mutate(cDOY = as.numeric(cDOY)) %>%
  # Convert Date column to date
  mutate(posDT = as.POSIXct(Date,
                            format = "%Y-%m-%dT%H:%M:%S")) %>%
  # Convert datetime to date 
  mutate(TSymd = lubridate::as_date(posDT)) %>%
  select(-Date, -cDOY_raw, -posDT) %>%
  rename(Date = TSymd) 

## Look at it
head(modNDVI)
head(modCDOY)

## Join datasets 
mydata = full_join(modNDVI,
                   modCDOY,
                   by = c("Date", "Position"))

# Use composite day of year to adjust date
mydata = mydata %>%
  mutate(Year = lubridate::year(Date),
         cDOY = ifelse(is.na(cDOY),
                       lubridate::yday(Date),
                       cDOY)) %>%
  mutate(yCDOY = paste0(Year, "-", cDOY)) %>%
  mutate(CompositeDate = as.Date(yCDOY, format = "%Y-%j")) %>%
  separate(Position,
           into = c("Site","Plot"),
           sep = "-") %>%
  mutate(Plot = str_pad(Plot,
                        width = 2,
                        side = "left",
                        pad = "0")) %>%
  mutate(Site = toupper(Site)) %>%
  unite(Site, Plot, col = Position, sep = "-") %>%
  mutate(NDVI = NDVI/10000) %>%
  select(-cDOY, -Date, -yCDOY, -Year) %>%
  rename(Date = CompositeDate)


#### Save ----
write.table(x = mydata,
            file = "data/modis/timelapseCams_NDVI_tidy.csv",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)



#### Plots ----
ggplot(mydata,
       aes(x = Date,
           y = NDVI,
           group = Position,
           color = Position)) +
  geom_line() +
  theme(legend.position = "none")
