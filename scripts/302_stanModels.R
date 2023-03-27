#### Examining snowmelt and green-up variation in the eastern Sierra using Bayesian mixed effects models 

#### Import packages ----
library(brms)
library(bayesplot)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(emmeans)
library(cowplot)

#### Establish STAN options ----
## Create output directory
if(!dir.exists("stan_models")){
  dir.create("stan_models")
}

options(mc.cores = parallel::detectCores())
## NUTS algorithm hyperparameters
control. <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.95,
  stepsize = 0.1,
  max_treedepth = 10
)

mySeed = 16270125 # Robert Boyle born
set.seed(mySeed)

#### Define parameters ----
minYr = 2019
maxYr = 2020

#### Import data ----
snowPhenology = read_csv("data/stanPrep/snowPhenology.csv")
gccPhenology = read_csv("data/stanPrep/gccPhenology.csv")
positionData = read_csv("data/stanPrep/positionData.csv")
modPhenology = read_csv("data/modis/modisPhenology.csv")

#### Data summaries ----
## camera-years
gccPhenology %>% 
  filter(!is.na(xmidS)) %>%
  group_by(Position,Year) %>%
  summarize(n = n()) %>%
  select(Year, Position) %>%
  pivot_wider(names_from="Year", values_from = "Position")
## 78 and 16 for 2020 and 2019 respectively

## camera-years
gccPhenology %>% 
  filter(!is.na(xmidS)) %>%
  group_by(paste0(Position,Year)) %>%
  summarize(n = n())
# A tibble: 94 × 2

## region-years
gccPhenology %>% 
  filter(!is.na(xmidS)) %>%
  group_by(paste0(Position,regionID)) %>%
  summarize(n = n())
# A tibble: 1,136 × 2

gccPhenology %>% 
  filter(!is.na(xmidS)) %>%
  group_by(Position, regionID, Year) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(posReg = paste0(Position, regionID)) %>%
  select(Year, posReg) %>%
  pivot_wider(names_from="Year", values_from = "posReg")
# 1,104 and 204 for 2020 and 2019 respectively

#### Data wrangling ----
#### > Snow ----
## Convert snow phenology data to a useful long format.
snowWide = snowPhenology %>%
  rename(snowFall = xmidS,
         fallRate = scalS,
         snowMelt = xmidA,
         meltRate = scalA) %>%
  separate(snowYear, into = c("FallYear","MeltYear"))
snowFallData = snowWide %>%
  select(Position, FallYear, snowFall, fallRate) %>%
  filter(!is.na(snowFall)) %>%
  rename(Year = FallYear)
snowMeltData = snowWide %>%
  select(Position, MeltYear, snowMelt, meltRate) %>%
  filter(!is.na(snowMelt)) %>%
  rename(Year = MeltYear)


snowWrangled = full_join(snowFallData,
                         snowMeltData,
                         by = c("Position", "Year")) %>%
  arrange(Position, Year) %>%
  select(Position, Year, snowMelt, meltRate, snowFall, fallRate) %>%
  filter(Year >= minYr,
         Year <= maxYr)
snowWrangled

## Compile snowmelt and positionData for snow questions
## Deals with autumn snowmelt by identifying when the date of snowmelt fell in the previous calendar year
## and modifies snowMeltDOY to be a negative value (i.e. super early) rather than an erroneously large
## positive value (i.e. snow melt happens super late the following fall)
snowMeta = left_join(snowWrangled,
                     positionData,
                     by = c("Position")) %>%
  mutate(cosAsp = cos((Aspect*pi)/180),
         elCent = scale(Elevation,center = TRUE, scale = F)[,1],
         latCent = scale(Latitude,center = TRUE, scale = F)[,1],
         elScal = scale(Elevation,center = TRUE, scale = T)[,1],
         latScal = scale(Latitude,center = TRUE, scale = T)[,1],
         yearID = as.integer(factor(Year)),
         posID = as.integer(factor(Position))) %>%
  mutate(snowMeltRAWYEAR = as.numeric(format(snowMelt, "%Y")),
         snowSeasonYear = as.numeric(Year),
         snowMeltDOY = case_when(snowSeasonYear == snowMeltRAWYEAR ~
                                   as.numeric(format(snowMelt, "%j")),
                                 snowSeasonYear == snowMeltRAWYEAR + 1 ~
                                   (as.numeric(format(snowMelt, "%j"))-365)))
min(snowMeta$Elevation)
max(snowMeta$Elevation)
mean(snowMeta$Elevation)
# 2889.847
sd(snowMeta$Elevation)
# 634.3736
# 1 standard deviation in elevation = 626 m elevation

#### > Green-up ----
gccWrangled = gccPhenology %>%
  mutate(fullRegionInfo = paste0(Position,"_",
                                 regionID),
         yearID = as.integer(factor(Year)),
         posID = as.integer(factor(Position)),
         roiID = as.integer(factor(fullRegionInfo)),
         sppID = as.integer(factor(species)),
         Year = as.character(Year),
         species = ifelse(species == "VEG VEG",
                          "Unknown/mixed",
                          species))

## Generate list of species observed in at least 3 camera locations;
## identify corresponding growth forms based on calflora (using graminoid
## rather than grass)
majorSpp = gccWrangled %>%
  group_by(species) %>%
  summarize(n = length(unique(posID))) %>%
  filter(n > 2) %>%
  pull(species)
# write_csv(x = data.frame("species" = majorSpp),
#           file = "data/majorSpecies.csv")
growthFormLookupTable = read_csv("data/majorSpecies_growthform.csv",
                                 skip = 1)

## Compile data for stan
y_data = left_join(gccWrangled,
                   positionData,
                   by = c("Position")) %>%
  left_join(snowWrangled,
            by = c("Position", "Year")) %>%
  left_join(growthFormLookupTable,
            by = c("species")) %>%
  filter(species %in% majorSpp) %>%
  mutate(cosAsp = cos((Aspect*pi)/180),
         elCent = scale(Elevation,center = TRUE, scale = F)[,1],
         latCent = scale(Latitude,center = TRUE, scale = F)[,1],
         elScal = scale(Elevation,center = TRUE, scale = T)[,1],
         latScal = scale(Latitude,center = TRUE, scale = T)[,1],
         species = case_when(
           species == "Packera cana\n" ~ "Packera cana",
           species == "Unknown forb\n" ~ "Unknown forb",
           species == "zlinanthus pungens" ~ "Linanthus pungens",
           species == "ROCK ROCK\n" ~ "ROCK ROCK",
           TRUE ~ species)) %>%
  mutate(snowMeltRAWYEAR = as.numeric(format(snowMelt, "%Y")),
         snowSeasonYear = as.numeric(Year),
         snowMeltDOY = case_when(
           snowSeasonYear == snowMeltRAWYEAR ~
             as.numeric(format(snowMelt, "%j")),
           snowSeasonYear == snowMeltRAWYEAR + 1 ~
             (as.numeric(format(snowMelt, "%j"))-365))) %>%
  filter(!is.na(xmidS),
         species != "ROCK ROCK",
         species != "TALUS TALUS",
         species != "DG DG",
         species != "Unknown/mixed",
         species != "Pinus sp") %>%
  filter(Year >= minYr,
         Year <= maxYr)
## Inspect
head(y_data)
# View(y_data)
summary(y_data)

y_data %>%
  group_by(Position, Elevation, Year) %>%
  summarize(snowmelt = mean(snowMeltDOY)) %>%
  ggplot(aes(x = Year,
             y = snowmelt,
             col = Elevation)) +
  geom_jitter(width = 0.1) +
  scale_color_viridis_c() +
  CJsBasics::BasicTheme


#### > Satellite ----
modPhenology = modPhenology %>%
  rename(modmidS = xmidS,
         modmidA = xmidA,
         modscalS = scalS,
         modscaclA = scalA) %>%
  mutate(Year = as.character(Year))

#### Data analysis ----
## > Green-up vs Snowmelt timing ----
gfFit = brms::brm(xmidS ~ snowMeltDOY + growthform +
                    (1 | Year + Position),
                  data = y_data, family = gaussian(),
                  warmup = 2500, iter = 5000, chains = 4,
                  control = control.,
                  seed = mySeed)
gfFit
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: xmidS ~ snowMeltDOY + growthform + (1 | Year + Position) 
# Data: y_data (Number of observations: 810) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    13.01      1.31    10.70    16.03 1.00     1046      640
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     9.63      8.78     0.76    34.56 1.00     1449      563
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                   99.63     10.86    77.91   121.21 1.00     1821     1369
# snowMeltDOY                  0.46      0.03     0.40     0.51 1.00     1751     3275
# growthformPerennialgrass    -2.81      7.62   -17.98    12.09 1.00     2862     4516
# growthformPerennialherb      0.18      7.68   -15.01    14.98 1.00     2793     4601
# growthformShrub              7.59      7.62    -7.20    22.39 1.00     2610     4732
# growthformTree               5.04     13.00   -20.32    30.94 1.00     4541     5829
# growthformUnknownDmixed      0.30      7.69   -14.81    15.19 1.00     2798     4751
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    11.66      0.31    11.08    12.28 1.00     8637     6528
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
saveRDS(gfFit,"stan_models/gfFit.RDS")
mcmc_plot(gfFit,type = "trace")
mcmc_plot(gfFit,type = "rhat_hist")
bayes_R2(gfFit)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.8537269 0.004399868 0.8444582 0.8616645

# How does the delay between snowmelt and green-up vary along a gradient of snowmelt timing? (Reviewer request)
brms::brm((xmidS-snowMeltDOY) ~ snowMeltDOY + growthform +
            (1 | Year + Position),
          data = y_data, family = gaussian(),
          warmup = 2500, iter = 5000, chains = 4,
          control = control.,
          seed = mySeed)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: (xmidS - snowMeltDOY) ~ snowMeltDOY + growthform + (1 | Year + Position) 
# Data: y_data (Number of observations: 810) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    12.94      1.29    10.62    15.68 1.01     1069     1433
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     8.47      7.95     0.75    30.66 1.00      648      270
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                   97.92     10.09    77.17   117.59 1.01      871      163
# snowMeltDOY                 -0.54      0.03    -0.60    -0.49 1.00     1299     3165
# growthformPerennialgrass    -2.89      7.52   -17.67    11.87 1.00     2351     4388
# growthformPerennialherb      0.15      7.56   -14.81    14.97 1.00     2348     4320
# growthformShrub              7.51      7.51    -7.36    22.26 1.00     2314     4270
# growthformTree               4.68     13.14   -20.74    30.58 1.00     3445     6042
# growthformUnknownDmixed      0.25      7.58   -14.69    15.03 1.00     2270     4110
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    11.67      0.31    11.09    12.28 1.00     2730     7376
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

## Longest intra-group green-up range
y_data %>%
  group_by(Year, growthform) %>%
  summarize(yrMin = min(xmidS),
            yrMax = max(xmidS)) %>%
  mutate(totalDiff = yrMax-yrMin) %>%
  group_by(Year) %>%
  filter(totalDiff == max(totalDiff))


y_data %>%
  ggplot(aes(x = snowMeltDOY,
             y = (xmidS-snowMeltDOY))) +
  geom_point() +
  stat_smooth(method = "lm")

## Figure 3
gfReg = y_data %>%
  filter(growthform != "Tree") %>%
  ggplot(aes(x = snowMeltDOY,
             y = xmidS)) +
  geom_point(col = "black",
             size = 1.5) +
  geom_point(aes(col = growthform),
             size = 1) +
  scale_color_manual("Growth form",
                     values = CJsBasics::KellyCols[8:20]) +
  stat_smooth(method = "lm",
              aes(col = growthform)) +
  geom_point(data = y_data[y_data$growthform=="Tree",],
             aes(x = snowMeltDOY,
                 y = xmidS),
             col = "black") +
  facet_grid(~ Year) +
  ylab("Green-up timing (day of year)") +
  xlab("Snowmelt timing (day of year)") +
  geom_abline()  +
  CJsBasics::BasicTheme +
  scale_x_continuous(breaks=c(0,100,200),
                     limits = c(min(c(y_data$xmidS, y_data$snowMeltDOY),
                                    na.rm = TRUE),
                                max(c(y_data$xmidS, y_data$snowMeltDOY),
                                    na.rm = TRUE))) +
  scale_y_continuous(breaks=c(0,100,200),
                     limits = c(min(c(y_data$xmidS, y_data$snowMeltDOY),
                                    na.rm = TRUE),
                                max(c(y_data$xmidS, y_data$snowMeltDOY),
                                    na.rm = TRUE)))+
  coord_equal()
ggsave(filename = "plots/fig3.jpg",
       gfReg,
       width = 7, height = 3, units = "in", dpi = 300)


## > Latitudinal variation in lapse rate ----
cor(y_data$Latitude, y_data$Elevation, method = "pearson")
cor(y_data$elScal, y_data$latScal, method = "pearson")

## Analyzed at the ROI level
## Fit an additive model
lapsePlusLat = brms::brm(xmidS ~ elScal + cosAsp + latScal +
                           (1 | Year + Position),
                         data = y_data, family = gaussian(),
                         prior = c(prior(normal(100, 100), 
                                         class = Intercept),
                                   prior(cauchy(0, 10), 
                                         class = sd),
                                   prior(cauchy(0, 10), 
                                         class = sigma)),
                         warmup = 2500, iter = 5000, chains = 4,
                         control = control.,
                         seed = mySeed,
                         save_pars = save_pars(all = TRUE))
lapsePlusLat
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: xmidS ~ elScal + cosAsp + latScal + (1 | Year + Position) 
# Data: y_data (Number of observations: 810) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    12.51      1.20    10.39    15.10 1.00     2241     4124
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    15.83     11.20     4.78    48.81 1.00     4025     3493
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept   157.65     13.10   125.58   181.72 1.00     2844     2191
# elScal       24.50      1.64    21.23    27.64 1.00     1518     2853
# cosAsp        8.63      2.32     4.00    13.23 1.00     1459     2835
# latScal       0.60      1.63    -2.59     3.80 1.00     1407     3083
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    11.99      0.31    11.40    12.62 1.00    11715     7505
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

## Interpreting cosAsp coefficient: Note that the estimate here is for 1 unit difference, which is north vs a slope with no north-south bias. So, to compare north- (cosAsp = 1) and south-facing slopes (cosAsp = -1), double the coefficients.
coefficients(lapsePlusLat)$Year[,,"cosAsp"][1,1]*2
# [1] 17.25764
coefficients(lapsePlusLat)$Year[,,"cosAsp"][1,2]*2
# [1] 4.644436
coefficients(lapsePlusLat)$Year[,,"cosAsp"][1,3]*2
# [1] 7.997995
coefficients(lapsePlusLat)$Year[,,"cosAsp"][1,4]*2
# [1] 26.46837


## Repeat but with unscaled covariates
lapsePlusLatUnsc = brms::brm(
  xmidS ~ Elevation + cosAsp + latScal +
    (1 | Year + Position),
  data = y_data, family = gaussian(),
  prior = c(prior(normal(100, 100), 
                  class = Intercept),
            prior(cauchy(0, 10), 
                  class = sd),
            prior(cauchy(0, 10), 
                  class = sigma)),
  warmup = 2500, iter = 5000, chains = 4,
  control = control.,
  seed = mySeed,
  save_pars = save_pars(all = TRUE))
lapsePlusLatUnsc
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: xmidS ~ Elevation + cosAsp + latScal + (1 | Year + Position) 
# Data: y_data (Number of observations: 810) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    12.47      1.21    10.31    15.07 1.00     1696     1212
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    16.47     12.39     4.82    50.69 1.00     2854     2217
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    32.68     16.91    -2.58    65.92 1.00     1126      769
# Elevation     0.04      0.00     0.04     0.05 1.00     1129      660
# cosAsp        8.67      2.31     4.02    13.16 1.00     1513     2834
# latScal       0.59      1.63    -2.59     3.88 1.00     1492     2704
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    11.99      0.32    11.38    12.63 1.00     8573     7362
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

## Convert elevational lapse rate to terms of lowest and highest elevation cameras
coefficients(
  lapsePlusLatUnsc
  )$Year[,,"Elevation"][1,1]*max(
    y_data$Elevation) - 
  coefficients(
    lapsePlusLatUnsc
    )$Year[,,"Elevation"][1,1]*min(
      y_data$Elevation
      )
# [1] 98.03672

saveRDS(lapsePlusLat, "stan_models/lapsePlusLat.RDS")
mcmc_plot(lapsePlusLat,type = "trace")
mcmc_plot(lapsePlusLat,type = "rhat_hist")
pp_check(lapsePlusLat, draws = 100)

## Fit model with an interaction term
lapseIntxLat = brms::brm(xmidS ~ elScal + cosAsp + latScal + 
                           elScal*latScal +
                           (1 | Year + Position),
                         data = y_data, family = gaussian(),
                         warmup = 2500, iter = 5000, chains = 4,
                         prior = c(prior(normal(100, 100), 
                                         class = Intercept),
                                   prior(cauchy(0, 10), 
                                         class = sd),
                                   prior(cauchy(0, 10), 
                                         class = sigma)),
                         control = control.,
                         seed = mySeed,
                         save_pars = save_pars(all = TRUE))
lapseIntxLat
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: xmidS ~ elScal + cosAsp + latScal + elScal * latScal + (1 | Year + Position) 
# Data: y_data (Number of observations: 810) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    12.59      1.23    10.44    15.21 1.01     1283     3530
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    16.64     12.31     4.91    51.30 1.01      331       71
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept        161.04     15.37   133.23   201.39 1.02      252       56
# elScal            23.78      1.81    20.23    27.41 1.01     1401     3702
# cosAsp             8.58      2.37     3.56    13.29 1.00     1151      637
# latScal            0.65      1.63    -2.63     3.88 1.01     1735     3115
# elScal:latScal    -1.13      1.59    -4.28     2.03 1.00     1775     3456
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    11.99      0.30    11.40    12.61 1.00     6489     7472
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
saveRDS(lapseIntxLat, "stan_models/lapseIntxLat.RDS")
mcmc_plot(lapseIntxLat,type = "trace")
mcmc_plot(lapseIntxLat,type = "rhat_hist")
pp_check(lapseIntxLat, draws = 100)

## Interpreting cosAsp coefficient: Note that the estimate here is for 1 unit difference, which is north vs a slope with no north-south bias. So, to compare north- (cosAsp = 1) and south-facing slopes (cosAsp = -1), double the coefficients.
coefficients(lapseIntxLat)$Year[,,"cosAsp"][1,1]*2
# [1] 17.15113
coefficients(lapseIntxLat)$Year[,,"cosAsp"][1,2]*2
# [1] 4.741412
coefficients(lapseIntxLat)$Year[,,"cosAsp"][1,3]*2
# [1] 7.123434
coefficients(lapseIntxLat)$Year[,,"cosAsp"][1,4]*2
# [1] 26.57502


## Repeat the unscaled covariate model but with interaction term
lapseIntxLatUnsc = brms::brm(
  xmidS ~ Elevation + cosAsp + latScal +
    Elevation*latScal +
    (1 | Year + Position),
  data = y_data, family = gaussian(),
  prior = c(prior(normal(100, 100), 
                  class = Intercept),
            prior(cauchy(0, 10), 
                  class = sd),
            prior(cauchy(0, 10), 
                  class = sigma)),
  warmup = 2500, iter = 5000, chains = 4,
  control = control.,
  seed = mySeed,
  save_pars = save_pars(all = TRUE))
lapseIntxLatUnsc
## Convert elevational lapse rate to terms of lowest and highest elevation cameras
coefficients(
  lapseIntxLatUnsc
)$Year[,,"Elevation"][1,1]*max(
  y_data$Elevation) - 
  coefficients(
    lapseIntxLatUnsc
  )$Year[,,"Elevation"][1,1]*min(
    y_data$Elevation
  )
# [1] 95.95115


## See https://psyc-bayes-notes.netlify.app/model-comparison-and-regularization.html
## for brms model comparison and interpretation
lapsePlusLat <- add_criterion(lapsePlusLat, "loo", 
                              reloo = TRUE, moment_match = TRUE)
lapseIntxLat <- add_criterion(lapseIntxLat, "loo", 
                              reloo = TRUE, moment_match = TRUE)
loo_compare(lapsePlusLat, lapseIntxLat, criterion = "loo")
# elpd_diff se_diff
# lapseIntxLat 0.0       0.0    
# lapsePlusLat 0.0       0.6    
bayes_R2(lapsePlusLat)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.8450638 0.004632916 0.8351615 0.8533049


## Figure 4
## Aggregate by camera position
z_data = y_data %>%
  group_by(Year, Position, elScal, Elevation, latScal, cosAsp) %>%
  summarize(xmidS = mean(xmidS, na.rm = T)) %>%
  ungroup() %>%
  mutate(latScalMean = Hmisc::cut2(latScal, 
                                   g = 3),
         latScalGroup = as.numeric(latScalMean),
         latScalGroup = ifelse(latScalGroup == 1,
                               "South",
                               ifelse(latScalGroup == 2,
                                      "Central",
                                      "North")),
         latScalGroup = factor(latScalGroup,
                               levels = c("South","Central","North"))) 

## elevation and latitude (which were not correlated; Pearson’s R = 0.39)
cor(z_data$elScal,z_data$latScal, method="pearson")
# [1] 0.385985

greenupLatFacets2020 = ggplot(data=z_data[z_data$Year==2020,],
                              aes(x=Elevation,y=xmidS)) +
  geom_point(col = "black",
             size = 2.5) +
  geom_point(aes(col = cosAsp),
             size = 2) +
  stat_smooth(method = "lm",
              col = "lightblue4") + 
  facet_wrap(~latScalGroup) +
  scale_color_viridis_c("Terrain\nAspect",
                        option = "A",limits = c(-1,1),
                        breaks = c(-1,0,1),
                        labels = c("South","","North")) +
  scale_x_continuous(breaks = c(1500,2500,3500)) +
  xlab("Elevation (m)") +
  ylab("Green-up timing (day of year)") +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.93,0.15),
        legend.title.align = 0,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.background = element_blank(),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.045, "in"),
        panel.spacing = unit(0.1,"in"))
greenupLatFacets2020
ggsave("plots/fig4.jpg",
       greenupLatFacets2020,
       width = 5, height = 3, units = "in", dpi = 300)


## > Species variation in lapse rate ----
## Retain only identifiable species for which there were at least 3 camera positions with observations within at least one year. 
q_data_sppSelection = y_data %>%
  separate(species, 
           into=c("Genus", "SpecificEpithet"), remove=FALSE) %>%
  filter(Genus != "Unknown",
         SpecificEpithet != "sp") %>%
  group_by(species,Year) %>%
  summarize(n = n_distinct(Position)) %>%
  arrange(desc(n)) %>%
  filter(n>=5) %>%
  pull(species)
q_data = y_data %>%
  separate(species, 
           into=c("Genus","SpecificEpithet"),
           remove=FALSE) %>%
  filter(species %in% q_data_sppSelection)

## Fit simple model
lapsePlusSpp = brms::brm(xmidS ~ elScal + cosAsp + species + 
                           (1 | Position + Year),
                         data = q_data, family = gaussian(),
                         warmup = 2500, iter = 5000, chains = 4,
                         prior = c(prior(normal(100, 100), 
                                         class = Intercept),
                                   prior(cauchy(0, 10), 
                                         class = sd),
                                   prior(cauchy(0, 10), 
                                         class = sigma)),
                         control = control.,
                         seed = mySeed,
                         save_pars = save_pars(all = TRUE))
lapsePlusSpp
saveRDS(lapsePlusSpp, "stan_models/lapsePlusSpp.RDS")
mcmc_plot(lapsePlusSpp,type = "rhat_hist")
pp_check(lapsePlusSpp, draws = 100)

## Fit model with interaction term
lapseIntxSpp = brms::brm(xmidS ~ elScal + cosAsp + 
                           species + elScal*species + 
                           (1 | Position + Year),
                         data = q_data, family = gaussian(),
                         prior = c(prior(normal(100, 100), 
                                         class = Intercept),
                                   prior(cauchy(0, 10),
                                         class = sd),
                                   prior(cauchy(0, 10), 
                                         class = sigma)),
                         warmup = 2500, iter = 5000, chains = 4,
                         control = control.,
                         seed = mySeed,
                         save_pars = save_pars(all = TRUE))
lapseIntxSpp
saveRDS(lapseIntxSpp, "stan_models/lapseIntxSpp.RDS")
mcmc_plot(lapseIntxSpp,type = "rhat_hist")
pp_check(lapseIntxSpp, draws = 100)

## Add cross-validation criteria and compare models
lapsePlusSpp <- add_criterion(lapsePlusSpp, "loo", 
                              reloo = TRUE, moment_match = TRUE)
lapseIntxSpp <- add_criterion(lapseIntxSpp, "loo", 
                              reloo = TRUE, moment_match = TRUE)
loo_compare(lapsePlusSpp, lapseIntxSpp, criterion = "loo")
# elpd_diff se_diff
# lapsePlusSpp  0.0       0.0   
# lapseIntxSpp -3.0       3.8   
bayes_R2(lapsePlusSpp)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.8543174 0.008511865 0.8356235 0.8688061
bayes_R2(lapseIntxSpp)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.8571558 0.008438126 0.8388261 0.8718611

## Generate plots
## Raw elevation trend (2020)
spp_elevTrendPlot = q_data %>%
  filter(Year == 2020) %>%
  group_by(Position, species, Elevation) %>%
  summarize(xmidS=mean(xmidS)) %>%
  ggplot(aes(x = Elevation,
             y = xmidS,
             group = species,
             col = species)) +
  geom_point(col = "black",
             size = 1.5) +
  geom_point(size = 1) +
  stat_smooth(method = "lm",alpha=0.15) + 
  ylab("Green-up timing (day of year)") +
  xlab("Elevation (m)") +
  scale_color_manual("Species",
                     values = CJsBasics::KellyCols[2:20]) +
  CJsBasics::BasicTheme +
  theme(legend.text = element_text(face = "italic"))
spp_elevTrendPlot

## Contrast plot - re-run brms model with unscaled elevation
lapseIntxSppUnsc =  brms::brm(
  xmidS ~ Elevation + cosAsp + species + Elevation*species + 
    (1 | Position + Year),
  data = q_data, family = gaussian(),
  prior = c(prior(normal(100, 100), class = Intercept),
            prior(cauchy(0, 10), class = sd),
            prior(cauchy(0, 10), class = sigma)),
  warmup = 2500, iter = 5000, chains = 4,
  control = control.,
  seed = mySeed,
  save_pars = save_pars(all = TRUE))

newQLabels = unique(q_data$species)[order(unique(q_data$species))][1:length(unique(q_data$species))]
trendPlotDF = plot(emtrends(lapseIntxSppUnsc, ~species, 
                            var="Elevation"),
                   plotit = F)


elevContrast = ggplot(trendPlotDF) + 
  geom_segment(aes(x = lower.HPD*100,
                   xend = upper.HPD*100,
                   y = species,
                   yend = species,
                   col = species),
               lwd = 1.5, alpha = .45) +
  geom_point(aes(x = the.emmean*100,
                 y = species,
                 col = species),
             size = 0.75) +
  scale_y_discrete(labels = rev(newQLabels), limits = rev) +
  scale_color_manual("Species",
                     values = CJsBasics::KellyCols[2:20]) +
  ylab("") +
  xlab("Elevational lapse rate\n(days/100m)") +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 8))

## Combine and save
sppSlopes = ggdraw() +
  draw_plot(spp_elevTrendPlot,
            x = 0, y = 0, 
            width = 1,
            height = 1) +
  draw_plot(elevContrast,
            x = 0.3275, y = 0.125,
            width = 0.3, height = 0.275)
sppSlopes
ggsave("plots/fig5.jpg",
       sppSlopes,
       width = 6, height = 4, units = "in", dpi = 300)  


## > Satellite comparison ----
## Join the cam phenology and satellite phenology data
s_data = y_data %>%
  group_by(Year, Position, 
           elScal, Elevation, cosAsp,
           LandCoverHigherOrder) %>%
  summarize(xmidS = mean(xmidS, na.rm = T)) %>%
  left_join(modPhenology,
            by = c("Position","Year"))
modBoxPrep = y_data %>%
  group_by(Position, Year) %>%
  summarize(cammidS = mean(xmidS, na.rm=T),
            cammidA = mean(xmidA, na.rm=T))  %>%
  pivot_longer(cols = c(cammidS,cammidA), 
               names_to = "Season", values_to = "Timing") %>%
  add_column("Sensor" = "Cameras")

modBoxplot = s_data %>%
  ungroup() %>%
  select(Position, modmidS, modmidA, Year) %>%
  pivot_longer(cols = c(modmidS,modmidA), 
               names_to = "Season", values_to = "Timing") %>%
  add_column("Sensor" = "MODIS") %>%
  rbind(modBoxPrep) %>%
  mutate(Season = case_when(Season == "modmidS" ~ "Green-up",
                            Season == "modmidA" ~ "Senescence",
                            Season == "cammidS" ~ "Green-up",
                            Season == "cammidA" ~ "Senescence")) %>%
  ggplot(aes(x = Year,
             y = Timing)) +
  geom_boxplot(aes(fill = Season, linetype = Sensor)) + 
  # facet_wrap(~Sensor) +
  scale_linetype_manual(values=c(1,5)) +
  scale_color_manual(values = alpha(CJsBasics::KellyCols[c(6,5)],
                                    0.5))  +
  scale_fill_manual(values = alpha(CJsBasics::KellyCols[c(6,5)],
                                   0.5)) +
  ylab("Day of year") +
  CJsBasics::BasicTheme
modBoxplot

# Finding correlations along sensors
lapply(unique(s_data$LandCoverHigherOrder),
       FUN = function(x){
         message(x)
         mydf = s_data %>%
           filter(LandCoverHigherOrder == x)
         message("Pearson's R:")
         cor(mydf$xmidS, mydf$modmidS, method = "pearson")
       })
# [[1]]
# [1] 0.8026863
# 
# [[2]]
# [1] 0.3324534
# 
# [[3]]
# [1] 0.05170887
# 
# [[4]]
# [1] 0.9239132


satComp = s_data %>%
  filter(Year == 2020) %>%
  filter(!is.na(xmidS)) %>%
  ggplot(aes(x = xmidS,
             y = modmidS)) +
  geom_point(col = "black",
             size = 2.5) +
  geom_point(aes(col = Elevation),
             size = 2) +
  geom_abline() +
  facet_wrap(~LandCoverHigherOrder) +
  stat_smooth(method = "lm",
              col = "lightblue4") +
  scale_color_viridis_c("Elevation (m)") +
  coord_equal() +
  xlim(-20,
       255) +
  ylim(-20,
       255) +
  xlab("Green-up timing (Camera)") +
  ylab("Green-up timing (MOD13Q1)") +
  CJsBasics::BasicTheme

modFig_compiled = cowplot::plot_grid(modBoxplot,
                                     satComp,
                                     nrow = 1,
                                     ncol = 2,
                                     rel_widths = c(0.4,0.6),
                                     labels = "AUTO")
modFig_compiled
ggsave("plots/fig6.jpg",
       modFig_compiled,
       width = 8, height = 4, units = "in", dpi = 300)

mod_elevLapseRate = brms::brm(
  modmidS ~ elScal + cosAsp + 
    (1 | Position + Year),
  data = s_data, family = gaussian(),
  warmup = 2500, iter = 5000, chains = 4,
  prior = c(prior(normal(100, 100), 
                  class = Intercept),
            prior(cauchy(0, 10), 
                  class = sd),
            prior(cauchy(0, 10), 
                  class = sigma)),
  control = control.,
  seed = mySeed,
  save_pars = save_pars(all = TRUE))
mod_elevLapseRate
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: modmidS ~ elScal + cosAsp + (1 | Position + Year) 
# Data: s_data (Number of observations: 90) 
# Draws: 4 chains, each with iter = 5000; warmup = 2500; thin = 1;
# total post-warmup draws = 10000
# 
# Group-Level Effects: 
#   ~Position (Number of levels: 76) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    15.57      4.68     2.48    22.25 1.01      535      887
# 
# ~Year (Number of levels: 2) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    29.20     22.10     9.36    85.85 1.00     3873     3534
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept   158.66     23.26   106.80   201.58 1.00     3159     2672
# elScal       32.58      2.36    28.00    37.34 1.00     4782     6089
# cosAsp        8.19      3.73     0.95    15.69 1.00     4956     5778
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    14.62      3.40     9.41    22.29 1.01      511     1092
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
plot(mod_elevLapseRate)
bayes_R2(mod_elevLapseRate)
# Estimate  Est.Error      Q2.5    Q97.5
# R2 0.8826786 0.05402613 0.7495244 0.949949
saveRDS(mod_elevLapseRate,file="stan_models/mod_elevLapseRate.RDS")


#### Purely terrain model, for Figure 7 (fig7; Fig 7; Fig. 7), see 303_aerialCoverage.R.
terrainModel = brms::brm(xmidS ~ 
                           Elevation + 
                           cosAsp +
                           (1 | Year + Position),
                         data = y_data, 
                         family = gaussian(),
                         prior = c(prior(normal(100, 100), 
                                         class = Intercept),
                                   prior(cauchy(0, 10), 
                                         class = sd),
                                   prior(cauchy(0, 10), 
                                         class = sigma)),
                         warmup = 2500, 
                         iter = 5000, 
                         chains = 4,
                         control = control.,
                         seed = mySeed,
                         save_pars = save_pars(all = TRUE))
plot(terrainModel)
bayes_R2(terrainModel)
# Estimate   Est.Error      Q2.5     Q97.5
# R2 0.8449318 0.004646443 0.8351117 0.8532258
saveRDS(terrainModel,file="stan_models/terrainModel.RDS")
