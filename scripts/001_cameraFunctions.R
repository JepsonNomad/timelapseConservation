## Functions for smoothing and fitting greening curve to phenocam GCC data
## Based on Bischof et al 2012 approach to NDVI time series

#### Data smoothing ----
# 2. Lower 0.1 quantile GCC for each pixel time series is used as winter GCC
# Here, winterGCC is a function that operates on a time series and sets all values that fall below the 0.1 quantile to the wGCC value
# Note that last 4 and first 3 values of the year can be set to the wGCC value.
winterGCC = function(x){
  wGCC = quantile(x, probs = 0.1, na.rm = T)
  x[x<wGCC] = wGCC
  ts_length = length(x)
  # Set NA pixels to wGCC value - median rolling window can't deal with all the NA's that are introduced by poor-quality pixels and poor-quality pixels are associated with non-growing season measurements
  # x[is.na(x)] = wGCC
  return(x)
}

# 3a. Use a moving median filter with search window = 7 images to smooth the time series.
library(zoo)
medianwindow = function(x, windowsize = 7L){
  w = zoo::rollmedian(x, k = windowsize, align = "center",
                      na.pad = TRUE)
  return(w)
}
quantilewindow = function(x, windowsize = 7L){
  w = zoo::rollapply(x, width=windowsize, quantile, 
                     probs=0.9, fill=NA,
                     na.rm = T)
  return(w)
}

#### Data scaling ----

# 3b. Only interested in relative timing of green-up and senescence, not absolute values, so scale values for each pixel from 0 to 1 based on the 0.925 quantile for each pixel.
library(scales)
ts_scale = function(x){
  if(sum(x!=0, na.rm = T) > 1){ 
    highval = quantile(x, probs = 0.925,na.rm=T)
    s = rescale(x, to=c(0,1), from=c(min(x,na.rm=T), highval))
    s[s > 2] = 1 # High-value outliers to facilitate curve fit
    s[s < -1] = 0 # Low-value outliers raised to 0 (See Beck et al)
  } else(s = rep(NA, length(x)))
  return(s)
}


#### Phenology derivation ----
# x = GCC time series as vector
# y = DOY time series as vector
# NOTE: yearphenology is only designed to deal with one year of data.
# Also note that because xmidS and xmidA are the inflection points of the double logistic function, they correspond with the timing of peak IRG (greenup) and peak IRS (senescence).
library(phenex)
library(DEoptim)
gccPheno = function(x,
                    lpy = 365){
  # Extract GCC and day of year values
  x.GCC = x$GCC[!is.na(x$GCC)]
  y.DOY = x$DOY[!is.na(x$GCC)]
  yDiffs = diff(y.DOY)
  
  if(sum(x.GCC > 0, na.rm=T) >= 1 &&
     sum(abs(diff(x.GCC)), na.rm=T) > 0 && 
     sum(!(c(1:365) %in% y.DOY)) < 90 &&
     all(yDiffs[y.DOY > 59 & y.DOY < 274] < 21, na.rm = TRUE)){
    # Remove points where the value is only ever 0;
    # Remove points where the GCC value never changes
    # Only analyze data that have less than 3 months missing data on the year
    # Remove points with large data gaps (i.e. 3 weeks or more) during the growing season
    # Create data.frame
    GCCDF = data.frame(GCC = x.GCC,
                        DOY = y.DOY)
    # DEoptim only works with complete cases
    GCCDF = GCCDF[complete.cases(GCCDF),]
    # The model we want to fit
    dlPredictor = function(xmidS,xmidA,scalS,scalA,DOY){
      GCC = 
        ((1/(1+exp((xmidS-DOY)/scalS))) -
           (1/(1+exp((xmidA-DOY)/scalA))))
      return(GCC)
    }
    # Calculate difference from fitted object and observed values (dlObs)
    dlDelta = function(x, DOY, dlObs){
      GCC = sum((dlPredictor(xmidS=x[1], xmidA=x[2], 
                              scalS=x[3], scalA=x[4], DOY)-dlObs)^2)
      return(ifelse((is.infinite(GCC)||is.nan(GCC)),1e50,GCC))
    }
    # Optimize for sum of squared errors
    dlModel = DEoptim(fn=dlDelta, DOY=GCCDF$DOY, dlObs=GCCDF$GCC,
                      lower = c(1,1,0,0),
                      upper = c(366,366,15,15), 
                      control=list(VTR=0, strategy=1, 
                                   NP=200, itermax = 200, 
                                   trace=FALSE, CR=0.9))
    
    # Access the midpoints of the spring and autumn seasons
    xmidS = dlModel$optim$bestmem[1]
    xmidA = dlModel$optim$bestmem[2]
    scalS = dlModel$optim$bestmem[3]
    scalA = dlModel$optim$bestmem[4]
    names(xmidS) = NULL
    names(xmidA) = NULL
    names(scalS) = NULL
    names(scalA) = NULL
    
    # Find midpoint of spring and autumn seasons
    GCCPheno = data.frame("xmidS" = xmidS, "scalS" = scalS, 
                           "xmidA" = xmidA, "scalA" = scalA)
    
    # Return two-layer array with melt season and snowfall onset
    return(GCCPheno)
  } else{return(data.frame("xmidS" = NA, "scalS" = NA, "xmidA" = NA, "scalA" = NA))}
}


## Create a function to calculate timing of onset and end of snow season. Expects a data.frame() containing 1 snowYear of data Snow data. Should contain a Date column with a formatted Date, and a Snow column with values between 0 and 1 inclusive. Only operates on cams that have n(unique days) > 300 and n(snowy days) > 10. SnowPheno1 returns a complete set of information, snowPheno2 supplies only a summary
snowPheno1 = function(x){
  myCam = x
  if(length(unique(myCam$Date)) < 200 |
     sum(myCam$Snow, na.rm = T) < 5){
    myCam$snowStart = NA
    myCam$snowStop = NA
    # Generate phenology curve from DEoptim function
    myCam$snowPred = NA
    myCam$obsDay = NA
  } else{
    startDate = min(myCam$Date)
    # Calculate time elapsed since first observation
    myCam$obsDay = as.numeric(myCam$Date - startDate)
    ## Create data.frame() to work through deOptim with
    snowDF = myCam %>%
      select(Snow, obsDay)
    # The model we want to fit - based on Beck 2006
    snowPredictor = function(xmidS, xmidA, scalS, scalA, DOY){
      SNOW = 
        ((1/(1+exp((xmidS-DOY)/scalS))) -
           (1/(1+exp((xmidA-DOY)/scalA))))
      return(SNOW)
    }
    # Calculate difference from fitted object and observed values (snowObs)
    snowDelta = function(x, DOY, snowObs){
      SNOW = sum((snowPredictor(xmidS=x[1], xmidA=x[2], 
                                scalS=x[3], scalA=x[4], DOY)-snowObs)^2)
      return(ifelse((is.infinite(SNOW)||is.nan(SNOW)),1e50,SNOW))
    }
    ## Optimize for sum of squared errors
    snowModel = DEoptim(fn=snowDelta, DOY=snowDF$obsDay, snowObs=snowDF$Snow,
                        lower = c(1,1,0,0),
                        upper = c(366,366,15,15), 
                        control=list(VTR=0, strategy=1, 
                                     NP=200, itermax = 200, 
                                     trace=FALSE, CR=0.9))
    # Access the midpoints of the spring and autumn seasons
    xmidS <- snowModel$optim$bestmem[1]
    xmidA <- snowModel$optim$bestmem[2]
    scalS <- snowModel$optim$bestmem[3]
    scalA <- snowModel$optim$bestmem[4]
    names(xmidS) = NULL
    names(xmidA) = NULL
    names(scalS) = NULL
    names(scalA) = NULL
    # Identify date of snow onset and end
    snowySeasonStart = xmidS + startDate + 1
    snowySeasonStop = xmidA + startDate + 1
    # Append dataframe with snow season info
    myCam$snowStart = snowySeasonStart
    myCam$snowStop = snowySeasonStop
    # Generate phenology curve from DEoptim function
    myCam$snowPred = ((1/(1+exp((xmidS-snowDF$obsDay)/scalS))) -
                        (1/(1+exp((xmidA-snowDF$obsDay)/scalA))))
  }
  return(myCam)
}


snowPheno2 = function(x){
  myCam = x
  if(length(unique(myCam$Date)) < 200 |
     sum(myCam$Snow, na.rm = T) < 5){
    return(data.frame("xmidS" = NA, "scalS" = NA, "xmidA" = NA, "scalA" = NA))
  } else{
    startDate = min(myCam$Date)
    # Calculate time elapsed since first observation
    myCam$obsDay = as.numeric(myCam$Date - startDate)
    ## Create data.frame() to work through deOptim with
    snowDF = myCam %>%
      select(Snow, obsDay)
    # The model we want to fit - based on Beck 2006
    snowPredictor = function(xmidS, xmidA, scalS, scalA, DOY){
      SNOW = 
        ((1/(1+exp((xmidS-DOY)/scalS))) -
           (1/(1+exp((xmidA-DOY)/scalA))))
      return(SNOW)
    }
    # Calculate difference from fitted object and observed values (snowObs)
    snowDelta = function(x, DOY, snowObs){
      SNOW = sum((snowPredictor(xmidS=x[1], xmidA=x[2], 
                                scalS=x[3], scalA=x[4], DOY)-snowObs)^2)
      return(ifelse((is.infinite(SNOW)||is.nan(SNOW)),1e50,SNOW))
    }
    ## Optimize for sum of squared errors
    snowModel = DEoptim(fn=snowDelta, DOY=snowDF$obsDay, snowObs=snowDF$Snow,
                        lower = c(1,1,0,0),
                        upper = c(366,366,15,15), 
                        control=list(VTR=0, strategy=1, 
                                     NP=200, itermax = 200, 
                                     trace=FALSE, CR=0.9))
    # Access the midpoints of the spring and autumn seasons
    xmidS <- snowModel$optim$bestmem[1]
    xmidA <- snowModel$optim$bestmem[2]
    scalS <- snowModel$optim$bestmem[3]
    scalA <- snowModel$optim$bestmem[4]
    names(xmidS) = NULL
    names(xmidA) = NULL
    names(scalS) = NULL
    names(scalA) = NULL
    # Identify date of snow onset and end
    snowySeasonStart = xmidS + startDate + 1
    snowySeasonStop = xmidA + startDate + 1
  }
  # Find midpoint of spring and autumn seasons
  snowPheno = data.frame("xmidS" = snowySeasonStart, "scalS" = scalS, 
                         "xmidA" = snowySeasonStop, "scalA" = scalA)
  # Return two-layer array with melt season and snowfall onset
  return(snowPheno)
}
