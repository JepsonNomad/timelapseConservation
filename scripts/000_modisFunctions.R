## Functions for following methods described by Bischof et al 2012
## See appendix A for more
#### Data filtering ----
# 1. NDVI pixels with QA == 2 or 3 already set to NA in Google Earth Engine

#### Data smoothing ----
# 2. Lower 0.025 quantile NDVI for each pixel time series is used as winter NDVI
# Here, winterNDVI is a function that operates on a time series and sets all values that fall below the 0.025 quantile to the wNDVI value
# Note that last 4 and first 3 values of the year can be set to the wNDVI value.
winterNDVI = function(x, lpy = 23){
  x[x==-99999] <- NA # Remove -99999 values (GEE Mask flag)
  x[x==0] <- NA # Remove exactly 0 values (Other mask)
  x[x==5000] <- NA # Some pixels stay stuck at 5000 all year; remove these
  wNDVI = quantile(x, probs = 0.025, na.rm = T)
  x[x<wNDVI] = wNDVI
  x[is.na(x)] = wNDVI
  ts_length = length(x)
  nyears = ts_length/lpy
  
  # Set measurements during unambiguous wintertime to wNDVI value
  for(i in 1:nyears){
    x[(i-1)*lpy + 1] = wNDVI # First datum of year
    x[(i-1)*lpy + 2] = wNDVI # Second datum of year
    x[(i-1)*lpy + 3] = wNDVI # Third datum of year
    x[i*lpy-3] = wNDVI # Fourth-last datum of year
    x[i*lpy-2] = wNDVI # Third-last datum of year
    x[i*lpy-1] = wNDVI # Second-last datum of year
    x[i*lpy] = wNDVI # Last datum of year
  }
  
  # Set NA pixels to wNDVI value - median rolling window can't deal with all the NA's that are introduced by poor-quality pixels and poor-quality pixels are associated with non-growing season measurements
  x[is.na(x)] = wNDVI
  
  return(x)
}

# 3a. Use a moving median filter with search window = 3 to smooth the time series. 
library(zoo)
medianwindow = function(x, windowsize = 3L){
  w = zoo::rollmedian(x, k = windowsize, align = "center",
                      na.pad = TRUE)
  return(w)
}

#### Data scaling ----

# 3b. Only interested in relative timing of green-up and senescence, not absolute values, so scale values for each pixel from 0 to 1 based on the 0.925 quantile for each pixel.
library(scales)
ts_scale = function(x, phenex = FALSE){
  if(phenex == TRUE){
    s = x/10000
    # Only look at pixels with at least 2 nonzero observations
  } else if(sum(x!=0,na.rm = T) > 1){ 
    highval = quantile(x, probs = 0.925,na.rm=T)
    s = rescale(x, to=c(0,1), from=c(min(x,na.rm=T), highval))
    s[s>2] = 1 # High-value outliers to facilitate model fit
    s[s< -1] = 0 # Low-value outliers raised to 0 (See Beck et al)
  } else(s = rep(NA, length(x)))
  return(s)
}

#### Curve fitting ----
# This function works on a dataframe (ts.df) with columns 
# "NDVI" = NDVI values
# "DOY" = Day of Year values

# Since there are 23 MODIS observations/year set default to 23 layers per year (lpy)
# x = NDVI time series as vector
# y = DOY time series as vector
# NOTE: yearphenology is only designed to deal with one year of data.
# Also note that because xmidS and xmidA are the inflection points of the double logistic function, they correspond with the timing of peak IRG (greenup) and peak IRS (senescence).
library(phenex)
library(DEoptim)
yearphenology = function(x, yearLayers,
                         lpy = 23){
  # Extract NDVI and day of year values
  x.NDVI = x$NDVI
  y.DOY = x$DOY
  
  if(sum(x.NDVI > 0, na.rm=T) >= 1 && sum(abs(diff(x.NDVI)), na.rm=T) > 0){
    # Only want pixels with at least 1 nonzero observation; pixel must have some unique observations
    # Create data.frame
    ndviDF = data.frame(NDVI = x.NDVI,
                        DOY = y.DOY)
    # DEoptim only works with complete cases
    ndviDF = ndviDF[complete.cases(ndviDF),]
    # The model we want to fit
    dlPredictor = function(xmidS,xmidA,scalS,scalA,DOY){
      NDVI = 
        ((1/(1+exp((xmidS-DOY)/scalS))) -
           (1/(1+exp((xmidA-DOY)/scalA))))
      return(NDVI)
    }
    # Calculate difference from fitted object and observed values (dlObs)
    dlDelta = function(x, DOY, dlObs){
      NDVI = sum((dlPredictor(xmidS=x[1], xmidA=x[2], 
                              scalS=x[3], scalA=x[4], DOY)-dlObs)^2)
      return(ifelse((is.infinite(NDVI)||is.nan(NDVI)),1e50,NDVI))
    }
    # Optimize for sum of squared errors
    dlModel = DEoptim(fn=dlDelta, DOY=ndviDF$DOY, dlObs=ndviDF$NDVI,
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
    ndviPheno = data.frame("xmidS" = xmidS, "scalS" = scalS, 
                  "xmidA" = xmidA, "scalA" = scalA)
    
    # Return two-layer array with melt season and snowfall onset
    return(ndviPheno)
  } else{return(data.frame("xmidS" = NA, "scalS" = NA, "xmidA" = NA, "scalA" = NA))}
}
