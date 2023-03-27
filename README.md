# Time-lapse Conservation

<img src="/plots/fig1.png" alt="Map (A), example image from time-lapse camera (B), and photo of Mt. Williamson (C)" height="210" width="425">

### Time-lapse camera network reveals intra- and inter-species variation in green-up timing across the range of Sierra Nevada bighorn sheep

Code and data to produce results presented in John et al. 2023, Remote Sensing in Ecology and Conservation. Please cite our article (will update throughout production process):

* John, C., Kerby, J., Stephenson, T.R., and Post, E. _In press._ Fine-scale landscape phenology revealed through time-lapse imagery: implications for conservation and management of an endangered migratory herbivore. Remote Sensing in Ecology and Conservation.

Basic folder structure should match this repository. Using RStudio Project, working directory and relative filepaths should be automatically assigned.

### Workflow

Follow the contents of the `scripts/` folder in numerical order. Descriptions follow:

000_modisFunctions.R: Curve fitting functions for MODIS data
001_cameraFunctions.R: Curve fitting functions for time-lapse GCC and snow cover data
101_GEE_tidy.R: Convert wide-format Earth Engine NDVI output to long-format data.frame
102_modPhenology.R: Derive phenology indices for MODIS data
201_organizeGCC.R: Compile and clean raw greenness data
202_camPhenology.R: Derive phenology indices for time-lapse camera data
301_networkSummary.R: Plot mean camera-level time series (Fig 2)
302_stanModels.R: Fit brms models and plot elevational lapse rates (Fig 3-6)
303_aerialCoverage.R: Measure temporal variation in availability of space at peak-greenup (Fig 7)


### Notes

Be sure to uncompress the exif.zip and extract.zip subdirectories. Note that NDVI (MODIS), elevation (DEP) and imagery (NAIP) were accessed via Google Earth Engine; code for reconstructing these data are available in their corresponding `data/` subdirectory but are required only for maps; wrangled data to reproduce analyses are in `data/` folders