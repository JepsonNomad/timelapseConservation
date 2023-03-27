Data downloaded from Earth Engine 4 Jan 2023

where huBuff is the 10km buffered region of SNBS herd units
and NED is
This is the seamless 3DEP DEM dataset for the U.S. with full coverage of the 48 conterminous states, Hawaii, and U.S. territories. Alaska coverage is partially available now and is being expanded to statewide coverage as part of the Alaska Mapping Initiative. Ground spacing is approximately 10 meters north/south, but variable east/west due to convergence of meridians with latitude.

Spatial metadata dataset is ingested as a separate asset USGS_3DEP_10m_metadata.

The 1m dataset is ingested as USGS_3DEP_1m.

Dataset uploaded by Farmers Business Network.

U.S. Geological Survey, 3D Elevation Program 10-Meter Resolution Digital Elevation Model.


//var ROI = ee.FeatureCollection('EPA/Ecoregions/2013/L4')
//  .filter(ee.Filter.eq("us_l3name", "Sierra Nevada"));
var ROI = ee.FeatureCollection("TIGER/2016/States")
.filter(ee.Filter.eq("NAME", "California"));
// Create clipping function so we're just dealing 
// with the ROI
var clipper = function(image){
  return image.clip(huBuff);
};

var DEM = clipper(DEP);
print("Elevation", DEM);

Map.addLayer(DEM, {min:-500,max:4000}, "Elevation");
Map.addLayer(herdunits, {color: '00FFFF10', opacity: 0}, "Herd units");
Map.addLayer(cams, {color: 'ff000099'}, "Camera positions");

// The following is based loosely on:
// https://gis.stackexchange.com/q/265392/67264
var fill = function(img, ini) {
  // gets the values for the points in the current img
  var ft2 = img.reduceRegions(ini, ee.Reducer.first(),30).set("name", "Elevation")
  // merges the FeatureCollections
  return ft2
};

// Iterates over the ImageCollection
var newft = ee.FeatureCollection(fill(DEP, cams));
print("new feature", newft);

var myCrs = "EPSG:32611";

// Export
Export.image.toDrive({
  image: DEM,
  fileNamePrefix: 'DEP_30m',
  description: "DEP_30m",
  folder: 'timelapseConservation',
  scale: 30,
  crs: myCrs,
  region: huBuff.geometry().bounds(),
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});
