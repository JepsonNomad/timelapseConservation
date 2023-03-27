NAIP imagery access available via Google Earth Engine. Code below:
Contact CDFW for Sierra bighorn ROI layer, alternatively define your own.


// Create clipping function so we're just dealing 
// with the ROI
var clipper = function(image){
  return image.clip(ROI);
};

var naip = ee.ImageCollection('USDA/NAIP/DOQQ')
  .filterDate("2016-01-01", "2020-07-01")
  .map(clipper)
  .reduce(ee.Reducer.median());
print("Raw image collection", naip);

Map.setCenter(-119,37,8);
var naipVis = {bands: ['R_median', 'G_median', 'B_median'], min: 0, max: 255};
var nV = naip.visualize(naipVis);
Map.addLayer(naip, naipVis, "NAIP");
Map.addLayer(cams, {color: 'ff000099'}, "Camera positions");

// Export
// Export the image, specifying scale and region.
// Scale = number of meters on each side of pixel
// be sure to set maxPixels to facilitate full export
Export.image.toDrive({
  image: nV.clip(ROI),
  description: 'NAIP',
  folder: 'SNBS_HU_NAIP',
  scale: 5,
  region: ROI.geometry().bounds(),
  maxPixels: 1e13});
