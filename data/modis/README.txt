README

Data downloaded from Google Earth Engine 23 Sept 2021 using SNBS_timelapseCams/MOD_NDVI
Note that .geo column has been obfuscated to mask precise camera locations in output tables.


CODE (for use in Earth Engine code editor):


// Based on Rodrigo Principe's work:
// https://gis.stackexchange.com/questions/265392/extracting-pixel-values-by-points-and-converting-to-table-in-google-earth-engine

// Set ROI
var ROI = ee.FeatureCollection(timelapseCams);
Map.centerObject(timelapseCams);

// Create a QA mask + clipping function
var masker = function(image){ 
  var mask = image.select('SummaryQA').lte(1);
  var maskedImage = image.updateMask(mask);
  return maskedImage.unmask(-99999).clip(ROI);
  
};

// Compile the data
var dataset = ee.ImageCollection('MODIS/006/MOD13Q1')
  .filterBounds(timelapseCams.geometry())
  .filter(ee.Filter.date('2018-01-01', '2021-12-31'))
  .map(masker);
print("MOD13Q1 Image Collection", dataset);


// Visualize ROI and imagery
Map.addLayer(dataset.reduce(ee.Reducer.median()), {bands: ['NDVI_median'], min: 0, max: 100}, "MOD13Q1");
Map.addLayer(ROI, {color: 'FF0000'}, "timelapseCams");

// function to map over the FeatureCollection
// For getting NDVI
var mapfunc1 = function(feat) {
  // Find feature name
  var id = ee.String(feat.get('m[,\"N\"]'));

  // get feature geometry
  var geom = feat.geometry();
  // function to iterate over the ImageCollection
  // the initial object for the iteration is the feature
  var addProp = function(img, f) {
    // cast Feature
    var newf = ee.Feature(f);
    // get date as string
    var date = img.date().format();
    // extract the NDVI value
    var NDVI = img.select('NDVI').rename("ndvi");
    var value = NDVI.reduceRegion(ee.Reducer.first(), geom, 30).get("ndvi");

    // Return values as a property of the feature. 
    // Return vals if present, "No data" if value = null
    // The name of the property will be the date.
    return ee.Feature(ee.Algorithms.If(value,
                                       newf.set(date, ee.String(value)),
                                       newf.set(date, ee.String('No data'))));
  };
  var newfeat = ee.Feature(dataset.iterate(addProp, feat), {name: id});
  return newfeat;
};

// For getting the composite day of year
var mapfunc2 = function(feat) {
  // Find feature name
  var id = ee.String(feat.get('m[,\"N\"]'));

  // get feature geometry
  var geom = feat.geometry();
  // function to iterate over the ImageCollection
  // the initial object for the iteration is the feature
  var addProp = function(img, f) {
    // cast Feature
    var newf = ee.Feature(f);
    // get date as string
    var date = img.date().format();
    // extract the NDVI value
    var DayOfYear = img.select('DayOfYear').rename("DayOfYear");
    var value = DayOfYear.reduceRegion(ee.Reducer.first(), geom, 30).get("DayOfYear");

    // Return values as a property of the feature. 
    // Return vals if present, "No data" if value = null
    // The name of the property will be the date.
    return ee.Feature(ee.Algorithms.If(value,
                                       newf.set(date, ee.String(value)),
                                       newf.set(date, ee.String('No data'))));
  };
  var newfeat = ee.Feature(dataset.iterate(addProp, feat), {name: id});
  return newfeat;
};

var newft1 = ROI.map(mapfunc1);
print("NDVI output", newft1);

var newft2 = ROI.map(mapfunc2);
print("DayOfYear output", newft2);

// Export
Export.table.toDrive({
  collection: newft1, 
  description: "timelapseCams_MOD13Q1_NDVI", 
  fileFormat: "CSV",
  folder: "timelapseConservation"});
Export.table.toDrive({
  collection: newft2, 
  description: "timelapseCams_MOD13Q1_cDOY", 
  fileFormat: "CSV",
  folder: "timelapseConservation"});
