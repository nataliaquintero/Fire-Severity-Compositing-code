//////////////////////////////////////////////////////////////
//       Fire Severity Compositing Code Repository          //
//                                                          //
//  This code calculates Landsat-based fire severity        //
//  extended composite. The pre-fire was set to one year    //
//  before the fire and the post-fire was set one year after//
//  the fire. The compositing criterion is the medoid.      //
//                                                          //
//  This code uses a Landsat composite function from the    //
//  external module `LandsatComposites`:                    //
//                                                          //
// var LandsatComposites =                                  //
// require('users/estro/Busqueda_imagenes:LandsatComposites')//
//                                                          //
//  Change the fire_perimeters and the study_area as needed.//
//  Fire perimeters should have a column named 'Year',      //
//  indicating the year of the fire. Dissolve your fire     //
//  perimeters based on the 'Year' column.                  //
//                                                          //
//  The time lapse to create the composite is set during    //
//  summer. Change it if your fire perimeters are from      //
//  another season of the year.                             //
//  GEE link: https://code.earthengine.google.com/fdf38e837b3c3d34303f987da4f45d42
//                                                          //
//  License: MIT License                                    //
//  Authors: Natalia Quintero, Olga Viedma,                 //
//           Sander Veraverbeke,                            //
//           José Manuel Moreno                             //
//  Affiliation: Universidad de Castilla La Mancha,         //
//               Vrije Universiteit Amsterdam,              //
//               University of East Anglia.                 //
//                                                          //
//  For more information, see the associated research paper://
//  "Optimising Regional Fire Severity Mapping using        //
//  Pixel-Based Image Compositing."                         //
// Database containing data generated using this code is    //
// accessible at:                                           //
//     Mendeley Data, V1, doi: 10.17632/dxp7p66gv3.1        //
//                                                          //
//////////////////////////////////////////////////////////////

// Convert the "year" property of each feature in the fires_summer collection to an integer format.
// This ensures that the "year" property is correctly handled as an integer for filtering and calculations.
var fires_summer = ee.FeatureCollection('users/estro/Spain_fires_summer_100'); 

fires_summer = fires_summer.map(function(feature) {
  var yearInt = ee.Number(feature.get('year')).toInt();
  return feature.set('year', yearInt);
});

///////////////////////////////////////////////////////////////
//  Function to calculate fire severity and related variables //
///////////////////////////////////////////////////////////////
var get_severity = function(fire) {
  // Parse the year from the fire feature properties and create a date object
  var Fireyear = ee.Date.parse('YYYY', fire.get('year'));

  // Create an image representing the year of the fire, useful for tracking the fire year in the output
  var year = ee.Image(1).multiply(ee.Number.parse(fire.get('year'))).rename('fire_year');

  // Function to rename bands and select them from the image
  var bands = function(image) {
    var bandNames = PreMedoid.first().bandNames(); // Retrieve band names from the pre-fire medoid
    var bandPositions = ee.List.sequence(1, bandNames.length().subtract(1)); // Get positions of the bands
    return image.reduce(ee.Reducer.min(bandNames.length()))
                 .select(bandPositions, bandNames.slice(1));
  };

  // Filtering Landsat collection for the pre-fire period (one year before the fire)
  var preFireNBR1_col = lsCol.filterBounds(study_site)
                              .filterDate(Fireyear.advance(-1, 'year'), Fireyear)
                              .filter(ee.Filter.dayOfYear(153, 274));

  // Calculating the medoid for the pre-fire period
  var PreMedoid = preFireNBR1_col.map(function(image) {
    var diff = ee.Image(image).subtract(preFireNBR1_col.median()).pow(ee.Image.constant(2)); 
    return diff.reduce('sum').addBands(image);  
  });

  // Get the pre-fire NBR by selecting the band positions and renaming the bands
  var preNBR = bands(PreMedoid);                      

  // Filtering Landsat collection for the post-fire period (one year after the fire)
  var postFireNBR1_col = lsCol.filterBounds(study_site)
                          .filterDate(Fireyear.advance(1, 'year'), Fireyear.advance(2, 'year'))
                          .filter(ee.Filter.dayOfYear(153, 274));

  // Calculating the medoid for the post-fire period
  var PostMedoid = postFireNBR1_col.map(function(image) {
    var diff = ee.Image(image).subtract(postFireNBR1_col.median()).pow(ee.Image.constant(2)); 
    return diff.reduce('sum').addBands(image);  
  });

  // Get the post-fire NBR by selecting the band positions and renaming the bands
  var postNBR = bands(PostMedoid);  
  
  // Calculate severity indices:
  // dNBR: Differenced Normalized Burn Ratio
  var dNBR = preNBR.subtract(postNBR).multiply(1000).select(['nbr'], ['dNBR']).toFloat();
  
  // RBR: Relativized Burn Ratio
  var rbr = dNBR.divide(preNBR.add(1.001)).select(['nbr'], ['rbr']).toFloat();
  
  // RdNBR: Relative differenced Normalized Burn Ratio
  var rdnbr = dNBR.divide(preNBR.divide(1000).abs().sqrt()).select(['nbr'], ['rdnbr']).toFloat();

  // Select and rename pre-fire and post-fire bands
  var preBands = preNBR.select('nbr', 'ndvi', 'nir', 'swir1', 'swir2', 'time', 'doy')
      .rename('nbr_pre', 'ndvi_pre', 'nir_pre', 'swir1_pre', 'swir2_pre', 'time_pre', 'doy_pre');
  var postBands = postNBR.select('nbr', 'ndvi', 'nir', 'swir1', 'swir2', 'time', 'doy')
      .rename('nbr_post', 'ndvi_post', 'nir_post', 'swir1_post', 'swir2_post', 'time_post', 'doy_post');

  // Combine all severity-related bands into a single image
  var severity = rbr.addBands([preBands]).addBands([postBands]).addBands([year]);

  // Clip the severity image to the fire perimeter and return it
  return ee.Image(severity).clip(fire).toFloat().set('id', fire.get('year'));
};

//------------------- INPUTS  -------------//

// Define the study area (Spain) using the FAO GAUL dataset
var study_site = ee.FeatureCollection('FAO/GAUL/2015/level0')
      .filter(ee.Filter.eq('ADM0_NAME', 'Spain')).geometry(); 

// Import the Landsat composites module
var LandsatComposites = require('users/estro/Busqueda_imagenes:LandsatComposites');
var lsCol = LandsatComposites.lsCol;

// Define the start and end years for the analysis
var year_start = 2000;
var year_end = 2022;

// Define the function to extract centroids of each image in the severity collection
var getcentroids = function(image) {
  var centroids = image.sample({
    region: study_site,
    geometries: true,
    scale: 30, // Landsat resolution
    tileScale: 8
  });
  return centroids.set("id", image.get("id"));
};

//------------------- Apply severity function and EXPORT images for each year  -------------//
for (var year = year_start; year <= year_end; year++) {
  // Filter the fire data for the current year
  var fires_year = fires_summer.filter(ee.Filter.eq('year', year));

  // Apply the severity calculation function to each fire in the current year
  var severity = ee.ImageCollection(fires_year.map(get_severity));

  // Create centroids for the severity ImageCollection
  var centroids = severity.map(getcentroids);

  // Uncomment the following code if you want to export the centroid points as a GeoJSON file:
  /*
  Export.table.toDrive({
    collection: centroids.flatten(),
    description: 'extended_medoid_' + year,
    folder: 'YourFolder',
    fileFormat: 'GEO_JSON'
  });
  */

  // Export the first severity image of the current year to Google Drive
  Export.image.toDrive({
    image: severity.first(), // Exporting the first image in the collection
    description: 'extended_medoid_' + year,
    folder: 'YourFolder',
    maxPixels: 1e13, 
    scale: 30
  });
}
