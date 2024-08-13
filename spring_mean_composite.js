//////////////////////////////////////////////////////////////
//       Fire Severity Compositing Code Repository          //
//                                                          //
//  This code calculates Landsat-based fire severity        //
//  spring composite. The pre-fire was set to the           //
//  spring season of the year of the fire and               //
//  the post-fire was set as the same year as the fire.     //
//  The compositing criterion is the mean                   //
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
//  GEE link: https://code.earthengine.google.com/6b92a52bf0d85a4f9d61f9f852e26d2e
//                                                          //
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
  var year = ee.Image(1).multiply(ee.Number.parse(fire.get('year'))).toUint16().rename('fire_year');

  // Pre-fire NBR: May 1 - Jun 1 of the Year Of Fire (YOF)
  var preFireNBR1 = lsCol.filterBounds(study_site)
                         .filterDate(Fireyear, Fireyear.advance(1, 'year'))
                         .filter(ee.Filter.dayOfYear(120, 152))
                         .mean();

  var preNBR = ee.Image(preFireNBR1);

  // Post-fire NBR: Jun 2 - Oct 1 of the Year Of Fire (YOF)
  var postFireNBR1 = lsCol.filterBounds(study_site)
                          .filterDate(Fireyear, Fireyear.advance(1, 'year'))
                          .filter(ee.Filter.dayOfYear(153, 274))
                          .mean();

  var postNBR = ee.Image(postFireNBR1);

  // Calculate fire severity indices:
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
  var severity = rbr.addBands([preBands, postBands, year]);

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
var year_end = 2023;

// Define the function to extract centroids of each image in the severity collection
var getcentroids = function(image) {
  var centroids = image.sample({
    region: study_site,
    geometries: true,
    scale: 30, // Landsat resolution
    tileScale: 16
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

  // Export the first severity image of the current year to Google Drive
  Export.image.toDrive({
    image: severity.first(), // Exporting the first image in the collection
    description: 'Spring_mean_' + year,
    folder: 'YourFolder',
    maxPixels: 1e13,
    scale: 30
  });

  // Uncomment the following code if you want to export the centroid points as a GeoJSON file:
  /*
  Export.table.toDrive({
    collection: centroids.flatten(),
    description: 'Spring_mean_' + year,
    folder: 'YourFolder',
    fileFormat: 'GEO_JSON'
  });
  */
}
