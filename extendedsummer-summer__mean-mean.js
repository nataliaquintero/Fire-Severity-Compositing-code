//////////////////////////////////////////////////////////////
//       Fire Severity Compositing Code Repository          //
//                                                          //
//  This code calculates Landsat-based fire severity        //
//  extended composite. The pre-fire was set to one year    //
//  before the fire and the post-fire was set one year after//
//  the fire. The compositing criterion is the mean.      //
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

var get_severity = function(fire){
  var Fireyear = ee.Date.parse('YYYY', fire.get('year'));
  var year =  ee.Image(1).multiply(ee.Number.parse(fire.get('year'))).rename('fire_year');

  //Severity indices
  var preFireNBR1 =  lsCol.filterBounds(study_site)
                          .filterDate(Fireyear.advance(-1, 'year'), Fireyear)
                          .filter(ee.Filter.dayOfYear(153, 274))
                          .mean();

  var preNBR = ee.Image(preFireNBR1);    

  var postFireNBR1 = lsCol.filterBounds(study_site)
                          .filterDate(Fireyear.advance(1, 'year'), Fireyear.advance(2, 'year'))
                          .filter(ee.Filter.dayOfYear(153, 274))
                          .mean();
                          
  var postNBR = ee.Image(postFireNBR1);
  
  //Severity indices
  var dNBR = preNBR.subtract(postNBR).multiply(1000).select(['nbr'], ['dNBR']).toFloat();
  var rbr = dNBR.divide(preNBR.add(1.001)).select(['nbr'], ['rbr']).toFloat();
  var rdnbr = dNBR.divide(preNBR.divide(1000).abs().sqrt()).select(['nbr'], ['rdnbr']).toFloat();
  
  var preBands = preNBR.select('nbr','ndvi', 'nir','swir1', 'swir2' ,'time', 'doy')
      .rename('nbr_pre', 'ndvi_pre', 'nir_pre','swir1_pre', 'swir2_pre', 'time_pre', 'doy_pre');
  var postBands = postNBR.select('nbr','ndvi', 'nir','swir1', 'swir2','time', 'doy')
      .rename('nbr_post', 'ndvi_post', 'nir_post','swir1_post', 'swir2_post', 'time_post', 'doy_post');
      
  var severity = rbr.addBands([preBands]).addBands([postBands]).addBands([year]);
                     
  return ee.Image(severity).clip(fire).toFloat().set('id', fire.get('year'));
};

//STUDY_SITE:
var study_site = ee.FeatureCollection('FAO/GAUL/2015/level0') 
      .filter(ee.Filter.eq('ADM0_NAME', 'Spain')).geometry(); 

//Landsat
var LandsatComposites = require('users/estro/Busqueda_imagenes:LandsatComposites');
var lsCol = LandsatComposites.lsCol;

//Fires
var year_start = 2000;
var year_end = 2022;

// Define the getcentroids function outside the loop
var getcentroids = function (image){
  var centroids = image.sample({
    region: study_site,
    geometries: true,
    scale: 30, //Landsat
    tileScale: 8
  });
  return centroids.set("id", image.get("id"));
};

//------------------- Apply severity function and EXPORT points for each year  -------------//
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
    description: 'Extended_mean_' + year,
    folder: 'YourFolder',
    maxPixels: 1e13,
    scale: 30
  });

  // Uncomment the following code if you want to export the centroid points as a GeoJSON file:
  /*
  Export.table.toDrive({
    collection: centroids.flatten(),
    description: 'Extended_mean_' + year,
    folder: 'YourFolder',
    fileFormat: 'GEO_JSON'
  });
  */
}
