//       Fire Severity Compositing Code Repository
//  This code  pre-process Landsat imagery
// Script to create landsat composites
//  L8 TC coefficients: 
//     -Derivation of a tasselled cap transformation based on Landsat 8 at- satellite reflectance
//  L4-7 TC coefficients: 
//     -Derivation of a Tasseled Cap Transformation Based On Landsat 7 At-Satellite Reflectance
//  Bit Values for masking clouds and cloud shadow
//     - L8-L9 Landsat 8-9 Collection 2 (C2) Level 2 Science Product (L2SP) Guide.
//     - L4-L7 Landsat 4-7 Collection 2 (C2) Level 2 Science Product (L2SP) Guide.

//////////////
  //  Inputs  //
  /////////////
  //------------------- Image Processing  -------------//
  
  var ls9SR = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
  ls8SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
  ls7SR = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
  ls5SR = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
  ls4SR = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2');
  
  ////////////////
  //  Funtions //
  //////////////
    
    // Returns vegetation indices for Landsat OLI sensor
  var ls8_9_Indices = function(lsImage){
    
    var opticalBands = lsImage.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = lsImage.select('ST_B.*').multiply(0.00341802).add(149.0);
    
    //Indices Calculation:
      //B1: ultraBlue, B2: Blue, B3: Green, B4: Red, B5: NIR, B6: SWIR1, B7: SWIR2;
    var lsImageBands = opticalBands.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',  'SR_B6', 'SR_B7']);
    
    //TC coefficients L8 
    
    var brightness_ = ee.Image.constant([0.3029,   0.2786,  0.4733, 0.5599,  0.508,   0.1872]);
    var greenness_ = ee.Image.constant([-0.2941,  -0.243,  -0.5424, 0.7276,  0.0713, -0.1608]);
    var wetness_ = ee.Image.constant([0.1511,   0.1973,  0.3283, 0.3407, -0.7117, -0.4559]);
    var sum = ee.call("Reducer.sum");
    
    var brightness = lsImageBands.multiply(brightness_).reduce(sum);
    var greenness = lsImageBands.multiply(greenness_).reduce(sum);
    var wetness = lsImageBands.multiply(wetness_).reduce(sum);
    
    //NBR: (NIR - SWIR2) / (NIR + SWIR2)
    var nbr = lsImageBands.normalizedDifference(['SR_B5', 'SR_B7']).toFloat();
    //NBR2: (SWIR1 – SWIR2) / (SWIR1 + SWIR2)
    var nbr2 = lsImageBands.normalizedDifference(['SR_B6', 'SR_B7']).toFloat();
    //NDVI: (NIR - R) / (NIR + R)
    var ndvi = lsImageBands.normalizedDifference(['SR_B5', 'SR_B4']).toFloat();
    //NDMI: (NIR – SWIR1) / (NIR + SWIR1)
    var ndmi = lsImageBands.normalizedDifference(['SR_B5', 'SR_B6']).toFloat();
    
    //Metadata bands
    var qa = lsImage.select(['QA_PIXEL']);
    var time = ee.Image(1).multiply(ee.Number(lsImage.get('system:time_start'))).toFloat();
    var doy = ee.Image(1).multiply(ee.Number.parse(ee.Date(lsImage.get('system:time_start')).format('D'))).toUint16();
    var nbr_sort = nbr.multiply(-1);
    var nir = lsImageBands.select('SR_B5');
    var swir1 = lsImageBands.select('SR_B6');
    var swir2 = lsImageBands.select('SR_B7');
    
    return  brightness.addBands([greenness]).addBands([wetness]).addBands([nbr])
    .addBands([nbr2]).addBands([ndvi]).addBands([ndmi])
    .addBands([qa]).addBands([time]).addBands([doy]).addBands([nbr_sort])
    .addBands([nir]).addBands([swir1]).addBands([swir2])
    .select([0,1,2,3,4,5,6,7,8,9,10, 11, 12, 13], 
            ['brightness','greenness','wetness','nbr','nbr2','ndvi',
             'ndmi','pixel_qa', 'time','doy','nbr_sort', 'nir', 'swir1', 'swir2'])
    .copyProperties(lsImage, 
                    ["LANDSAT_PRODUCT_ID", 'system:time_start', 'CLOUD_COVER',
                     "EARTH_SUN_DISTANCE", "DATE_ACQUIRED", "SUN_AZIMUTH", 
                     "SUN_ELEVATION", "SENSOR_ID"]);
  };
  
  // Returns vegetation indices for Landsat TM and ETM+ sensors
  var ls4_7_Indices = function(lsImage){
    
    var opticalBands = lsImage.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = lsImage.select('ST_B6').multiply(0.00341802).add(149.0);
    
    var lsImageBands = opticalBands.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']);
    
    //TC coefficients  
    //B1: Blue, B2: Green, B3: Red, B4: NIR, B5: SWIR1, B7: SWIR2; 
    var brightness_= ee.Image([0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596]);
    var greenness_= ee.Image([-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630]);
    var wetness_= ee.Image([0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388]);
    var sum = ee.call("Reducer.sum");
    
    //Indices Calculation:  
      var brightness = lsImageBands.multiply(brightness_).reduce(sum);
      var greenness = lsImageBands.multiply(greenness_).reduce(sum);
      var wetness = lsImageBands.multiply(wetness_).reduce(sum);
      
      //NBR: (NIR - SWIR2) / (NIR + SWIR2)
      var nbr = lsImageBands.normalizedDifference(['SR_B4', 'SR_B7']).toFloat();
      //NBR2: (SWIR1 – SWIR2) / (SWIR1 + SWIR2)  
      var nbr2 = lsImageBands.normalizedDifference(['SR_B5', 'SR_B7']).toFloat();
      //NDVI: (NIR - R) / (NIR + R)  
      var ndvi = lsImageBands.normalizedDifference(['SR_B4', 'SR_B3']).toFloat();
      //NDMI: (NIR – SWIR1) / (NIR + SWIR1)  
      var ndmi = lsImageBands.normalizedDifference(['SR_B4', 'SR_B5']).toFloat();
      
      //Metadata bands
      var qa = lsImage.select(['QA_PIXEL']);
      var time = ee.Image(1).multiply(ee.Number(lsImage.get('system:time_start'))).toFloat();
      var doy = ee.Image(1).multiply(ee.Number.parse(ee.Date(lsImage.get('system:time_start')).format('D'))).toUint16();
      var nbr_sort = nbr.multiply(-1);
      var nir = lsImageBands.select('SR_B4');
      var swir1 = lsImageBands.select('SR_B5');
      var swir2 = lsImageBands.select('SR_B7');
      
      return  brightness.addBands([greenness]).addBands([wetness]).addBands([nbr])
      .addBands([nbr2]).addBands([ndvi]).addBands([ndmi])
      .addBands([qa]).addBands([time]).addBands([doy]).addBands([nbr_sort])
      .addBands([nir]).addBands([swir1]).addBands([swir2])
      .select([0,1,2,3,4,5,6,7,8,9,10,11, 12, 13], 
              ['brightness','greenness','wetness','nbr','nbr2','ndvi',
               'ndmi','pixel_qa', 'time','doy','nbr_sort', 'nir', 'swir1', 'swir2'])
      .copyProperties(lsImage, 
                      ["LANDSAT_PRODUCT_ID", 'system:time_start', 'CLOUD_COVER',
                       "EARTH_SUN_DISTANCE", "DATE_ACQUIRED", "SUN_AZIMUTH", 
                       "SUN_ELEVATION", "SENSOR_ID"]);
  };
  
  // Mask Landsat surface reflectance images
  // Creates a mask for clear pixels Bit Value: 
    var lsCfmask = function(lsImg){
      var quality =lsImg.select(['pixel_qa']);
      var clear = quality.bitwiseAnd(8).eq(0) // cloud shadow
      .and(quality.bitwiseAnd(32).eq(0) // cloud
           .and(quality.bitwiseAnd(4).eq(0) // water
                .and(quality.bitwiseAnd(16).eq(0)))); // snow
      return lsImg.updateMask(clear).select([0,1,2,3,4,5,6,7,8,9, 10, 11, 12, 13])                                    
      .copyProperties(lsImg, ["LANDSAT_PRODUCT_ID", 'system:time_start', 'CLOUD_COVER',
                              "EARTH_SUN_DISTANCE", "DATE_ACQUIRED", "SUN_AZIMUTH", 
                              "SUN_ELEVATION", "SENSOR_ID"]);  
    };
  
  
  ////////////////
  // Processing  //
  ///////////////
    
//ls9SR = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
// Map functions across Landsat Collections
//.map iterate over all images into imagecollection
  
  var ls9 = ls9SR.map(ls8_9_Indices)
  .map(lsCfmask);
  var ls8 = ls8SR.map(ls8_9_Indices)
  .map(lsCfmask);
  var ls7 = ls7SR.map(ls4_7_Indices)
  .map(lsCfmask); 
  var ls5 = ls5SR.map(ls4_7_Indices)
  .map(lsCfmask); 
  var ls4 = ls4SR.map(ls4_7_Indices)
  .map(lsCfmask); 
  
  
  // Merge Landsat Collections
  var lsCol = ee.ImageCollection(ls9.merge(ls8).merge(ls7).merge(ls5).merge(ls4));
  exports.lsCol = ee.ImageCollection(ls9.merge(ls8).merge(ls7).merge(ls5).merge(ls4));
  