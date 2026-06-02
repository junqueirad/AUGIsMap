// sort stratified spatialPoints to use as training samples

// define string to use as metadata
var version = 2;

// define ibges' carta id
var file_name = 'depressaoPeriferica';

// input metadata
var version_output = 2;

// output directory
var output_dir = 'users/deisejunqueira/AUGI-MAP';

// define sample size
var sampleSize = 7000;     // number of samples to be sorted 
var nSamplesMin = 700;     // minimum sample size by class
var nSamplesAUI = 700;     // number of samples of AUI 

// training samples
var gtPolys = ee.FeatureCollection("projects/ee-deisejunqueira/assets/ManuscriptPolygons_metrics")
  .map(function(feature) {
    return feature.set('reference', 9);
  });

// rasterize
var AUI = ee.Image().paint(gtPolys, 'reference').rename('reference'); 

// read training mask
var trainingMask = ee.Image('projects/ee-ipam-cerrado/assets/Collection_11/masks/cerrado_trainingMask_1985_2024_v4')
  .rename('reference')
  .blend(AUI);  // Overlay the AUGIs on the stable samples

// read areas reference 
var referenceAreas = ee.FeatureCollection('projects/ee-deisejunqueira/assets/AUGI-MAP/depressaoPeriferica_area_v2');

// read study area
var carta = ee.FeatureCollection('projects/ee-deisejunqueira/assets/DepressaoPeriferica');

// import mapbiomas module
var vis = {
    'min': 0,
    'max': 62,
    'palette': require('users/mapbiomas/modules:Palettes.js').get('classification8')
};

// plot stable pixels
Map.addLayer(trainingMask, vis, 'trainingMask', false);

// define function to get trainng samples
var getTrainingSamples = function (feature) {
  // read the area for each class
  var forest = ee.Number(feature.get('3'));
  var wetland = ee.Number(feature.get('11'));
  var pasture = ee.Number(feature.get('15'));
  var agriculture = ee.Number(feature.get('18'));
  var non_vegetated = ee.Number(feature.get('25'));
  var water = ee.Number(feature.get('33'));

  // compute the total area 
  var total = forest.add(wetland).add(pasture).add(agriculture)
                  .add(non_vegetated).add(water);
              
  // define the equation to compute the n of samples per class
  var computeSize = function (number) {
    return number.divide(total).multiply(sampleSize).round().int16().max(nSamplesMin);
  };
  
  // apply the equation to compute the number of samples
  var n_forest = computeSize(ee.Number(forest));
  var n_wetland = computeSize(ee.Number(wetland));
  var n_pasture = computeSize(ee.Number(pasture));
  var n_agriculture = computeSize(ee.Number(agriculture));
  var n_non_vegetated = computeSize(ee.Number(non_vegetated));
  var n_water = computeSize(ee.Number(water));
  var n_AUI = ee.Number(nSamplesAUI);

  // get the geometry of the region
  var region_i_geometry = ee.Feature(feature).geometry();
  // clip stablePixels only to the region 
  var referenceMap =  trainingMask.clip(region_i_geometry);
                      
  // generate the sample points
  var training = referenceMap.stratifiedSample(
                                {'scale': 10,
                                 'classBand': 'reference', 
                                 'numPoints': 0,
                                 'region': feature.geometry(),
                                 'seed': 1,
                                 'geometries': true,
                                 'classValues': [3, 9, 11, 15, 18, 25, 33],
                                 'classPoints': [n_forest, n_AUI, n_wetland, n_pasture, n_agriculture, n_non_vegetated, n_water]
                                  }
                                );
  return training;
 };

// apply function and get sample points
var samplePoints = referenceAreas.map(getTrainingSamples).flatten(); 

// apply style over the points
var paletteMapBiomas = require('users/mapbiomas/modules:Palettes.js').get('classification8');
var newSamplesStyled = samplePoints.map(
    function (feature) {
        return feature.set('style', {
            'color': ee.List(paletteMapBiomas)
                .get(feature.get('reference')),
            'width': 1,
        });
    }
).style(
    {
        'styleProperty': 'style'
    }
);

// plot points
Map.addLayer(newSamplesStyled, {}, 'samplePoints');

// export as GEE asset
Export.table.toAsset({'collection': samplePoints,
                      'description': file_name + '_samplePoints_v' + version,
                      'assetId':  output_dir + '/' + file_name + '_samplePoints_v' + version});
