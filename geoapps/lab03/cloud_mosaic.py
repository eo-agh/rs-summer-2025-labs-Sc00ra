import ee
import geemap

def create_cloud_free_mosaic(aoi, start_date, end_date):
    def ESAcloudMask(img):
        qa = img.select('QA60')
        cloud_bit = int(2 ** 10)
        cirrus_bit = int(2 ** 11)
        clear = qa.bitwiseAnd(cloud_bit).eq(0).And(qa.bitwiseAnd(cirrus_bit).eq(0))
        cloud = clear.Not().rename(['ESA_clouds'])
        return img.addBands(cloud)
    ic = ee.ImageCollection('COPERNICUS/S2_HARMONIZED') \
        .filterBounds(aoi) \
        .filterDate(start_date, end_date)

    img = ee.Image(ic.first())

    toa = img.select(
        ['B1', 'B2', 'B3', 'B4', 'B6', 'B8A', 'B9', 'B10', 'B11', 'B12'],
        ['aerosol', 'blue', 'green', 'red', 'red2', 'red4', 'h2o', 'cirrus', 'swir1', 'swir2']
    ).divide(10000).addBands(img.select('QA60')) \
     .set('solar_azimuth', img.get('MEAN_SOLAR_AZIMUTH_ANGLE')) \
     .set('solar_zenith', img.get('MEAN_SOLAR_ZENITH_ANGLE'))

    toa = ESAcloudMask(toa)

    Map = geemap.Map(center = aoi.centroid().coordinates().getInfo()[::-1], zoom=10)
    Map.addLayer(toa.select(['red', 'green', 'blue']), {'min': 0, 'max': 0.3}, 'RGB')
    Map.addLayer(toa.select('ESA_clouds'), {'min': 0, 'max': 1, 'palette': ['red']}, 'Chmury (QA60)')
    Map.addLayerControl()

    return Map