import ee
import folium

# Authenticate and initialize Earth Engine
ee.Authenticate()
ee.Initialize(project='ee-skoraremotesensing')

# EE config
aoi = ee.Geometry.Rectangle([19.8, 50.0, 20.0, 50.2])
START_DATE = '2023-06-01'
END_DATE = '2023-07-01'
CLOUD_FILTER = 60
CLD_PRB_THRESH = 50
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 1
BUFFER = 50

def get_s2_sr_cld_col(aoi, start_date, end_date):
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))
    
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))
    
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))

def add_cloud_bands(img):
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def add_shadow_bands(img):
    not_water = img.select('SCL').neq(6)
    dark_pixels = img.select('B8').lt(NIR_DRK_THRESH * 1e4).multiply(not_water).rename('dark_pixels')
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST * 10)
                .reproject(crs=img.select(0).projection(), scale=100)
                .select('distance')
                .mask()
                .rename('cloud_transform'))
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shdw_mask(img):
    img_cloud = add_cloud_bands(img)
    img_cloud_shadow = add_shadow_bands(img_cloud)
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
        .reproject(crs=img.select([0]).projection(), scale=20)
        .rename('cloudmask'))
    return img_cloud_shadow.addBands(is_cld_shdw)

def apply_cld_shdw_mask(img):
    not_cld_shdw = img.select('cloudmask').Not()
    return img.select('B.*').updateMask(not_cld_shdw)

# Folium helper
def add_ee_layer(self, ee_image_object, vis_params, name, show=True, opacity=1, min_zoom=0):
    map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
    folium.raster_layers.TileLayer(
        tiles=map_id_dict['tile_fetcher'].url_format,
        attr='Map Data © Google Earth Engine',
        name=name,
        show=show,
        opacity=opacity,
        min_zoom=min_zoom,
        overlay=True,
        control=True
    ).add_to(self)

folium.Map.add_ee_layer = add_ee_layer

# Main functions for external use
def display_cloud_map():
    col = get_s2_sr_cld_col(aoi, START_DATE, END_DATE).map(add_cld_shdw_mask)
    img = col.mosaic()

    center = aoi.centroid(10).coordinates().reverse().getInfo()
    m = folium.Map(location=center, zoom_start=12)

    m.add_ee_layer(img, {'bands': ['B4', 'B3', 'B2'], 'min': 0, 'max': 2500}, 'S2 image')
    m.add_ee_layer(img.select('probability'), {'min': 0, 'max': 100}, 'cloud probability', False)
    m.add_ee_layer(img.select('clouds').selfMask(), {'palette': 'e056fd'}, 'clouds', False)
    m.add_ee_layer(img.select('cloud_transform'), {'min': 0, 'max': 1, 'palette': ['white', 'black']}, 'cloud_transform', False)
    m.add_ee_layer(img.select('dark_pixels').selfMask(), {'palette': 'orange'}, 'dark_pixels', False)
    m.add_ee_layer(img.select('shadows').selfMask(), {'palette': 'yellow'}, 'shadows', False)
    m.add_ee_layer(img.select('cloudmask').selfMask(), {'palette': 'orange'}, 'cloudmask', True, 0.5)

    m.add_child(folium.LayerControl())
    display(m)

def display_cloud_free_mosaic():
    col = get_s2_sr_cld_col(aoi, START_DATE, END_DATE)
    mosaic = (col.map(add_cld_shdw_mask)
                 .map(apply_cld_shdw_mask)
                 .median())

    center = aoi.centroid(10).coordinates().reverse().getInfo()
    m = folium.Map(location=center, zoom_start=12)
    m.add_ee_layer(mosaic,
                   {'bands': ['B4', 'B3', 'B2'], 'min': 0, 'max': 2500, 'gamma': 1.1},
                   'S2 cloud-free mosaic')
    m.add_child(folium.LayerControl())
    display(m)