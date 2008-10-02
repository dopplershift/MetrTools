from numpy import cos, sin, hypot
from numpy import arctan2 as atan2
from constants import *

class polar_stereo:
  def __init__(self, std_lat, std_lon, map_scale = 1.0, radius = 6.371e6):
    self.std_lat = std_lat
    self.std_lon = std_lon
    self.origin_lat = 90.0
    self.map_scale = map_scale
    self.radius = radius
    self.__img_factor = 1 + sin(self.std_lat * rad_per_deg)
  def image_scale(self, lat):
    '''Returns the image scale factor for the given latitude, lat, given in
       degrees'''
    return self.__img_factor / (1 + sin(lat * rad_per_deg))
  def geog_to_map(self, lat, lon):
    '''Converts the latitude and longitude coordinates, in degrees, to
       coordinates on the map'''
    lon_dev = (self.std_lon - lon) * rad_per_deg
    scale_factor = self.map_scale * self.image_scale(lat) * self.radius \
        * cos(lat * rad_per_deg)
    x = scale_factor * cos(lon_dev)
    y = scale_factor * sin(lon_dev)
    return (x,y)
  def map_to_geog(self, x, y):
    '''Converts x,y coordinates on the map to latitude and longitude in
       degrees'''
    lon = self.std_lon - atan2(y,x) / rad_per_deg
    lat = (pi/2.0  - 2 * atan2(hypot(x,y)/self.map_scale,\
        self.radius*self.__img_factor)) / rad_per_deg
    return (lat,lon)
  def orient_wind(self, u, v, lon):
    '''Given a u,v in normal coords, generate U,V relative to the map coords'''
    l = (self.std_lon - lon) * rad_per_deg
    U = -u * sin(l) - v * cos(l)
    V = u * cos(l) - v * sin(l)
    return U,V
  def to_basemap_xy(self, bm_proj, x, y):
    '''Given x,y in funky coordinates, convert to x,y in basemap space'''
    offset = bm_proj(-self.std_lon, self.origin_lat)
    return offset[0] + y / self.map_scale, offset[1] - x / self.map_scale
