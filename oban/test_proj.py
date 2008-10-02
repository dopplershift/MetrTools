from matplotlib.toolkits.basemap import Basemap
import numpy as N
import pylab as P
import matplotlib as M
import projection, oban_text, oban
from constants import *

##ps = projection.polar_stereo(60.0, 110.0)
##w = 1500 * m_per_km
##h = 1000 * m_per_km
##m = Basemap(projection='stere',lat_0=90.0,lon_0=-110.0,lat_ts=60.0,height=h,
##  width=w, rsphere=ps.radius)
##
##test_lon = -97
##test_lat = 31.25
##
##my_xy = ps.geog_to_map(test_lat, -test_lon)
##bm_xy = m(test_lon, test_lat)
##print "Mine -- \tX: %f Y: %f" % my_xy
##print "Base Conv --\tX: %f Y: %f" % (-(bm_xy[1] - h/2.0), bm_xy[0] - w/2.0)
##print "Basemap --\tX: %f Y: %f" % bm_xy

ps2 = projection.polar_stereo(60.0, 110.0, 1.0/15000000.0)
m2 = Basemap(projection='stere',lat_0=90.0,lon_0=-110.0,lat_ts=60.0,
  rsphere=ps2.radius,llcrnrlon=-120.5,llcrnrlat=25.1,urcrnrlon=-61.8,urcrnrlat=43.4,
  resolution='l',area_thresh=5000)
##m2 = Basemap(projection='stere',lat_0=90.0,lon_0=-110.0,lat_ts=60.0,
##  rsphere=ps2.radius,llcrnrlon=-131.3,llcrnrlat=19.8,urcrnrlon=-82.8,urcrnrlat=55.6,
##  resolution='l',area_thresh=5000)
##my_xy2 = N.array(ps2.geog_to_map(test_lat, -test_lon))
###my_xy2 /= ps2.map_scale
##offset = N.array(m2(-110,90))
##convy = offset[1] - my_xy2[0] / ps2.map_scale
##convx = offset[0] + my_xy2[1] / ps2.map_scale
##bm_xy2 = m2(test_lon, test_lat)
##print "Mine -- \tX: %f Y: %f" % tuple(my_xy2)
##print "Conv -- \tX: %f Y: %f" % (convx, convy)
##print "Basemap --\tX: %f Y: %f" % bm_xy2
##print ps2.to_basemap_xy(m2, *my_xy2)

data = oban_text.get_arrays('obj01.dat')
data.x, data.y = ps2.geog_to_map(data.lats, data.lons)
data.x *= cm_per_m
data.y *= cm_per_m

##testx, testy = m(-data.lons, data.lats)
##newx = -(testy - h/2.0) * ps.map_scale * cm_per_m
##newy = (testx - w/2.0) * ps.map_scale * cm_per_m

delta = 1.27
x0 = 22.86
y0 = -8.89
x = N.arange(22) * delta + x0
y = N.arange(28) * delta + y0
x_grid,y_grid = N.meshgrid(x,y)

x_bm,y_bm = ps2.to_basemap_xy(m2, x_grid / cm_per_m, y_grid / cm_per_m)
#Radius factor used in all of the analyses
R = 4.32

heights_cress = oban.analyze_grid(data.heights, x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, R)
contours = N.arange(5300., 6000., 60.0)
parallels = N.arange(25., 60.0, 10.0)
meridians = N.arange(-120., -60.0, 10.0)

##m2.drawstates()
##m2.drawcountries()
##m2.drawcoastlines()
##m2.drawparallels(parallels, labels = [1,1,0,0])
##m2.drawmeridians(meridians, labels = [0,0,1,1])
ob_x, ob_y = ps2.to_basemap_xy(m2, data.x / cm_per_m, data.y / cm_per_m)
#m2.plot(ob_x, ob_y, 'bx')
#m2.plot(x_bm, y_bm, 'g.')
for name in M.cm.cmapnames:
  f = P.figure()
  m2.drawstates()
  m2.drawcountries()
  m2.drawcoastlines()
  m2.drawparallels(parallels, labels = [1,1,0,0])
  m2.drawmeridians(meridians, labels = [0,0,1,1])
  c = m2.contour(x_bm, y_bm, heights_cress, contours, cmap=M.cm.get_cmap(name))
  P.title(name)
  #P.clabel(c)
  P.show()
