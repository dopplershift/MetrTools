import sys
import numpy as N
import pylab as P
import matplotlib as M
from matplotlib.toolkits.basemap import Basemap
import oban_text, projection, oban
from constants import *

data = oban_text.get_arrays('obj01.dat')
data.wind_speeds *= ms_per_kt
f = 2 * earth_omega * N.sin(data.lats * rad_per_deg)
geos_factor = f / g

ps = projection.polar_stereo(60.0, 110.0, 1.0/15000000.0)
data.sigma = ps.image_scale(data.lats)

#Convert the ob location to map coordinates in cm
data.x, data.y = ps.geog_to_map(data.lats, data.lons)
data.x *= cm_per_m
data.y *= cm_per_m

#Calculate wind components on the map
data.u, data.v = oban.get_wind_comps(data.wind_speeds, data.wind_dirs)
data.u_map, data.v_map = ps.orient_wind(data.u, data.v, data.lons)

#Generate grid of x,y positions
delta = 1.27
x0 = 22.86
y0 = -8.89
x = N.arange(22) * delta + x0
y = N.arange(28) * delta + y0
x_grid,y_grid = N.meshgrid(x,y)

#Radius factor used in all of the analyses
R = 4.32

#Perform analysis of height obs using uniform weighting function
heights_uniform = oban.analyze_grid(data.heights, x_grid, y_grid, data.x,\
  data.y, oban.uniform_weights, R)

#Perform analysis of height obs using Cressman weights
heights_cress = oban.analyze_grid(data.heights, x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, R)

#Generate adjusted height field using the (assumed) geostrophic wind on the map 
#to estimate height gradients
map_factor = geos_factor / (ps.map_scale * data.sigma * cm_per_m)
adj_height = oban.adjust_field(data.heights, map_factor * data.v_map,
  -map_factor * data.u_map, x_grid, y_grid, data.x, data.y)
  
#Perform analysis of height obs using winds to estimate height gradients and
#a uniform weighting function
heights_uniforbm_ps = oban.analyze_grid(adj_height, x_grid, y_grid, data.x, data.y,
  oban.uniform_weights, R)

#Perform analysis of height obs using winds to estimate height gradients and
#a uniform weighting function
heights_cress2 = oban.analyze_grid(adj_height, x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, R)

#Generate a grid for basemap plotting
bm_ps = Basemap(projection='stere',lat_0=90.0,lon_0=-110.0,lat_ts=60.0,
  rsphere=ps.radius,llcrnrlon=-120.5,llcrnrlat=25.1,urcrnrlon=-61.8,
  urcrnrlat=43.4,resolution='l',area_thresh=5000)

#Transform the grid to basemap space for plotting
x_bm,y_bm = ps.to_basemap_xy(bm_ps, x_grid / cm_per_m, y_grid / cm_per_m)

#Check if we want to view or save output
if len(sys.argv) > 1 and sys.argv[1].startswith('silent'):
  save_work = True
  colormap = M.cm.get_cmap('gist_gray')
  fontsize = 10
else:
  save_work = False
  colormap = M.cm.get_cmap('jet')
  fontsize = 14

#Do the necessary processing for each analysis
contours = N.arange(5300., 6000., 60.0)
parallels = N.arange(25., 60.0, 10.0)
meridians = N.arange(-120., -60.0, 10.0)
fields = (heights_uniform, heights_cress, heights_uniforbm_ps, heights_cress2)
names = ('Uniform Weights', 'Cressman Weights',\
  'Uniform weights with wind correction',\
  'Cressman weights with wind correction')
filenames = ('uniform1','cress1','uniform2','cress2')
for field,name,filename in zip(fields, names, filenames):
  f = P.figure()
  bm_ps.drawstates()
  bm_ps.drawcountries()
  bm_ps.drawcoastlines()
  bm_ps.drawparallels(parallels, labels = [1,1,0,0])
  bm_ps.drawmeridians(meridians, labels = [0,0,1,1])
  cp = bm_ps.contour(x_bm, y_bm, field, contours, cmap=colormap)
  P.clabel(cp, fontsize=fontsize, fmt='%.1f')
  f.text(0.5,0.95,name,horizontalalignment='center',fontsize=16)
  if save_work:
    P.savefig(filename + '.png', dpi = 300)
  else:
    P.show()
  height_diff = list()
  if save_work:
    sys.stdout = open(filename + '.txt', 'w')
  print name
  print "\nAnalysis Array:"
  for n in N.arange(y.size):
    print "\t%d" % (n + 1),
  print
  for (num, row) in enumerate(field.T):
    print "%d\t" % (num + 1),
    for col in row:
      print "%7.2f" % col,
    print
  print "\nk\tSTID\tZ_o\t\tZ_a\t\tZ_o - Z_a"
  for i in range(len(data.x)):
    try:
      interp_height = oban.bilinear(x, y, field, data.x[i], data.y[i])
      hd = data.heights[i] - interp_height
      height_diff.append(hd)
      print "%d\t%.0f\t%.4f\t%.4f\t%.4f" % (i, data.stids[i], data.heights[i],\
        interp_height, hd)
    except ValueError:
      print "%d\t%.0f\t%.4f\tOff Grid\tOff Grid" % (i, data.stids[i],\
        data.heights[i])
  height_diff = N.array(height_diff)
  rms = N.sqrt(N.average(height_diff * height_diff))
  print "RMS: %.4f m" % rms
  if save_work:
    sys.stdout.close()
