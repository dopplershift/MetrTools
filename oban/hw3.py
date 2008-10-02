import sys
import numpy as N
import pylab as P
import matplotlib as M
from matplotlib.toolkits.basemap import Basemap
import oban_text, projection, oban, misc
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

#Mean station separation
d = 2.30564

#Cressman first pass
cress1 = oban.analyze_grid(data.heights, x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, 4 * d)

#Get observation increments for this pass
cress1_incs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, cress1,
  cressman_radius = 4 * d)
cress1_rms = misc.rms(cress1_incs)
print "Cressman 1st Pass rms: %f" % cress1_rms

#Cressman 2nd pass
cress2 = oban.analyze_grid(N.array(cress1_incs), x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, 2.5 * d) + cress1

#Get observation increments for this pass
cress2_incs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, cress2,
  cressman_radius = 2.5 * d)
cress2_rms = misc.rms(cress2_incs)
print "Cressman 2nd Pass rms: %f" % cress2_rms

#Cressman 3rd pass
cress3 = oban.analyze_grid(N.array(cress2_incs), x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, 1.5 * d) + cress2

#Get observation increments for this pass
cress3_incs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, cress3,
  cressman_radius = 1.5 * d)
cress3_rms = misc.rms(cress3_incs)
print "Cressman 3rd Pass rms: %f" % cress3_rms

#Get RMS of differences
print "Cressman 1st - 2nd rms: %f" % misc.rms(cress1 - cress2)
print "Cressman 1st - 3rd rms: %f" % misc.rms(cress1 - cress3)
print "Cressman 2nd - 3rd rms: %f" % misc.rms(cress2 - cress3)

#Cressman single pass for comparison with Barnes passes
cress_187 = oban.analyze_grid(data.heights, x_grid, y_grid, data.x, data.y,
  oban.cressman_weights, 1.87 * d)
diffs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, cress_187,
  cressman_radius = 1.87 * d)
cress_187_rms = misc.rms(diffs)
print "Cressman 1.87d rms: %f" % cress_187_rms

#Do the list of 2-pass barnes analyses
gammas = N.array([1.0,0.4,0.2])
barnes_analyses = []
barnes_rms = []
kappa0 = oban.calc_barnes_param(d)
for g in gammas:
  field = oban.analyze_grid_multipass(data.heights, x_grid, y_grid, data.x, 
    data.y, 2, oban.barnes_weights, (kappa0, g))
  barnes_analyses.append(field)
  diffs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, field)
  barnes_rms.append(misc.rms(diffs))
  print "Barnes 2-pass gamma=%.1f rms: %f" % (g, barnes_rms[-1])

#Do the barnes 3-pass
barnes_3pass = oban.analyze_grid_multipass(data.heights, x_grid, y_grid,
  data.x, data.y, 3, oban.barnes_weights, (kappa0, 1.0))
diffs = oban.get_ob_incs(data.x, data.y, data.heights, x, y, barnes_3pass)
barnes_3pass_rms = misc.rms(diffs)
print "Barnes 3-pass rms: %f" % barnes_3pass_rms

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
fields = (cress1, cress2, cress3, cress_187, barnes_analyses[0],
  barnes_analyses[1], barnes_analyses[2], barnes_3pass)
names = ('Cressman 1st Pass', 'Cressman 2nd Pass', 'Cressman 3rd Pass',
  'Cressman R=1.87d', r'$\rm{Barnes} \gamma=1.0$',r'$\rm{Barnes} \gamma=0.4$',
  r'$\rm{Barnes} \gamma=0.2$', r'$\rm{Barnes 3-pass} \gamma=1.0$')
filenames = ('cress1', 'cress2', 'cress3', 'cress187', 'barnes1', 'barnes2',
  'barnes3', 'barnes3pass')
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

filename = 'cress12diff'
f = P.figure()
bm_ps.drawstates()
bm_ps.drawcountries()
bm_ps.drawcoastlines()
bm_ps.drawparallels(parallels, labels = [1,1,0,0])
bm_ps.drawmeridians(meridians, labels = [0,0,1,1])
cp = bm_ps.contour(x_bm, y_bm, cress1 - cress2, cmap=colormap)
P.clabel(cp, fontsize=fontsize, fmt='%.1f')
f.text(0.5,0.95,'Cressman 1st - 2nd', horizontalalignment='center',
  fontsize=16)
if save_work:
  P.savefig(filename + '.png', dpi = 300)
else:
  P.show()

filename = 'cress13diff'
f = P.figure()
bm_ps.drawstates()
bm_ps.drawcountries()
bm_ps.drawcoastlines()
bm_ps.drawparallels(parallels, labels = [1,1,0,0])
bm_ps.drawmeridians(meridians, labels = [0,0,1,1])
cp = bm_ps.contour(x_bm, y_bm, cress1 - cress3, cmap=colormap)
P.clabel(cp, fontsize=fontsize, fmt='%.1f')
f.text(0.5,0.95,'Cressman 1st - 3rd', horizontalalignment='center',
  fontsize=16)
if save_work:
  P.savefig(filename + '.png', dpi = 300)
else:
  P.show()

filename = 'cress23diff'
f = P.figure()
bm_ps.drawstates()
bm_ps.drawcountries()
bm_ps.drawcoastlines()
bm_ps.drawparallels(parallels, labels = [1,1,0,0])
bm_ps.drawmeridians(meridians, labels = [0,0,1,1])
cp = bm_ps.contour(x_bm, y_bm, cress2 - cress3, cmap=colormap)
P.clabel(cp, fontsize=fontsize, fmt='%.1f')
f.text(0.5,0.95,'Cressman 2nd - 3rd',horizontalalignment='center',
  fontsize=16)
if save_work:
  P.savefig(filename + '.png', dpi = 300)
else:
  P.show()

##  if save_work:
##    sys.stdout = open(filename + '.txt', 'w')
##  print name
##  print "\nAnalysis Array:"
##  for n in N.arange(y.size):
##    print "\t%d" % (n + 1),
##  print
##  for (num, row) in enumerate(field.T):
##    print "%d\t" % (num + 1),
##    for col in row:
##      print "%7.2f" % col,
##    print
##  print "\nk\tSTID\tZ_o\t\tZ_a\t\tZ_o - Z_a"
##  for i in range(len(data.x)):
##    try:
##      interp_height = oban.bilinear(x, y, field, data.x[i], data.y[i])
##      hd = data.heights[i] - interp_height
##      height_diff.append(hd)
##      print "%d\t%.0f\t%.4f\t%.4f\t%.4f" % (i, data.stids[i], data.heights[i],\
##        interp_height, hd)
##    except ValueError:
##      print "%d\t%.0f\t%.4f\tOff Grid\tOff Grid" % (i, data.stids[i],\
##        data.heights[i])
##  height_diff = N.array(height_diff)
##  rms = N.sqrt(N.average(height_diff * height_diff))
##  print "RMS: %.4f m" % rms
##  if save_work:
##    sys.stdout.close()
