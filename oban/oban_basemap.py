import numpy as N
import pylab as P
import matplotlib as M
from matplotlib.toolkits.basemap import Basemap

def to_basemap_xy(bm_proj, x, y, origin_lon, origin_lat, map_scale):
  '''Given x,y in funky coordinates, convert to x,y in basemap space'''
  offset = bm_proj(origin_lon, origin_lat)
  return offset[0] + y / map_scale, offset[1] - x / map_scale

# Read data from text file where values are were written separated by spaces
# and newlines, iterating over y first, then x -- for instance, a text output
# of x in columns and y in rows. (x and y here are in the normal sense)
data = N.fromfile('data.txt', sep=' ').reshape(28,22) 

#Generate grid of x,y positions
delta = 1.27
x0 = 22.86
y0 = -8.89
x = N.arange(22) * delta + x0
y = N.arange(28) * delta + y0

#Makes correct 2D grids of x and y -- x in columns, y in rows
x_grid,y_grid = N.meshgrid(x,y) 

# Generate a grid for basemap plotting
# The llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat are the lower-left and
# upper right corners of the desired plot area, and can be set to anything
# Resolution specifies what resolution map outlines to use, and area_thresh
# specifies how large a map feature must be to be plotted
bm_ps = Basemap(projection='stere',lat_0=90.0,lon_0=-110.0,lat_ts=60.0,
  rsphere=6371000,llcrnrlon=-120.5,llcrnrlat=25.1,urcrnrlon=-61.8,
  urcrnrlat=43.4,resolution='l',area_thresh=5000)

# Transform the grid to basemap space for plotting
cm_per_m = 100.0
x_bm,y_bm = to_basemap_xy(bm_ps, x_grid / cm_per_m, y_grid / cm_per_m,
  -110.0, 90.0, 1.0/15000000.0)

#Specify contour intervals, as well as what parallels and meridians to show
contours = N.arange(5300., 6000., 60.0)
parallels = N.arange(25., 60.0, 10.0)
meridians = N.arange(-120., -60.0, 10.0)

#These draw the various maps (self explanatory)
bm_ps.drawstates()
bm_ps.drawcountries()
bm_ps.drawcoastlines()

#Draw grid lines, labels specify where to label the grid lines
#[left, right, top, bottom]
bm_ps.drawparallels(parallels, labels = [1,1,0,0])
bm_ps.drawmeridians(meridians, labels = [0,0,1,1])

#Remove the cmap argument below to get the normal, colored colortable
cp = bm_ps.contour(x_bm, y_bm, data, contours, cmap=M.cm.get_cmap('gist_gray'))

#Label contours
P.clabel(cp)
#Normally just P.title('My plot'), but the longitude labels get overlapped
P.gcf().text(0.5,0.95,'My Plot',horizontalalignment='center',fontsize=16)

#P.savefig('plot.png', dpi = 300) # This saves to a file
P.show()
