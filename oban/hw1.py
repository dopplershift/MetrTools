#!/usr/bin/python

import math, sys
import numpy as N
import pylab as pl
import oban_text, projection

deg_to_rad = math.pi / 180.0
m_to_km = 1.0 / 1000.0
kt_to_ms = 0.5144444

data = oban_text.get_obs('obj01.dat')
ps = projection.polar_stereo(60.0, 110.0, 1.0/15000000.0)

#Convert the ob location to map coordinates
for ob in data:
  ob.x,ob.y = ps.geog_to_map(ob.lat, ob.lon)

#Find distance to nearest station for each station
for ob in data:
  dist_stid = [(math.hypot(ob.x - n_ob.x, ob.y - n_ob.y), n_ob.stid) \
    for n_ob in data if n_ob is not ob]
  min_loc = N.argmin([ds[0] for ds in dist_stid])
  ob.min_dist,ob.nearest = dist_stid[min_loc]

#Compute average minimum distance to nearest station
avg_min_dist = N.average([ob.min_dist for ob in data])
avg_min_dist_km = avg_min_dist / ps.map_scale * m_to_km

#Compute earth-relative velocity components in m s^-1
for ob in data:
  ob.u = -ob.speed * kt_to_ms * math.sin(ob.dir * deg_to_rad)
  ob.v = -ob.speed * kt_to_ms * math.cos(ob.dir * deg_to_rad)

#Generate grid of x,y positions
delta = 1.27
x0 = 22.86
y0 = -8.89
x = N.arange(22) * delta + x0
y = N.arange(28) * delta + y0
x_grid,y_grid = N.meshgrid(x,y)

obx = N.array([ob.x for ob in data]) * 100.0
oby = N.array([ob.y for ob in data]) * 100.0

#Calculate array of distance of each grid point to each ob station
dist_grid = N.hypot(x_grid[...,N.newaxis] - obx[N.newaxis,N.newaxis,...],
  y_grid[...,N.newaxis] - oby[N.newaxis,N.newaxis,...])

#Get the minimum along the axis of ob stations
min_grid_dists = dist_grid.min(axis=2)

#Calculate histogram
hist_delta = 0.2
num_bins = 20
freq,left_edge = N.histogram(min_grid_dists, bins = num_bins,\
  range = [0.0, hist_delta*num_bins])

#And cumulative distribution
cum_dist = freq.cumsum()

#sys.stdout = open('hw1_out1.txt', 'w')
print "STID\tX(cm)\tY(cm)\tNearest STID\tNearest Dist (cm)\tu (m/s)\tv (m/s)"
for ob in data:
  print "%5.0f\t%7.4f\t%7.4f\t%5.0f\t\t%7.4f\t\t\t%7.4f\t%7.4f" % \
    (ob.stid, ob.x*100.0, ob.y*100.0, ob.nearest, ob.min_dist*100.0, ob.u, ob.v)
print "Average distance to closest station -- %f cm %f km" % \
  (avg_min_dist * 100.0, avg_min_dist_km)

#sys.stdout = open('hw1_out2.txt', 'w')
print "Array of distance from each grid point to nearest station (cm), where"
print "each [...] is a row:"
print N.round(min_grid_dists.T, 4)
print "\nFor distance from grid point to nearest station:"
print "\tAverage:\t%.4f cm" % min_grid_dists.mean()
print "\tMaximum:\t%.4f cm" % min_grid_dists.max()
print "\tMinimum:\t%.4f cm" % min_grid_dists.min()
print "\tMode:\t\t%.4f cm" % 1.4
print "\tMedian:\t\t%.4f cm" % N.median(min_grid_dists.flat)
print "\nDistribution information:"
print "Bin Start (cm)\tBin End (cm)\tFreq.\tCum. Freq."
for i in range(len(left_edge)):
  print "%4.2f\t\t%4.2f\t\t%d\t%d" % \
    (left_edge[i], left_edge[i] + hist_delta, freq[i], cum_dist[i])

pl.subplot(2,1,1)
pl.bar(left_edge, freq, hist_delta, color='white')
pl.xlabel('Distance to nearest station (cm)')
pl.ylabel('# of stations')
pl.title('Frequency distribution of distances from grid point to nearest station')
pl.subplot(2,1,2)
pl.bar(left_edge, cum_dist, hist_delta, color='white')
pl.xlabel('Distance to nearest station (cm)')
pl.ylabel('Cumulative # of stations')
pl.title('Cumulative frequency distribution of distances from grid points to nearest station')
pl.subplots_adjust(hspace=0.4)
#pl.show()
#pl.savefig('histograms.png', dpi=300)
