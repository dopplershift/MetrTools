from scipy.interpolate import interpolate
from numpy.random import randn
from numpy import meshgrid,arange,array,dot
import oban

data = randn(16).reshape(4,4)
x = y = arange(4)
X,Y = meshgrid(x,y)

x_loc = 1.4
y_loc = 2.7

full_bilin = interpolate.interp2d(X,Y,data,kind='linear')
partial_bilin = interpolate.interp2d(X[1:,1:],Y[1:,1:],data[1:,1:],kind='linear')
single_bilin = interpolate.interp2d(X[2:4,1:3],Y[2:4,1:3],data[2:4,1:3],kind='linear')

print full_bilin(x_loc,y_loc)
print partial_bilin(x_loc,y_loc)
print single_bilin(x_loc,y_loc)
top = data[2,1] + (data[2,2] - data[2,1])*0.4
bottom = data[3,1] + (data[3,2] - data[3,1])*0.4
print top + (bottom - top)*0.7
xw = 0.4
yw = 0.7
xws = array([1-xw, xw])
yws = array([1-yw, yw])
print dot(yws,dot(data[2:4,1:3],xws))
print oban.bilinear(x,y,data,x_loc,y_loc)

print data[2:4,1:3]
