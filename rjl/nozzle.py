from pylab import *

xcell = linspace(0,32,200)
rad1 = 1-0.3*(1 + cos(pi*(xcell-16)/8.))
radius = where(abs(xcell-16)<8, rad1, 1.)
clf()
plot(xcell, radius, 'b')
plot(xcell, -radius, 'b')
