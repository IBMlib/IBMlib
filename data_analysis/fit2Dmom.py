#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################################################
#    Fit 2D distribution to rotated 2D normal via second moments (mxx,myy,mxy)
#    expressed as ellipse with axes a,b and rotation phi of a axe from x-axis
#    Co-plot data and fit; scale ellipse so it contains 68 % of PDF
#
#                                     1
#    vcfac2d = Sqrt[-2 Log[1 - Erf[-------]]]
#                                  Sqrt[2]
#
#############################################################################################
from numpy import *
import scipy
import matplotlib.pyplot as plt

a   = 3.05  # example a axis
b   = 1.77  # example b axis
phi = 0.43  # example rotation, in radian
n   = 300   # example set size

# --- in rotated coordinates distribution factorizes
xp   = scipy.stats.norm.rvs(scale=a, size=n)
yp   = scipy.stats.norm.rvs(scale=b, size=n)

# --- rotate to base Cartesian
x = xp*cos(phi) - yp*sin(phi)
y = xp*sin(phi) + yp*cos(phi)

mxx = mean(x*x)
myy = mean(y*y)
mxy = mean(x*y)

# ------------ fit data ------------

phifit = 0.5*arctan(2*mxy/(mxx-myy))
s      = sin(phifit)
c      = cos(phifit)
denom  = c**4 - s**4
afit   = sqrt(( mxx*c**2 - myy*s**2)/denom)
bfit   = sqrt((-mxx*s**2 + myy*c**2)/denom)

print("fit            = ", afit, bfit, phifit)
print("data generator = ", a,    b,    phi)

# ------------ plot data + fit ------------

vcfac2d = 1.51517 # scale ellipse so it contains 68 % of PDF
plt.scatter(x,y)
nplt          = 200
unitc         = linspace(0,2*pi,nplt)
xye           = zeros((2,nplt),float)
xv            = vcfac2d*afit*cos(unitc)
yv            = vcfac2d*bfit*sin(unitc)
xye[0,:]      = xv*c - yv*s
xye[1,:]      = xv*s + yv*c
plt.plot(xye[0,:] ,xye[1,:], c='k')      # principal ellipsis
plt.show()
