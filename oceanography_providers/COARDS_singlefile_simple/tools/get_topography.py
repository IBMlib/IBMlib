#!/usr/bin/env python
#############################################################################################################
#
#   Project GEBCO topography onto grid of hydrographic cells corresponding to COARDS compliant output
#
#   Usage :
#       get_topography.py   <GEBCOfile>  <hydrofile.nc>  <varname>  <sub_cell_sampling>  <topography_file> 
#
#   Rely on masked variables in hydrofile.nc to pick up topography:
#   variable.mask[iy,ix] == True => (ix,iy) is a dry point
#
#   get_topography.py  gebco_08.nc  /home/data/CMEMS/NWS/MetO-NWS-PHYS-hi-TEM.nc  votemper  10   wd.ixiy
#############################################################################################################
from   numpy   import *
import netCDF4 as NetCDF
import sys
from water_depth_GEBCO import GEBCO_GridInterpolator
# in this case apply probing variable :
#    short votemper(time, depth, lat, lon)

# -------- parse input --------
assert len(sys.argv) == 6

print "reading topography database from           %s" % sys.argv[1]
gebco = GEBCO_GridInterpolator(sys.argv[1])

ncfile = NetCDF.Dataset(sys.argv[2], 'r')
varname = sys.argv[3]
print "reading hydrographic grid from %s using variable %s" % (sys.argv[2],sys.argv[3])

nsampl = int(sys.argv[4])
print "sub cell divisions of hydrographic cells = %d" % nsampl

fout = open(sys.argv[5], "w")
print "writing sampled topography to              %s " % sys.argv[5]

# -------- process input --------

avar   = ncfile.variables[varname][0,0,:,:] # dims=(time, depth, lat, lon)  pick surface layer, first time frame
lon    = ncfile.variables['lon'][:]
lat    = ncfile.variables['lat'][:]
# (lon0, lat0) is cell center of SW corner
lon0   = lon[0]           
dlon   = lon[1]-lon[0]
lat0   = lat[0]
dlat   = lat[1]-lat[0]
assert dlon > 0
assert dlat > 0

npt  = 0               # processed wet points
wmin =  1.0e20
wmax = -1.0e20
wmin_acceptable = 0.1  # minimum acceptable water depth for a wet point
ny,nx  = avar.shape
for ix in range(nx):
    for iy in range(ny):
        if avar.mask[iy,ix] == False: # (ix,iy) is a wet point
            x0 = lon0 + dlon*(ix - 0.5)
            x1 = lon0 + dlon*(ix + 0.5)
            y0 = lat0 + dlat*(iy - 0.5)
            y1 = lat0 + dlat*(iy + 0.5)
            wd = gebco.get_average_depth(x0, x1, y0, y1, nsampl, nsampl) # assume grid close to isotropic
            wd = max(wmin_acceptable, wd)
            wmin = min(wmin, wd)
            wmax = max(wmax, wd)
            fout.write("%d  %d  %12.7f\n" % (ix+1,iy+1,wd)) # print cell with fortran indexing
            npt += 1           
fout.close()
print "processed %d wet points" % npt
print "wet fraction = %12.7f"   % (1.0*npt/nx/ny)
