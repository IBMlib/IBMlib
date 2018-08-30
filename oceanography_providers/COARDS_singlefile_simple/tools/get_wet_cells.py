#!/usr/bin/env python
#################################################################################################
#
#   Extract a wet cell file in format compartible with that expected by water_depth_GEBCO.py
#
#   Usage :
#             get_wet_cells.py   < hydrofile.nc >   < wetpt_file >
#  
#   Adapt this script in relation to actual variable names in hydrography file and wet/dry marking
#################################################################################################
from   numpy   import *
import netCDF4 as NetCDF
import sys

# in this case apply probing variable :
#    short votemper(time, depth, lat, lon)

ncfile = NetCDF.Dataset(sys.argv[1], 'r')
avar   = ncfile.variables['votemper'][0,0,:,:] # pick surface layer
lon    = ncfile.variables['lon'][:]
lat    = ncfile.variables['lat'][:]
fout = open(sys.argv[2], "w")
fout.write("%12.5f %12.5f\n" % (lon[0], lon[1]-lon[0]))
fout.write("%12.5f %12.5f\n" % (lat[0], lat[1]-lat[0]))
ny,nx  = avar.shape
for ix in range(nx):
    for iy in range(ny):
        if avar[iy,ix] > -32768: # wet test for this variable
            #print ix+1,iy+1      # fortran indexing
            fout.write("%d %d\n" % (ix+1,iy+1)) # fortran indexing
fout.close()
