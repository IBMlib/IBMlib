#!/usr/bin/env python
#############################################################################################################
#
#   Aggregate native isotropic grid into sub grid of larger cells with corresponding to wet fraction > minwetfrac
#
#   get_wet_boxes.py  ncfile                                        varname   lonmin latmin lonmax latmax  nsp  minwetfrac
#   get_wet_boxes.py  /home/data/CMEMS/NWS/MetO-NWS-PHYS-hi-TEM.nc  votemper    -999  -999   999     999    1     0.99   > wet_allpts.lonlat
#   get_wet_boxes.py  /home/data/CMEMS/NWS/MetO-NWS-PHYS-hi-TEM.nc  votemper    -2.5   51     9      57     2     0.49   > wet_src.lonlat
#   get_wet_boxes.py  /home/data/CMEMS/NWS/MetO-NWS-PHYS-hi-TEM.nc  votemper    -4     49     12     62     2     0.49   > wet_dest.lonlat
#
#############################################################################################################
from   numpy   import *
import netCDF4 as NetCDF
import sys

def get_grid_coor(x,y): return (x-lon0)/dlon,  (y-lat0)/dlat  # lon,lat -> C grid coor; SW @ (0,0) + cell centers == integers
def get_lonlat(sx,sy):  return lon0 + sx*dlon, lat0 + sy*dlat #  C grid coor -> lon,lat
def get_subgrid_coor(x,y): return (x-lon0_sub)/dlon_sub,  (y-lat0_sub)/dlat_sub  # lon,lat -> C grid coor; SW @ (0,0) + cell centers == integers
def nint(x): return int(0.5+x)
    
# -------- parse input --------

ncfile     = NetCDF.Dataset(sys.argv[1], 'r')
varname    = sys.argv[2]
lonmin     = float(sys.argv[3])
latmin     = float(sys.argv[4])
lonmax     = float(sys.argv[5])
latmax     = float(sys.argv[6])
nsp        = int(sys.argv[7])
minwetfrac = float(sys.argv[8])

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

wetmask = where(avar.mask == False, 1, 0) # False -> wet
wetmask = transpose(wetmask)              # now (lon,lat)

#nsp = 1
nx,ny   = wetmask.shape
nx0,ny0 = map(nint, get_grid_coor(lonmin,latmin))  # SW corner of scan
nx1,ny1 = map(nint, get_grid_coor(lonmax,latmax))  # NE corner of scan

nx0 = max(0, nx0)  # confine to full grid
ny0 = max(0, ny0)  # confine to full grid
nx1 = min(nx, nx1) # confine to full grid
ny1 = min(ny, ny1) # confine to full grid

# ---- grid descriptor of sub grid
lon0_sub, lat0_sub = get_lonlat(nx0+0.5*(nsp-1), ny0+0.5*(nsp-1)) # SW corner of sub grid
dlon_sub = dlon*nsp
dlat_sub = dlat*nsp
nx_sub = len(range(nx0,nx1,nsp))
ny_sub = len(range(ny0,ny1,nsp))
print "# nxsub,nysub= %d %d" % (nx_sub,   ny_sub)
print "# SW-corner=   %f %f" % (lon0_sub, lat0_sub)
print "# dlonlat=     %f %f" % (dlon_sub, dlat_sub)
print "# x0         y0         x1       y1       ix(offset0) iy(offset0)" 
# ---- sub grid scan on native grid ----
for ix in range(nx0,nx1,nsp):     # fill grid scannx0,ny0 = map(nint, get_grid_coor(lonmin,latmin))
    for iy in range(ny0,ny1,nsp): # fill grid scan
        #if wetmask[ix,iy]: print ix,iy
        #if wetmask[ix,iy]: print "%f %f" % get_lonlat(ix,iy)
        x0,y0   = get_lonlat(ix - 0.5,       iy - 0.5)       # SW corner
        x1,y1   = get_lonlat(ix - 0.5 + nsp, iy - 0.5 + nsp) # NE corner
        xc      = 0.5*(x0+x1) # x center of aggregated cell
        yc      = 0.5*(y0+y1) # y center of aggregated cell
        ix_sub, iy_sub = map(nint, get_subgrid_coor(xc,yc))
        wetfrac = 1.0*sum(wetmask[ix:ix+nsp, iy:iy+nsp])/nsp/nsp
        #if wetfrac>minwetfrac: print xc,yc
        if wetfrac>minwetfrac: print x0,y0,x1,y1, ix_sub, iy_sub
        #if wetfrac>minwetfrac: print ix_sub, iy_sub
        #fout.write("%d  %d  %12.7f\n" % (ix+1,iy+1,wd)) # print cell with fortran indexing
        #npt += 1           
#fout.close()
#print "processed %d wet points" % npt
#print "wet fraction = %12.7f"   % (1.0*npt/nx/ny)
