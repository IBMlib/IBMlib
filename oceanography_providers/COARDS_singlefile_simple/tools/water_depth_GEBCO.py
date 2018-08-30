#!/usr/bin/env python
#################################################################################################
#
#   Interface to GEBCO global topography:
#          0.5 minute (GEBCO/gebco_08.nc) + 1 minute (GEBCO/gridone.nc) resolution
#          documentation: GEBCO/gebco_08.pdf  GEBCO/gridone.pdf
#          http://www.gebco.net/data_and_products/gridded_bathymetry_data/
#          data sets: 'GEBCO/gridone.nc', 'GEBCO/gebco_08.nc'
#          vertical resolution of data is 1 meter steps (data represented as integer in data sets)
#
#   script usage:
#        water_depth_GEBCO.py   <GEBCOfile>  <input_cells>  <outputfile>  <sub_cell_sampling>
#
#   Trimmed interface derived from ~/DTU/JomfruhummerGUDP/Data/DybdeKort/water_depth.py
#   Scientific.IO.NetCDF replaced with netCDF4, fixed API differences
#   removed FFT kit
#   default load == full Earth
#################################################################################################
from   numpy   import *
import math
#from Scientific.IO import NetCDF
import netCDF4 as NetCDF
import os
import sys

earth_radius      = 6371.0    # in kilometers
deg2rad           = pi/180

pad_depth_invalid = 0.0       # vertex value assigned to invalid/land points

NWEurope = (  -8,  30,  49,  62)   # GEBCO data exceeds mem, allow default subgrid (E,W,S,N)
Earth    = (-180, 180, -90,  90)   # full Earth

verbose  = False 

def nint(x): return int(x+0.5)

class LonLatGrid:
    # ---------------------------------------------------------
    # Vertex oriented lon-lat grid (data on vertices)
    # 
    # Grid coordinates: (x,y) = (0,0) ->  (lon,lat) = (lambda0,phi0)
    # integer (x,y) corresponds to vertices
    # (lon,lat) in degrees
    # grid axes are oriented (W -> E)  and  (S -> N)
    # ---------------------------------------------------------
    def __init__(self,nx,ny,lambda0,phi0,dlambda,dphi):
        self.nx      = nx
        self.ny      = ny
        self.lambda0 = lambda0
        self.phi0    = phi0
        self.dlambda = dlambda
        self.dphi    = dphi
        self.drdx    = lambda phi, dl=dlambda: earth_radius*dl*deg2rad*cos(phi*deg2rad) # kilometers/grid_step
        self.drdy    = earth_radius*dphi*deg2rad                                        # kilometers/grid_step
   
    def inside_grid(self, lam, phi):
        x,y = self.get_xy(lam, phi)
        return (0 <= x <= self.nx-1) and (0 <= y <= self.ny-1)
    
    def grid_range(self):
        #   --- interior grid range ---
        return ((self.lambda0, self.lambda0 + (self.nx-1)*self.dlambda),
                (self.phi0,    self.phi0    + (self.ny-1)*self.dphi))
    
    # convert grid coordinate (x,y) to (lambda,phi). All angles in degrees 
    def get_lon(self, x):       return self.lambda0 + x*self.dlambda    
    def get_lat(self, y):       return self.phi0    + y*self.dphi      
    def get_lonlat(self, x,y):  return (self.get_lon(x), self.get_lat(y))
    
    def get_x(self, lam):       return (lam - self.lambda0) / self.dlambda 
    def get_y(self, phi):       return (phi - self.phi0   ) / self.dphi   
    def get_xy(self, lam, phi): return (self.get_x(lam), self.get_y(phi))
    #
    def get_interpolation_vars(self, lam, phi):
        # --- projects outliers onto rim ---
        #     ix0,iy0,ix1,iy1 should correspondsto valid vertices
        x,y = self.get_xy(lam, phi)
        ix0 = min(max(nint(math.floor(x)), 0), self.nx-1)
        iy0 = min(max(nint(math.floor(y)), 0), self.ny-1)
        ix1 = min(ix0+1, self.nx-1)
        iy1 = min(iy0+1, self.ny-1)
        sx  = x-ix0  
        sy  = y-iy0
        return ix0, ix1, sx, iy0, iy1, sy
    def trilinear_interpolation(self, lam, phi, data, deriv=False):
        ix0, ix1, sx, iy0, iy1, sy = self.get_interpolation_vars(lam, phi)
        v00 = data[ix0,iy0]
        v01 = data[ix0,iy1]
        v10 = data[ix1,iy0]
        v11 = data[ix1,iy1]
        ax  = v10 - v00
        ay  = v01 - v00
        axy = v11 + v00 - v01 - v10
        if deriv:   # gradient in units value/kilometer
            dfdx  = ax + axy*sy
            dfdy  = ay + axy*sx
            dxdr  = 1.0/self.drdx(phi)
            dydr  = 1.0/self.drdy
            return array([dfdx*dxdr, dfdy*dydr])
        else:       # value 
            return v00 + ax*sx + ay*sy + axy*sx*sy
 
    

class GridInterpolator:
    # --------------------------------------------------------------------
    #    Super class for trilinear interpolation vertex-oriented data set
    #    on a lon-lat grid.
    #    Sub classes should set (at least) the following attributes at instantiation:
    #       grid     = LonLatGrid like instance
    #       data     = array like of data on vertices. Water depths positive in water
    #       wet      = array like of 1.0/0.0 on vertices indicating wet/dry point respectively
    #       pad      = padvalue (for vertices with no data)
    #       filename = name of file associated with data set
    #
    #    gradient returned as value/kilometer
    # --------------------------------------------------------------------
    def _print_wet_points(self, fout):
        nx,ny = self.data.shape
        for ix in range(nx):
            for iy in range(ny):
                if self.wet[ix,iy] < 1e-3:
                    fout.write("%12.7f %12.7f\n" % (self.grid.get_lon(ix), self.grid.get_lat(iy)))
        
    def is_wet(self, lon, lat):
        wetness = self.grid.trilinear_interpolation(lon, lat, self.wet)
        return wetness>0.5
    
    def is_inside_grid(self, lon, lat):
        return self.grid.inside_grid(lon, lat)
    
    def __call__(self, lon, lat, deriv=False):
        # ---- trilinear interpolation
        if self.grid.inside_grid(lon, lat):
            return self.grid.trilinear_interpolation(lon, lat, self.data, deriv)
        else:
            return None
    def get_average_depth(self, x0, x1, y0, y1, nsx, nsy):
        # ------------------------------------------------
        # evaluate the average water depth within
        #   x0 < lon < x1   and   y0 < lat < y1 
        # for by sampling (nsx, nsy) along (x,y) cell sides
        # Assign water depth = 0 for dry points, so average
        # is over full cell, not wet fraction
        # Use HP offset, so cell lines does not get biased
        # ------------------------------------------------
        xsgrid = x0 + (x1-x0)*linspace(1.0/2.0/nsx, 1.0-1.0/2.0/nsx, nsx)
        ysgrid = y0 + (y1-y0)*linspace(1.0/2.0/nsy, 1.0-1.0/2.0/nsy, nsy)
        wdsum = 0.
        for x in xsgrid:
            for y in ysgrid:
                wdsum += max(0.0, self(x,y)) # clip dry points
        return wdsum/nsx/nsy
                    
class GEBCO_GridInterpolator(GridInterpolator):
    # ---------------------------------------------------------------
    #    setup GridInterpolator for GEBCO format data sets
    #
    #    GEBCO: depths negative, heights positive, all grid points defined
    #           on a regular lon-lat grid
    #           The data start at the Northwest corner of the files and
    #           are arranged in latitudinal bands
    #    High resolution GEBCO data exceeds memory, so allow subgrid, given
    #    as (E,W,S,N) grid bounds in degrees
    # ---------------------------------------------------------------
    def __init__(self, filename, padvalue = pad_depth_invalid, (lonmin,lonmax,latmin,latmax) = Earth):
        if verbose: print "loading %s -> GEBCO_GridInterpolator" % filename
        ncfile    = NetCDF.Dataset(filename, 'r')
        xa0,xa1   = ncfile.variables['x_range'][:]    # full grid dim
        ya0,ya1   = ncfile.variables['y_range'][:]    # full grid dim
        dlon,dlat = ncfile.variables['spacing'][:]
        assert dlon > 0
        assert dlat > 0
        nxa,nya   = ncfile.variables['dimension'][:]  # full grid dim
        data = transpose(reshape(ncfile.variables['z'][:] , (nya,nxa)))
        # flip y-axis to S->N scan
        for iy in range(nya/2):
            tmp              = 1.0*data[:,iy]
            data[:,iy]       = data[:,nya-1-iy]
            data[:,nya-1-iy] = tmp
        # extract sub grid
        ix0 =     int((lonmin-xa0)/dlon)  # define sub grid indices in full grid
        ix1 = 1 + int((lonmax-xa0)/dlon)  # define sub grid indices in full grid
        iy0 =     int((latmin-ya0)/dlat)  # define sub grid indices in full grid
        iy1 = 1 + int((latmax-ya0)/dlat)  # define sub grid indices in full grid
        lon0 = xa0 + ix0*dlon             # sub grid offset
        lat0 = ya0 + iy0*dlat             # sub grid offset
        # configure object
        self.data     = -1.0*data[ix0:ix1, iy0:iy1] # truncate array to window and recast to float64 after truncation, with depths positive (minus)
        nx,ny         = data.shape
        self.grid     = LonLatGrid(nx,ny, lon0, lat0, dlon, dlat)
        self.wet      = where(self.data>0, 1.0, 0.0)
        self.pad      = padvalue
        self.filename = filename # for reference
        ncfile.close()
        



###############################################################
#                  module test section                        #
###############################################################
if __name__ == "__main__":
    ### -------------------- basic raster scan ----------------------------
    ##f = open("jj", "w")
    ##(lonmin,lonmax,latmin,latmax) = NWEurope
    ##gebco = GEBCO_GridInterpolator('gebco_08.nc')
    #### gebco = GEBCO_GridInterpolator('./GEBCO/gridone.nc')
    ##for lon in arange(lonmin, lonmax, 0.005):
    ##    for lat in arange(latmin, latmax, 0.005):
    ##         if gebco(lon, lat)<0: # dry
    ##             f.write("%12.7f  %12.7f\n" % (lon,lat))
    ##f.close()
    ### --------------- generate a water depth file -----------------------
    ##
    assert len(sys.argv) == 5
    gebco = GEBCO_GridInterpolator(sys.argv[1])
    print "reading topography from     %s" % sys.argv[1]
    fin  = open(sys.argv[2], "r")
    print "reading wet cells from      %s" % sys.argv[2]
    fout = open(sys.argv[3], "w")
    print "printing cell bathymetry to %s" % sys.argv[3]
    nsampl = int(sys.argv[4])
    #
    lon1, dlon = map(float, fin.readline().split())
    lat1, dlat = map(float, fin.readline().split())
    print "SW corner (lon1,lat1)       = %12.5f %12.5f" % (lon1, lat1)
    print "longitude grid step dlon    = %12.5f" % dlon
    print "latitude  grid step dlat    = %12.5f" % dlat
    print "sub grid sampling divisions = %d" % nsampl
    #
    npt = 0
    wmin =  1.0e20
    wmax = -1.0e20
    wmin_acceptable = 0.1  # minimum acceptable water depth for a wet point
    # rely on fortran indexing starting at 1
    for line in fin.readlines():
        ix,iy = map(int, line.split())[:2]  # ignore possible prepended information
        x0 = lon1 + dlon*(ix-1.5)
        x1 = lon1 + dlon*(ix-0.5)
        y0 = lat1 + dlat*(iy-1.5)
        y1 = lat1 + dlat*(iy-0.5)
        wd = gebco.get_average_depth(x0, x1, y0, y1, nsampl, nsampl) # assume grid close to isotropic
        wd = max(wmin_acceptable, wd)
        wmin = min(wmin, wd)
        wmax = max(wmax, wd)
        fout.write("%d  %d  %12.7f\n" % (ix,iy,wd))
        npt += 1
    print "evaluated %d grid cells" % npt
    print "probed depth range: %12.7f < wd < %12.7f" % (wmin,wmax)
    fout.close()
    fin.close()
    #
