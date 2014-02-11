#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#
#  Provide tools to check correctness of implementation and facilitate access to
#  Atlantis input data
# 
##########################################################################################
import sys
import os

from Scientific.IO.NetCDF import *
from numpy import *


class HydroInput:
    ################################################################################
    # 
    #   Represent hydro input file and provide tools to check data
    #
    ################################################################################
    def __init__(self, ncfilename):
        self.ncfile    = NetCDFFile(ncfilename, 'r')
        self.nframes   = self.ncfile.dimensions['t']
        self.nboxes    = self.ncfile.dimensions['b']
        self.nlayers   = self.ncfile.dimensions['z']
        self.ndest     = self.ncfile.dimensions['dest']
        self.exchange  = self.ncfile.variables['exchange'].getValue()
        self.time      = self.ncfile.variables['t'].getValue()
        self.dest_b    = self.ncfile.variables['dest_b'].getValue()
	self.dest_k    = self.ncfile.variables['dest_k'].getValue()
        self.dest_fill = self.ncfile.variables['dest_k']._FillValue
        
        
    def get_exchange(self, sb, sl, tb, tl, frame=-1):
        # ----------------------------------------------------
        # extract transport (sb, sl) -> (tb, tl) counted
        # positive in this direction (source -> target)
        # ----------------------------------------------------
        xc = 0.
        # scan (sb, sl) as source
        for ide in range(self.ndest):
            if self.dest_b[frame,sb,sl,ide] == self.dest_fill: break  # end of connections, in sync with dest_k
            # --- count s->t transport positive
            if (self.dest_b[frame,sb,sl,ide] == tb) and (self.dest_k[frame,sb,sl,ide] == tl):
                xc += self.exchange[frame,sb,sl,ide]
        # scan (tb, tl) as source
        for ide in range(self.ndest):
            if self.dest_b[frame,tb,tl,ide] == self.dest_fill: break  # end of connections, in sync with dest_k
            # --- count t->s transport negative
            if (self.dest_b[frame,tb,tl,ide] == sb) and (self.dest_k[frame,tb,tl,ide] == sl):
                xc -= self.exchange[frame,sb,sl,ide]    
        return xc

    
    def get_outflow(self, sb, sl,frame=-1):
        # ----------------------------------------------------
        # extract outflow from source (sb, sl)
        # inflows counted negative
        # ----------------------------------------------------
        xc = 0.
        # --- add any outflows from source (sb, sl)
        for ide in range(self.ndest):
            if self.dest_b[frame,sb,sl,ide] == self.dest_fill: break  # end of connections, in sync with dest_k
            xc += self.exchange[frame,sb,sl,ide]
        # --- subtract any inflows to (sb, sl)
        for ibox in range(self.nboxes):
            for ilay  in range(self.nlayers):
                if (ibox == sb) and (ilay == sl): break # exclude source arg
                for ide in range(self.ndest):
                    if self.dest_b[frame,ibox,ilay,ide] == self.dest_fill:
                        break  # in sync with dest_k
                    if (self.dest_b[frame,ibox,ilay,ide] == sb) and (self.dest_k[frame,ibox,ilay,ide] == sl):
                        xc -= self.exchange[frame,ibox,ilay,ide]
        return xc


class BoxAverage:
    ################################################################################
    # 
    #   Super class for box averaged quantities like salinity and temperature
    #
    ################################################################################
    def __init__(self, ncfilename, varname):
        self.ncfile   = NetCDFFile(ncfilename, 'r')
        self.nframes  = self.ncfile.dimensions['t']
        self.nboxes   = self.ncfile.dimensions['b']
        self.nlayerss = self.ncfile.dimensions['z']  # including sediment == wet layers + 1
        self.data     = self.ncfile.variables[varname].getValue()  # shape = (t, b, z)
        self.time     = self.ncfile.variables['t'].getValue()
    def __get_item__(self, what):
        return self.data[what]


class TemperatureInput(BoxAverage):
    def __init__(self, ncfilename):
        BoxAverage.__init__(self, ncfilename, "temperature")

class SalinityInput(BoxAverage):
    def __init__(self, ncfilename):
        BoxAverage.__init__(self, ncfilename, "salinity")


    
#################### self test #########################
if __name__ == "__main__":
    hydro = HydroInput("jjh.nc")
    print hydro.get_exchange(21, 0, 22, 2)
    print hydro.get_exchange(21, 1, 22, 3)
    print hydro.get_exchange(21, 2, 22, 4)
    print hydro.get_exchange(21, 3, 22, 5)
    print
    print hydro.get_exchange(21, 3, 21, 4)
    print hydro.get_exchange(21, 2, 21, 3)
    print hydro.get_exchange(21, 1, 21, 2)
    print hydro.get_exchange(21, 0, 21, 1)
    print hydro.get_exchange(21, 1, 21, 0)
    #
    temp = TemperatureInput("jjt.nc")
    salt = SalinityInput("jjs.nc")
    print temp.data[-1][21,0:] # shape = (t, b, z)
    print salt.data[-1][21,0:] # shape = (t, b, z)
    #
