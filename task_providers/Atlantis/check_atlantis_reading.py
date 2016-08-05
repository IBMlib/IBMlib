#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#
#  Provide tools to check correctness of implementation and facilitate access to
#  Atlantis input data
# 
##########################################################################################
from atlantis_input import *
from exceptions import *

class ReadExchanges:
    # -----------------------------------------------------------
    # Read from Atlantis text logging
    # add up exchanges if interface has multiple segments
    # -----------------------------------------------------------
    def __init__(self, fname):
        data = open(fname, "r").readlines()
        self.data = {}
        # b - source box
        # k - source layer
        # d - destination layer # (from 1-10)
        # bb - sink box
        # kk - sink layer
        # -1 for lack of exchange
        # 0                   1                  2     3     4     5       6 
        # Time: 0.000000e+00, amt: 2.930815e+08, b: 0, k: 0, d: 0, bb:  0, kk: 1
        for line in data[1:]: # skip header
            items = line.split(',')
            #print items
            b     = int(items[2].split(':')[-1])
            k     = int(items[3].split(':')[-1])
            bb    = int(items[5].split(':')[-1])
            kk    = int(items[6].split(':')[-1])
            ext   = float(items[1].split(':')[-1])
            if bb != -1:                       # only set exchange for existing connection
                xkey = (b,k,bb,kk)
                if self.data.has_key(xkey):
                    self.data[xkey] += ext # add up exchanges if interface has multiple segments
                else:
                    self.data[xkey] = ext  # possibly zero
        
    def get_exchange(self,i0,i1,i2,i3,ifr):
        if ifr>0: raise InputError
        key = (i0,i1,i2,i3)
        if self.data.has_key(key):
            return self.data[key]
        else:
            return 0.0
        
rex = ReadExchanges("/home/asbjorn/ModellingPackages/Atlantis/hydrodynamics_from_IBMlib/flux_problem/initial_exchanges_read_by_atlantis.txt")
hy = HydroInput("~/ModellingPackages/Atlantis/hydrodynamics_from_IBMlib/trunk/results/hydro_14Feb2014.nc")
nkeys = len(rex.data.keys())
miss  = 0
for key in rex.data.keys():
    x1 = rex.get_exchange(*(key+(0,)))
    x2 = hy.get_exchange(*(key+(0,)))
    diff = abs(2*(x1-x2)/(1e-4 + x1 + x2))
    if diff>1e-4:
        miss += 1
        #print miss, x1/x2, x1,x2
        print key, diff
    #print key, diff
#print nkeys, miss
## ifr = 0
## print ifr, hy.get_exchange(2, 0, 2, 2, ifr), hy.get_exchange(2, 0, 2, 1, ifr)
## #
## print hy.exchange.shape

## for ifr in range(hy.exchange.shape[0]):
##     print ifr, hy.get_exchange(2, 0, 2, 2, ifr), hy.get_exchange(2, 0, 2, 1, ifr)


