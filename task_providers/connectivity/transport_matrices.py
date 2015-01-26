#!/usr/bin/env python
# -----------------------------------------------------------------
# Basic uptake of transport matrices dumped by module connectivity
# TODO def plot_sink_index(self, landboxes=None, bbox=None): prior ...
#      def plot_dist_tmat(self ...)
#      def cluster analysis
# ----------------------------------------------------------------

verbose = False
from Scientific.IO.NetCDF import *
from numpy  import *
from string import * # aqquire split
import os
import sys

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

thisdir = os.path.dirname(__file__)  # allow remote import

earth_radius  = 6371.0   # Earth average radius in km


def nice_matrix(m, fmt="%12.7f"):
    # -------------------------------------------
    # convert matrix m to a readable string
    # formatted by fmt
    # -------------------------------------------
    res = ""
    for row in m:
        for item in row: res += fmt % item
        res += "\n"
    return res


# ------------ bounding box functionality ------------

def get_common_bbox(bbox1, bbox2):
        xmin = min(bbox1[0], bbox2[0])
        xmax = max(bbox1[2], bbox2[2])
        ymin = min(bbox1[1], bbox2[1])
        ymax = max(bbox1[3], bbox2[3])
        return (xmin, ymin, xmax, ymax)

def get_center_point(bbox):
	return 0.5*(bbox[0]+bbox[2]), 0.5*(bbox[1]+bbox[3])

_init_bbox =  (sys.float_info.max, sys.float_info.max, -sys.float_info.max, -sys.float_info.max) # bbox initialization constant


def lonlat_dist(lon0,lat0, lon1,lat1):
	# -------------------------------------------------------
	# Compute distance between (lon0, lat0) and (lon1, lat1)
	# point by point is slow ...
	# -------------------------------------------------------
	lam0 = lon0*pi/180.
	lam1 = lon1*pi/180.
	phi0 = lat0*pi/180.
	phi1 = lat1*pi/180.
	rho  = cos(lam0)*cos(phi0)*cos(lam1)*cos(phi1) + sin(lam0)*cos(phi0)*sin(lam1)*cos(phi1) + sin(phi0)*sin(phi1)
        rho = min(max(rho, 0.0), 1.0)
	return earth_radius*arccos(rho)



pie_colors = ['red','blue','green','yellow','magenta','purple']
 
def draw_pie(ax,ratios,X,Y,size):
	# -------------------------------------------------------
	# Draw pie chart from scratch
	# Code from
	#   http://www.geophysique.be/2010/11/15/matplotlib-basemap-tutorial-05-adding-some-pie-charts/
	# Adaptations: Added assertion that sum(ratios) ~ 1
	#              Applied trigonometry/pi from numpy corresponding to all-import 
	# Arguments:
	#   ax:     subplot, where pie chart should be added
	#   ratios: pie ratios
	#   X,Y:    pie center in plot
	#   size:   pie size
	# -------------------------------------------------------
	N = len(ratios)
	assert abs(sum(ratios)-1) < 1e-4
	xy = []
	start = 0.
	for ratio in ratios:
		x = [0] + cos(linspace(2*pi*start,2*pi*(start+ratio), 30)).tolist()
		y = [0] + sin(linspace(2*pi*start,2*pi*(start+ratio), 30)).tolist()
		xy1 = zip(x,y)
		xy.append(xy1)
		start += ratio
	for i, xyi in enumerate(xy):
		ax.scatter([X],[Y] , marker=(xyi,0), s=size, facecolor=pie_colors[i], zorder=100)



class TransportMatrix:
    # ==========================================================
    # Interface to a time series of transport matrices
    # Main attributes:
    #    tmat[itime,idest,isrc]: transport probability isrc -> idest
    #                            at itime
    #    time[itime]           : time value associated with tmat[itime]
    #    [destinations]        : list of destination objects (attached after instantiation)
    #    [sources]             : list of sources objects  (attached after instantiation)
    #    [parent]              : if instance is a aggregation of an underlying set
    #                            of transport matrices
    # ==========================================================
    def __init__(self, fname=None):
        # ---------------------------------------------
        # Compartible with fortran module connectivity 
        # ---------------------------------------------
        if fname is not None:
            f = NetCDFFile( os.path.join(thisdir, fname), "r")
            self.tmat = f.variables['tmat'].getValue()                   # [ntime,ndest,nsource]
            self.time = f.variables['time'].getValue().tolist()          # [ntime] - make index() available by list cast
            f.close()
            assert len(self.time) == len(self.tmat)
        
    def __get_item__(self, time):
        try:
            it = self.time.index(time)
        except ValueError:
            raise IndexError
        return self.tmat[it]


    def aggregate(self, dest_groups, src_groups):
        # ------------------------------------------------------
        # Create a new TransportMatrix instance corresponding to
        # aggregated areas. destinations/sources are not set, but 
        # time vector is inherited
        # ------------------------------------------------------
        # --- create aggregated transport matrix ---
        nt = len(self.time)
        nd = len(dest_groups)
        ns = len(src_groups)
        tagg = zeros((nt,nd,ns), float)
        for it in range(nt):
            for (ide, degr) in enumerate(dest_groups):
                for (isc, scgr) in enumerate(src_groups):
                    tagg[it, ide,isc] = sum(take(take(self.tmat[it], degr, axis=0), scgr, axis=1))/len(scgr)
        # --- configure child instance ---
        child = TransportMatrix() # clean object
        child.tmat   = tagg
        child.time   = self.time # no copy, shared object
        child.parent = self
        return child

    def get_survival_by_source(self, frame=-1):
        return sum(self.tmat[frame], axis=0)
    def get_sinkiness_by_destination(self, frame=-1):
        return sum(self.tmat[frame], axis=1)

    
    def plot_survival_by_source(self, landboxes=None, bbox=None, colmap='jet', frame=-1):
        # ------------------------------------------------------
        # Plot survival probability by source
        #
        # If landboxes is provided, they will be plottet along
        # If bbox is provided, this will be used, otherwise extracted
        # colmap is the color map
        # frame is the time frame to plot
        #
        # Attribute sources must have been set
        # ------------------------------------------------------
        fig      = plt.figure() 
        ax       = fig.add_subplot(111)
        #
        if landboxes is not None:
            actual_bbox = self._add_landboxes(ax, landboxes) 
        else:
            actual_bbox = _init_bbox # out of range constant
        #
        cmap     = cm.get_cmap(colmap) # or any other one
        nti, nde, nsc = self.tmat.shape
        ts = self.get_survival_by_source(-1) # last frame
	tsmin = ts.min() 
	tsmax = ts.max()
	norm = matplotlib.colors.Normalize(tsmin, tsmax) # min_val, max_val
        for isc in range(nsc):
            xsw, ysw, xne, yne = self.sources[isc]
            actual_bbox = get_common_bbox(actual_bbox, (xsw, ysw, xne, yne))
            dx   = xne - xsw
            dy   = yne - ysw
            rect = matplotlib.patches.Rectangle((xsw, ysw), dx, dy, color=cmap(norm(ts[isc])))
            ax.add_patch(rect)
        if bbox is None: # apply generated bbox
            plt.xlim([actual_bbox[0],actual_bbox[2]])
            plt.ylim([actual_bbox[1],actual_bbox[3]])
        else:
            plt.xlim(bbox[0], bbox[2])
            plt.ylim(bbox[1], bbox[3])
	#
	# --- set colorbar ---
	ncolors = 20 # number of color in color map, could be an input parameter
	cmmapable = cm.ScalarMappable(norm, cmap)
	cmmapable.set_array( linspace(tsmin, tsmax, ncolors) )
	plt.colorbar(cmmapable)
	
	# --- render plot ---
        plt.show()


    def plot_sinkiness_by_destination(self, landboxes=None, bbox=None, colmap='jet', frame=-1):
        # ------------------------------------------------------
        # Plot sinkiness probability by destination
        #
        # If landboxes is provided, they will be plottet along
        # If bbox is provided, this will be used, otherwise extracted
        # colmap is the color map
        # frame is the time frame to plot
        #
        # Attribute destinations must have been set
        # ------------------------------------------------------
        fig      = plt.figure() 
        ax       = fig.add_subplot(111)
        #
        if landboxes is not None:
            actual_bbox = self._add_landboxes(ax, landboxes) 
        else:
            actual_bbox = _init_bbox # out of range constant
        #
        cmap     = cm.get_cmap(colmap) # or any other one
        nti, nde, nsc = self.tmat.shape
        sink = self.get_sinkiness_by_destination(-1) # last frame
	sinkmin = sink.min() 
	sinkmax = sink.max()
	norm = matplotlib.colors.Normalize(sinkmin, sinkmax) # min_val, max_val
        for ide in range(nde):
            xsw, ysw, xne, yne = self.destinations[ide]
            actual_bbox = get_common_bbox(actual_bbox, (xsw, ysw, xne, yne))
            dx   = xne - xsw
            dy   = yne - ysw
            rect = matplotlib.patches.Rectangle((xsw, ysw), dx, dy, color=cmap(norm(sink[ide])))
            ax.add_patch(rect)
        if bbox is None: # apply generated bbox
            plt.xlim([actual_bbox[0],actual_bbox[2]])
            plt.ylim([actual_bbox[1],actual_bbox[3]])
        else:
            plt.xlim(bbox[0], bbox[2])
            plt.ylim(bbox[1], bbox[3])
	#
	# --- set colorbar ---
	ncolors = 20 # number of color in color map, could be an input parameter
	cmmapable = cm.ScalarMappable(norm, cmap)
	cmmapable.set_array( linspace(sinkmin, sinkmax, ncolors) )
	plt.colorbar(cmmapable)
	
	# --- render plot ---
        plt.show()
	
    def _add_landboxes(self, ax, landboxes, landcol='yellow'):    
        # ------------------------------------------------------
        # Plot add landboxes (BoxSet like instance) to
        # subplot handle ax with color landcol
        # ------------------------------------------------------
        for isc in range(len(landboxes)):  
            xsw, ysw, xne, yne = landboxes[isc]
            dx   = xne - xsw
            dy   = yne - ysw
            rect = matplotlib.patches.Rectangle((xsw, ysw), dx, dy, color=landcol)
            ax.add_patch(rect)
        return landboxes.bbox


    def _get_dist_tmat_lists(self, tmin, frame=-1, dest="all", src="all"):
        # ------------------------------------------------------
	# Resolve lists of tmat[i,j], dist[i,j] for all
	# elements where tmat > tmin
	# dest/src are optional indices of groups (default all destinations/sources)
	# ------------------------------------------------------
	nti, nde, nsc = self.tmat.shape
	tlist = []
	dlist = []
	if dest   == "all":
		destset = range(nde)
	else:
		destset = dest
	if src == "all":
		srcset = range(nsc)
	else:
		srcset = src
	for ide in destset:
	        for isc in srcset:
			t = self.tmat[frame,ide,isc]
			if t>tmin:
				tlist.append(t)
				lon0,lat0 = get_center_point(self.sources[isc])
				lon1,lat1 = get_center_point(self.destinations[ide])
				dlist.append( lonlat_dist(lon0,lat0, lon1,lat1) )
	return dlist, tlist

    def plot_dist_tmat(self, tmin, frame=-1, dest="all", src="all"):
        # ------------------------------------------------------
        # Make scatter plot of (distances, transport_probability)
	# Attributes destinations and sources must have been set
	# dest/src are optional indices of groups (default all destinations/sources)
        # ------------------------------------------------------
        dlist, tlist = self._get_dist_tmat_lists(tmin, frame, dest, src)
	plt.scatter(dlist, tlist)
	plt.xlabel('distance[km]')
	plt.ylabel('transport probability')
	plt.show()

    def plot_dist_histogram(self, tmin, frame=-1, dest="all", src="all"):
        # ------------------------------------------------------
        # Make scatter plot of (distances, transport_probability)
	# Attributes destinations and sources must have been set
	# dest/src are optional indices of groups (default all destinations/sources)
        # ------------------------------------------------------
        dlist, tlist = self._get_dist_tmat_lists(tmin, frame,  dest, src)
	plt.hist(dlist)
	plt.xlabel('distance[km]')
	plt.ylabel('frequency')
	plt.show()



class BoxSet(list):
    # ==========================================================
    # Represent and ordered list of lon-lat boxes
    # Attributes:
    #    self[ibx]:  box ibx as (SWlon, SWlat, NElon, NElat)
    #    tags[ibx]:  optional tag for box ibx (default is "")
    #    bbox:       minimal bounding box for box set as (SWlon, SWlat, NElon, NElat)
    # ========================================================== 
    def __init__(self, fname):
        # ---------------------------------------------------
        # File fname contains lines of lon-lat box corners
        #   SWlon, SWlat, NElon, NElat [,tag]
        # where order of appearence determines sequence indices of the
        # boxes; tag is optional, any elements after tag are ignored
        # Minimal bounding box is generated during initialization
        # ---------------------------------------------------
        f = open(fname, "r")
        self.tags = []
        self.bbox = _init_bbox
        for line in f.readlines():
            items = line.split()
            self.append(map(float, items[:4]))
            if len(items)>4: # ignore 5+ elements
                self.tags.append(items[4]) # keep as string
            else:
                self.tags.append("")
            self.bbox = get_common_bbox(self.bbox, self[-1])
        f.close()

	
    def get_groups(self, tagmap = lambda x:x):
	# ---------------------------------------------------
        # Return a mapping of lon-lat boxes
	#    {tag1: indices_of_tag1,
	#     tag2: indices_of_tag2, ...}
	# where tags corresponds to items in tagmap(self.tags)
	# By default transformation tagmap is the identity map
	# indices (in self) are those boxes whose tagmap(fulltag) maps to tag key
	# It is not assumed that tags appear consequetively in self.tags
	# tag tested identity by == operator on tagmap(fulltag)
        # ---------------------------------------------------    
	groupmap = {}
	for it,fulltag in enumerate(self.tags):
            tag = tagmap(fulltag) # possibly same
	    if not tag in groupmap.keys():
	         groupmap[tag] = []
	    groupmap[tag].append(it)
	return groupmap

    def get_group(self, tag, tagmap = lambda x:x):
	# ---------------------------------------------------
        # Return indices of boxes with tag==tagmap(fulltag) in self.tags
	# tag tested identity by == operator on tagmap(fulltag)
	# By default transformation tagmap is the identity map
        # ---------------------------------------------------
	indices = []
	for it in range(len(self.tags)):
            if tag == tagmap(self.tags[it]):
	         indices.append(it)
	return indices # probably empty
	 
        

if __name__ == "__main__":
    tmat         = TransportMatrix("../../runs/connect.nc")
    #print min(tmat.tmat.flat), max(tmat.tmat.flat)
    #stop
    landboxes         = BoxSet("../../../Habitats/land.bbx")
    tmat.sources      = BoxSet("../../../Habitats/spawning_habitats.bbx")
    tmat.destinations = BoxSet("../../../Habitats/settlement_habitats.bbx")
    #ts = tmat.get_survival_by_source(-1) # last frame   
    #for i,t in enumerate(ts):
    #    print i,t
    # --- plot variants ---
    #tmat.plot_survival_by_source(landboxes, bbox=(-2,53,15,59))
    tmat.plot_sinkiness_by_destination(landboxes, bbox=(8,55,15,58))
    #tmat.plot_dist_tmat(1e-5)
    #tmat.plot_dist_histogram(1e-5)
    #jutland zeeland sweden jammerbugt
    #tmat.plot_dist_histogram(1e-5, dest = tmat.destinations.get_group("jammerbugt", lambda x: x.split("_")[0]))   
    #
