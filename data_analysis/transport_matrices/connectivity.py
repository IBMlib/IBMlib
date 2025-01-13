#!/usr/bin/env python
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
from   numpy import *
import netCDF4 as netcdf
import sys
from   mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt



# check minimal consistency of aggregates:
# 1) all groups contain at least one mmeber
# 2) all indices are integers
# 3) that all indices are positive and less than upper limit imax
# it is not required that all parent nodes should be included, and
# nodes can appear in several groups
#
# return True - or - False, message
def check_grouping(groups, imax):
    # flatten groups ()
    indices = []
    for subgr in groups:
        if len(subgr)==0:
            return False, "groups contain at least one empty subgroup"
        indices.append( list(subgr) )
    n = len(indices)
    # check index bounds
    if min(indices) < 0:
        return False, "offset < 0"
    if max(indices) >= imax:
        return False, "highest index invalid"
    # check all is integers
    for i in range(n):
        if not isinstance(indices[i],int):
            return False, "groups contain non integer"
    return True # all tests passed

# avoid repetition for static default option
def get_option_or_set_default(dict, key, default):
    if key in dict:
        return dict[key]
    else:
        return default
        #


# ===========================================================================
# Represent, transform, analyze and plot a spatial connectivity matrix
# with elements bounded in [0,1]
# 
# connectivity matrix without spatial association can be plotted as
# matshow(cmat) or imshow(cmat)
# ===========================================================================
class connectivity_matrix:
    # ----------------------------------------------
    # lon: longitude of nodes
    # lat: latitude  of nodes
    # cmat: connectivity matrix as (ndest, nsource)
    # ----------------------------------------------
    def __init__(self, lon, lat, cmat):
        self.lon  = 1.0*array(lon, float)  # force copy+array
        self.lat  = 1.0*array(lat, float)  # force copy+array
        self.cmat = 1.0*array(cmat, float) # force copy+array (ndest, nsource)
        assert cmat.shape == (len(lon), len(lon))
    # ----------------------------------------------
    # create a new connectivity matrix based on probability interpretation
    # aggmat = avg_sources( sum_destina( cmat(destina, sources) ))
    # groups is a list of lists containing indices of aggregates
    # number of aggregates in groups can be different
    # ----------------------------------------------
    def aggregate(self, groups, centers="avg", check=True):
        n       = len(groups) # number of aggregations
        aggcmat = zeros((n,n), float)
        agglon  = zeros(n, float)
        agglat  = zeros(n, float)
        if check:
            assert check_grouping(groups, n)
        # --- create aggregation centers
        if centers != "avg":
            raise NotImplementedError("centers option %s" % str(centers))
        for igr in range(n):
            agglon[igr] = mean(take(self.lon, groups[igr]))
            agglat[igr] = mean(take(self.lat, groups[igr]))
        # --- create aggregated connectivity
        for ides in range(n):
            for isrc in range(n):
                submat = take(take(self.cmat, groups[ides], axis=0), groups[isrc], axis=1)
                aggcmat[ides,isrc] = sum(submat.flat)/len(groups[isrc]) # avg over sources
        #
        return connectivity_matrix(agglon, agglat, aggcmat)
    
    # ====================== plotters ======================
    
    # common initializations for Basemap plot
    #   pltwindow = (W,S,E,N) - default lon/lat limits +/- 1
    def init_basemap_plot(self, **args):  
        if "pltwindow" in args:
            (W,S,E,N) = args["pltwindow"]
        else:
            (W,S,E,N) = (amin(self.lon)-1, amin(self.lat)-1,
                         amax(self.lon)+1, amax(self.lat)+1) # default
        #
        self.bmap = Basemap(projection='cyl',
                            llcrnrlon=W,
                            llcrnrlat=S,
                            urcrnrlon=E,
                            urcrnrlat=N,
                            resolution='h')
        # draw coastlines, country boundaries, fill continents.
        meridians = arange(0,360,2)
        parallels = arange(-90,90,2)
        self.bmap.drawmeridians(meridians, labels=[True,False,False,True])  # labels = [left,right,top,bottom]
        self.bmap.drawparallels(parallels, labels=[True,False,False,True])  # labels = [left,right,top,bottom]
        self.bmap.fillcontinents(color='burlywood',lake_color='grey')
        self.bmap.drawcoastlines(linewidth=0.25)
        #
        if "add_topo" in args and args["add_topo"]: # left-to-right evaluation
            self.bmap.etopo()
        

    # common wrapups for Basemap plot
    def wrapup_basemap_plot(self, **args):
        if "title" in args:
            plt.title(args["title"], fontsize = 30)
        #
        if "colorbar" in args and args["colorbar"]:
            self.bmap.colorbar(self.pobj)     
        #
        plt.tight_layout()
        if "saveas" in args:
            plt.savefig(args["saveas"], dpi=300)
        else:
            plt.show()
        #

    # -------------------------------------------------------    
    # Plot connectivity as lines between nodes
    # with relative thickness set by symmetrized connectivity
    # -------------------------------------------------------    
    def bmapplot_symmetrized_connectivity(self, **args):
        self.init_basemap_plot(**args)  # create self.bmap
        lth = get_option_or_set_default(args, "line_thickness", 3)
        #
        csym  = 0.5*(self.cmat + transpose(self.cmat))
        csym *= lth/amax(triu(csym, k=1)) # normalize range to [0;lth]
        n     = len(self.lon)
        for ides in range(n):
            x2,y2 = self.lon[ides],self.lat[ides]
            for isrc in range(ides+1,n):
                x1,y1 = self.lon[isrc],self.lat[isrc]
                self.bmap.plot([x1,x2],[y1,y2], 'k-',linewidth=csym[ides,isrc])
        #
        # --- plot nodes ---
        node_size  = get_option_or_set_default(args, "node_size", 3)
        node_color = get_option_or_set_default(args, "node_color", "k")
        #
        for ides in range(n):
            self.bmap.scatter(self.lon[ides],self.lat[ides], s=node_size, c=node_color)
        #
        self.wrapup_basemap_plot(**args)
       
    # -------------------------------------------------------    
    # Plot survival as colored nodes
    # -------------------------------------------------------    
    def bmapplot_survival(self, **args):
        self.init_basemap_plot(**args) # create self.bmap
        w    = sum(self.cmat, axis=0)  # survival/sourciness
        node_size  = get_option_or_set_default(args, "node_size", 3)
        self.pobj = self.bmap.scatter(self.lon, self.lat, c=w, cmap='coolwarm', s=node_size)
        self.wrapup_basemap_plot(**args)
    # -------------------------------------------------------    
    # Plot sinkiness as colored nodes
    # -------------------------------------------------------    
    def bmapplot_sinkiness(self, **args):
        self.init_basemap_plot(**args) # create self.bmap
        w    = sum(self.cmat, axis=1)  # sinkiness
        node_size  = get_option_or_set_default(args, "node_size", 3)
        self.pobj = self.bmap.scatter(self.lon, self.lat, c=w, cmap='coolwarm', s=node_size)
        self.wrapup_basemap_plot(**args)
    # -------------------------------------------------------    
    # Plot retention as colored nodes
    # -------------------------------------------------------    
    def bmapplot_retention(self, **args):
        self.init_basemap_plot(**args) # create self.bmap
        w    = diag(self.cmat)         # retention
        node_size  = get_option_or_set_default(args, "node_size", 3)
        self.pobj = self.bmap.scatter(self.lon, self.lat, c=w, cmap='coolwarm', s=node_size)
        self.wrapup_basemap_plot(**args)
            


##############################################################
if __name__ == "__main__":
    #cmat = random.random((3,3))
    #lon  = arange(8,  11)
    #lat  = arange(54, 57)
    #lon  = array([ 8,  8.5,  9.5])
    #lat  = array([54, 55.5, 54.5])
    #c1 = connectivity_matrix(lon,lat,cmat)
    #c1.bmapplot()
    nc   = netcdf.Dataset("/home/asbjorn/People/OleHenriksen/PELA/IBMlib/run_14Jun2022/tmat_2016_np1000.nc", "r")
    cmat = nc.variables["tmat"][0,:,:]  # (ndest, nsource)
    lon = []
    lat = []
    for line in open("/home/asbjorn/People/OleHenriksen/PELA/IBMlib/run_14Jun2022/NSammodyt_habitats_13Jun2022.lpn", "r").readlines():
        (x,y,no) = map(float, line.split())
        lon.append(x)
        lat.append(y)
    lon = array(lon)
    lat = array(lat)
    cob = connectivity_matrix(lon,lat,cmat)
    cob.bmapplot_symmetrized_connectivity()
    #cob.bmapplot_survival(colorbar=1)
    #cob.bmapplot_sinkiness(colorbar=1)
    #
    #tmat = nc.variables["tmat"][0,:,:]  # (ndest, nsource)
    #ndest, nsource = tmat.shape
    #for p in sum(tmat, axis=0):
    #    print(p)
    #print(sum(tmat)/len(tmat))
    #print(amax(tmat))
    #print("til  fra  P(til,fra)")
    #for ide in range(ndest):
    #    for isr in range(nsource):
    #        print(ide+1, isr+1, tmat[ide,isr])
    #nc.close()

