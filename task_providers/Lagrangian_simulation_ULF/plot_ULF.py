#!/usr/bin/env python
# -------------------------------------------------------------------
# 
# -------------------------------------------------------------------
from mpl_toolkits.basemap import Basemap
import numpy as np 
import matplotlib.pyplot as plt
from   scipy.io import netcdf
# --- load topography data ---
ncfile = netcdf.NetCDFFile("result.nc", "r")
lon = ncfile.variables["topolon"][:]
print lon; stop
lat = ncfile.variables["topolat"][:]
z   = ncfile.variables["topography"][:,:]

# --- load particle data ---
x  = ncfile.variables["lon"][:,:]   # (nframes, nparticles)
y  = ncfile.variables["lat"][:,:]   # (nframes, nparticles
x0 = x[2,:]
y0 = y[2,:]
xe = x[-1,:]
ye = y[-1,:]

	
fig = plt.Figure()
# ------ basic regional map ------
bmap = Basemap(projection='cyl',
               llcrnrlon=lon[0],
               llcrnrlat=lat[0],
               urcrnrlon=lon[-1],
               urcrnrlat=lat[-1],
               resolution='h')
# draw coastlines, country boundaries, fill continents.
bmap.drawcoastlines(linewidth=0.25)
bmap.drawcountries(linewidth=0.25)
bmap.drawmeridians(np.arange(0,360,1))
bmap.drawparallels(np.arange(-90,90,0.5))
bmap.fillcontinents(color='coral',lake_color='aqua')


# add color plot
xt, yt = np.meshgrid(lon,lat,indexing="ij")
plt.pcolor(xt, yt, z)
plt.colorbar()
plt.scatter(x0,y0,color='black')
plt.scatter(xe,ye,color='red')

plt.show()
#plt.savefig("avgzoo.pdf")
