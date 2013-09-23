#
# NetCDF4-Python test to read gridded data via OPeNDAP and create Arc Raster
#
# Rich Signell, original coding
# Curtis Price, small tweaks and doc
# Rich Signell, added plotting, checking for uniform grid spacing and
#               generalized to handle any bathy/topo data at:
#               http://geoport.whoi.edu/thredds/bathy_catalog.html
#
# Prerequisite: Christoph Gohlke's netCDF4-Python library 
#               specially compiled for ArcGIS 10 at:
#   http://www.lfd.uci.edu/~gohlke/pythonlibs/
#     netCDF4-0.9.7-ArcGIS-10.0.win32-py2.6.exe (ArcGIS 10.0)
#     netCDF4-0.9.7-ArcGIS-10.1.win32-py2.7.exe (ArcGIS 10.1 beta)

import os
import netCDF4
import numpy as np


# extract data from OPeNDAP or NetCDF URL 
# url='c:/rps/python/traster_RasterToNetCDF.nc'
# Test Latitude coordinate decreasing:
url='http://geoport.whoi.edu/thredds/dodsC/bathy/crm_vol1.nc'
# Test Latitude coordinate increasing:
#url='http://geoport.whoi.edu/thredds/dodsC/bathy/etopo1_bed_g2'
# Test failure because lat is not uniformly spaced:
#url='http://geoport.whoi.edu/thredds/dodsC/bathy/smith_sandwell_v11'
box = [-71.4,41,-70.2,42]
#box=[-71.0,41.0,-67.0,46.0]

nc=netCDF4.Dataset(url)
print "Source name: %s" % nc.title
lon=nc.variables['lon'][:]
lat=nc.variables['lat'][:]
bi=(lon>=box[0])&(lon<=box[2])
bj=(lat>=box[1])&(lat<=box[3])
z=nc.variables['topo'][bj,bi]
nc.close()

lonmin=np.min(lon[bi])
latmin=np.min(lat[bj])
dx=np.diff(lon)
dy=np.diff(lat)
# check if dx or dy vary by more than one percent
assert np.abs(np.ptp(dx)/np.mean(dx))<=0.01,'longitude spacing is not uniform'
assert np.abs(np.ptp(dy)/np.mean(dy))<=0.01,'latitude spacing is not uniform'
if dy[0]>0:  # lat increasing
    z=np.array(z[::-1,:])
if dx[0]<0:  # lon decreasing
    z=np.array(z[:,::-1])    
dx=np.abs(np.mean(dx))
dy=np.abs(np.mean(dy))

import matplotlib.pylab as plt
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.1, 0.15, 0.8, 0.8])
pc = ax.pcolormesh(lon[bi], lat[bj], z,
                   vmin=np.min(z), vmax=0,
                   cmap=plt.cm.RdBu_r)
cbax = fig.add_axes([0.1, 0.08, 0.8, 0.02])
cb = plt.colorbar(pc, cax=cbax, 
                  orientation='horizontal')
cb.set_label('Depth [m]')
ax.set_aspect(1.0/np.cos(latmin * np.pi / 180.0))
plt.show()

try:
    import arcpy
    xyOrig = arcpy.Point(float(lonmin),float(latmin))

    # create Arc Raster

    arcpy.workspace  = "C:\\workspace"
    arcpy.env.overwriteOutput = True
    rasterName = "traster"
    outRaster = os.path.normpath(os.path.join(arcpy.workspace,rasterName))
    grid1=arcpy.NumPyArrayToRaster(z,xyOrig,dx,dy)
    grid1.save(os.path.join(arcpy.workspace,outRaster))
    strPrj = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID"\
             "['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],"\
             "UNIT['Degree',0.0174532925199433]]"
    arcpy.DefineProjection_management(outRaster,strPrj)
    print "Written: %s" % grid1
    arcpy.AddMessage("Written: %s" % grid1)
    # plot using Matplotlib
except:
    print "arcpy not found"


    

