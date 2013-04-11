''' extract historical climate data from PRISM based on 
    bounding box and lat/lon range
'''
import os
import netCDF4
import numpy as np
import arcpy

# mean precip over specified bounding box, start and stop time
def mean_precip(nc,bbox=None,start=None,stop=None):
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]
    tindex0=netCDF4.date2index(start,nc.variables['time'],select='nearest')
    tindex1=netCDF4.date2index(stop,nc.variables['time'],select='nearest')
    bi=(lon>=box[0])&(lon<=box[2])
    bj=(lat>=box[1])&(lat<=box[3])
    p=nc.variables['precip_mean'][tindex0:tindex1,bj,bi]
    latmin=np.min(lat[bj])
    p=np.mean(p,axis=0)
    lon=lon[bi]
    lat=lat[bj]
    return p,lon,lat

url='http://geoport.whoi.edu/thredds/dodsC/prism3/monthly'
box = [-102,36.5,-101,37]  # Bounding box for Texas County, Oklahoma
nc=netCDF4.Dataset(url)
z,lon,lat=mean_precip(nc,bbox=box,start=datetime.datetime(1936,11,1,0,0),stop=datetime.datetime(1937,4,1,0,0))
latmin=np.min(lat)
lonmin=np.min(lon)
latmin=np.min(lat)
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
xyOrig = arcpy.Point(float(lonmin),float(latmin))

# create Arc Raster
arcpy.workspace  = "C:\\workspace"
arcpy.env.overwriteOutput = True
rasterName = "precip"
outRaster = os.path.normpath(os.path.join(arcpy.workspace,rasterName))
grid1=arcpy.NumPyArrayToRaster(z,xyOrig,dx,dy)
grid1.save(os.path.join(arcpy.workspace,outRaster))
strPrj = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID"\
         "['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],"\
         "UNIT['Degree',0.0174532925199433]]"
arcpy.DefineProjection_management(outRaster,strPrj)
print "Written: %s" % grid1
arcpy.AddMessage("Written: %s" % grid1)