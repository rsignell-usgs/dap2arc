# ugrid2arcShape.py 
# Write unstructured grid model output to ESRI Shapefile

# Uses NetCDF4-python to access remote data via OPeNDAP
# Christoph Gohlke has Windows installers for ArcGIS 10 users at
# http://www.lfd.uci.edu/~gohlke/pythonlibs/#netcdf4

# Uses Pyshp from http://code.google.com/p/pyshp/

# Instructions: 
# 1. In script below: Modify the URL, the output shapefile name, the variables you want, and (if appropriate) the time step
#    If unsure of valid variable names and time steps, you can add an ".html" on the 
#    URLs below to see the dataset info in your browser.

# Rich Signell (rsignell@usgs.gov)

import numpy as np
import netCDF4
import datetime
import shapefile
#######################################
# MAKE EDITS HERE

file_point='c:/rps/python/shapefiles/vdatum_gom_mean_point'

# Specify the URL


# adcirc vdatum tides (steady component)
url='http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/vdatum_gulf_of_maine/adcirc54_38_orig.nc'

# processing starts here

# Open OpenDAP URL to get remote triangular mesh + data

nc=netCDF4.Dataset(url)

# read lon,lat 
lon=nc.variables['lon'][:]
lat=nc.variables['lat'][:]
nnodes=len(lon)

# read connectivity array
nv=nc.variables['ele'][:,:]
[ncells,three]=np.shape(nv)
nv=nv-1  # python is 0-based


# now read data at nodes:
# read water depth at nodes
h=nc.variables['depth'][:]

# find time index to read
#hours_from_now=0   # Examples: 0=>nowcast, 3 => forecast 3 hours from now, etc. 
#date=datetime.datetime.utcnow()+datetime.timedelta(0,3600*hours_from_now)  
#date=datetime.datetime(2011,9,9,17,00)  # specific time (UTC)
#tindex=netCDF4.date2index(date,nc.variables['time'],select='nearest')
tindex=0
# read water level at nodes at a specific time
#z=nc.variables['zeta'][-1,:]   # -1 is the last time step of the forecast
#z=nc.variables['zeta'][tindex,:]  # index for date specified above

# read significant wave height at nodes at specified time step 
#z=nc.variables['hs'][tindex,:]   # Note: 'hs' is only in GOM3 wave model

# write Shapefile using pyshp

# Test 1. Write polygons for each triangle, and create record
# values for each triangle that are the average of the 3 nodal values
# (since depth and water level are defined at nodes)


ilev=-1  # [0] surface, [-1] bottom
u=nc.variables['u_amp'][0,:,0]
v=nc.variables['v_amp'][0,:,0]
ang=np.angle(u + v*1j)*180/np.pi
spd=np.abs(u + v*1j)

# Test 2. Write points for the center of each node, and create
# record values with u,v velocity components (since u,v defined at cell centers)

w = shapefile.Writer(shapefile.POINT)
w.field('ang(math)','F',8,3)
w.field('speed(m/s)','F',8,3)
#w.field('u(m/s)','F',8,3)
#w.field('v(m/s)','F',8,3)
for i in range(nnodes):
   w.point(lon[i],lat[i])
   w.record(ang[i],spd[i])
#   w.record(ang[i],spd[i],u[i],v[i])
 
w.save(file_point)

# create the PRJ file
prj = open("%s.prj" % file_point, "w")
epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
prj.write(epsg)
prj.close()