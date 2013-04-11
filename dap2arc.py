"""
Functions for demonstrating OPeNDAP data access in ArcMAP10+

This module provides a number of functions to show ArcMAP users
how to access OPeNDAP datasets using the NetCDF4-Python module.

Functions
--------------
- `dem2raster`      -- Extract data from Digital Elevation Model as Raster
- `ugrid_bathy2tin` -- Extract bathymetry from the FVCOM unstructured 
                       grid model as TIN, via LandXML[1]
                       
Notes
-----
 To run these routine, you need:
 - ArcGIS 10.0 or 10.1
 - NetCDF4-Python using installers at
   http://www.lfd.uci.edu/~gohlke/pythonlibs/#netcdf4
 - In ArcMap: Geoprocessing=>Python, and in the Python Command Window, type
       >>> sys.path.append('.')
       >>> import dap2arc
     
       >>> dap2arc.ugrid_bathy2tin()  # specify your paths
       Then in ArcMap, go to Add Data, select your TIN and add it!

References
----------
.. [1] LandXML Group, http://landxml.org/

Examples
--------
>>> import dap2arc
>>> dap2arc.ugrid_bathy2tin(url = 'http://geoport.whoi.edu/thredds/dodsC/usgs/data1/rsignell/models/fvcom/GOM2_2008/gom2_200804.nc',
    landXMLfile = 'c:/rps/python/esri_demos/landXML.xml',
    outputFolder = 'c:/rps/python/esri_demos/tins',
    outputBase = 'bathy',
    tinPosition = "1")

"""

import numpy as np
import netCDF4
import time
import shutil
import arcpy
import os
import datetime
import shapefile

def ugrid_bathy2tin(url = 'http://geoport.whoi.edu/thredds/dodsC/usgs/data1/rsignell/models/fvcom/GOM2_2008/gom2_200804.nc',
    landXMLfile = 'c:/rps/python/esri_demos/landXML.xml',
    outputFolder = 'c:/rps/python/esri_demos/tins',
    outputBase = 'bathy',
    tinPosition = "1"):
    # make Tin directory if it doesn't exist 
    try:
        shutil.os.mkdir(outputFolder)
    except:
        print('Output Folder {0} exists'.format(outputFolder))
    
    nc=netCDF4.Dataset(url)
    x=nc.variables['lon'][:]
    y=nc.variables['lat'][:]
    # read water depth at nodes
    h=-nc.variables['h'][:]
    # read connectivity array
    nv=nc.variables['nv'][:,:]
    nc.close()
    nv=nv.T
    nnodes=len(h)
    nele,three=np.shape(nv)

    # remove existing TIN before writing (don't know how to overwrite)
    try:
        shutil.rmtree(outputFolder+'/'+outputBase)
    except:
        print('Creating {0}'.format(outputBase))
    
    xml_header = '''<?xml version="1.0"?>
    <LandXML version="1.1" date="2008-06-09" time="08:33:4"
     xmlns="http://www.landxml.org/schema/LandXML-1.1"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.1 http://www.landxml.org/schema/LandXML-1.1/LandXML-1.1.xsd">
     <Units>
      <Metric areaUnit="squareMeter" linearUnit="meter" volumeUnit="cubicMeter"
       temperatureUnit="celsius" pressureUnit="milliBars"/>
     </Units>
     <Application name="FVCOM" version="10.2.0 (Build 400)" manufacturer="SMAST"
      manufacturerURL="http://fvcom.smast.umassd.edu/FVCOM/index.html"><Author createdBy="rsignell"/></Application>
     
     
     <Surfaces>
      <Surface name="FVCOM_GOM2_Grid">
    '''   
    time1=time.clock()
    f=open(landXMLfile, 'w')
    f.write(xml_header)
    f.write('  <Definition surfType=\"TIN\" elevMin=\"%f\" elevMax=\"%f\">\n' % (min(h[:]), max(h[:])))
    f.write('    <Pnts>\n')

    Lines = [\
        '      <P id=\"%d\"> %f %f %f </P>\n' % (i+1, y[i], x[i], h[i]) \
        for i in range(nnodes)]
    f.writelines(Lines)


    f.write('   </Pnts>\n')
    f.write('   <Faces>\n')

    Lines = [\
        '      <F> %d %d %d </F>\n' % (nv[i,0], nv[i,1], nv[i,2]) \
        for i in range(nele)]
    f.writelines(Lines)

    xml_trailer = '''    </Faces>
       </Definition>
      </Surface>
     </Surfaces>
    </LandXML>
    '''
    f.write(xml_trailer)
    f.close()
    time2=time.clock()
    print 'elapsed time %s' % (time2-time1)

    try:
        arcpy.env.overwriteOutput = True
        arcpy.LandXMLToTin_3d(landXMLfile,outputFolder,outputBase,tinPosition)
    except:
        print "arcpy not found, or problem using arcpy.LandXMLToTin_3d"

def dem2raster(url ='http://geoport.whoi.edu/thredds/dodsC/bathy/crm_vol1.nc',
    box = [-71.4,41,-70.2,42]) :
    '''
    NetCDF4-Python test to read DEM data via OPeNDAP and create Arc Raster
    also tests out writing a small plot using Matplotlib
    
    Rich Signell, original coding
    Curtis Price, small tweaks and doc
    Rich Signell, added plotting, checking for uniform grid spacing and
                generalized to handle any bathy/topo data at:
                http://geoport.whoi.edu/thredds/bathy_catalog.html
    

     #extract data from OPeNDAP or NetCDF URL 
     url='c:/rps/python/traster_RasterToNetCDF.nc'
     #Test Latitude coordinate decreasing:
     url='http://geoport.whoi.edu/thredds/dodsC/bathy/crm_vol1.nc'
     #Test Latitude coordinate increasing:
     url='http://geoport.whoi.edu/thredds/dodsC/bathy/etopo1_bed_g2'
     #Test failure because lat is not uniformly spaced:
     url='http://geoport.whoi.edu/thredds/dodsC/bathy/smith_sandwell_v11
     box=[-71.0,41.0,-67.0,46.0]
    '''
    nc=netCDF4.Dataset(url)
    print "Source name: %s" % nc.title
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]
    bi=(lon>=box[0])&(lon<=box[2])
    bj=(lat>=box[1])&(lat<=box[3])
    z=nc.variables['topo'][bj,bi]
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

    try:
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
    except:
        print "matplotlib not found"
        
    nc.close()
    
def ww2raster(url ='http://motherlode.ucar.edu/thredds/dodsC/fmrc/NCEP/WW3/Regional_US_West_Coast/NCEP-WW3-Regional_US_West_Coast_best.ncd',
    box = [-132.95925,35.442,-117.279,51.12225],
    var = 'Significant_height_of_combined_wind_waves_and_swell'):
    
    '''
    NetCDF4-Python test to read DEM data via OPeNDAP and create Arc Raster
    also tests out writing a small plot using Matplotlib
    Global: http://motherlode.ucar.edu/thredds/dodsC/fmrc/NCEP/WW3/Global/NCEP-WW3-Global_best.ncd
    West Coast: http://motherlode.ucar.edu/thredds/dodsC/fmrc/NCEP/WW3/Regional_US_West_Coast/NCEP-WW3-Regional_US_West_Coast_best.ncd
    '''
    nc = netCDF4.Dataset(url)
    print "Source name: %s" % nc.title
    lon = nc.variables['lon'][:]-360.0
    lat = nc.variables['lat'][:]
    bi = (lon>=box[0]) & (lon<=box[2])
    bj = (lat>=box[1]) & (lat<=box[3])
    
    # find time index to read
    hours_from_now = 0   # Examples: 0=>nowcast, 3 => forecast 3 hours from now, etc. 
    date = datetime.datetime.utcnow()+datetime.timedelta(0,3600*hours_from_now)  
    #date=datetime.datetime(2011,9,9,17,00)  # specific time (UTC)
    
    tindex = netCDF4.date2index(date,nc.variables['time'],select='nearest')
    z = nc.variables[var][tindex,bj,bi]
    lonmin = np.min(lon[bi])
    latmin = np.min(lat[bj])
    dx = np.diff(lon)
    dy = np.diff(lat)
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
    rasterName = "sig_ht"
    outRaster = os.path.normpath(os.path.join(arcpy.workspace,rasterName))
    print outRaster
    grid1=arcpy.NumPyArrayToRaster(z,xyOrig,dx,dy)
    grid1.save(os.path.join(arcpy.workspace,outRaster))
    strPrj = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID"\
             "['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],"\
             "UNIT['Degree',0.0174532925199433]]"
    arcpy.DefineProjection_management(outRaster,strPrj)
    print "Written: %s" % grid1
    arcpy.AddMessage("Written: %s" % grid1)
    nc.close()
    
    
 
def ugrid_tidecur2shape(url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc',
    file_poly='c:/rps/python/shapefiles/tidal_current_poly',
    file_point='c:/rps/python/shapefiles/tidal_current_point'):

    '''

    # Specify the URL

    # IOOS Northeast Coastal Ocean Forecast System (NECOFS) using triangle-based FVCOM model

    # GOM2: lower resolution Gulf of Maine 
    #url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM2_FORECAST.nc'

    # GOM3: higher resolution along the coast and larger domain (GOM3):
    #url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'

    # MASSBAY: mass bay model
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'

    # processing starts here

    # Open OpenDAP URL to get remote triangular mesh + data
    '''


    nc=netCDF4.Dataset(url)

    # read lon,lat 
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]

    # read connectivity array
    nv=nc.variables['nv'][:,:].T  # transpose to get (ncells,3)
    [ncells,three]=np.shape(nv)
    nv=nv-1  # python is 0-based


    # now read data at nodes:
    # read water depth at nodes
    h=nc.variables['h'][:]

    # find time index to read
    hours_from_now=0   # Examples: 0=>nowcast, 3 => forecast 3 hours from now, etc. 
    date=datetime.datetime.utcnow()+datetime.timedelta(0,3600*hours_from_now)  
    #date=datetime.datetime(2011,9,9,17,00)  # specific time (UTC)
    tindex=netCDF4.date2index(date,nc.variables['time'],select='nearest')

    # read water level at nodes at a specific time
    #z=nc.variables['zeta'][-1,:]   # -1 is the last time step of the forecast
    #z=nc.variables['zeta'][tindex,:]  # index for date specified above

    # read significant wave height at nodes at specified time step 
    #z=nc.variables[var][tindex,:]   # Note: 'hs' is only in GOM3 wave model

    # write Shapefile using pyshp

    # Test 1. Write polygons for each triangle, and create record
    # values for each triangle that are the average of the 3 nodal values
    # (since depth and water level are defined at nodes)


    ilev=0  # [0] surface, [-1] bottom
    u=nc.variables['u'][tindex,ilev,:]
    v=nc.variables['v'][tindex,ilev,:]
    ang=np.angle(u + v*1j)*180/np.pi
    spd=np.abs(u + v*1j)
    lonc=nc.variables['lonc'][:]
    latc=nc.variables['latc'][:]
    nc.close()

    w = shapefile.Writer(shapefile.POLYGON)
    w.field('Depth(m)','F',8,2)
    w.field('Speed(m/s)','F',8,3)
    for i in range(ncells):
       w.poly(parts=[[ [lon[nv[i,0]],lat[nv[i,0]]], [lon[nv[i,1]],lat[nv[i,1]]], [lon[nv[i,2]],lat[nv[i,2]]] ]])
       w.record(np.mean(h[nv[i,:]]),spd[i])

    w.save(file_poly)
    # create the PRJ file
    prj = open("%s.prj" % file_poly, "w")
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()

    # Test 2. Write points for the center of each element, and create
    # record values with u,v velocity components (since u,v defined at cell centers)

    w = shapefile.Writer(shapefile.POINT)
    w.field('ang(math)','F',8,3)
    w.field('speed(m/s)','F',8,3)
    #w.field('u(m/s)','F',8,3)
    #w.field('v(m/s)','F',8,3)
    for i in range(ncells):
       w.point(lonc[i],latc[i])
       w.record(ang[i],spd[i])
    #   w.record(ang[i],spd[i],u[i],v[i])
     
    w.save(file_point)

    # create the PRJ file
    prj = open("%s.prj" % file_point, "w")
    epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()