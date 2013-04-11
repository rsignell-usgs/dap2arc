import numpy as np
import netCDF4
import time

url='http://geoport.whoi.edu/thredds/dodsC/usgs/data1/rsignell/models/fvcom/GOM2_2008/gom2_200804.nc'
nc=netCDF4.Dataset(url)
x=nc.variables['lon'][:]
y=nc.variables['lat'][:]
# read water depth at nodes
h=-nc.variables['h'][:]
# read connectivity array
nv=nc.variables['nv'][:,:]
nv=nv.T

nnodes=len(h)
nele,three=np.shape(nv)

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
f=open('c:/rps/python/test.xml', 'w')
f.write(xml_header)
f.write('  <Definition surfType=\"TIN\" elevMin=\"%f\" elevMax=\"%f\">\n' % (min(h[:]), max(h[:])))
f.write('    <Pnts>\n')

Lines = [\
    '      <P id=\"%d\"> %f %f %f </P>\n' % (i+1, y[i], x[i], h[i]) \
    for i in range(nnodes)]
f.writelines(Lines)

#for i in arange(nnodes):
#    f.write('      <P id=\"%d\"> %f %f %f </P>\n' % (i+1, x[i], y[i], h[i]))

f.write('   </Pnts>\n')
f.write('   <Faces>\n')

Lines = [\
    '      <F> %d %d %d </F>\n' % (nv[i,0], nv[i,1], nv[i,2]) \
    for i in range(nele)]
f.writelines(Lines)
#for i in arange(nele):
#    f.write('      <F> %d %d %d </F>\n' % (nv[i,0], nv[i,1], nv[i,2]))

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
    from arcpy import env
    outputFolder="tins"
    outputBase = "test"
    grab = "1"
    env.workspace="c:/rps/python/"
    arcpy.LandXMLToTin_3d("c:/rps/python/test.xml",outputFolder,outputBase,grab)
except:
    print "arcpy not found"
