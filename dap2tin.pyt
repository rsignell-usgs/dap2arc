import arcpy
import os
import netCDF4
import numpy as np
import pyproj
import datetime as dt
import shutil

def DefineProjectionForTin(tin,prj):
    '''Workaround for Arc bug - define a projection for a TIN
    1. Create a temporary raster
    2. Define its projection
    3. Copy the prj.adf file to the TIN

    '''
    wks = os.path.dirname(tin)
    tempGrid = arcpy.CreateScratchName("xx","","RasterDataset",wks)
    arcpy.CreateRasterDataset_management(wks,os.path.basename(tempGrid))
    arcpy.DefineProjection_management(tempGrid,prj)
    shutil.copyfile(os.path.join(tempGrid,"prj.adf"),
                    os.path.join(tin,"prj.adf"))
    arcpy.Delete_management(tempGrid)

def writeLandXML(nv,x,y,h,xmlTemp):
    '''Write TIN components to a LandXML file.

    Arguments:
    nv = 3 column connectivity array
    x = x locations of triangle nodes
    y = y locations of triangle nodes
    h = the data value at the triangle nodes
    xmlTemp = the full path of the XML file to be written

    '''
    nnodes = len(h)
    nele,three = np.shape(nv)
    # Write the LandXML file
    f = open(xmlTemp, 'w')
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


    f.write(xml_header)
    f.write('<Definition surfType=\"TIN\" elevMin=\"%f\" elevMax=\"%f\">\n' % \
            (min(h[:]), max(h[:])))
    f.write('<Pnts>\n')
    Lines = ['<P id=\"%d\"> %f %f %f </P>\n' % \
            (i+1, y[i], x[i], h[i]) for i in range(nnodes)]
    f.writelines(Lines)
    f.write('</Pnts>\n')
    f.write('<Faces>\n')
    Lines = ['<F> %d %d %d </F>\n' % \
            (nv[i,0], nv[i,1], nv[i,2]) for i in range(nele)]
    f.writelines(Lines)
    xml_trailer = '''    </Faces>
       </Definition>
      </Surface>
     </Surfaces>
    </LandXML>
    '''
    f.write(xml_trailer)
    f.close()


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Dap2tin"
        self.alias = "dap2tin"

        # List of tool classes associated with this toolbox
        self.tools = [Dap2tin]


class Dap2tin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Dap2tin"
        self.description = "Access data from the FVCOM Ocean Model via DAP and return a TIN"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Parameter: URL
        url = arcpy.Parameter(
            displayName="OPeNDAP URL",
            name="url",
            datatype="String",
            parameterType="Required",
            direction="Input")
        # set default value to 30 year hindcast monthly mean dataset
        url.value = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean'
        url.filter.type = "ValueList"
        url.filter.list = ['http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean',
        'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/wave_gom3',
        'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_gom3',
        'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_mb',
        'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc',
        'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc']

        # Parameter: Variable name
        dataset_var = arcpy.Parameter(
            displayName="Variable",
            name="dataset_var",
            datatype="String",
            parameterType="Required",
            direction="Input")
        # set default value to temperature
        dataset_var.value = 'temp'
        dataset_var.filter.type = "ValueList"
        dataset_var.filter.list = ["temp","salinity","hs","h"]

        # Parameter: Year
        iyear = arcpy.Parameter(
            displayName="Year",
            name="iyear",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        # set default value to 2nd time step
        iyear.value = 2010

        # Parameter: Month
        imonth = arcpy.Parameter(
            displayName="Month",
            name="imonth",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        # set default value to 2nd time step
        imonth.value = 1

        # Parameter: Day
        iday = arcpy.Parameter(
            displayName="Day",
            name="iday",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        # set default value to 2nd time step
        iday.value = 1

        # Parameter: Hour
        ihour = arcpy.Parameter(
            displayName="Hour",
            name="ihour",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        # set default value to 2nd time step
        ihour.value = 0

        # Parameter: Level
        klev = arcpy.Parameter(
            displayName="Vertical layer (0 is surface, -1 is bottom)",
            name="klev",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        # set default value to surface layer
        klev.value = 0

        # Parameter: Output Tin
        outTin = arcpy.Parameter(
            displayName="Output TIN",
            name="outTin",
            datatype="TIN",
            parameterType="Required",
            direction="Output")
	    # set default name of output TIN
        outTin.value = 'c:/rps/python/tins/necofs'


        # set location of layer files
        layer_dir = 'c:/users/rsignell/documents/github/dap2arc/'
        if dataset_var.value == "temp":
            outTin.symbology = layer_dir + 'temperature.lyr'
        elif dataset_var.value == "salinity":
            outTin.symbology = layer_dir + 'salinity.lyr'
        elif dataset_var.value == "hs":
            outTin.symbology = layer_dir + 'wave_height.lyr'
        else:
            outTin.symbology = layer_dir + 'temperature.lyr'

        return [url, dataset_var, iyear,imonth,iday,ihour, klev, outTin]

	def isLicensed(self):
		"""LandXMLToTin_3d used in this routine requires the ArcGIS 3D Analyst extension
		to be available."""
		try:
			if arcpy.CheckExtension("3D") == "Available":
				arcpy.CheckOutExtension("3D")
			else:
				raise Exception
		except:
			return False # tool cannot be executed

		return True # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        url = parameters[0].valueAsText
        dataset_var = parameters[1].valueAsText
        iyear = int(parameters[2].valueAsText)
        imonth = int(parameters[3].valueAsText)
        iday = int(parameters[4].valueAsText)
        ihour = int(parameters[5].valueAsText)

        klev = int(parameters[6].valueAsText)
        outTin = parameters[7].valueAsText

        # create dictionary for WKT strings
        prj={}
        #ESRI WKT for Geographic, WGS84 (EPSG:4326)
        prj['4326'] = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',"\
             "SPHEROID['WGS_1984',6378137.0,298.257223563]],"\
             "PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

        #ESRI WKT for Massachusetts state plane, mainland (EPSG:26986)
        prj['26986'] = "PROJCS['NAD_1983_StatePlane_Massachusetts_Mainland_FIPS_2001',"\
             "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',"\
             "SPHEROID['GRS_1980',6378137.0,298.257222101]],"\
             "PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],"\
             "PROJECTION['Lambert_Conformal_Conic'],"\
             "PARAMETER['False_Easting',200000.0],"\
             "PARAMETER['False_Northing',750000.0],"\
             "PARAMETER['Central_Meridian',-71.5],"\
             "PARAMETER['Standard_Parallel_1',41.71666666666667],"\
             "PARAMETER['Standard_Parallel_2',42.68333333333333],"\
             "PARAMETER['Latitude_Of_Origin',41.0],UNIT['Meter',1.0]]"

        arcpy.env.outputCoordinateSystem = prj['26986']

        # output location
        outputFolder= os.path.dirname(outTin)
        outputBase = os.path.basename(outTin)
        arcpy.env.workspace = outputFolder
        os.chdir(arcpy.env.workspace)

        arcpy.AddMessage("Reading OPeNDAP data from %s" % url)
        arcpy.AddMessage((repr(dataset_var)))

        nc = netCDF4.Dataset(url)
        x = nc.variables['lon'][:]
        y = nc.variables['lat'][:]

        # determine time step from [year,month,day,hour]
        times = nc.variables['time']
        start = dt.datetime(iyear,imonth,iday,ihour)
        itime = netCDF4.date2index(start,times,select='nearest')
        arcpy.AddMessage("Reading time step %d" % itime)

        # read connectivity array
        nv = nc.variables['nv'][:,:]
        nv=nv.T

        # convert Lon/Lat to Mass State Plane using PyProj(Proj4)
        p1 = pyproj.Proj(init='epsg:4326')   # geographic WGS84
        p2 = pyproj.Proj(init='epsg:26986') # Mass State Plane
        x,y = pyproj.transform(p1,p2,x,y)

        # read data at nodes for selected variable
        # handle 1D, 2D, and 3D variables
        hvar = nc.variables[dataset_var]
        numdims = len(hvar.dimensions)
        if numdims==3:   # time, level, node (e.g. temperature 'temp')
            h = nc.variables[dataset_var][itime,klev,:]
        elif numdims==2: # time, node (e.g. water level 'zeta')
            h = nc.variables[dataset_var][itime,:]
        elif numdims==1: # node (e.g. water depth 'h')
            h = nc.variables[dataset_var][:]
        else:
            arcpy.AddMessage("Variable %s has %d dimensions" % (dataset_var,numdims))
            raise Exception
        nc.close()

        # Name of temporary LandXML file
        xmlTemp = os.path.normpath(os.path.join(arcpy.env.workspace,"foo.xml"))
        arcpy.AddMessage("Writing LandXML file %s" % xmlTemp)

        # Write the LandXML file
        writeLandXML(nv,x,y,h,xmlTemp)

        arcpy.AddMessage("Converting LandXML to TIN...")
        grab = "1"  # grab the 1st TIN (there is only 1 TIN)
        arcpy.LandXMLToTin_3d(xmlTemp,outputFolder,outputBase,grab)

        #params = arcpy.GetParameterInfo()
        #arcpy.AddMessage("parameter info %s" % params)

        DefineProjectionForTin(outTin,prj['26986'])
        arcpy.AddMessage(arcpy.GetMessages())
        #arcpy.DefineProjection_management(outTin,dataPrj) # not working!
        return




