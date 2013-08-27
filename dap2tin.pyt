import arcpy
import os
import netCDF4
import numpy as np
import datetime as dt
import shutil
try:
    import pyproj
    pyproj_enabled = True
except ImportError:
    pyproj_enabled = False

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
        self.canRunInBackground = False
        self.url = None
        self.dataset = None
        self.cols = {
            'url': 0,
            'dataset_var': 1,
            'iyear': 2,
            'imonth': 3,
            'iday': 4,
            'ihour': 5,
            'klev': 6,
            'sr' : 7,
            'out_tin': 8
        }
        self.default_urls = ['http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/hindcasts/30yr_gom3/mean',
        'http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/hindcasts/wave_gom3',
        'http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/archives/necofs_gom3',
        'http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/archives/necofs_mb',
        'http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/necofs/forecasts/necofs_gom3_forecast.nc',
        'http://www.smast.umassd.edu:8080/thredds/dodsc/fvcom/necofs/forecasts/necofs_fvcom_ocean_massbay_forecast.nc']

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
        url.filter.list = self.default_urls

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

        sr = arcpy.Parameter(
            displayName="Spatial Reference",
            name="sr",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")
        # default to USGS Mass Plane for now; could also use actual data projection.
        sr.value = arcpy.SpatialReference(26986).factoryCode

        # Parameter: Output Tin
        out_tin = arcpy.Parameter(
            displayName="Output TIN",
            name="outTin",
            datatype="TIN",
            parameterType="Required",
            direction="Output")
        # set default name of output TIN
        out_tin.value = 'c:/rps/python/tins/necofs'

        # set location of layer files
        current_dir = os.path.dirname(__file__)
        layer_dir = os.path.join(current_dir, 'layers')
        layers = {"temp": 'temperature.lyr',
                  "salinity": 'salinity.lyr',
                  "hs" : 'wave_height.lyr'}
        if dataset_var.value in layers:
            symbology_file = layers[dataset_var.value]
        else:
            symbology_file = 'temperature.lyr'
        out_tin.symbology = os.path.join(layer_dir, symbology_file) 
        
        return [url, dataset_var, iyear, imonth, iday, ihour, klev, sr, out_tin]

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
        url = parameters[0]
        variable = parameters[1]

        if url.value is not None and self.url != url.value:
            self.url = url.value
            try:
                self.dataset = netCDF4.Dataset(self.url, filter_out_nd_coordinates=True)
            except:
                arcpy.AddError("Unable to open NetCDF source: {}".format(self.url))
                return
               
            var_names = list(self.dataset.variables.keys())
            parameters[1].filter.list = var_names
            # add back the new url to to filter list
            parameters[0].filter.list = [self.url] + self.default_urls
            if variable.value not in var_names:
                parameters[1].value = var_names[0]
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
        sr_string = parameters[7].valueAsText
        out_tin = parameters[8].valueAsText

        sr = arcpy.SpatialReference()
        sr.loadFromString(sr_string)
        arcpy.env.outputCoordinateSystem = sr

        # output location
        outputFolder = os.path.dirname(out_tin)
        outputBase = os.path.basename(out_tin)
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

        # convert Lon/Lat to Mass State Plane
        if pyproj_enabled:
            # use PyProj(Proj4)
            p1 = pyproj.Proj(init='epsg:4326')   # geographic WGS84
            try:
                # factory code / EPSG relationship documented at:
                #  http://gis.stackexchange.com/q/18651/143
                epsg_code = sr.factoryCode
                p2 = pyproj.Proj(init='epsg:{0}'.format(epsg_code))
            except:
                arcpy.AddError("Unable to use PyProj to reproject to EPSG {0}".format(epsg_code))
                sys.exit()
            x,y = pyproj.transform(p1,p2,x,y)
        else:
            for i in range(len(x)):
                # create a point; cast to float as Point doesn't know about numpy types
                point = arcpy.Point(float(x[i]), float(y[i]))
                point_geom = arcpy.PointGeometry(point, arcpy.SpatialReference(4326))
                proj_point = point_geom.projectAs(sr)
                (x[i], y[i]) = (proj_point.firstPoint.X, proj_point.firstPoint.Y)

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




