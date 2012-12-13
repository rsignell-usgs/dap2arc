import arcpy
import os
import netCDF4
import numpy as np

def DefineProjectionForTin(tin,prj):
# workaround for a bug - define a projection for a TIN
# 1. Create a temporary raster
# 2. Define its projection
# 3. Copy the prj.adf file to the TIN
    import os
    import shutil
    import arcpy
    wks = os.path.dirname(tin)
    tempGrid = arcpy.CreateScratchName("xx","","RasterDataset",wks)
    arcpy.CreateRasterDataset_management(wks,os.path.basename(tempGrid))
    arcpy.DefineProjection_management(tempGrid,prj)
    shutil.copyfile(os.path.join(tempGrid,"prj.adf"),
                    os.path.join(tin,"prj.adf"))
    arcpy.Delete_management(tempGrid)

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
        url = arcpy.Parameter(
            displayName="OPeNDAP URL",
            name="url",
            datatype="String",
            parameterType="Required",
            direction="Input")	
        # set default value to 30 year hindcast monthly mean dataset
        url.value = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean'

        dataset_var = arcpy.Parameter(
            displayName="Variable",
            name="dataset_var",
            datatype="String",
            parameterType="Required",
            direction="Input")	
        # set default value to temperature
        dataset_var.value = 'temp'
        dataset_var.filter.type = "ValueList"
        dataset_var.filter.list = ["temp","salinity"]
        
        itime = arcpy.Parameter(
            displayName="Time Step (First time step is 0)",
            name="itime",
            datatype="Long",
            parameterType="Required",
            direction="Input")	
        # set default value to 2nd time step
        itime.value = 1
        klev = arcpy.Parameter(
            displayName="Vertical layer (0 is surface, -1 is bottom)",
            name="klev",
            datatype="Long",
            parameterType="Required",
            direction="Input")	
        # set default value to surface layer
        klev.value = 0
		
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
        else:
            outTin.symbology = layer_dir + 'temperature.lyr'
        return [url, dataset_var, itime, klev, outTin]

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
        itime = int(parameters[2].valueAsText)
        klev = int(parameters[3].valueAsText)
        outTin = parameters[4].valueAsText
        prj = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',"\
             "SPHEROID['WGS_1984',6378137.0,298.257223563]],"\
             "PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
        arcpy.env.outputCoordinateSystem = prj

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
        # read water depth at nodes
        h = nc.variables[dataset_var][itime,klev,:]
        # read connectivity array
        nv = nc.variables['nv'][:,:]
        nv=nv.T
        nc.close()

        nnodes = len(h)
        nele,three = np.shape(nv)


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
        # name of temporary LandXML file
        xmlTemp = os.path.normpath(os.path.join(arcpy.env.workspace,"foo.xml"))
        arcpy.AddMessage("Writing LandXML file %s" % xmlTemp)
        f = open(xmlTemp, 'w')
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

        arcpy.AddMessage("Temporary LandXML file written")

        arcpy.AddMessage("Converting to TIN...")
        grab = "1"  # grab the 1st TIN (there is only 1 TIN)
        arcpy.LandXMLToTin_3d(xmlTemp,outputFolder,outputBase,grab)

        #params = arcpy.GetParameterInfo()
        #arcpy.AddMessage("parameter info %s" % params)

        DefineProjectionForTin(outTin,prj)
        arcpy.AddMessage(arcpy.GetMessages())
        #arcpy.DefineProjection_management(outTin,dataPrj) # not working!        
        return
        



