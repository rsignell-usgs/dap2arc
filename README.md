DAP2ARC
=======

Python tools for bridging the gap between OPeNDAP and ArcGIS. 

This project is developing tools and techniques for ArcGIS users to remotely access ocean and atmospheric model data and create customized data products.  Model output is typically stored GRIB, HDF or NetCDF files, and made available via the OPeNDAP web service, which allows extraction of specific variables and geospatial subsets.   Managers often use ArcGIS, which has no native ability to read data from OPeNDAP.   We bridge this gap here using the OPeNDAP support in NetCDF4-Python, which leverages the Unidata NetCDF C-library for OPeNDAP support. With OPeNDAP access via Python, we can build ArcGIS tools that allow managers to access ocean model data on both rectilinear and unstructured triangular native grids, avoiding loss of information through decimation or interpolation, and allowing the creation of custom products to meet specific needs. 

##Requirements

ArcGIS 10.1 (The Python Toolbox was introduced with 10.1).  

NetCDF4-Python:  This module is now distributed by ESRI as part of the Multidimensional Supplemental Tools package:
http://esriurl.com/MultidimensionSupplementalTools

##Usage

The dap2arc tools are either python scripts (ending in .py) or python toolboxes (ending in .pyt)

Python scripts can be executed by going to Geoprocessing=>Python and cut-and-pasting the commands in the scripts into the Python command window. 

Python toolboxes can be executed by using ArcToolbox to navigate to the appropriate directory and double clicking the python toolbox file, which will open an interface. 

##Tutorial
I made a 12 minute video tutorial at:
https://www.youtube.com/watch?v=OccLjc4M2Tg

best to watch it in HD:
https://www.youtube.com/v/OccLjc4M2Tg

Good Luck,
Rich
