DAP2ARC
=======

Python tools for bridging the gap between OPeNDAP and ArcGIS. 

This project is developing tools and techniques for ArcGIS users to remotely access ocean and atmospheric model data and create customized data products.  Model output is typically stored GRIB, HDF or NetCDF files, and made available via the OPeNDAP web service, which allows extraction of specific variables and geospatial subsets.   Managers often use ArcGIS, which has no native ability to read data from OPeNDAP.   This gap can be bridged using Python, available with ArcGIS 9.3 and higher.  We use here the OPeNDAP Python module NetCDF4-Python, which leverages the Unidata NetCDF C-library for OPeNDAP support. With OPeNDAP access via Python, we can build ArcGIS tools that allow managers to access ocean model data on both rectilinear and unstructured triangular native grids, avoiding loss of information through decimation or interpolation, and allowing the creation of custom products to meet specific needs.  

##Requirements

ArcGIS 10.1 (The Python Toolbox was introduced with 10.1).  

NetCDF4-Python, enabled with OPeNDAP.  We also use the Enthought Python Distribution, which includes the NetCDF4-Python binary for both 32 and 64 bit windows, enabling foreground 32-bit processing as well as 64-bit background processing. 

If you don't have the EPD, you can download the 32-bit NetCDF4-Python Windows Binary for ArcGIS 10.1 at:
http://www.lfd.uci.edu/~gohlke/pythonlibs/#netcdf4