# SUMMARY:      wshdslope.py (QGIS version)
# USAGE:        Computes slope and aspect from a DEM and saves them as GeoTIFF files
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted for QGIS by Y. Song
# E-MAIL:       yiyun.song@duke.edu
# ORIG-DATE:    Apr-2017
# LAST CHANGE:  2025-09-09
# DESCRIPTION:  Converts ArcGIS `arcpy` script to work with QGIS using GDAL and NumPy
#               Computes terrain slope and aspect based on DEM derivatives.
# REQUIREMENTS: QGIS with Python 3, GDAL, and NumPy

import os
import numpy as np
from osgeo import gdal, osr

def read_raster(file_path):
    """ Reads a raster file and returns a NumPy array, affine transformation parameters, and projection information. """
    dataset = gdal.Open(file_path)
    if dataset is None:
        raise FileNotFoundError(f"Unable to open file: {file_path}")
    
    band = dataset.GetRasterBand(1)
    array = band.ReadAsArray().astype(np.float32)
    
    # Handle NoData values
    nodata_value = band.GetNoDataValue()
    if nodata_value is not None:
        array[array == nodata_value] = np.nan
    
    return array, dataset.GetGeoTransform(), dataset.GetProjection()

def write_raster(output_path, array, transform, projection, nodata_value=np.nan):
    """ Saves a NumPy array as a GeoTIFF raster file. """
    driver = gdal.GetDriverByName("GTiff")
    rows, cols = array.shape
    dataset = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    
    dataset.SetGeoTransform(transform)
    dataset.SetProjection(projection)
    
    band = dataset.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(nodata_value)
    
    dataset.FlushCache()
    dataset = None  # Close the file
    print(f"Output file saved: {output_path}")

def compute_slope_aspect(dem_path, slope_output, aspect_output):
    """ Computes slope and aspect from a DEM and saves the results as GeoTIFF files. """
    dem, transform, projection = read_raster(dem_path)
    
    # Get DEM resolution
    x_res = transform[1]  # Resolution in the X direction
    y_res = -transform[5]  # Resolution in the Y direction (negative to correct sign)

    # Compute elevation change (dzdx, dzdy)
    dzdx = np.gradient(dem, axis=1) / (2 * x_res)
    dzdy = np.gradient(dem, axis=0) / (2 * y_res)
    
    # Compute slope
    slope = np.sqrt(dzdx**2 + dzdy**2)
    
    # Compute aspect using atan2 and convert to degrees
    aspect = np.arctan2(dzdx, dzdy) * (180.0 / np.pi)
    aspect = (aspect + 360) % 360  # Ensure values are within [0, 360] degrees
    
    # Handle NaN values
    slope[np.isnan(dem)] = np.nan
    aspect[np.isnan(dem)] = np.nan

    # Save slope and aspect rasters
    write_raster(slope_output, slope, transform, projection)
    write_raster(aspect_output, aspect, transform, projection)

# Define input DEM file and output files
dem_file = "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CAAR/Reprojected_DEM/USGS_1_n36w084_20220725_UTM17.tif"
slope_output = "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CAAR/Processed/slope_UTM17.tif"
aspect_output = "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CAAR/Processed/aspect_UTM17.tif"

# Run slope and aspect computation
compute_slope_aspect(dem_file, slope_output, aspect_output)
