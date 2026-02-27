# SUMMARY:      soildepthscript.py (QGIS version)
# USAGE:        Computes soil depth based on slope, source area, and elevation.
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted for QGIS by Y.Song
# E-MAIL:       yiyun.song@duke.edu
# ORIG-DATE:    Apr-2017
# LAST CHANGE:  2025-09-09
# DESCRIPTION:  Converts ArcGIS `arcpy` script to work with QGIS `QgsRasterLayer`

from qgis.core import QgsRasterLayer
import processing

def soildepthfun(flowacc_path, elev_path, mindepth, maxdepth, soildepth_output):
    """
    Generates a soil depth raster based on DEM slope, flow accumulation, and elevation.

    Parameters:
    flowacc_path (str): Path to the flow accumulation raster (.tif)
    elev_path (str): Path to the DEM raster (.tif)
    mindepth (float): Minimum soil depth (m)
    maxdepth (float): Maximum soil depth (m)
    soildepth_output (str): Output path for soil depth raster (.tif)
    """

    # Define weighting factors
    wtslope = 0.7
    wtsource = 0.0
    wtelev = 0.3
    maxslope = 30.0
    maxsource = 100000.0
    maxelev = 1500
    powslope = 0.25
    powsource = 1.0
    powelev = 0.75

    print(f"Processing soil depth with min: {mindepth}m, max: {maxdepth}m")
    print(f"Slope weight: {wtslope}, Source weight: {wtsource}, Elevation weight: {wtelev}")

    # Compute slope from DEM
    slope_raster = elev_path.replace(".tif", "_slope.tif")
    processing.run("gdal:slope", {
        'INPUT': elev_path,
        'BAND': 1,
        'OUTPUT': slope_raster
    })
    print(f"Slope raster generated: {slope_raster}")

    # Apply threshold conditions
    flowacc_limited = flowacc_path.replace(".tif", "_flowacc_limited.tif")
    elev_limited = elev_path.replace(".tif", "_elev_limited.tif")
    slope_limited = slope_raster.replace(".tif", "_slope_limited.tif")

    processing.run("gdal:rastercalculator", {
        'INPUT_A': flowacc_path, 'BAND_A': 1,
        'FORMULA': f"A*(A<={maxsource}) + {maxsource}*(A>{maxsource})",
        'OUTPUT': flowacc_limited
    })

    processing.run("gdal:rastercalculator", {
        'INPUT_A': elev_path, 'BAND_A': 1,
        'FORMULA': f"A*(A<={maxelev}) + {maxelev}*(A>{maxelev})",
        'OUTPUT': elev_limited
    })

    processing.run("gdal:rastercalculator", {
        'INPUT_A': slope_raster, 'BAND_A': 1,
        'FORMULA': f"A*(A<={maxslope}) + {maxslope}*(A>{maxslope})",
        'OUTPUT': slope_limited
    })

    print("Threshold constraints applied.")

    # Compute soil depth using weighted function
    soildepth_formula = (
        f"{mindepth} + ({maxdepth} - {mindepth}) * ("
        f"{wtslope} * (1 - (A / {maxslope}) ** {powslope}) + "
        f"{wtsource} * ((B / {maxsource}) ** {powsource}) + "
        f"{wtelev} * (1 - (C / {maxelev}) ** {powelev})"
        f")"
    )

    processing.run("gdal:rastercalculator", {
        'INPUT_A': slope_limited, 'BAND_A': 1,
        'INPUT_B': flowacc_limited, 'BAND_B': 1,
        'INPUT_C': elev_limited, 'BAND_C': 1,
        'FORMULA': soildepth_formula,
        'OUTPUT': soildepth_output
    })

    print(f"Soil depth raster created: {soildepth_output}")

