# ----------------------------------------------------------------------
# Script: dhsvm_area_config_from_dem.py
#
# Purpose:
#   Read a DEM raster (e.g., elev_clipped.tif) in QGIS and print a
#   DHSVM [AREA] configuration block to the console, including:
#     - Extreme North (northern edge of the grid, in map units)
#     - Extreme West (western edge of the grid, in map units)
#     - Number of Rows / Columns
#     - Grid spacing (assumed square cells, in meters)
#     - Center Latitude / Longitude (in decimal degrees, WGS84)
#     - Time Zone Meridian (approximate, derived from center longitude)
#     - Watershed Area (km^2), computed from non-NoData DEM cells
#
# Usage:
#   1. Open QGIS.
#   2. Open the Python Console.
#   3. Paste this script and run it.
#   4. Copy the printed [AREA] block into your DHSVM INPUT file.
#
# Notes:
#   - The script assumes a projected DEM in meters (e.g., UTM).
#   - It transforms the DEM center to EPSG:4326 to obtain lat/lon.
#   - Watershed area is computed from non-NoData cells of band 1.
# ----------------------------------------------------------------------

import math
from qgis.core import (
    QgsRasterLayer,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsProject,
    QgsRasterBandStats,
)

# ----------------------------------------------------------------------
# 1. Set the path to your DEM raster
# ----------------------------------------------------------------------
DEM_PATH = "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_1201/Reprojected_DEM/elev_clipped.tif"

# If you prefer to use the currently selected raster layer in QGIS,
# you can comment out DEM_PATH above and uncomment the following lines:
#
# layer = iface.activeLayer()
# if not isinstance(layer, QgsRasterLayer):
#     raise RuntimeError("Active layer is not a raster. Please select your DEM layer.")
# dem_layer = layer

# ----------------------------------------------------------------------
# 2. Load the DEM as a QgsRasterLayer
# ----------------------------------------------------------------------
dem_layer = QgsRasterLayer(DEM_PATH, "dem_for_dhsvm")

if not dem_layer.isValid():
    raise RuntimeError(f"Could not load DEM from: {DEM_PATH}")

# ----------------------------------------------------------------------
# 3. Extract raster geometry: extent, grid size, cell size
# ----------------------------------------------------------------------
extent = dem_layer.extent()
xmin = extent.xMinimum()
xmax = extent.xMaximum()
ymin = extent.yMinimum()
ymax = extent.yMaximum()

ncols = dem_layer.width()
nrows = dem_layer.height()

# Pixel size in X/Y (may be negative depending on geotransform orientation)
px = dem_layer.rasterUnitsPerPixelX()
py = dem_layer.rasterUnitsPerPixelY()

# Use absolute value of cell size
cell_size_x = abs(px)
cell_size_y = abs(py)

# Check if cells are (approximately) square
if not math.isclose(cell_size_x, cell_size_y, rel_tol=1e-6):
    print("WARNING: Pixel size in X and Y differ. DHSVM assumes square cells.")
grid_spacing = cell_size_x  # for DHSVM [AREA] block

# Extreme North / West in map units (e.g., meters for UTM)
extreme_north = ymax
extreme_west = xmin

# ----------------------------------------------------------------------
# 4. Compute center point and transform to WGS84 (lat/lon)
# ----------------------------------------------------------------------
cx = 0.5 * (xmin + xmax)
cy = 0.5 * (ymin + ymax)

src_crs = dem_layer.crs()
wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")

coord_transform = QgsCoordinateTransform(src_crs, wgs84, QgsProject.instance())
center_point = coord_transform.transform(cx, cy)

center_lon = center_point.x()  # longitude (degrees)
center_lat = center_point.y()  # latitude (degrees)

# ----------------------------------------------------------------------
# 5. Approximate time zone meridian from center longitude
#    (round to nearest multiple of 15 degrees)
# ----------------------------------------------------------------------
time_zone_meridian = 15 * round(center_lon / 15.0)

# ----------------------------------------------------------------------
# 6. Compute watershed area from DEM cell count (non-NoData)
# ----------------------------------------------------------------------
provider = dem_layer.dataProvider()

# Use band 1 for statistics over the full extent
stats = provider.bandStatistics(
    1,
    QgsRasterBandStats.All,
    extent,
    0  # no special statistics flags
)

valid_pixels = stats.elementCount  # number of non-NoData pixels
cell_area_m2 = cell_size_x * cell_size_y
watershed_area_m2 = valid_pixels * cell_area_m2
watershed_area_km2 = watershed_area_m2 / 1.0e6

# Also compute a simple "full extent" area for comparison (includes NoData)
full_extent_pixels = nrows * ncols
full_extent_area_km2 = full_extent_pixels * cell_area_m2 / 1.0e6

# ----------------------------------------------------------------------
# 7. Print DHSVM [AREA] block
# ----------------------------------------------------------------------
print("################################################################################")
print("# MODEL AREA SECTION")
print("################################################################################")
print("")
print("[AREA]                                   # Model area")
print("Coordinate System    =  UTM              # UTM or USER_DEFINED")
print(f"Extreme North        =  {extreme_north:.8f} # Coordinate for northern edge of grid")
print(f"Extreme West         =  {extreme_west:.8f} # Coordinate for western edge of grid")
print(f"Center Latitude      =  {center_lat:.6f}   # Central parallel of basin (deg, WGS84)")
print(f"Center Longitude     =  {center_lon:.6f}   # Central meridian of basin (deg, WGS84)")
print(f"Time Zone Meridian   =  {time_zone_meridian:.1f}         # Standard time meridian (deg)")
print(f"Number of Rows       =  {nrows:4d}              # Number of rows")
print(f"Number of Columns    =  {ncols:4d}              # Number of columns")
print(f"Grid spacing         =  {grid_spacing:.3f}         # Grid resolution (m)")
print("Point North          =                    # For POINT mode only")
print("Point East           =                    # For POINT mode only")
print("")
print("################################################################################")
print("# Additional metadata (for reference, not required in [AREA])")
print("################################################################################")
print(f"# DEM path: {DEM_PATH}")
print(f"# DEM CRS:  {src_crs.authid()} ({src_crs.description()})")
print(f"# Center (projected) X, Y: {cx:.3f}, {cy:.3f}")
print(f"# Center (lat, lon) in decimal degrees: {center_lat:.6f}, {center_lon:.6f}")
print(f"# Cell size (X, Y) in meters: {cell_size_x:.3f}, {cell_size_y:.3f}")
print(f"# Watershed area (non-NoData cells): {watershed_area_km2:.3f} km^2")
print(f"# Full extent area (rows*cols):      {full_extent_area_km2:.3f} km^2")
print(f"# Valid pixels (non-NoData):         {valid_pixels}")
print(f"# Full extent pixels (rows*cols):    {full_extent_pixels}")
