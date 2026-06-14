"""Compute a DHSVM [AREA] block from a clipped DEM (rasterio + pyproj).

Standalone, no-QGIS port of qgis_CA/dhsvm_area_config_from_dem.py. Reads the DEM
grid geometry and returns the [AREA] field values, so make_dhsvm_config.py fills
them into the template at config time instead of the template freezing one
resolution. The same template then produces a correct [AREA] for any DEM (28 m,
10 m, or another watershed).

The time zone meridian is NOT derived from the basin longitude. The qgis script
rounds longitude to the nearest 15 deg, which returns -90 (Central) for the
Southern Appalachians and has to be hand-corrected to -75 every time. The basins
are Eastern Time, whose standard meridian is -75 (UTC-5 x 15) by definition, so
the value comes from paths.TIME_ZONE_MERIDIAN as a declared constant. Correct,
and no hand edit.

Geometry mirrors the qgis script: extreme north and west are the outer grid
edges, rows, cols and grid spacing come from the raster, and the grid center is
transformed to WGS84 for the central parallel and meridian.

CLI:
    python area.py [DEM]      # DEM defaults to paths.ELEV_CLIPPED; prints the block
"""
import sys
from pathlib import Path

import rasterio
from pyproj import Transformer

# paths.py is the single source of truth (sibling pipeline/ dir).
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pipeline"))
import paths  # noqa: E402


def area_fields(dem_path, meridian=None):
    """[AREA] field values computed from dem_path.

    meridian defaults to paths.TIME_ZONE_MERIDIAN (the declared constant).
    """
    if meridian is None:
        meridian = paths.TIME_ZONE_MERIDIAN

    with rasterio.open(str(dem_path)) as src:
        bounds = src.bounds              # left, bottom, right, top: outer grid edges
        nrows, ncols = src.height, src.width
        res_x, res_y = src.res           # cell size, both positive
        epsg = src.crs.to_epsg()

    if abs(res_x - res_y) > 1e-6:
        print(f"[area] WARNING non-square cells {res_x} x {res_y}; DHSVM assumes square",
              file=sys.stderr)

    cx = 0.5 * (bounds.left + bounds.right)
    cy = 0.5 * (bounds.bottom + bounds.top)
    to_wgs84 = Transformer.from_crs(epsg, 4326, always_xy=True)
    center_lon, center_lat = to_wgs84.transform(cx, cy)

    return {
        "extreme_north": bounds.top,
        "extreme_west":  bounds.left,
        "center_lat":    center_lat,
        "center_lon":    center_lon,
        "meridian":      float(meridian),
        "nrows":         nrows,
        "ncols":         ncols,
        "grid_spacing":  res_x,
    }


def format_area_block(fields):
    """Render the [AREA] section text from area_fields(). Same field formatting as
    the qgis script, so a 28 m DEM reproduces the calibrated block to displayed
    precision. DHSVM reads the section free-format, so spacing is not significant.
    """
    f = fields
    return (
        "[AREA]                                   # Model area\n"
        "Coordinate System    =  UTM              # UTM or USER_DEFINED\n"
        f"Extreme North        =  {f['extreme_north']:.8f} # Coordinate for northern edge of grid\n"
        f"Extreme West         =  {f['extreme_west']:.8f} # Coordinate for western edge of grid\n"
        f"Center Latitude      =  {f['center_lat']:.6f}   # Central parallel of basin (deg, WGS84)\n"
        f"Center Longitude     =  {f['center_lon']:.6f}   # Central meridian of basin (deg, WGS84)\n"
        f"Time Zone Meridian   =  {f['meridian']:.1f}         # Standard time meridian (deg)\n"
        f"Number of Rows       =  {f['nrows']}              # Number of rows\n"
        f"Number of Columns    =  {f['ncols']}              # Number of columns\n"
        f"Grid spacing         =  {f['grid_spacing']:.3f}         # Grid resolution (m)\n"
        "Point North          =                    # For POINT mode only\n"
        "Point East           =                    # For POINT mode only\n"
    )


if __name__ == "__main__":
    dem = Path(sys.argv[1]) if len(sys.argv) > 1 else paths.ELEV_CLIPPED
    sys.stdout.write(format_area_block(area_fields(dem)))
