"""Stage 2 (Python part): finish the slope raster and validate against qgis_CA.

Takes the raw r.slope.aspect output (slope_raw.tif, degrees, produced by
run_slope_grass.sh), stamps the CRS back (r.out.gdal cannot write it on this
install, the gcs.csv issue), then conditions it with slope_fill using the
clipped DEM — exactly as qgis_CA does (fill_slope_nodata overwrites in place
there; here we write slope_filled.tif). Finally compares against the qgis_CA
reference stream_slope_filled.tif. Both sides use the same slope_fill code, so
any residual difference comes from r.slope.aspect itself, not conditioning.
"""
import sys
sys.path.insert(0, ".")

import rasterio
from rasterio.crs import CRS

from paths import ELEV_CLIPPED, SLOPE_RAW, SLOPE_FILLED, REF_SLOPE, EPSG
from slope_fill import fill_slope_nodata
from compare import compare_rasters


def stamp_crs(path, epsg):
    with rasterio.open(path, "r+") as s:
        if s.crs is None:
            s.crs = CRS.from_epsg(epsg)
    with rasterio.open(path) as s:
        print(f"  CRS after stamp: {s.crs}")


def main():
    # 1. Stamp CRS onto the raw GRASS export.
    print("[slope] stamping CRS on raw r.slope.aspect output")
    stamp_crs(SLOPE_RAW, EPSG)

    # 2. Condition: fill r.slope.aspect's boundary/void NoData ring using the
    #    clipped DEM as the valid-extent reference. Writes SLOPE_FILLED (does
    #    not overwrite the raw output, so we can inspect both).
    print("[slope] conditioning with slope_fill (fill NoData by neighbour mean)")
    stats = fill_slope_nodata(str(SLOPE_RAW), str(ELEV_CLIPPED),
                              out_path=str(SLOPE_FILLED), verbose=True)
    print(f"  fill stats: {stats}")

    # 3. Validate against qgis_CA reference.
    print("\n[slope] comparing against qgis_CA reference")
    ok = compare_rasters(str(SLOPE_FILLED), str(REF_SLOPE), atol=0.0)
    print("\nSLOPE MATCH (atol=0):", ok)
    if not ok:
        # slope is an algorithm-vs-algorithm comparison; show tolerance bands
        # so we can judge whether differences are small+boundary (acceptable)
        # or large+interior (a real parameter problem).
        print("\n[slope] tolerance sweep (same r.slope.aspect both sides, so")
        print("        differences should be ~floating-point, not structural):")
        for tol in (1e-6, 1e-4, 1e-2, 0.1, 0.5):
            ok_t = compare_rasters(str(SLOPE_FILLED), str(REF_SLOPE), atol=tol)
            print(f"  atol={tol}: MATCH={ok_t}")
            print()


if __name__ == "__main__":
    main()
