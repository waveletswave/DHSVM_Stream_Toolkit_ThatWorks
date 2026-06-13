# -*- coding: utf-8 -*-
# =====================================================================
# prep_dem_10m.py  -  build the 10 m clipped DEM for the resolution test
#
# Reprojects a USGS 3DEP 1/3 arc-second DEM (~10 m native, usually EPSG:4326)
# to the pipeline CRS EPSG:32617 (UTM 17N) on a clean 10 m grid, then masks to
# the watershed polygon. The output is the 10 m equivalent of elev_clipped.tif:
# feed it straight to drop_analysis.py, and (for the network density look) place
# it as <10m_OUT>/elev_clipped.tif and run the hydrology stage, which derives
# the cell threshold from the same STREAM_SOURCE_AREA_M2 (47571.5 m2). At 10 m
# that is ~476 cells, so drainage density is held fixed across resolutions.
#
# Resampling is BILINEAR, not nearest: 1/3 arc-second is not exactly 10 m, so
# the reprojection to a clean 10 m UTM grid interpolates and mildly smooths the
# surface, which is acceptable for routing.
#
# Masking to the watershed is required: drop_analysis routes the whole raster,
# so an unmasked rectangle would pull in neighbouring drainage and give the
# wrong A_c.
#
# Download the source tile first (this script does not fetch it; no network):
#   - py3dep (HyRiver): py3dep.get_dem(watershed_geometry, resolution=10)
#   - or The National Map: USGS 1/3 arc-second product, tile n36w084
#
# Run:  python3 prep_dem_10m.py <src_3dep.tif> [<out_elev_10m.tif>]
# =====================================================================

import sys
import math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.features import geometry_mask
import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import WATERSHED, EPSG

DST_CRS = "EPSG:%d" % EPSG       # 32617
RES = 10.0                       # target cell size, metres
PAD_M = 200.0                    # small extent pad so the basin is not flush to the edge
NODATA = -9999.0


def target_grid(watershed_utm):
    """10 m grid in the target CRS over the watershed bounds + a small pad,
    snapped so cell edges fall on multiples of RES."""
    minx, miny, maxx, maxy = watershed_utm.total_bounds
    minx = math.floor((minx - PAD_M) / RES) * RES
    miny = math.floor((miny - PAD_M) / RES) * RES
    maxx = math.ceil((maxx + PAD_M) / RES) * RES
    maxy = math.ceil((maxy + PAD_M) / RES) * RES
    width = int(round((maxx - minx) / RES))
    height = int(round((maxy - miny) / RES))
    transform = from_origin(minx, maxy, RES, RES)
    return transform, width, height, (minx, miny, maxx, maxy)


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: prep_dem_10m.py <src_3dep.tif> [<out_elev_10m.tif>]")
    src_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "elev_clipped_10m.tif"

    wsf = gpd.read_file(WATERSHED).to_crs(DST_CRS)
    transform, width, height, bounds = target_grid(wsf)
    print("target: %s  %.0f m  %d rows x %d cols" % (DST_CRS, RES, height, width))
    print("extent (UTM17): %s" % (tuple(round(b, 1) for b in bounds),))

    dst = np.full((height, width), NODATA, dtype="float32")
    with rasterio.open(src_path) as src:
        print("source: %s  %dx%d  res~%s  nodata=%s"
              % (src.crs, src.width, src.height, tuple(round(r, 6) for r in src.res), src.nodata))
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=transform,
            dst_crs=DST_CRS,
            dst_nodata=NODATA,
            resampling=Resampling.bilinear,
        )

    # mask to the watershed polygon: True outside -> nodata
    outside = geometry_mask(wsf.geometry, out_shape=(height, width),
                            transform=transform, invert=False)
    dst[outside] = NODATA

    profile = dict(driver="GTiff", dtype="float32", count=1, nodata=NODATA,
                   width=width, height=height, crs=DST_CRS, transform=transform,
                   compress="lzw")
    with rasterio.open(out_path, "w", **profile) as f:
        f.write(dst, 1)

    valid = dst[dst != NODATA]
    print("wrote %s" % out_path)
    print("  valid cells %d / %d  (nodata %d)"
          % (valid.size, width * height, width * height - valid.size))
    if valid.size:
        print("  elev range %.1f .. %.1f m" % (float(valid.min()), float(valid.max())))


if __name__ == "__main__":
    main()
