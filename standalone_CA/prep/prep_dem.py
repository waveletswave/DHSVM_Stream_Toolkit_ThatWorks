# -*- coding: utf-8 -*-
# =====================================================================
# prep_dem.py  -  build the clipped DEM for an arbitrary watershed
#
# General front door to the pipeline. Reprojects a source DEM (for example a
# USGS 3DEP tile, in any CRS) to the pipeline CRS (paths.EPSG, default
# EPSG:32617) on a clean target-resolution grid, then masks to the watershed
# polygon. The output is the elev_clipped.tif every downstream stage consumes:
# place it at <DHSVM_OUT>/elev_clipped.tif and run slope, bins, hydrology, and
# the rest. The hydrology stage derives the stream cell threshold from
# STREAM_SOURCE_AREA_M2, so drainage density is held fixed as resolution changes.
#
# clip.py stays the CA 28 m byte-match reproducer for the regression test. Use
# clip.py to reproduce the calibrated CA grid against the QGIS reference; use
# this script for any other basin or resolution.
#
# Resampling is BILINEAR, not nearest: a source tile is rarely exactly on the
# target grid, so the reprojection interpolates and mildly smooths the surface,
# which is acceptable for routing.
#
# Masking to the watershed is required: the hydrology stage routes the whole
# raster, so an unmasked rectangle would pull in neighbouring drainage and give
# the wrong contributing area.
#
# Download the source tile first (this script does not fetch it; no network):
#   - py3dep (HyRiver): py3dep.get_dem(watershed_geometry, resolution=<res>)
#   - or The National Map: USGS 1/3 arc-second product (tile n36w084 for CA)
#
# Run:
#   python3 prep_dem.py <src_dem.tif> [<out_elev.tif>] [--res M] \
#                       [--watershed SHP] [--epsg N]
#
# Defaults: out is paths.ELEV_CLIPPED (<DHSVM_OUT>/elev_clipped.tif), res 10 m
# (env DHSVM_DEM_RES), watershed paths.WATERSHED (env DHSVM_WATERSHED), epsg
# paths.EPSG (env DHSVM_EPSG).
# =====================================================================

import os
import sys
import math
import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.features import geometry_mask
import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import WATERSHED, EPSG, ELEV_CLIPPED

DEFAULT_RES = 10.0    # target cell size in metres; override with --res or DHSVM_DEM_RES
PAD_M = 200.0         # extent pad so the basin is not flush to the grid edge
NODATA = -9999.0


def target_grid(watershed_utm, res, pad_m=PAD_M):
    """Target-CRS grid over the watershed bounds plus a small pad, snapped so
    cell edges fall on multiples of res."""
    minx, miny, maxx, maxy = watershed_utm.total_bounds
    minx = math.floor((minx - pad_m) / res) * res
    miny = math.floor((miny - pad_m) / res) * res
    maxx = math.ceil((maxx + pad_m) / res) * res
    maxy = math.ceil((maxy + pad_m) / res) * res
    width = int(round((maxx - minx) / res))
    height = int(round((maxy - miny) / res))
    transform = from_origin(minx, maxy, res, res)
    return transform, width, height, (minx, miny, maxx, maxy)


def prep_dem(src_path, out_path, watershed_path, res, epsg, pad_m=PAD_M):
    """Reproject src_path to EPSG:epsg on a clean res-metre grid and mask it to
    the watershed polygon. Write out_path. Return (valid_cells, total_cells)."""
    dst_crs = "EPSG:%d" % epsg
    wsf = gpd.read_file(watershed_path).to_crs(dst_crs)
    transform, width, height, bounds = target_grid(wsf, res, pad_m)
    print("target: %s  %.3f m  %d rows x %d cols" % (dst_crs, res, height, width))
    print("extent (%s): %s" % (dst_crs, tuple(round(b, 1) for b in bounds)))

    dst = np.full((height, width), NODATA, dtype="float32")
    with rasterio.open(src_path) as src:
        print("source: %s  %dx%d  res~%s  nodata=%s"
              % (src.crs, src.width, src.height,
                 tuple(round(r, 6) for r in src.res), src.nodata))
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=transform,
            dst_crs=dst_crs,
            dst_nodata=NODATA,
            resampling=Resampling.bilinear,
        )

    # mask to the watershed polygon: True outside the polygon -> nodata
    outside = geometry_mask(wsf.geometry, out_shape=(height, width),
                            transform=transform, invert=False)
    dst[outside] = NODATA

    profile = dict(driver="GTiff", dtype="float32", count=1, nodata=NODATA,
                   width=width, height=height, crs=dst_crs, transform=transform,
                   compress="lzw")
    with rasterio.open(out_path, "w", **profile) as f:
        f.write(dst, 1)

    valid = dst[dst != NODATA]
    print("wrote %s" % out_path)
    print("  valid cells %d / %d  (nodata %d)"
          % (valid.size, width * height, width * height - valid.size))
    if valid.size:
        print("  elev range %.1f .. %.1f m" % (float(valid.min()), float(valid.max())))
    return int(valid.size), int(width * height)


def main():
    ap = argparse.ArgumentParser(
        description="Build elev_clipped.tif for an arbitrary watershed: reproject "
                    "a source DEM to the pipeline CRS on a clean target-resolution "
                    "grid, then mask to the watershed polygon. General front door; "
                    "clip.py is the CA 28 m byte-match reproducer.")
    ap.add_argument("src_dem", help="source DEM (e.g. a USGS 3DEP tile), any CRS")
    ap.add_argument("out_dem", nargs="?", default=str(ELEV_CLIPPED),
                    help="output clipped DEM (default paths.ELEV_CLIPPED)")
    ap.add_argument("--res", type=float,
                    default=float(os.environ.get("DHSVM_DEM_RES", DEFAULT_RES)),
                    help="target cell size in metres (default %.1f, env DHSVM_DEM_RES)"
                         % DEFAULT_RES)
    ap.add_argument("--watershed", default=str(WATERSHED),
                    help="watershed polygon (default paths.WATERSHED, env DHSVM_WATERSHED)")
    ap.add_argument("--epsg", type=int, default=EPSG,
                    help="target EPSG (default paths.EPSG = %d, env DHSVM_EPSG)" % EPSG)
    args = ap.parse_args()

    prep_dem(args.src_dem, args.out_dem, args.watershed, args.res, args.epsg)


if __name__ == "__main__":
    main()
