# -*- coding: utf-8 -*-
# =====================================================================
# fetch_dem.py  -  fetch a source DEM for a watershed from USGS 3DEP
#
# Networked step 0 for the general pipeline. Pulls a DEM over the watershed
# bounding box from the 3DEP service via py3dep and writes a GeoTIFF. That
# GeoTIFF is the source DEM for prep_dem.py, which reprojects it to the pipeline
# grid (paths.EPSG) at the target resolution and masks to the polygon.
#
# Kept separate from prep_dem.py on purpose: prep_dem and every downstream stage
# stay no-network and run on a compute node. fetch_dem is the one networked step.
# Run it on a node with outbound internet (a login node, not a compute node) to
# produce the source tile, then run the pipeline offline. py3dep can live in its
# own env, since its only output is this GeoTIFF.
#
# 3DEP staged resolutions are 10, 30, and 60 m (10 m is the finest over CONUS).
# Fetch the finest sensible source (default 10 m); prep_dem then resamples to
# whatever target grid resolution you want, so fetch res and pipeline res are
# independent. The bbox is padded by a margin so prep_dem's own pad plus the
# bilinear halo has source coverage to the polygon edge.
#
# py3dep must be installed in the env and reachable from the node. HyRiver now
# also ships Seamless3DEP, a lighter alternative with the same data source.
#
# Run (on a networked node):
#   python3 fetch_dem.py [<out_src_dem.tif>] [--res 10|30|60] \
#                        [--watershed SHP] [--margin M]
# Then feed the output to prep_dem.py as the source DEM.
# =====================================================================

import os
import sys
import argparse
from pathlib import Path

import geopandas as gpd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import WATERSHED

DEFAULT_RES = 10        # 3DEP staged resolutions: 10, 30, 60 m
MARGIN_M = 400.0        # pad beyond the polygon so prep_dem pad + bilinear halo has coverage


def fetch_dem(watershed_path, out_path, res, margin_m=MARGIN_M):
    """Fetch a 3DEP DEM over the (padded) watershed bbox and write a GeoTIFF.
    Source CRS is whatever 3DEP returns; prep_dem.py reprojects it later."""
    import py3dep
    import rioxarray  # noqa: F401  registers the .rio accessor on the DataArray

    wsf = gpd.read_file(watershed_path)
    if wsf.crs is None:
        sys.exit("watershed has no CRS; set one before fetching")
    minx, miny, maxx, maxy = wsf.total_bounds
    bbox = (minx - margin_m, miny - margin_m, maxx + margin_m, maxy + margin_m)
    print("fetch: 3DEP res=%d m  geom_crs=%s" % (res, wsf.crs))
    print("  bbox (+%.0f m): %s" % (margin_m, tuple(round(v, 1) for v in bbox)))

    dem = py3dep.get_dem(bbox, resolution=res, crs=wsf.crs)
    dem.rio.to_raster(out_path)

    with rasterio.open(out_path) as s:
        print("wrote %s" % out_path)
        print("  crs: %s  size: %dx%d  res: %s"
              % (s.crs, s.width, s.height, tuple(round(r, 6) for r in s.res)))
        a = s.read(1, masked=True)
        print("  elev range %.1f .. %.1f m" % (float(a.min()), float(a.max())))


def main():
    ap = argparse.ArgumentParser(
        description="Fetch a source DEM for a watershed from USGS 3DEP and write a "
                    "GeoTIFF for prep_dem.py. Networked step 0; run on a node with "
                    "internet. The pipeline itself stays no-network.")
    ap.add_argument("out_dem", nargs="?", default="src_dem_3dep.tif",
                    help="output source DEM GeoTIFF (default src_dem_3dep.tif)")
    ap.add_argument("--res", type=int,
                    default=int(os.environ.get("DHSVM_FETCH_RES", DEFAULT_RES)),
                    help="3DEP source resolution m, one of 10/30/60 "
                         "(default %d, env DHSVM_FETCH_RES)" % DEFAULT_RES)
    ap.add_argument("--watershed", default=str(WATERSHED),
                    help="watershed polygon (default paths.WATERSHED, env DHSVM_WATERSHED)")
    ap.add_argument("--margin", type=float, default=MARGIN_M,
                    help="pad beyond the polygon in metres (default %.0f)" % MARGIN_M)
    args = ap.parse_args()
    fetch_dem(args.watershed, args.out_dem, args.res, args.margin)


if __name__ == "__main__":
    main()
