# -*- coding: utf-8 -*-
# =====================================================================
# hydrology.py  -  standalone hydrology core (DCC, path a, option i)
#
# Stage: flow accumulation + drainage + stream raster + stream vector.
# Reproduces the three grass7: processing.run calls in prep_dhsvm_inputs.py:
#   r.watershed  ->  r.stream.extract  ->  r.to.vect (fallback)
# all chained in a single GRASS location (see run_hydrology_grass.sh).
#
# This wrapper:
#   1. runs the GRASS chain via the shim (the .sh does the GRASS work),
#   2. stamps EPSG back onto the exported GeoTIFFs with rasterio in r+ mode --
#      r.out.gdal cannot write the CRS tag on this install (gcs.csv missing),
#      exactly as in the slope stage.
#
# Paths via paths.py. DEM_TIF is the clipped DEM from the clip stage; it drives
# the region and the 74x82 grid. The GRASS shim path is passed to the .sh.
#
# Validation (separate, against REF_INTERMED):
#   flow_acc.tif / flow_dir.tif / stream_raster.tif  -> raster data-equivalence
#   streamfile.shp (41 line features)                -> feature count + geometry
# =====================================================================

import subprocess
from pathlib import Path
from rasterio.crs import CRS
import rasterio

from paths import ELEV_CLIPPED, OUT, EPSG, FLOW_ACC, FLOW_DIR, STREAM_RASTER

DEM_TIF  = ELEV_CLIPPED
OUT_DIR  = OUT
HYDRO_SH = Path(__file__).resolve().parent / "run_hydrology_grass.sh"

# GeoTIFF exports that need a CRS stamp (vector shp carries its own .prj).
RASTER_OUTPUTS = [FLOW_ACC, FLOW_DIR, STREAM_RASTER]


def run_grass_chain():
    if not DEM_TIF.exists():
        raise FileNotFoundError(f"[ERROR] clipped DEM not found: {DEM_TIF}")
    print(f"[step] GRASS hydrology chain (r.watershed -> r.stream.extract -> r.to.vect)")
    print(f"       DEM={DEM_TIF.name}  out={OUT_DIR}")
    # NOTE: run_hydrology_grass.sh still hardcodes the GRASS shim path internally.
    # Parameterizing the .sh (read shim from an env var / arg) is a separate
    # follow-up; the .sh contract here remains (dem_tif, out_dir).
    subprocess.run(
        ["bash", str(HYDRO_SH), str(DEM_TIF), str(OUT_DIR)],
        check=True,
    )


def stamp_crs():
    """r.out.gdal leaves crs=None on this install; rewrite the header only."""
    crs = CRS.from_epsg(EPSG)
    for p in RASTER_OUTPUTS:
        if not p.exists():
            print(f"  [warn] missing export, skip CRS stamp: {p.name}")
            continue
        with rasterio.open(p, "r+") as s:
            if s.crs is None:
                s.crs = crs
                print(f"  -> stamped EPSG:{EPSG} on {p.name}")
            else:
                print(f"  -> {p.name} already has CRS {s.crs}, left as is")


def run():
    print("\n=======================================================")
    print("  DHSVM HYDROLOGY CORE (standalone)")
    print("=======================================================")
    run_grass_chain()
    stamp_crs()
    print("=======================================================\n")


if __name__ == "__main__":
    run()
