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
#   2. stamps EPSG:32617 back onto the exported GeoTIFFs with rasterio in
#      r+ mode -- r.out.gdal cannot write the CRS tag on this install
#      (gcs.csv missing), exactly as in the slope stage.
#
# Paths use a local block for now; fold into paths.py in the parameterization
# pass (step 9). DEM_TIF is the clipped DEM from the clip stage; it drives the
# region and the 74x82 grid.
#
# Validation (separate, against qgis_CA_ref/Intermediate_GIS):
#   flow_acc.tif / flow_dir.tif / stream_raster.tif  -> raster data-equivalence
#   streamfile.shp (41 line features)                -> feature count + geometry
# =====================================================================

import subprocess
from pathlib import Path
from rasterio.crs import CRS
import rasterio

# ----------------------------- paths (fold into paths.py later) -------------
OUT_DIR  = Path("/work/ys451/dhsvm_ca/standalone_dev/outputs")
DEM_TIF  = OUT_DIR / "elev_clipped.tif"          # from the clip stage
HYDRO_SH = Path(__file__).resolve().parent / "run_hydrology_grass.sh"

EPSG = 32617

# GeoTIFF exports that need a CRS stamp (vector shp carries its own .prj).
RASTER_OUTPUTS = ["flow_acc.tif", "flow_dir.tif", "stream_raster.tif"]


def run_grass_chain():
    if not DEM_TIF.exists():
        raise FileNotFoundError(f"[ERROR] clipped DEM not found: {DEM_TIF}")
    print(f"[step] GRASS hydrology chain (r.watershed -> r.stream.extract -> r.to.vect)")
    print(f"       DEM={DEM_TIF.name}  out={OUT_DIR}")
    subprocess.run(
        ["bash", str(HYDRO_SH), str(DEM_TIF), str(OUT_DIR)],
        check=True,
    )


def stamp_crs():
    """r.out.gdal leaves crs=None on this install; rewrite the header only."""
    crs = CRS.from_epsg(EPSG)
    for name in RASTER_OUTPUTS:
        p = OUT_DIR / name
        if not p.exists():
            print(f"  [warn] missing export, skip CRS stamp: {name}")
            continue
        with rasterio.open(p, "r+") as s:
            if s.crs is None:
                s.crs = crs
                print(f"  -> stamped EPSG:{EPSG} on {name}")
            else:
                print(f"  -> {name} already has CRS {s.crs}, left as is")


def run():
    print("\n=======================================================")
    print("  DHSVM HYDROLOGY CORE (standalone)")
    print("=======================================================")
    run_grass_chain()
    stamp_crs()
    print("=======================================================\n")


if __name__ == "__main__":
    run()
