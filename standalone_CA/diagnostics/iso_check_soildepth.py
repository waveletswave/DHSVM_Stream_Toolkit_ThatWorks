# -*- coding: utf-8 -*-
# iso_check_soildepth.py
# One-off diagnostic for the soildepth float32 diff (412 cells, max ~4.77e-7).
#
# Two questions, isolated:
#   RESULT : recompute soildepth from the REFERENCE inputs (ref elev + ref
#            filled slope + ref flow_acc) and compare to the reference bin.
#            diff=0 means the formula + nodata handling reproduce the ref bin
#            exactly, so any standalone diff comes only from the inputs.
#   SLOPE  : compare the reference filled slope to the standalone filled slope
#            directly. >0 here means the slope inputs already differ at the
#            float32 last bit, which the soildepth power terms amplify to 1 ULP.
#
# Run:  python3 iso_check_soildepth.py

import sys
from pathlib import Path
import numpy as np
from osgeo import gdal

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
import paths

# Reference (QGIS) and standalone output roots, from paths.py. Subpaths are
# joined below; both roots honour DHSVM_REF / DHSVM_OUT.
REF = str(paths.REF)
STD = str(paths.OUT)


def read(fp):
    ds = gdal.Open(fp)
    if ds is None:
        raise FileNotFoundError(fp)
    b = ds.GetRasterBand(1)
    return b.ReadAsArray().astype(np.float32), b.GetNoDataValue()


# ----- reference inputs -----
elev, end = read(f"{REF}/elev_clipped.tif")
slope, snd = read(f"{REF}/stream_slope_filled.tif")        # filled slope, ref root
fac, fnd = read(f"{REF}/Intermediate_GIS/flow_acc.tif")

# ----- PNNL formula (verbatim) -----
MIN_DEPTH, MAX_DEPTH = 2.0, 6.0
WT_SLOPE, WT_SOURCE, WT_ELEV = 0.7, 0.0, 0.3
MAX_SLOPE, MAX_SOURCE, MAX_ELEV = 30.0, 100000.0, 1500.0
POW_SLOPE, POW_SOURCE, POW_ELEV = 0.25, 1.0, 0.75
ND = -9999.0

vm = (elev != end) & (~np.isnan(elev))
ec = np.where(vm & (elev >= 0), elev, 0.0)
sc = np.where(vm & (slope != snd) & (~np.isnan(slope)) & (slope >= 0), slope, 0.0)
fc = np.where(vm & (fac != fnd) & (~np.isnan(fac)), np.abs(fac), 0.0)
el = np.clip(ec, 0, MAX_ELEV)
sl = np.clip(sc, 0, MAX_SLOPE)
fl = np.clip(fc, 0, MAX_SOURCE)
ts = WT_SLOPE * (1 - (sl / MAX_SLOPE) ** POW_SLOPE)
to = WT_SOURCE * ((fl / MAX_SOURCE) ** POW_SOURCE)
te = WT_ELEV * (1 - (el / MAX_ELEV) ** POW_ELEV)
dc = MIN_DEPTH + (MAX_DEPTH - MIN_DEPTH) * (ts + to + te)
fd = np.full(elev.shape, ND, dtype=np.float32)
fd[vm] = np.clip(dc[vm], MIN_DEPTH, MAX_DEPTH)

ref = np.fromfile(f"{REF}/DHSVM_input_binaries/soildepth.bin",
                  dtype=np.float32).reshape(elev.shape)
d = np.abs(fd.astype("f8") - ref.astype("f8"))
print("RESULT cells_differ =", int(np.sum(d > 0)), " max_abs_diff =", d.max())

# ----- slope ref vs standalone -----
std_slope, ssnd = read(f"{STD}/slope_filled.tif")
sd = np.abs(slope.astype("f8") - std_slope.astype("f8"))
print("SLOPE ref-vs-standalone cells_differ =", int(np.sum(sd > 0)),
      " max_abs_diff =", sd.max())
