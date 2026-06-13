# -*- coding: utf-8 -*-
# iso_check2_soildepth.py
# Follow-up diagnostic. Two clean checks, nodata/nan handled properly.
#
#   A) ref filled slope vs standalone filled slope, masked to valid cells only.
#      Tells whether the slope INPUTS actually differ (the prior check was
#      polluted by nan in unmasked subtraction).
#   B) recompute soildepth from REF inputs vs ref bin, then characterize the
#      412 differing cells: are they interior or boundary, what is the diff
#      distribution, and is the ref bin itself float32-granular there.
#
# Run:  python3 iso_check2_soildepth.py

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


# ---------- A) slope ref vs standalone, masked ----------
ref_slope, rsnd = read(f"{REF}/stream_slope_filled.tif")
std_slope, ssnd = read(f"{STD}/slope_filled.tif")

print("=== A) filled slope: ref vs standalone ===")
print(f"  ref  nodata={rsnd}  nan_count={int(np.isnan(ref_slope).sum())}")
print(f"  std  nodata={ssnd}  nan_count={int(np.isnan(std_slope).sum())}")

# valid where both are finite and not their respective nodata
rvalid = np.isfinite(ref_slope) & (ref_slope != rsnd if rsnd is not None else True)
svalid = np.isfinite(std_slope) & (std_slope != ssnd if ssnd is not None else True)
both = rvalid & svalid
print(f"  ref valid cells={int(rvalid.sum())}  std valid cells={int(svalid.sum())}  both={int(both.sum())}")
print(f"  valid-mask mismatch cells (one valid, other not) = {int((rvalid ^ svalid).sum())}")

sd = np.abs(ref_slope[both].astype('f8') - std_slope[both].astype('f8'))
if sd.size:
    print(f"  slope diff on both-valid: cells_differ={int((sd>0).sum())}  max={sd.max():.3e}")
    if (sd > 0).any():
        print(f"    p50={np.percentile(sd[sd>0],50):.3e}  p99={np.percentile(sd[sd>0],99):.3e}")

# ---------- B) recompute from REF inputs, characterize the 412 ----------
elev, end = read(f"{REF}/elev_clipped.tif")
fac, fnd = read(f"{REF}/Intermediate_GIS/flow_acc.tif")
slope, snd = ref_slope, rsnd

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
d = np.abs(fd.astype('f8') - ref.astype('f8'))
diffmask = d > 0

print("\n=== B) soildepth recompute (REF inputs) vs ref bin ===")
print(f"  cells_differ={int(diffmask.sum())}  max={d.max():.3e}")

# boundary vs interior: a cell is boundary if any of its 4-neighbours is nodata/off-grid
valid2d = vm
nb = np.zeros_like(valid2d)
nb[1:, :]  |= ~valid2d[:-1, :]
nb[:-1, :] |= ~valid2d[1:, :]
nb[:, 1:]  |= ~valid2d[:, :-1]
nb[:, :-1] |= ~valid2d[:, 1:]
boundary = valid2d & nb
nd_on_boundary = int((diffmask & boundary).sum())
nd_interior = int((diffmask & valid2d & ~boundary).sum())
print(f"  of differing cells: on boundary ring={nd_on_boundary}  interior={nd_interior}")

# are the differences exactly 1 ULP of float32 at those magnitudes?
vals = ref[diffmask].astype('f8')
ulp = np.spacing(vals.astype('f4')).astype('f8')   # float32 ULP at each value
within_1ulp = int((d[diffmask] <= ulp + 1e-12).sum())
print(f"  differing cells within 1 float32 ULP = {within_1ulp} / {int(diffmask.sum())}")
print(f"  example differing cells (flat idx, ref, recompute, diff):")
idx = np.where(diffmask.ravel())[0][:5]
fr = fd.ravel(); rr = ref.ravel()
for i in idx:
    print(f"    {i}: ref={rr[i]:.7f} recompute={fr[i]:.7f} diff={abs(float(rr[i])-float(fr[i])):.3e}")
