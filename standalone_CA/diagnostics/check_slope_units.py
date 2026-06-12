"""Determine the slope unit empirically by comparing slope_filled against the
slope recomputed from the DEM in three candidate units. Smallest median
absolute difference (over finite cells only) wins. Independent of any
generator code.
"""
import sys
import numpy as np
import rasterio

run = sys.argv[1] if len(sys.argv) > 1 else "."
with rasterio.open(f"{run}/elev_clipped.tif") as d:
    z = d.read(1, masked=True).astype(float)
    res = d.res[0]
with rasterio.open(f"{run}/slope_filled.tif") as s:
    sl = s.read(1, masked=True).astype(float)

zf = np.ma.filled(z, np.nan)
gy, gx = np.gradient(zf, res)
grad = np.hypot(gx, gy)                       # rise/run = tan(slope)
cand = {
    "degrees":  np.degrees(np.arctan(grad)),
    "percent":  grad * 100.0,
    "fraction": grad,
}
slf = np.ma.filled(sl, np.nan)
valid = np.isfinite(slf)
print(f"slope_filled range: {np.nanmin(slf):.3g} .. {np.nanmax(slf):.3g}"
      f"  median {np.nanmedian(slf):.3g}")
print("median |slope_filled - recomputed|  (finite cells only):")
best, bestdiff = None, np.inf
for name, ref in cand.items():
    m = valid & np.isfinite(ref)
    if not m.any():
        print(f"  {name:9s} no overlap")
        continue
    diff = float(np.median(np.abs(slf[m] - ref[m])))
    print(f"  {name:9s} {diff:.4g}    (n={int(m.sum())})")
    if diff < bestdiff:
        best, bestdiff = name, diff
print(f"=> slope_filled is most consistent with: {best.upper() if best else 'UNKNOWN'}")
