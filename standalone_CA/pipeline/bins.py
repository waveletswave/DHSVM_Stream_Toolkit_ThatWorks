"""Stage 3: DHSVM base-map binaries from the clipped DEM.

Standalone port of qgis_CA's dem_to_dhsvm_bins.py. Reads elev_clipped.tif and
writes dem.bin (Float32), mask/soil/veg.bin (Int8, 1 inside basin / 0 outside),
and the uniform soil-depth baselines (Float32). These are raw .tofile() dumps:
C-order, native (little-endian on DCC) byte order, no header. The valid-pixel
mask is taken from the DEM nodata, which clip already aligned to the qgis_CA
reference byte-for-byte, so these binaries are expected to match byte-identical.

Logic mirrors the qgis_CA script exactly; only paths differ. soil.bin and
veg.bin here are the uniform placeholders (type 1 inside basin) that the qgis_CA
base-map script writes; the per-class veg_bs binaries come from the veg-class
scripts, a separate stage.
"""
import sys
sys.path.insert(0, ".")

import numpy as np
from osgeo import gdal

from paths import ELEV_CLIPPED, OUT

DHSVM_NODATA = -9999.0
UNIFORM_SOIL_DEPTHS = [2.0, 2.5, 3.0, 3.5, 4.0]

BIN_OUT = OUT / "bins"
BIN_OUT.mkdir(parents=True, exist_ok=True)


def generate_basemaps():
    print(f"[bins] reading {ELEV_CLIPPED}")
    ds = gdal.Open(str(ELEV_CLIPPED))
    band = ds.GetRasterBand(1)
    dem_arr = band.ReadAsArray().astype(np.float32)
    dem_nd = band.GetNoDataValue()

    # Master valid mask from DEM nodata — identical logic to qgis_CA.
    valid_mask = (dem_arr != dem_nd) & (~np.isnan(dem_arr))
    print(f"[bins] valid pixels: {int(valid_mask.sum())} / {dem_arr.size}")

    # DEM: Float32, nodata-filled outside.
    dem_out = np.full(dem_arr.shape, DHSVM_NODATA, dtype=np.float32)
    dem_out[valid_mask] = dem_arr[valid_mask]

    # MASK / SOIL / VEG: Int8, 1 inside basin, 0 outside.
    mask_out = np.zeros(dem_arr.shape, dtype=np.int8)
    mask_out[valid_mask] = 1
    soil_out = np.zeros(dem_arr.shape, dtype=np.int8)
    soil_out[valid_mask] = 1
    veg_out = np.zeros(dem_arr.shape, dtype=np.int8)
    veg_out[valid_mask] = 1

    dem_out.tofile(BIN_OUT / "dem.bin")
    mask_out.tofile(BIN_OUT / "mask.bin")
    soil_out.tofile(BIN_OUT / "soil.bin")
    veg_out.tofile(BIN_OUT / "veg.bin")
    print("[bins] wrote dem/mask/soil/veg .bin")

    for depth in UNIFORM_SOIL_DEPTHS:
        depth_out = np.full(dem_arr.shape, DHSVM_NODATA, dtype=np.float32)
        depth_out[valid_mask] = depth
        fn = f"soildepth_uniform_{depth:.1f}m.bin"
        depth_out.tofile(BIN_OUT / fn)
        print(f"[bins] wrote {fn}")

    ds = None


if __name__ == "__main__":
    generate_basemaps()
