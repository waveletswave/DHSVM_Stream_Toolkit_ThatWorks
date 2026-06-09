"""Stage 1: clip the source DEM to the watershed.

Standalone replacement for qgis_CA's gdal:cliprasterbymasklayer. The output
grid is locked to the qgis_CA reference window (its bounds, on the shared
source grid) so the standalone produces the exact same raster dimensions as
qgis_CA. This avoids the crop-to-cutline rounding-direction difference between
rasterio.mask(crop=True) and GDAL warp, which otherwise leaves a one-pixel
boundary ring (76x84 vs the reference 74x82). Pixels inside the window but
outside the watershed polygon are set to nodata=-9999, matching the reference.
"""
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.features import geometry_mask

from paths import SRC_DEM, WATERSHED, ELEV_CLIPPED, REF_ELEV_CLIPPED

NODATA = -9999.0


def clip_dem():
    # Target window = the reference raster's footprint, snapped to the source grid.
    with rasterio.open(REF_ELEV_CLIPPED) as ref:
        ref_bounds = ref.bounds
        ref_shape = ref.shape  # (rows, cols) we must reproduce

    with rasterio.open(SRC_DEM) as src:
        win = from_bounds(*ref_bounds, transform=src.transform)
        # round to whole pixels so the window lands exactly on the source grid
        win = win.round_offsets().round_lengths()
        data = src.read(1, window=win).astype("float32")
        win_transform = src.window_transform(win)

    # Mask out pixels whose center falls outside the watershed polygon.
    shapes = gpd.read_file(WATERSHED)
    geoms = [g.__geo_interface__ for g in shapes.geometry]
    outside = geometry_mask(
        geoms,
        out_shape=data.shape,
        transform=win_transform,
        all_touched=False,
        invert=False,  # True where OUTSIDE the geometry
    )
    data[outside] = NODATA

    meta = {
        "driver": "GTiff",
        "height": data.shape[0],
        "width": data.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": rasterio.crs.CRS.from_epsg(32617),
        "transform": win_transform,
        "nodata": NODATA,
    }
    with rasterio.open(ELEV_CLIPPED, "w", **meta) as dst:
        dst.write(data, 1)

    print(f"wrote {ELEV_CLIPPED}")
    with rasterio.open(ELEV_CLIPPED) as s:
        a = s.read(1)
        print("  shape:", s.shape, "| (ref shape:", ref_shape, ")")
        print("  dtype:", a.dtype, "| nodata:", s.nodata)
        print("  bounds:", tuple(round(b, 3) for b in s.bounds))
        finite = a[a != NODATA]
        if finite.size:
            print("  data min/max:", float(finite.min()), float(finite.max()))


if __name__ == "__main__":
    clip_dem()
