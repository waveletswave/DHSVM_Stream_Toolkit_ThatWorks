# -*- coding: utf-8 -*-
# =====================================================================
# vector_attrs.py  -  standalone port of prep's per-segment attribute steps
#
# Sub-step B (part 1) of the #6 vector-IO rebuild, option (ii): read the
# hydrology streamfile.shp (geometry only: cat/value/label) and populate the
# per-segment attributes that prep_dhsvm_inputs.py computes in-place, writing a
# NEW file streamfile_attr.shp (the geometry file is left untouched).
#
# Ports three prep blocks, sampling semantics copied verbatim:
#   _ensure_fields_rowcol_len  (prep line 275): arcid, Shape_Leng, Row, Col
#       Row/Col use TOP-LEFT origin int((ymax-c.y())/py)  (Tier D2 fix).
#       arcid = 1-based iteration index. centroid-based row/col.
#   _sample_mean_slope_deg     (prep line 334): slope_deg
#       BASE_SLOPE_SAMPLES=12 points at i/(12+1)*L, max(0,val), MEAN.
#   meanmsq                     (prep line 375):
#       15 points at k/(15-1)*L (endpoints included), take values >0, then
#       the MIDPOINT-position value vals[len//2] (NOT the mean); empty -> cell_area.
#
# geopandas/shapely/rasterio replace QgsVectorLayer/QgsGeometry/provider.sample.
# Paths use a local block; fold into paths.py in the parameterization pass.
#
# Run:  python3 vector_attrs.py

import numpy as np
import geopandas as gpd
import rasterio
from shapely.geometry import LineString, MultiLineString
from pathlib import Path

from paths import (STREAMFILE, STREAMFILE_ATTR, ELEV_CLIPPED, SLOPE_FILLED, FLOW_ACC)

STREAM_IN  = STREAMFILE            # hydrology output (geometry only)
STREAM_OUT = STREAMFILE_ATTR       # this stage's output (with attrs)
ELEV_TIF   = ELEV_CLIPPED          # DEM: grid origin + cell size
SLOPE_TIF  = SLOPE_FILLED          # filled slope (degrees)
FACC_TIF   = FLOW_ACC              # flow accumulation

BASE_SLOPE_SAMPLES = 12                            # prep constant


def _longest_line(geom):
    """Mirror prep _line_coords: for MultiLineString take the longest part."""
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, MultiLineString):
        parts = list(geom.geoms)
        return max(parts, key=lambda g: g.length) if parts else None
    return geom


def _sample_point(ds, band_arr, gt_inv, x, y, default=0.0):
    """Sample raster value at map coord (x,y). Returns default if off-grid/nodata."""
    # gt_inv maps (x,y)->(col,row) fractional; floor to cell (GDAL pixel-is-area).
    col = int(gt_inv[0] + gt_inv[1] * x + gt_inv[2] * y)
    row = int(gt_inv[3] + gt_inv[4] * x + gt_inv[5] * y)
    h, w = band_arr.shape
    if 0 <= row < h and 0 <= col < w:
        v = band_arr[row, col]
        if v is not None and np.isfinite(v):
            return float(v)
    return default


def main():
    for tag, p in [("streamfile", STREAM_IN), ("elev", ELEV_TIF),
                   ("slope", SLOPE_TIF), ("flow_acc", FACC_TIF)]:
        if not p.exists():
            raise FileNotFoundError(f"[error] missing {tag}: {p}")

    gdf = gpd.read_file(STREAM_IN)
    print(f"[vector_attrs] read {len(gdf)} features from {STREAM_IN.name}")

    # --- DEM grid geometry (top-left origin), matching prep ---
    dem = rasterio.open(ELEV_TIF)
    px = abs(dem.transform.a)
    py = abs(dem.transform.e)
    xmin = dem.bounds.left
    ymax = dem.bounds.top
    cell_area = px * py                                  # prep cell_area (same cell size)

    # --- raster arrays + inverse geotransforms for sampling (GDAL convention) ---
    from osgeo import gdal
    def _open(path):
        ds = gdal.Open(str(path))
        arr = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
        gt = ds.GetGeoTransform()
        inv = gdal.InvGeoTransform(gt)
        return ds, arr, inv
    sds, sarr, sinv = _open(SLOPE_TIF)
    fds, farr, finv = _open(FACC_TIF)

    arcid = []
    shape_leng = []
    row_l = []
    col_l = []
    slope_deg_l = []
    meanmsq_l = []

    for i, geom in enumerate(gdf.geometry, start=1):
        g = _longest_line(geom)
        L = g.length if g is not None else 0.0

        # --- rowcol/len: centroid, top-left origin (Tier D2) ---
        if g is not None:
            c = geom.centroid
            cx, cy = c.x, c.y
        else:
            cx, cy = xmin, ymax
        row = int((ymax - cy) / py) if py > 0 else 0
        col = int((cx - xmin) / px) if px > 0 else 0
        arcid.append(i)
        shape_leng.append(float(L))
        row_l.append(int(row))
        col_l.append(int(col))

        # --- slope_deg: 12 samples i/(12+1)*L, max(0,val), MEAN ---
        if g is None or L <= 0:
            slope_deg_l.append(0.0)
        else:
            vals = []
            for k in range(1, BASE_SLOPE_SAMPLES + 1):
                d = (k / (BASE_SLOPE_SAMPLES + 1.0)) * L
                p = g.interpolate(d)
                v = _sample_point(sds, sarr, sinv, p.x, p.y, default=None)
                if v is not None:
                    vals.append(max(0.0, float(v)))
            slope_deg_l.append(sum(vals) / len(vals) if vals else 0.0)

        # --- meanmsq: 15 samples k/(15-1)*L incl endpoints, >0 -> midpoint value ---
        if g is None or L <= 0:
            meanmsq_l.append(float(cell_area))
        else:
            n = 15
            vals = []
            for k in range(n):
                d = (k / (n - 1.0)) * L
                p = g.interpolate(d)
                v = _sample_point(fds, farr, finv, p.x, p.y, default=None)
                if v is not None:
                    vals.append(float(v))
            vals = [v for v in vals if v is not None and v > 0]
            val = vals[len(vals) // 2] if vals else cell_area
            meanmsq_l.append(float(val))

    gdf["arcid"] = arcid
    gdf["Shape_Leng"] = shape_leng
    gdf["Row"] = row_l
    gdf["Col"] = col_l
    gdf["slope_deg"] = slope_deg_l
    gdf["meanmsq"] = meanmsq_l

    gdf.to_file(STREAM_OUT)
    print(f"[vector_attrs] wrote {STREAM_OUT.name} with attrs: "
          f"arcid, Shape_Leng, Row, Col, slope_deg, meanmsq")
    # quick sanity echo
    print(f"  slope_deg range: {min(slope_deg_l):.3f}..{max(slope_deg_l):.3f}")
    print(f"  meanmsq   range: {min(meanmsq_l):.1f}..{max(meanmsq_l):.1f}")


if __name__ == "__main__":
    main()
