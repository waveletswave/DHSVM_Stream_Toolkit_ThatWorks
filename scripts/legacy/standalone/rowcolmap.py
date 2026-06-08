# SUMMARY:      rowcolmap.py (standalone / DCC-ready, rasterio + geopandas)
# USAGE:        (A) Raster -> Polygon + ROW/COL
#               (B) Write ROW/COL & SegID to stream line layer
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted by Y. Song
# LAST UPDATE:  2026-05-02
#
# Standalone (non-QGIS) sibling of qgis/roadaspect.py (robust rowcolmap).
# Uses rasterio.features.shapes() + geopandas instead of QgsRasterLayer /
# QgsVectorLayer / processing.run("gdal:polygonize").
#
# IMPORTANT — floating-point parity:
#   We deliberately use an explicit (y_top - y) / dy formula instead of
#   rasterio.transform.rowcol() to guarantee byte-for-byte parity with
#   qgis/roadaspect.py. When a centroid sits exactly on a cell boundary,
#   rasterio's internal affine inverse and the manual division differ by
#   ~1e-11 in opposite directions, producing off-by-one ROW values that
#   would break DHSVM stream.map.dat consistency. Using the explicit
#   formula keeps both pipelines on the same floating-point trajectory.
#
# CLI usage:
#   python rowcolmap.py polygonize <dem.tif> <out_polygon.shp>
#   python rowcolmap.py annotate <stream.shp> <dem.tif> [--seg-field SegID]
#
# Library usage:
#   from rowcolmap import rowcolmapfun, write_col_row_and_segid
#
# Dependencies: rasterio, geopandas, shapely, numpy, pandas

import os
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import shapes as rio_shapes
from shapely.geometry import shape as shp_shape


# ---------------------------
# Internal helpers
# ---------------------------
def _driver_for(path):
    """Pick the OGR driver based on file extension; None lets fiona guess."""
    ext = os.path.splitext(str(path))[1].lower()
    return {
        ".shp":     "ESRI Shapefile",
        ".gpkg":    "GPKG",
        ".geojson": "GeoJSON",
        ".json":    "GeoJSON",
    }.get(ext)


def _save_vector(gdf, path):
    """Save a GeoDataFrame, preserving the source format."""
    drv = _driver_for(path)
    if drv:
        gdf.to_file(path, driver=drv)
    else:
        gdf.to_file(path)


def _crs_match(crs_a, crs_b):
    """Permissive CRS equality (returns True if either is missing)."""
    if crs_a is None or crs_b is None:
        return True
    try:
        return crs_a.to_wkt() == crs_b.to_wkt()
    except Exception:
        return str(crs_a) == str(crs_b)


def _is_line_layer(gdf):
    types = set(g.geom_type for g in gdf.geometry if g is not None)
    return types.issubset({"LineString", "MultiLineString"})


def _xy_to_rowcol(transform, x, y, height, width):
    """
    Convert world (x, y) to integer (row, col) with TOP-LEFT origin
    using the explicit formula from qgis/roadaspect.py (robust rowcolmap).

    This intentionally avoids rasterio.transform.rowcol() — that function
    routes through an Affine inverse whose floating-point error sign can
    differ from the manual division at exact cell-boundary inputs (raw
    row = 36.0 exact may become 36.000000000003 vs 35.999999999985),
    producing off-by-one ROW assignments. The explicit formula matches
    QGIS bit-for-bit.

    Top-left origin: row 0 is the top of the raster, row increases
    downward. Matches DHSVM map-file convention.
    """
    dx = transform.a            # pixel width  (positive)
    dy = abs(transform.e)       # pixel height (transform.e is negative for north-up)
    x_origin = transform.c      # x of the upper-left pixel's left edge
    y_top    = transform.f      # y of the upper-left pixel's top edge (= y_max)

    col = int((x - x_origin) / dx)
    row = int((y_top - y)    / dy)

    # Clamp to valid index range
    col = max(0, min(col, width  - 1))
    row = max(0, min(row, height - 1))
    return row, col


# ===================================================================
# (A) Raster -> Polygon + ROW/COL
# ===================================================================
def rowcolmapfun(elev_path, output_polygon):
    """
    Convert DEM raster into a polygon vector. Each polygon corresponds to
    a connected region of equal-valued DEM cells (mostly per-pixel for
    continuous-valued DEMs). Each polygon is annotated with:

        VALUE   the underlying DEM value (typically elevation in m)
        ROW     row index of the polygon centroid (TOP-LEFT origin)
        COL     column index of the polygon centroid

    Parameters
    ----------
    elev_path : str | Path
        Path to the DEM raster file (.tif).
    output_polygon : str | Path
        Output vector path (e.g. .shp / .gpkg).
    """
    elev_path = str(elev_path)
    output_polygon = str(output_polygon)

    if not os.path.exists(elev_path):
        raise FileNotFoundError(f"Unable to load DEM: {elev_path}")

    out_dir = os.path.dirname(os.path.abspath(output_polygon))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Read DEM and capture metadata (transform / CRS / shape)
    with rasterio.open(elev_path) as src:
        arr = src.read(1)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
        height, width = src.shape

    # Build validity mask (True where data is valid)
    if nodata is not None:
        mask = arr != nodata
    elif np.issubdtype(arr.dtype, np.floating):
        mask = ~np.isnan(arr)
    else:
        mask = np.ones_like(arr, dtype=bool)

    # rasterio.features.shapes wants float32 or int32 for the source
    if arr.dtype not in (np.int32, np.float32):
        arr_for_shapes = arr.astype(np.float32, copy=False)
    else:
        arr_for_shapes = arr

    # 4-connectivity to mirror QGIS call (EIGHT_CONNECTEDNESS=False)
    geoms = []
    values = []
    for geom_dict, val in rio_shapes(arr_for_shapes,
                                     mask=mask,
                                     transform=transform,
                                     connectivity=4):
        geoms.append(shp_shape(geom_dict))
        values.append(val)

    # Compute ROW/COL from each polygon's centroid (top-left origin)
    rows = np.empty(len(geoms), dtype=np.int32)
    cols = np.empty(len(geoms), dtype=np.int32)
    for i, geom in enumerate(geoms):
        cen = geom.centroid
        r, c = _xy_to_rowcol(transform, cen.x, cen.y, height, width)
        rows[i] = r
        cols[i] = c

    gdf = gpd.GeoDataFrame(
        {
            "VALUE": values,
            "ROW":   rows,
            "COL":   cols,
            "geometry": geoms,
        },
        crs=crs,
    )

    _save_vector(gdf, output_polygon)
    print(f"[rowcolmap] Polygonized DEM -> {output_polygon}")
    print(f"[rowcolmap] {len(gdf)} polygons written; ROW/COL fields populated.")


# ===================================================================
# (B) Write ROW/COL & SegID to stream lines
# ===================================================================
def write_col_row_and_segid(vector_path, dem_path, seg_field="SegID"):
    """
    For a stream LINE layer:
      - Write integer Col / Row fields aligned to the DEM grid
        (top-left origin; row increases downward; matches DHSVM map files).
      - Write a continuous integer SegID (1..N). If SegID already exists,
        keep valid integer values and fill gaps with the next-available ID;
        otherwise assign 1..N in feature order.

    Parameters
    ----------
    vector_path : str | Path
        Path to a LINE layer of the stream network (.shp / .gpkg).
    dem_path : str | Path
        Path to the DEM raster used by DHSVM.
    seg_field : str, default "SegID"
        Name of the per-segment ID field to write/ensure.
    """
    vector_path = str(vector_path)
    dem_path = str(dem_path)

    if not os.path.exists(vector_path):
        raise FileNotFoundError(f"Cannot load vector layer: {vector_path}")
    if not os.path.exists(dem_path):
        raise FileNotFoundError(f"DEM not valid: {dem_path}")

    # Load DEM metadata only (no need to read the array)
    with rasterio.open(dem_path) as src:
        transform = src.transform
        dem_crs = src.crs
        height, width = src.shape

    # Load stream layer
    layer = gpd.read_file(vector_path)
    if len(layer) == 0:
        raise RuntimeError(f"Stream layer is empty: {vector_path}")
    if not _is_line_layer(layer):
        raise RuntimeError("Input vector is not a LINE geometry layer.")

    # CRS sanity check (warn but proceed)
    if not _crs_match(layer.crs, dem_crs):
        print("[rowcolmap] WARNING: CRS mismatch between stream layer and DEM. "
              "Consider reprojecting to the DEM CRS to avoid index drift.")

    # ---------------------------
    # Two-pass SegID assignment (mirrors QGIS logic)
    #   Pass 1: collect existing valid IDs and identify rows needing IDs
    #   Pass 2: fill gaps with the next available ID, starting from 1
    # ---------------------------
    n = len(layer)
    if seg_field not in layer.columns:
        # Field absent entirely -> assign 1..N in row order
        new_ids = np.arange(1, n + 1, dtype=np.int64)
    else:
        existing_ids = set()
        need_id_pos_set = set()
        seg_values = layer[seg_field].tolist()

        for pos, val in enumerate(seg_values):
            if pd.isna(val):
                need_id_pos_set.add(pos)
                continue
            try:
                existing_ids.add(int(val))
            except (TypeError, ValueError):
                need_id_pos_set.add(pos)

        new_ids = np.empty(n, dtype=np.int64)
        next_id = 1
        for pos in range(n):
            if pos in need_id_pos_set:
                while next_id in existing_ids:
                    next_id += 1
                new_ids[pos] = next_id
                existing_ids.add(next_id)
                next_id += 1
            else:
                try:
                    new_ids[pos] = int(seg_values[pos])
                except (TypeError, ValueError):
                    # Defensive fallback; should not be reached
                    new_ids[pos] = pos + 1

    # ---------------------------
    # Compute Col/Row from each line's centroid
    # ---------------------------
    rows = np.empty(n, dtype=np.int32)
    cols = np.empty(n, dtype=np.int32)
    for i, geom in enumerate(layer.geometry):
        cen = geom.centroid
        r, c = _xy_to_rowcol(transform, cen.x, cen.y, height, width)
        rows[i] = r
        cols[i] = c

    # ---------------------------
    # Write fields back
    # ---------------------------
    layer["Col"]      = cols
    layer["Row"]      = rows
    layer[seg_field]  = new_ids.astype(np.int32)

    _save_vector(layer, vector_path)

    # ---------------------------
    # Console summary + uniqueness check
    # ---------------------------
    print(f"[rowcolmap] Col/Row & {seg_field} written to "
          f"{os.path.basename(vector_path)}.")

    if layer[seg_field].duplicated().any():
        print("[rowcolmap] WARNING: Duplicate SegID detected. "
              "Consider reassigning SegID.")
    else:
        print(f"[rowcolmap] SegID uniqueness check: OK (N={n}).")


# ===================================================================
# CLI
# ===================================================================
def _build_parser():
    p = argparse.ArgumentParser(
        description="rowcolmap utilities for DHSVM stream-map preparation."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser(
        "polygonize",
        help="Convert a DEM raster to per-region polygons with ROW/COL/VALUE.",
    )
    pp.add_argument("dem", help="Path to DEM raster (.tif)")
    pp.add_argument("output", help="Output polygon vector (.shp / .gpkg)")

    ap = sub.add_parser(
        "annotate",
        help="Write Col/Row/SegID to an existing stream line layer.",
    )
    ap.add_argument("stream", help="Stream line layer (.shp / .gpkg)")
    ap.add_argument("dem",    help="DEM raster (.tif) used as the index reference")
    ap.add_argument("--seg-field", default="SegID",
                    help="Name of the segment ID field (default: SegID)")

    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()

    if args.cmd == "polygonize":
        rowcolmapfun(args.dem, args.output)
    elif args.cmd == "annotate":
        write_col_row_and_segid(args.stream, args.dem, seg_field=args.seg_field)