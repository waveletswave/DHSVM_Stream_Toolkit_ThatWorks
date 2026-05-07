# SUMMARY:      roadaspect.py (QGIS version)
# USAGE:        Computes aspect and length of stream arc within each DEM cell.
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted for QGIS by Y.Song
# LAST UPDATE:  2025-09-09 (field checking added)

# -*- coding: utf-8 -*-
# SUMMARY:      rowcolmap.py (QGIS version, robust)
# USAGE:        (A) Raster→Polygon + ROW/COL; (B) Write ROW/COL & SegID to stream lines
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted for QGIS by Y. Song (refined)
# LAST CHANGE:  2025-10-30

from qgis.core import (
    QgsRasterLayer, QgsVectorLayer, QgsField, QgsWkbTypes, QgsCoordinateTransformContext
)
from qgis.PyQt.QtCore import QVariant
import processing
import os
import math

# ---------------------------
# helpers
# ---------------------------
def _raster_metrics(dem: QgsRasterLayer):
    """
    Return (x_min, x_max, y_min, y_max, width, height, dx, dy)
    """
    ext = dem.extent()
    W, H = dem.width(), dem.height()
    if W <= 0 or H <= 0:
        raise RuntimeError("DEM has invalid width/height.")
    x0, x1 = ext.xMinimum(), ext.xMaximum()
    y0, y1 = ext.yMinimum(), ext.yMaximum()
    dx, dy = (x1 - x0) / float(W), (y1 - y0) / float(H)
    return x0, x1, y0, y1, W, H, dx, dy


def _ensure_fields(layer: QgsVectorLayer, fields):
    """
    Ensure given fields exist on the layer. `fields` is list of (name, QVariant.Type).
    """
    layer.startEditing()
    for name, qtype in fields:
        if layer.fields().indexOf(name) < 0:
            layer.addAttribute(QgsField(name, qtype))
    layer.commitChanges()


def _is_linestring(layer: QgsVectorLayer) -> bool:
    return layer.geometryType() == QgsWkbTypes.LineGeometry


def _same_crs(layer_a, layer_b) -> bool:
    try:
        return layer_a.crs() == layer_b.crs()
    except Exception:
        return True  # be permissive


# ---------------------------
# (A) Raster→Polygon + ROW/COL
# ---------------------------
def rowcolmapfun(elev_path: str, output_polygon: str):
    """
    Converts DEM raster into a polygon representation where each polygon represents a DEM cell.
    Each polygon is assigned row and column numbers (origin at TOP-LEFT; row increases downward).

    Parameters
    ----------
    elev_path : str
        Path to the DEM raster file (.tif)
    output_polygon : str
        Output vector path (e.g., .shp / .gpkg) for polygon representation
    """

    dem = QgsRasterLayer(elev_path, "dem_raster")
    if not dem.isValid():
        raise RuntimeError(f"Unable to load DEM: {elev_path}")

    out_dir = os.path.dirname(os.path.abspath(output_polygon))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Raster → polygon
    processing.run("gdal:polygonize", {
        "INPUT": elev_path,
        "BAND": 1,
        "FIELD": "VALUE",
        "EIGHT_CONNECTEDNESS": False,
        "EXTRA": "",
        "OUTPUT": output_polygon
    })
    print(f"[rowcolmap] Polygonized DEM → {output_polygon}")

    poly = QgsVectorLayer(output_polygon, "rowcolpoly", "ogr")
    if not poly.isValid():
        raise RuntimeError(f"Unable to load polygon layer: {output_polygon}")

    # Prepare fields
    _ensure_fields(poly, [("ROW", QVariant.Int), ("COL", QVariant.Int)])

    x0, x1, y0, y1, W, H, dx, dy = _raster_metrics(dem)

    # Assign ROW/COL (TOP-LEFT origin)
    poly.startEditing()
    for ft in poly.getFeatures():
        cen = ft.geometry().centroid().asPoint()
        col = int((cen.x() - x0) / dx)
        row = int((y1 - cen.y()) / dy)  # top-left origin: row increases downward
        # clamp
        col = max(0, min(col, W - 1))
        row = max(0, min(row, H - 1))
        ft["ROW"] = row
        ft["COL"] = col
        poly.updateFeature(ft)
    poly.commitChanges()
    print("[rowcolmap] ROW/COL written on polygonized grid.")


# ---------------------------
# (B) Write ROW/COL & SegID to stream lines
# ---------------------------
def write_col_row_and_segid(vector_path: str, dem_path: str, seg_field: str = "SegID"):
    """
    For a stream LINE layer:
      - Write integer ROW/COL fields aligned to the DEM grid
      - Write a continuous integer SegID (1..N) for stable IDs used by map.dat/network.dat
        (If SegID exists, keep/exhaustive-fill; otherwise assign new 1..N)
    Notes
    -----
    - ROW/COL origin at TOP-LEFT; row increases downward (consistent with DHSVM map files)
    - DEM and vector should be in the same CRS; a mismatch will warn but still proceed

    Parameters
    ----------
    vector_path : str
        Path to a LINE layer of the stream network (e.g., .shp/.gpkg)
    dem_path : str
        Path to the DEM raster used by DHSVM
    seg_field : str, default "SegID"
        Name of the per-segment ID field to write/ensure
    """

    dem = QgsRasterLayer(dem_path, "dem")
    if not dem.isValid():
        raise RuntimeError(f"DEM not valid: {dem_path}")

    layer = QgsVectorLayer(vector_path, "streams", "ogr")
    if not layer.isValid():
        raise RuntimeError(f"Cannot load vector layer: {vector_path}")
    if not _is_linestring(layer):
        raise RuntimeError("Input vector is not a LINE geometry layer.")

    if not _same_crs(layer, dem):
        print("[rowcolmap] WARNING: CRS mismatch between stream layer and DEM. "
              "Consider reprojecting to the DEM CRS to avoid index drift.")

    _ensure_fields(layer, [("Col", QVariant.Int), ("Row", QVariant.Int), (seg_field, QVariant.Int)])

    x0, x1, y0, y1, W, H, dx, dy = _raster_metrics(dem)

    # First pass: determine which features need new SegID
    need_id = []
    have_ids = set()
    for ft in layer.getFeatures():
        val = ft[seg_field]
        if val is None:
            need_id.append(ft.id())
        else:
            try:
                have_ids.add(int(val))
            except Exception:
                need_id.append(ft.id())

    # Build a continuous pool of available IDs starting from 1
    next_id = 1
    assigned = {}
    for fid in need_id:
        while next_id in have_ids:
            next_id += 1
        assigned[fid] = next_id
        have_ids.add(next_id)
        next_id += 1

    # Write back Col/Row and SegID
    layer.startEditing()
    for ft in layer.getFeatures():
        cen = ft.geometry().centroid().asPoint()
        col = int((cen.x() - x0) / dx)
        row = int((y1 - cen.y()) / dy)  # TOP-LEFT origin
        col = max(0, min(col, W - 1))
        row = max(0, min(row, H - 1))

        ft["Col"] = col
        ft["Row"] = row
        if ft.id() in assigned:
            ft[seg_field] = int(assigned[ft.id()])
        else:
            # sanitize existing value
            try:
                ft[seg_field] = int(ft[seg_field])
            except Exception:
                # in case of odd types
                ft[seg_field] = int(ft.id())  # fallback but still deterministic

        layer.updateFeature(ft)
    layer.commitChanges()

    print(f"[rowcolmap] Col/Row & {seg_field} written to {os.path.basename(vector_path)}.")
    # Quick sanity: SegID should be unique per feature
    ids = set()
    dup = False
    for ft in layer.getFeatures():
        sid = int(ft[seg_field])
        if sid in ids:
            dup = True
            break
        ids.add(sid)
    if dup:
        print("[rowcolmap] WARNING: Duplicate SegID detected. Consider reassigning SegID.")
    else:
        print(f"[rowcolmap] SegID uniqueness check: OK (N={len(ids)}).")

