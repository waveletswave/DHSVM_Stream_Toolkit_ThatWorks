# -*- coding: utf-8 -*-
# =====================================================================
# createstreamnetwork.py — QGIS + GRASS (portable, robust, DHSVM-ready)
# =====================================================================
# PURPOSE
#   Build a stream network and DHSVM inputs (stream.map.dat,
#   stream.network.dat, stream.class.dat) in QGIS using GRASS tools,
#   with clean fallbacks, portable paths, and stable raster sampling.
#
# RUN INSIDE: QGIS Python Console (Processing + GRASS provider enabled)
#
# INPUT (auto-discovered)
#   Workspace (WS)    = parent folder of this script (…/DEM_<BASIN>/)
#   DEM (projected)   = WS/Reprojected_DEM/* or WS/*
#   Watershed polygon = WS/Reprojected_Watersheds/*.shp or WS/*.shp
#
# OUTPUTS (all in WS)
#   stream_raster.tif
#   streamfile.shp
#   stream_slope.tif   (degrees)
#   stream.map.dat
#   stream.network.dat (6 cols: SegID Nin/Str Slope(tan) Length_m ClassID DownSegID)
#   stream.class.dat
#
# KEY CONVENTIONS
#   - Threshold is specified in *cells*; code prints the equivalent area (m²).
#   - Slope raster is in degrees; stream.network.dat stores tan(slope).
#   - ClassID must match stream.class.dat; read from 'chanclass' if present.
#
# Author: Y. Song (Duke / PNNL pipeline) — 2026-02-27
# Notes : Integrated raster-lifetime hotfix + FA-directed topology + single-outlet consolidation
# =====================================================================

import os, time, math
from math import tan, radians, hypot, cos, sin
from pathlib import Path
import importlib.util
from collections import defaultdict, deque

from qgis.core import (
    QgsApplication, QgsRasterLayer, QgsVectorLayer, QgsField, QgsWkbTypes,
    QgsPointXY, QgsFeature, QgsGeometry, QgsSpatialIndex
)
from qgis.PyQt.QtCore import QVariant
import processing

# ----------------------------- #
#   User-tunable configuration  #
# ----------------------------- #
MIN_SRC_CELLS   = 60          # r.stream.extract threshold in number of cells
NIN_MODE        = "indegree"   # "indegree" or "shreve" for Col(2)
BASE_SLOPE_SAMPLES = 12        # samples per line when computing tan(slope)
SNAP_TOL_LIST   = [0.75, 1.25, 1.75]  # unused in current FA approach; kept for reference
FA_TOL_MULT     = 1.10         # FA matching tolerance = 1.10 × pixel-diagonal
AHEAD_STEPS     = (1.0, 2.0)   # forward probing (in pixel units) to bridge tiny gaps

# ----------------------------- #
#  Portable path configuration  #
# ----------------------------- #
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

def p(*xs) -> str:
    return str(WS.joinpath(*xs))

# ----------------------------- #
#     Input auto-discovery      #
# ----------------------------- #

def _pick_first_existing(cands):
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None

# --- DEM auto-discovery (more robust across basins) ---
DEM_DIR = WS / "Reprojected_DEM"

# 1) Common names we have used before
primary_dem_candidates = [
    p("Reprojected_DEM", "Cropped_DEM.tif"),
    p("Reprojected_DEM", "dem.tif"),
    p("Reprojected_DEM", "elev.tif"),
    p("Reprojected_DEM", "elev_clipped.tif"),
]

# 2) Fallback: any other GeoTIFF in Reprojected_DEM (e.g. USGS_1_..._UTM17.tif)
generic_tifs = []
if DEM_DIR.exists():
    generic_tifs = [str(x) for x in sorted(DEM_DIR.glob("*.tif"))]

# 3) Last resort: dem.tif / elev.tif directly under WS
DEM_CANDIDATES = primary_dem_candidates + generic_tifs + [
    p("dem.tif"),
    p("elev.tif"),
]

def _find_watershed_shp():
    wd = WS/"Reprojected_Watersheds"
    if wd.exists():
        shp = sorted([str(x) for x in wd.glob("*.shp")],
                     key=lambda s: ("watershed" not in s.lower(), s.lower()))
        if shp:
            return shp[0]
    shp2 = sorted([str(x) for x in WS.glob("*.shp")],
                  key=lambda s: ("watershed" not in s.lower(), s.lower()))
    return shp2[0] if shp2 else None

elev_raw = _pick_first_existing(DEM_CANDIDATES)
wshed    = _find_watershed_shp()
if not elev_raw:
    raise FileNotFoundError("Projected DEM not found. Put it in Reprojected_DEM/.")

if wshed: print(f"[ok] Watershed polygon: {wshed}")
print(f"[ok] DEM: {elev_raw}")

# Derived paths
elev           = p("Reprojected_DEM","elev_clipped.tif") if wshed else elev_raw
flow_dir       = p("flow_dir.tif")
flow_acc       = p("flow_acc.tif")
stream_raster  = p("stream_raster.tif")
stream_vec     = p("stream_order.shp")   # r.stream.extract vector (if available)
streamfile_fbk = p("streamfile.shp")     # strict fallback from raster to vector
stream_slope   = p("stream_slope.tif")

# ----------------------------- #
#          Raster hotfix        #
# ----------------------------- #
# Keep QgsRasterLayer objects alive for the whole run (never cache providers alone).
_RASTER_LAYER_CACHE = {}

def _raster_layer_for(raster_path: str) -> QgsRasterLayer:
    """Return a live QgsRasterLayer and keep it cached to avoid GC issues."""
    key = os.path.abspath(raster_path)
    rl = _RASTER_LAYER_CACHE.get(key)
    if (rl is None) or (not isinstance(rl, QgsRasterLayer)) or (not rl.isValid()):
        rl = QgsRasterLayer(raster_path, f"ras::{os.path.basename(raster_path)}")
        if not rl.isValid():
            raise RuntimeError(f"Cannot load raster: {raster_path}")
        _RASTER_LAYER_CACHE[key] = rl
    return rl

def _provider_for(raster_path: str):
    """Always get the provider from a live layer (prevents 'provider deleted' crashes)."""
    return _raster_layer_for(raster_path).dataProvider()

def _sample_raster_val(ras_path: str, x: float, y: float, default=0.0):
    """Sample band-1 value at (x, y); robust against provider GC."""
    prov = _provider_for(ras_path)
    val, ok = prov.sample(QgsPointXY(x, y), 1)
    if ok and (val is not None):
        try:
            fv = float(val)
            if math.isfinite(fv):
                return fv
        except Exception:
            pass
    return default

def _cell_area_m2(raster_path: str) -> float:
    rl = _raster_layer_for(raster_path)
    return abs(rl.rasterUnitsPerPixelX() * rl.rasterUnitsPerPixelY())

def _pixel_size_and_diag(raster_path: str):
    rl = _raster_layer_for(raster_path)
    px = abs(rl.rasterUnitsPerPixelX()); py = abs(rl.rasterUnitsPerPixelY())
    diag = (px**2 + py**2) ** 0.5
    return px, py, diag

# ----------------------------- #
#         Helpers (vector)      #
# ----------------------------- #

def _geom_type(path: str) -> str:
    v = QgsVectorLayer(path,"chk","ogr")
    return QgsWkbTypes.displayString(v.wkbType()) if v.isValid() else "Unknown"

def _is_line_layer(path: str) -> bool:
    g = _geom_type(path)
    return ("LineString" in g)

def _feature_count(path: str) -> int:
    v = QgsVectorLayer(path, "chk", "ogr")
    return sum(1 for _ in v.getFeatures()) if v.isValid() else 0

def _is_null(v):
    return (v is None) or (isinstance(v, QVariant) and v.isNull())

def _to_float(v, default=None):
    try:
        return float(v) if (v is not None) else default
    except Exception:
        return default

# ----------------------------- #
#      Threshold & clipping     #
# ----------------------------- #
cell_area    = _cell_area_m2(elev_raw)
MIN_SRC_AREA = MIN_SRC_CELLS * cell_area
print(f"[info] Cell area ≈ {cell_area:.3f} m²")
print(f"[info] MIN_SRC_AREA (m²) = {MIN_SRC_AREA:.3f}  (from {MIN_SRC_CELLS} cells)")

if wshed:
    print("[step] Clipping DEM to watershed…")
    processing.run("gdal:cliprasterbymasklayer", {
        "INPUT": elev_raw, "MASK": wshed, "OUTPUT": elev, "NODATA": -9999
    })
else:
    elev = elev_raw
print(f"[ok] DEM ready: {elev}")

# ----------------------------- #
#  Flow dir / accumulation      #
# ----------------------------- #
print("[step] GRASS r.watershed → flow direction & accumulation…")
processing.run("grass7:r.watershed", {
    'elevation': elev, 'accumulation': flow_acc, 'drainage': flow_dir,
    'convergence': 5, 'memory': 300,
    'GRASS_REGION_PARAMETER': elev, 'GRASS_REGION_CELLSIZE_PARAMETER': 0,
    'GRASS_OUTPUT_TYPE_PARAMETER': 0
})
print("[ok] r.watershed completed.")

# ----------------------------- #
#  r.stream.extract raster+vec  #
# ----------------------------- #
print("[step] GRASS r.stream.extract → raster + vector (with topology)…")
processing.run("grass7:r.stream.extract", {
    'elevation': elev,
    'accumulation': flow_acc,
    'direction': flow_dir,          # request vector with topology when possible
    'threshold': MIN_SRC_CELLS,
    'stream_raster': stream_raster,
    'stream_vector': stream_vec,
    'memory': 300, '-m': True,
    'GRASS_REGION_PARAMETER': elev, 'GRASS_REGION_CELLSIZE_PARAMETER': 0,
    'GRASS_OUTPUT_TYPE_PARAMETER': 0
})
for _ in range(12):
    if os.path.exists(stream_raster) and os.path.getsize(stream_raster)>0: break
    time.sleep(0.25)
print(f"[ok] Stream raster: {stream_raster}")

vec_count = _feature_count(stream_vec)
print(f"[info] stream_vec feature count = {vec_count}")

def _force_lines_from_raster(raster_path: str, out_path: str) -> str:
    """Try several r.to.vect variants. Accept only proper line geometry."""
    trials = [
        {'type': 0},                  # often MultiLineString
        {'type': 'line'},
        {'type': 1},
        {'type': 'line', '-s': True}, # smoothed
    ]
    for opt in trials:
        try:
            processing.run("grass7:r.to.vect", {
                'input': raster_path,
                'type': opt.get('type', 1),
                'output': out_path,
                **({'-s': True} if opt.get('-s') else {}),
                'GRASS_REGION_PARAMETER': raster_path,
                'GRASS_REGION_CELLSIZE_PARAMETER': 0
            })
            g = _geom_type(out_path); n = _feature_count(out_path)
            print(f"[debug] r.to.vect{opt} → {g}, features={n}")
            if n > 0 and _is_line_layer(out_path):
                return out_path
        except Exception as e:
            print(f"[warn] r.to.vect{opt} failed: {e}")
    raise RuntimeError("Failed to obtain LineString geometry from stream raster.")

vector_path = stream_vec
if vec_count == 0 or not _is_line_layer(vector_path):
    print("[warn] r.stream.extract vector empty or not lines → fallback to r.to.vect")
    vector_path = _force_lines_from_raster(stream_raster, streamfile_fbk)

print(f"[ok] Streams ready → {vector_path} ({_geom_type(vector_path)}, features={_feature_count(vector_path)})")

# Load once and keep this layer alive for the rest of the pipeline
vl_lines = QgsVectorLayer(vector_path,"streams_live","ogr")
if not vl_lines.isValid():
    raise RuntimeError("Failed to load streams vector.")

# ----------------------------- #
#   Slope raster (DEGREES)      #
# ----------------------------- #
print("[step] GRASS r.slope.aspect → slope raster…")
processing.run("grass7:r.slope.aspect",{
    'elevation':elev,'slope':stream_slope,'format':1, # 1 = degrees
    'GRASS_REGION_PARAMETER':elev,'GRASS_REGION_CELLSIZE_PARAMETER':0
})
print(f"[ok] Slope raster: {stream_slope}")

# ----------------------------- #
#  Normalize fields for map.dat #
# ----------------------------- #
def _ensure_fields_rowcol_len(layer: QgsVectorLayer):
    if sum(1 for _ in layer.getFeatures()) == 0:
        raise RuntimeError("Stream line layer has 0 features after fallback.")
    if QgsWkbTypes.geometryType(layer.wkbType()) != QgsWkbTypes.LineGeometry:
        # allow MultiLineString; check name
        if "LineString" not in QgsWkbTypes.displayString(layer.wkbType()):
            raise RuntimeError("Stream vector is not LineString/MultiLineString.")
    layer.startEditing()
    need = [("arcid",QVariant.Int),("Shape_Leng",QVariant.Double),("Row",QVariant.Int),("Col",QVariant.Int)]
    names = layer.fields().names()
    for n,t in need:
        if n not in names: layer.addAttribute(QgsField(n,t))
    layer.updateFields()

    dem = _raster_layer_for(elev)
    px,py = abs(dem.rasterUnitsPerPixelX()), abs(dem.rasterUnitsPerPixelY())
    xmin,ymin = dem.extent().xMinimum(), dem.extent().yMinimum()

    idx_arc = layer.fields().indexFromName("arcid")
    idx_len = layer.fields().indexFromName("Shape_Leng")
    idx_row = layer.fields().indexFromName("Row")
    idx_col = layer.fields().indexFromName("Col")

    i=0
    for ft in layer.getFeatures():
        i+=1
        g=ft.geometry(); L=g.length() if g and not g.isEmpty() else 0.0
        c=g.centroid().asPoint() if g else QgsPointXY(xmin,ymin)
        row=int((c.y()-ymin)/py) if py>0 else 0
        col=int((c.x()-xmin)/px) if px>0 else 0
        ft[idx_arc]=i; ft[idx_len]=float(L); ft[idx_row]=row; ft[idx_col]=col
        layer.updateFeature(ft)
    layer.commitChanges()
    print(f"[ok] Attributes normalized. arcid=1..{i}; Shape_Leng>0 expected.")

_ensure_fields_rowcol_len(vl_lines)

# ----------------------------- #
#  Mean slope fields on lines   #
# ----------------------------- #
def _sample_mean_slope_deg(geom, slope_raster_path: str, n_samples=BASE_SLOPE_SAMPLES):
    if (geom is None) or geom.isEmpty() or QgsWkbTypes.geometryType(geom.wkbType())!=QgsWkbTypes.LineGeometry:
        return 0.0
    L = geom.length()
    if L<=0: return 0.0
    prov = _provider_for(slope_raster_path)
    vals=[]
    for i in range(1, n_samples+1):
        d=(i/(n_samples+1.0))*L
        p=geom.interpolate(d).asPoint()
        ok,val = prov.sample(p,1)
        if ok and val is not None:
            try: vals.append(max(0.0, float(val)))
            except: pass
    return sum(vals)/len(vals) if vals else 0.0

def _write_slope_deg_field(layer: QgsVectorLayer):
    layer.startEditing()
    if "slope_deg" not in layer.fields().names():
        layer.addAttribute(QgsField("slope_deg", QVariant.Double))
        layer.updateFields()
    idx = layer.fields().indexFromName("slope_deg")
    for ft in layer.getFeatures():
        ft[idx] = _sample_mean_slope_deg(ft.geometry(), stream_slope, BASE_SLOPE_SAMPLES)
        layer.updateFeature(ft)
    layer.commitChanges()
    print("[ok] slope_deg written (mean along segment).")

_write_slope_deg_field(vl_lines)

# ----------------------------- #
#  Accumulation mean (meanmsq)  #
# ----------------------------- #
def _find_mean_field(layer, prefix):
    names = {f.name().lower(): f.name() for f in layer.fields()}
    for key in ("mean","average","avg"):
        cands = [names[n] for n in names if n.startswith(prefix.lower()) and key in n]
        if cands: return cands[0]
    cands = [names[n] for n in names if n.startswith(prefix.lower())]
    return cands[0] if cands else None

def _qgis_zonal(layer_path, raster_path, prefix):
    algs = {a.id() for a in QgsApplication.processingRegistry().algorithms()}
    created = None
    if "qgis:zonalstatisticsfb" in algs:
        try:
            processing.run("qgis:zonalstatisticsfb", {
                "INPUT": layer_path, "INPUT_RASTER": raster_path,
                "RASTER_BAND": 1, "COLUMN_PREFIX": prefix,
                "STATISTICS": [3],  # mean
                "OUTPUT": layer_path
            }, is_child_algorithm=True)
            vl = QgsVectorLayer(layer_path, "tmp", "ogr")
            created = _find_mean_field(vl, prefix)
            if created: return created
        except Exception:
            pass
    if "qgis:zonalstatistics" in algs:
        try:
            processing.run("qgis:zonalstatistics", {
                "INPUT_RASTER": raster_path, "RASTER_BAND": 1,
                "INPUT_VECTOR": layer_path,
                "COLUMN_PREFIX": prefix, "STATISTICS": [2]  # mean
            }, is_child_algorithm=True)
            vl = QgsVectorLayer(layer_path, "tmp2", "ogr")
            created = _find_mean_field(vl, prefix)
            return created
        except Exception:
            pass
    return None

def _grass_stats_join(layer_path, raster_path, prefix, ws_dir):
    tmp_gpkg = os.path.join(ws_dir, f"{prefix}tmp.gpkg")
    processing.run("grass7:v.rast.stats", {
        "map": layer_path, "raster": raster_path,
        "column_prefix": prefix, "method": 4,  # average
        "output": tmp_gpkg, "overwrite": True,
        "GRASS_REGION_PARAMETER": raster_path,
        "GRASS_REGION_CELLSIZE_PARAMETER": 0
    }, is_child_algorithm=True)
    out_join = os.path.join(ws_dir, f"{prefix}joined.gpkg")
    processing.run("native:joinattributesbylocation", {
        "INPUT": layer_path, "JOIN": tmp_gpkg,
        "PREDICATE": [0], "JOIN_FIELDS": [],
        "METHOD": 0, "DISCARD_NONMATCHING": False,
        "OUTPUT": out_join
    }, is_child_algorithm=True)
    jl = QgsVectorLayer(out_join, "joined", "ogr")
    mean_name = _find_mean_field(jl, prefix)
    vl = QgsVectorLayer(layer_path, "orig", "ogr")
    vl.startEditing()
    if vl.fields().indexOf("mean_tmp") < 0:
        vl.addAttribute(QgsField("mean_tmp", QVariant.Double))
        vl.updateFields()
    vals = []
    for f in jl.getFeatures():
        vals.append(_to_float(f[mean_name], default=None))
    i = 0
    for ft in vl.getFeatures():
        val = vals[i] if i < len(vals) else None
        ft["mean_tmp"] = val
        vl.updateFeature(ft); i += 1
    vl.commitChanges()
    return "mean_tmp"

def _sample_line_raster_path(raster_path: str, geom: QgsGeometry, n=15):
    """Sample a raster along a line at n evenly spaced points (robust provider)."""
    if (geom is None) or geom.isEmpty(): return []
    L = geom.length() or 0.0
    if L <= 0: return []
    prov = _provider_for(raster_path)
    vals=[]
    for k in range(n):
        d = (k/(n-1.0))*L
        p = geom.interpolate(d).asPoint()
        v, ok = prov.sample(p, 1)
        if ok and (v is not None):
            try: vals.append(float(v))
            except: pass
    return vals

def ensure_meanmsq(stream_layer: QgsVectorLayer, acc_raster_path: str, ws_dir: str,
                   out_field="meanmsq", pixel_area_m2=1.0):
    created = _qgis_zonal(stream_layer.source(), acc_raster_path, prefix="acc_")
    if not created:
        try:
            created = _grass_stats_join(stream_layer.source(), acc_raster_path, "acc_", ws_dir)
        except Exception:
            created = None
    vl = stream_layer
    vl.startEditing()
    if vl.fields().indexOf(out_field) < 0:
        vl.addAttribute(QgsField(out_field, QVariant.Double)); vl.updateFields()
    for ft in vl.getFeatures():
        val = _to_float(ft[created], default=None) if created else None
        if (val is None) or (not math.isfinite(val)) or (val <= 0):
            vals = [v for v in _sample_line_raster_path(acc_raster_path, ft.geometry(), n=15)
                    if (v is not None) and math.isfinite(v) and (v > 0)]
            if vals:
                vals.sort(); val = vals[len(vals)//2]
            else:
                val = max(1.0, float(pixel_area_m2))
        ft[out_field] = float(val)
        vl.updateFeature(ft)
    vl.commitChanges()
    print("[ok] meanmsq written (per-segment contributing area, m²).")

# Compute meanmsq BEFORE channel classes
ensure_meanmsq(vl_lines, flow_acc, str(WS), out_field="meanmsq", pixel_area_m2=cell_area)

# ----------------------------- #
#   rd_aspect (optional)        #
# ----------------------------- #
def _calc_rd_aspect_fallback(layer: QgsVectorLayer):
    if "rd_aspect" in layer.fields().names(): return
    layer.startEditing(); layer.addAttribute(QgsField("rd_aspect",QVariant.Int))
    idx = layer.fields().indexFromName("rd_aspect")
    for ft in layer.getFeatures():
        g=ft.geometry()
        if (g is None) or g.isEmpty() or (QgsWkbTypes.geometryType(g.wkbType())!=QgsWkbTypes.LineGeometry):
            ft[idx]=0; layer.updateFeature(ft); continue
        line = g.asMultiPolyline()[0] if g.isMultipart() else g.asPolyline()
        if len(line)>=2:
            dx=line[-1].x()-line[0].x(); dy=line[-1].y()-line[0].y()
            ang=(math.degrees(math.atan2(dy,dx))+360)%360; ft[idx]=int(round(ang))
        else:
            ft[idx]=0
        layer.updateFeature(ft)
    layer.commitChanges(); print("[ok] rd_aspect computed by fallback.")

print("[step] Compute rd_aspect via roadaspect.py…")
try:
    spec = importlib.util.spec_from_file_location("roadaspect", str(SCRIPT_DIR/"roadaspect.py"))
    mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    if hasattr(mod,"roadaspectfun"):
        mod.roadaspectfun(vl_lines.source())
    else:
        print("[warn] roadaspectfun() not found; fallback."); _calc_rd_aspect_fallback(vl_lines)
except Exception as e:
    print(f"[warn] roadaspect step failed; fallback: {e}"); _calc_rd_aspect_fallback(vl_lines)

# ----------------------------- #
#   Channel classes (external)  #
# ----------------------------- #
print("[step] Classify channels via channelclass.py…")
try:
    spec = importlib.util.spec_from_file_location("channelclass", str(SCRIPT_DIR/"channelclass.py"))
    mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.channelclassfun(vl_lines.source(), str(WS))
    print("[ok] stream.class.dat written.")
except Exception as e:
    print(f"[warn] Channel classification skipped (continuing): {e}")

# ----------------------------- #
#  Width/Height for map.dat     #
# ----------------------------- #
# Mirror to the SAME in-memory layer used for topology (stable FIDs).
vl_lines.startEditing()
if vl_lines.fields().indexOf("Width") < 0:
    vl_lines.addAttribute(QgsField("Width", QVariant.Double))
if vl_lines.fields().indexOf("Height") < 0:
    vl_lines.addAttribute(QgsField("Height", QVariant.Double))
vl_lines.updateFields()

names = {f.name().lower(): f.name() for f in vl_lines.fields()}
w_name = next((names[nm] for nm in ("effwidth","hydwidth","w","width") if nm in names), None)
h_name = next((names[nm] for nm in ("effdepth","hyddepth","d","depth","bankht","bankheight") if nm in names), None)

for ft in vl_lines.getFeatures():
    try:
        ft["Width"] = float(ft[w_name]) if w_name else 0.50
    except Exception:
        ft["Width"] = 0.50
    try:
        ft["Height"] = float(ft[h_name]) if h_name else 0.10
    except Exception:
        ft["Height"] = 0.10
    vl_lines.updateFeature(ft)
vl_lines.commitChanges()

# ----------------------------- #
#  Field helpers                #
# ----------------------------- #
def _field_name_ci(vl: QgsVectorLayer, name: str):
    return {n.lower(): n for n in vl.fields().names()}.get(name.lower())

def _first_existing_field(vl: QgsVectorLayer, candidates):
    for cand in candidates:
        nm = _field_name_ci(vl, cand)
        if nm is not None: return nm
    return None

# ----------------------------- #
#  FA-directed network builder  #
# ----------------------------- #
def _line_coords(geom: QgsGeometry):
    if geom.isMultipart():
        parts = geom.asMultiPolyline()
        if not parts: return []
        return max(parts, key=lambda L: len(L))  # pick the densest part
    else:
        return geom.asPolyline()

def _azimuth_0_360_from(p0: QgsPointXY, p1: QgsPointXY) -> float:
    # East-based CCW to North-based CW
    theta = math.degrees(math.atan2(p1.y() - p0.y(), p1.x() - p0.x()))
    return (90.0 - theta) % 360.0

def _deflection(a: float, b: float) -> float:
    # minimal angular difference in [0..180]
    return abs((a - b + 180.0) % 360.0 - 180.0)

def _step_along_az(x: float, y: float, az_deg: float, step_m: float):
    theta = radians((90.0 - az_deg) % 360.0)
    return (x + step_m * cos(theta), y + step_m * sin(theta))

def _build_directed_by_FA(vl: QgsVectorLayer, acc_path: str, dem_path: str):
    feats = list(vl.getFeatures())
    fid_to_feat = {f.id(): f for f in feats}

    px, py, diag = _pixel_size_and_diag(dem_path)
    tol = FA_TOL_MULT * diag  # meters

    # endpoints and azimuths (set downstream toward larger FA)
    up_end = {}
    down_end = {}
    az = {}

    for f in feats:
        g = f.geometry()
        line = _line_coords(g)
        if len(line) < 2:
            p = g.centroid().asPoint()
            up_end[f.id()] = (p.x(), p.y())
            down_end[f.id()] = (p.x(), p.y())
            az[f.id()] = 0.0
            continue
        a = line[0]; b = line[-1]
        fa_a = _sample_raster_val(acc_path, a.x(), a.y(), 0.0)
        fa_b = _sample_raster_val(acc_path, b.x(), b.y(), 0.0)
        if fa_b > fa_a:
            up_end[f.id()] = (a.x(), a.y())
            down_end[f.id()] = (b.x(), b.y())
            az[f.id()] = _azimuth_0_360_from(a, b)
        else:
            up_end[f.id()] = (b.x(), b.y())
            down_end[f.id()] = (a.x(), a.y())
            az[f.id()] = _azimuth_0_360_from(b, a)

    # spatial index over upstream endpoints
    sidx = QgsSpatialIndex()
    for fid, (ux, uy) in up_end.items():
        feat = QgsFeature(fid)
        feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(ux, uy)))
        sidx.addFeature(feat)

    def _best_neighbor(cur_fid: int, x: float, y: float):
        cand_ids = sidx.nearestNeighbor(QgsPointXY(x, y), 8)
        best = None; best_d = tol
        fa_here = _sample_raster_val(acc_path, x, y, 0.0)
        for cid in cand_ids:
            if cid == cur_fid: continue
            ux, uy = up_end[cid]
            d = hypot(ux - x, uy - y)
            if d > tol: continue
            fa_cand = _sample_raster_val(acc_path, ux, uy, 0.0)
            # do not flow uphill in FA space
            if fa_cand + 1e-9 < fa_here:
                continue
            # enforce gentle turning (< 90°)
            if _deflection(az[cur_fid], az[cid]) >= 90.0:
                continue
            if d < best_d:
                best_d = d; best = cid
        return best

    down_map = {fid: -1 for fid in fid_to_feat.keys()}
    indeg    = defaultdict(int)

    for fid in fid_to_feat.keys():
        x, y = down_end[fid]
        j = _best_neighbor(fid, x, y)
        if j is None:
            for step_cells in AHEAD_STEPS:  # forward probe
                x2, y2 = _step_along_az(x, y, az[fid], step_cells * max(px, py))
                j = _best_neighbor(fid, x2, y2)
                if j is not None:
                    break
        if j is not None:
            down_map[fid] = j
            indeg[j] += 1

    # ---- Single-outlet consolidation (if multiple sinks remain) ----
    sinks = [fid for fid, dv in down_map.items() if dv == -1]
    if len(sinks) > 1:
        # choose main outlet by max FA at downstream endpoint
        best_sink = max(sinks, key=lambda s: _sample_raster_val(acc_path, *down_end[s], 0.0))
        for s in sinks:
            if s == best_sink: continue
            down_map[s] = best_sink
            indeg[best_sink] += 1
        sinks = [best_sink]

    # Shreve numbers (optional)
    upstream = defaultdict(list)
    for u, v in down_map.items():
        if v != -1:
            upstream[v].append(u)
    memo = {}
    def _shreve(u):
        if u in memo: return memo[u]
        ups = upstream.get(u, [])
        memo[u] = 1 if not ups else sum(_shreve(x) for x in ups)
        return memo[u]
    for k in set(list(upstream.keys()) + list(down_map.keys())):
        _shreve(k)

    # Topo order per connected component
    def _toposort(feats_ids, down):
        nbrs = defaultdict(set)
        for u in feats_ids:
            v = down.get(u, -1)
            if v != -1:
                nbrs[u].add(v); nbrs[v].add(u)
        seen=set(); comps=[]
        for u in feats_ids:
            if u in seen: continue
            q=[u]; seen.add(u); comp=[u]
            for x in q:
                for y in nbrs[x]:
                    if y not in seen: seen.add(y); q.append(y); comp.append(y)
            comps.append(comp)
        order=[]
        for comp in comps:
            indeg_local = defaultdict(int)
            for u in comp:
                v = down.get(u,-1)
                if v in comp and v!=-1: indeg_local[v]+=1
            dq=deque([u for u in comp if indeg_local[u]==0])
            od=[]
            while dq:
                u=dq.popleft(); od.append(u)
                v=down.get(u,-1)
                if v in comp and v!=-1:
                    indeg_local[v]-=1
                    if indeg_local[v]==0: dq.append(v)
            if len(od)!=len(comp):
                rest=[u for u in comp if u not in od]; od.extend(rest)
            order.extend(od)
        return order

    topo_fids = _toposort(list(fid_to_feat.keys()), down_map)
    new_id = {fid: i+1 for i, fid in enumerate(topo_fids)}
    n_outlets = sum(1 for fid in topo_fids if down_map.get(fid, -1) == -1)

    return fid_to_feat, down_map, indeg, memo, az, topo_fids, new_id, n_outlets

# ----------------------------- #
#  Slope sampler for network    #
# ----------------------------- #
def _sample_mean_slope_tan_along_line(geom, slope_raster_path: str, n_samples=BASE_SLOPE_SAMPLES):
    if (geom is None) or geom.isEmpty() or QgsWkbTypes.geometryType(geom.wkbType())!=QgsWkbTypes.LineGeometry:
        return 0.01
    L = geom.length()
    if L<=0: return 0.01
    prov = _provider_for(slope_raster_path)
    vals=[]
    for i in range(1, n_samples+1):
        d=(i/(n_samples+1.0))*L
        p=geom.interpolate(d).asPoint()
        val,ok = prov.sample(p,1)
        if ok and val is not None:
            try: vals.append(max(0.0, tan(radians(float(val)))))
            except: pass
    return sum(vals)/len(vals) if vals else 0.01

# ----------------------------- #
#  Write stream.network.dat     #
# ----------------------------- #
def _write_stream_network_FA(vl: QgsVectorLayer, slope_raster_path: str, out_dir: str, nin_mode=NIN_MODE):
    fid_to_feat, down_map, indeg_map, shreve_map, az, topo_fids, new_id, n_out = \
        _build_directed_by_FA(vl, flow_acc, elev)

    class_field  = _first_existing_field(vl, ["chanclass","class_id","class"])
    length_field = _first_existing_field(vl, ["Shape_Leng","length","len"])

    out_path = os.path.join(out_dir,"stream.network.dat")
    with open(out_path, "w") as f:
        for fid in topo_fids:
            ft   = fid_to_feat[fid]
            geom = ft.geometry()
            # length
            length_m = None
            if length_field and ft[length_field] not in (None,"",QVariant()):
                try: length_m = float(ft[length_field])
                except: pass
            if length_m is None:
                length_m = geom.length() if geom else 0.0
            # slope tan (mean along segment)
            slope_tan = _sample_mean_slope_tan_along_line(geom, slope_raster_path, n_samples=BASE_SLOPE_SAMPLES)
            # class id
            if class_field and ft[class_field] not in (None,"",QVariant()):
                try: class_id = int(ft[class_field])
                except: class_id = 1
            else:
                class_id = 1
            # Nin (enforce Nin >= 1 to avoid 0 in stream.network.dat)
            if nin_mode == "shreve":
                nin_raw = int(shreve_map.get(fid, 1))
            else:
                nin_raw = int(indeg_map.get(fid, 0))
            nin = max(1, nin_raw)
            # downstream remap
            dwn_fid = down_map.get(fid,-1)
            dwn_new = -1 if dwn_fid==-1 else int(new_id.get(dwn_fid,-1))
            # seg id
            segid = new_id[fid]
            f.write(f"{segid:d} {nin:d} {slope_tan:0.5f} {length_m:0.5f} {class_id:d} {dwn_new:d}\n")

    if n_out>1:
        print(f"[warn] {n_out} outlets detected (multiple -1).")
    else:
        print("[ok] Single outlet: only the last line has DownSegID = -1.")
    print(f"[ok] stream.network.dat written: {out_path}")

    return fid_to_feat, az, topo_fids, new_id

# ----------------------------- #
#  Write stream.map.dat (0–360) #
# ----------------------------- #
def _write_stream_map(vl: QgsVectorLayer, out_dir: str, topo_fids, new_id, azimuth_map):
    """
    Write a snake-style stream.map.dat file with per–grid-cell channel length,
    using a row/column convention consistent with the DHSVM mask/SoilDepth
    grids (row index increases from north to south, i.e., top to bottom).

    For each stream segment (line feature):

      - We walk along the polyline at approximately half-pixel spacing.
      - Each sample has a fixed arc-length (L_total / n_steps).
      - Each sample is assigned to the DEM grid cell that contains its
        midpoint.
      - For each (row, col) cell, we accumulate the total length contributed
        by all samples that fall inside that cell.

    As a result:
      - The same SegID may appear on multiple rows, tracing out a “snake-like”
        channel in grid space.
      - The sum of per-cell lengths for a given SegID is close to the
        polyline’s geometric length.
      - The length stored in stream.map.dat is per–grid-cell, which is
        consistent with the way CutBankGeometry expects the channel/road
        surface area to be defined (Area ≈ Width × Length_in_cell).
    """
    # 1) DEM geometry for row/col mapping
    rl_dem = _raster_layer_for(elev)
    px = abs(rl_dem.rasterUnitsPerPixelX())
    py = abs(rl_dem.rasterUnitsPerPixelY())
    ext = rl_dem.extent()
    xmin = ext.xMinimum()
    xmax = ext.xMaximum()
    ymin = ext.yMinimum()
    ymax = ext.yMaximum()
    ncols = rl_dem.width()
    nrows = rl_dem.height()

    # Pixel diagonal: physical upper bound on the path length inside a single cell
    diag = (px * px + py * py) ** 0.5

    # Base sampling step ≈ half a pixel, so each crossed cell is hit by several samples
    base_step = max(px, py) * 0.5

    len_field = _first_existing_field(vl, ["Shape_Leng", "length", "len"])
    w_field   = _first_existing_field(vl, ["Width", "effwidth", "hydwidth", "w", "width"])
    h_field   = _first_existing_field(vl, ["Height", "effdepth", "hyddepth", "d",
                                           "depth", "bankht", "bankheight"])

    out_path = os.path.join(out_dir, "stream.map.dat")
    with open(out_path, "w") as f:
        # Header (kept consistent with previous versions)
        f.write("###### This file has been automatically generated #####\n")
        f.write("######             EDIT WITH CARE!!!              #####\n")
        f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Workspace: {str(WS)}\n")
        f.write("#                   Segment  Cut/Bank     Cut     Segment\n")
        f.write("#  Col  Row  ID      Length   Height     Width     Aspect   SINK?\n")
        f.write("#                     (m)      (m)        (m)       (d)    (optional)\n#\n")

        # Cache features to avoid repeated iteration
        fid_to_ft = {ft.id(): ft for ft in vl.getFeatures()}

        # Loop over segments in topological order
        for fid in sorted(topo_fids, key=lambda x: new_id[x]):
            ft = fid_to_ft.get(fid)
            if ft is None:
                continue

            geom = ft.geometry()
            if geom is None or geom.isEmpty():
                continue

            segid = new_id[fid]

            # ---- Segment-level attributes ----
            # Total segment length (used as a consistency target, not written per cell)
            try:
                if len_field and ft[len_field] not in (None, "", QVariant()):
                    L_total = float(ft[len_field])
                else:
                    L_total = float(geom.length() or 0.0)
            except Exception:
                L_total = float(geom.length() or 0.0)

            # Width / height
            try:
                width = float(ft[w_field]) if w_field else 0.50
            except Exception:
                width = 0.50
            try:
                height = float(ft[h_field]) if h_field else 0.10
            except Exception:
                height = 0.10

            # Aspect: FA-directed azimuth, 0–360°, 0 = North, clockwise
            aspect = int(round(azimuth_map.get(fid, 0.0))) % 360

            # Degenerate case: zero or negative length
            if L_total <= 0.0:
                c = geom.centroid().asPoint()
                # Col: left → right
                col = int((c.x() - xmin) / px) if px > 0 else 0
                # Row: top → bottom (north → south), consistent with mask/SoilDepth
                row = int((ymax - c.y()) / py) if py > 0 else 0
                if 0 <= col < ncols and 0 <= row < nrows:
                    f.write(
                        f"{col:5d}{row:6d}{segid:6d}"
                        f"{L_total:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n"
                    )
                continue

            line = _line_coords(geom)
            if not line:
                # Fallback: centroid if asPolyline/asMultiPolyline failed
                c = geom.centroid().asPoint()
                col = int((c.x() - xmin) / px) if px > 0 else 0
                row = int((ymax - c.y()) / py) if py > 0 else 0
                if 0 <= col < ncols and 0 <= row < nrows:
                    f.write(
                        f"{col:5d}{row:6d}{segid:6d}"
                        f"{L_total:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n"
                    )
                continue

            # ---------------------------------------------------------
            # 1. Sample along the line with fixed arc-length steps
            # ---------------------------------------------------------
            n_steps = max(1, int(math.ceil(L_total / max(base_step, 1e-6))))
            step_len = L_total / n_steps

            # Accumulate length contribution per (row, col)
            cell_lengths = {}

            for k in range(n_steps):
                # Sample at step midpoints to avoid double-counting vertices
                d = ((k + 0.5) / n_steps) * L_total
                pt = geom.interpolate(d).asPoint()
                x = pt.x()
                y = pt.y()

                # Column: left → right
                col = int((x - xmin) / px) if px > 0 else 0
                # Row: top → bottom (north → south), consistent with mask/SoilDepth
                row = int((ymax - y) / py) if py > 0 else 0

                if not (0 <= col < ncols and 0 <= row < nrows):
                    continue

                key = (row, col)
                cell_lengths[key] = cell_lengths.get(key, 0.0) + step_len

            if not cell_lengths:
                # Safety fallback: centroid only
                c = geom.centroid().asPoint()
                col = int((c.x() - xmin) / px) if px > 0 else 0
                row = int((ymax - c.y()) / py) if py > 0 else 0
                if 0 <= col < ncols and 0 <= row < nrows:
                    f.write(
                        f"{col:5d}{row:6d}{segid:6d}"
                        f"{L_total:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n"
                    )
                continue

            # Optional normalization so that the sum of per-cell lengths
            # is close to the total segment length L_total
            sum_cells = sum(cell_lengths.values())
            scale = (L_total / sum_cells) if sum_cells > 0 else 1.0

            # ---------------------------------------------------------
            # 2. Emit one row per distinct (row, col) with per-cell length,
            #    and clamp length_in_cell <= pixel diagonal for safety.
            # ---------------------------------------------------------
            for (row, col), L_cell_raw in sorted(cell_lengths.items()):
                # Rescale to match the total segment length
                L_cell = L_cell_raw * scale
                # Physical upper bound: cannot exceed the pixel diagonal
                L_cell = min(L_cell, diag)

                f.write(
                    f"{col:5d}{row:6d}{segid:6d}"
                    f"{L_cell:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n"
                )

    print("[ok] stream.map.dat written (snake-style, per-cell length,"
          " aspect = 0–360°, North-bearing, clockwise).")
    return out_path

# ----------------------------- #
#  Run writers                  #
# ----------------------------- #
print("[step] Write stream.network.dat (FA-directed, 6 cols)…")
fid2f, az_map, topo_fids, new_id = _write_stream_network_FA(vl_lines, stream_slope, str(WS), nin_mode=NIN_MODE)

print("[step] Write stream.map.dat (0–360° aspect)…")
_ = _write_stream_map(vl_lines, str(WS), topo_fids, new_id, az_map)

# ----------------------------- #
#  Generate Soil Depth Binary   #
# ----------------------------- #
print("[step] Compute Soil Depth via soildepthscript.py…")
try:
    spec = importlib.util.spec_from_file_location("soildepthscript", str(SCRIPT_DIR/"soildepthscript.py"))
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    
    if hasattr(mod, "generate_soildepth"):
        bin_out = p("soildepth.bin")
        tif_out = p("soildepth.tif")
        
        mod.generate_soildepth(elev, stream_slope, flow_acc, bin_out, tif_out)
    else:
        print("[warn] generate_soildepth() function not found in soildepthscript.py.")
except Exception as e:
    print(f"[error] Soil depth generation failed: {e}")

# ----------------------------- #
#              Done             #
# ----------------------------- #
print("\nPipeline finished successfully.")
print("Artifacts:")
print(f"  - {stream_raster}")
print(f"  - {vl_lines.source()}")
print(f"  - {stream_slope}")
print(f"  - {p('stream.map.dat')}")
print(f"  - {p('stream.network.dat')}")
print(f"  - {p('stream.class.dat')}")
print(f"  - {p('soildepth.bin')}")