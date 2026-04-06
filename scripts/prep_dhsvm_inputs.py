# -*- coding: utf-8 -*-
# =====================================================================
# prep_dhsvm_inputs.py
# =====================================================================
# PURPOSE:
#   Main script to prepare all spatial inputs for DHSVM.
#   Builds the hydrologically routed stream network using QGIS/GRASS 
#   tools, and coordinates sub-modules (soildepth, channelclass, basemaps, 
#   states) to generate ready-to-use binary grids and stream files.
#
# FILE ORGANIZATION:
#   Keeps the workspace clean by separating intermediate GIS artifacts 
#   from final model inputs:
#     -> /DHSVM_input_binaries   : Flat binaries (DEM, Mask, Soil, Veg, Depth)
#     -> /DHSVM_input_streams    : Topological inputs (Map, Network, Class)
#     -> /DHSVM_input_modelstate : Python-generated initial conditions
#     -> /Intermediate_GIS       : Diagnostic rasters and shapefiles
#
# DEPENDENCIES:
#   Run natively inside the QGIS Python Console.
#
# AUTHOR: Yiyun Song
# DATE:   2026-04-06
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
MIN_SRC_CELLS      = 60           # r.stream.extract threshold in number of cells
NIN_MODE           = "indegree"   # "indegree" or "shreve" for Col(2)
BASE_SLOPE_SAMPLES = 12           # samples per line when computing tan(slope)
SNAP_TOL_LIST      = [0.75, 1.25, 1.75]  
FA_TOL_MULT        = 1.10         # FA matching tolerance = 1.10 × pixel-diagonal
AHEAD_STEPS        = (1.0, 2.0)   # forward probing (in pixel units) to bridge tiny gaps

# ----------------------------- #
#  Portable path & Directories  #
# ----------------------------- #
SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR.parent

# Define robust output directories
DIR_BINARIES = WS / "DHSVM_input_binaries"
DIR_STREAMS  = WS / "DHSVM_input_streams"
DIR_STATES   = WS / "modelstate"
DIR_INTERMED = WS / "Intermediate_GIS"

# Ensure directories exist
for d in [DIR_BINARIES, DIR_STREAMS, DIR_STATES, DIR_INTERMED]:
    d.mkdir(parents=True, exist_ok=True)

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

DEM_DIR = WS / "Reprojected_DEM"
primary_dem_candidates = [
    p("Reprojected_DEM", "Cropped_DEM.tif"),
    p("Reprojected_DEM", "dem.tif"),
    p("Reprojected_DEM", "elev.tif"),
    p("Reprojected_DEM", "elev_clipped.tif"),
]

generic_tifs = []
if DEM_DIR.exists():
    generic_tifs = [str(x) for x in sorted(DEM_DIR.glob("*.tif"))]

DEM_CANDIDATES = primary_dem_candidates + generic_tifs + [p("dem.tif"), p("elev.tif")]

def _find_watershed_shp():
    wd = WS/"Reprojected_Watersheds"
    if wd.exists():
        shp = sorted([str(x) for x in wd.glob("*.shp")], key=lambda s: ("watershed" not in s.lower(), s.lower()))
        if shp: return shp[0]
    shp2 = sorted([str(x) for x in WS.glob("*.shp")], key=lambda s: ("watershed" not in s.lower(), s.lower()))
    return shp2[0] if shp2 else None

elev_raw = _pick_first_existing(DEM_CANDIDATES)
wshed    = _find_watershed_shp()
if not elev_raw:
    raise FileNotFoundError("Projected DEM not found. Put it in Reprojected_DEM/.")

print("\n=======================================================")
print("  DHSVM SPATIAL PIPELINE INITIALIZING")
print("=======================================================")
if wshed: print(f"[ok] Watershed polygon: {Path(wshed).name}")
print(f"[ok] Master DEM: {Path(elev_raw).name}")

# Route derived files to their specific directories
elev           = str(DIR_INTERMED / "elev_clipped.tif") if wshed else elev_raw
flow_dir       = str(DIR_INTERMED / "flow_dir.tif")
flow_acc       = str(DIR_INTERMED / "flow_acc.tif")
stream_raster  = str(DIR_INTERMED / "stream_raster.tif")
stream_vec     = str(DIR_INTERMED / "stream_order.shp")   
streamfile_fbk = str(DIR_INTERMED / "streamfile.shp")     
stream_slope   = str(DIR_INTERMED / "stream_slope.tif")

# ----------------------------- #
#         Raster hotfix         #
# ----------------------------- #
_RASTER_LAYER_CACHE = {}

def _raster_layer_for(raster_path: str) -> QgsRasterLayer:
    key = os.path.abspath(raster_path)
    rl = _RASTER_LAYER_CACHE.get(key)
    if (rl is None) or (not isinstance(rl, QgsRasterLayer)) or (not rl.isValid()):
        rl = QgsRasterLayer(raster_path, f"ras::{os.path.basename(raster_path)}")
        if not rl.isValid():
            raise RuntimeError(f"Cannot load raster: {raster_path}")
        _RASTER_LAYER_CACHE[key] = rl
    return rl

def _provider_for(raster_path: str):
    return _raster_layer_for(raster_path).dataProvider()

def _sample_raster_val(ras_path: str, x: float, y: float, default=0.0):
    prov = _provider_for(ras_path)
    val, ok = prov.sample(QgsPointXY(x, y), 1)
    if ok and (val is not None):
        try:
            fv = float(val)
            if math.isfinite(fv): return fv
        except: pass
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
    return ("LineString" in _geom_type(path))

def _feature_count(path: str) -> int:
    v = QgsVectorLayer(path, "chk", "ogr")
    return sum(1 for _ in v.getFeatures()) if v.isValid() else 0

def _to_float(v, default=None):
    try: return float(v) if (v is not None) else default
    except: return default

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
print(f"[ok] Working DEM ready: {Path(elev).name}")

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

# ----------------------------- #
#  r.stream.extract raster+vec  #
# ----------------------------- #
print("[step] GRASS r.stream.extract → extracting streams…")
processing.run("grass7:r.stream.extract", {
    'elevation': elev, 'accumulation': flow_acc, 'direction': flow_dir,
    'threshold': MIN_SRC_CELLS, 'stream_raster': stream_raster,
    'stream_vector': stream_vec, 'memory': 300, '-m': True,
    'GRASS_REGION_PARAMETER': elev, 'GRASS_REGION_CELLSIZE_PARAMETER': 0,
    'GRASS_OUTPUT_TYPE_PARAMETER': 0
})
for _ in range(12):
    if os.path.exists(stream_raster) and os.path.getsize(stream_raster)>0: break
    time.sleep(0.25)

vec_count = _feature_count(stream_vec)

def _force_lines_from_raster(raster_path: str, out_path: str) -> str:
    trials = [{'type': 0}, {'type': 'line'}, {'type': 1}, {'type': 'line', '-s': True}]
    for opt in trials:
        try:
            processing.run("grass7:r.to.vect", {
                'input': raster_path, 'type': opt.get('type', 1), 'output': out_path,
                **({'-s': True} if opt.get('-s') else {}),
                'GRASS_REGION_PARAMETER': raster_path, 'GRASS_REGION_CELLSIZE_PARAMETER': 0
            })
            if _feature_count(out_path) > 0 and _is_line_layer(out_path): return out_path
        except: pass
    raise RuntimeError("Failed to obtain LineString geometry.")

vector_path = stream_vec
if vec_count == 0 or not _is_line_layer(vector_path):
    vector_path = _force_lines_from_raster(stream_raster, streamfile_fbk)

vl_lines = QgsVectorLayer(vector_path,"streams_live","ogr")

# ----------------------------- #
#    Slope raster (DEGREES)     #
# ----------------------------- #
print("[step] GRASS r.slope.aspect → calculating slope…")
processing.run("grass7:r.slope.aspect",{
    'elevation':elev,'slope':stream_slope,'format':1,
    'GRASS_REGION_PARAMETER':elev,'GRASS_REGION_CELLSIZE_PARAMETER':0
})

# ----------------------------- #
#  Normalize fields for map.dat #
# ----------------------------- #
def _ensure_fields_rowcol_len(layer: QgsVectorLayer):
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

_ensure_fields_rowcol_len(vl_lines)

# ----------------------------- #
#  Mean slope fields on lines   #
# ----------------------------- #
def _sample_mean_slope_deg(geom, slope_raster_path: str, n_samples=BASE_SLOPE_SAMPLES):
    if (geom is None) or geom.isEmpty(): return 0.0
    L = geom.length()
    if L<=0: return 0.0
    prov = _provider_for(slope_raster_path); vals=[]
    for i in range(1, n_samples+1):
        d=(i/(n_samples+1.0))*L
        p=geom.interpolate(d).asPoint()
        ok,val = prov.sample(p,1)
        if ok and val is not None:
            try: vals.append(max(0.0, float(val)))
            except: pass
    return sum(vals)/len(vals) if vals else 0.0

vl_lines.startEditing()
if "slope_deg" not in vl_lines.fields().names():
    vl_lines.addAttribute(QgsField("slope_deg", QVariant.Double))
    vl_lines.updateFields()
idx = vl_lines.fields().indexFromName("slope_deg")
for ft in vl_lines.getFeatures():
    ft[idx] = _sample_mean_slope_deg(ft.geometry(), stream_slope, BASE_SLOPE_SAMPLES)
    vl_lines.updateFeature(ft)
vl_lines.commitChanges()

# ----------------------------- #
#  Accumulation mean (meanmsq)  #
# ----------------------------- #
def _sample_line_raster_path(raster_path: str, geom: QgsGeometry, n=15):
    if (geom is None) or geom.isEmpty(): return []
    L = geom.length() or 0.0
    if L <= 0: return []
    prov = _provider_for(raster_path); vals=[]
    for k in range(n):
        d = (k/(n-1.0))*L
        p = geom.interpolate(d).asPoint()
        v, ok = prov.sample(p, 1)
        if ok and (v is not None):
            try: vals.append(float(v))
            except: pass
    return vals

vl_lines.startEditing()
if vl_lines.fields().indexOf("meanmsq") < 0:
    vl_lines.addAttribute(QgsField("meanmsq", QVariant.Double)); vl_lines.updateFields()
for ft in vl_lines.getFeatures():
    vals = [v for v in _sample_line_raster_path(flow_acc, ft.geometry(), n=15) if v is not None and v > 0]
    val = vals[len(vals)//2] if vals else cell_area
    ft["meanmsq"] = float(val)
    vl_lines.updateFeature(ft)
vl_lines.commitChanges()

# ----------------------------- #
#    Channel classes (external) #
# ----------------------------- #
print("[step] Routing classifications to DHSVM_input_streams/…")
try:
    spec = importlib.util.spec_from_file_location("channelclass", str(SCRIPT_DIR/"channelclass.py"))
    mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.channelclassfun(vl_lines.source(), str(DIR_STREAMS))
except Exception as e:
    print(f"[warn] Channel classification skipped: {e}")

# ----------------------------- #
#  Width/Height for map.dat     #
# ----------------------------- #
vl_lines.startEditing()
if vl_lines.fields().indexOf("Width") < 0: vl_lines.addAttribute(QgsField("Width", QVariant.Double))
if vl_lines.fields().indexOf("Height") < 0: vl_lines.addAttribute(QgsField("Height", QVariant.Double))
vl_lines.updateFields()

names = {f.name().lower(): f.name() for f in vl_lines.fields()}
w_name = next((names[nm] for nm in ("effwidth","hydwidth","w","width") if nm in names), None)
h_name = next((names[nm] for nm in ("effdepth","hyddepth","d","depth","bankht","bankheight") if nm in names), None)

for ft in vl_lines.getFeatures():
    ft["Width"] = float(ft[w_name]) if w_name else 0.50
    ft["Height"] = float(ft[h_name]) if h_name else 0.10
    vl_lines.updateFeature(ft)
vl_lines.commitChanges()

# ----------------------------- #
#  FA-directed network builder  #
# ----------------------------- #
def _line_coords(geom: QgsGeometry):
    if geom.isMultipart():
        parts = geom.asMultiPolyline()
        return max(parts, key=len) if parts else []
    else: return geom.asPolyline()

def _azimuth_0_360_from(p0: QgsPointXY, p1: QgsPointXY) -> float:
    theta = math.degrees(math.atan2(p1.y() - p0.y(), p1.x() - p0.x()))
    return (90.0 - theta) % 360.0

def _deflection(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)

def _step_along_az(x: float, y: float, az_deg: float, step_m: float):
    theta = radians((90.0 - az_deg) % 360.0)
    return (x + step_m * cos(theta), y + step_m * sin(theta))

def _build_directed_by_FA(vl: QgsVectorLayer, acc_path: str, dem_path: str):
    feats = list(vl.getFeatures())
    fid_to_feat = {f.id(): f for f in feats}

    px, py, diag = _pixel_size_and_diag(dem_path)
    tol = FA_TOL_MULT * diag 

    up_end, down_end, az = {}, {}, {}
    for f in feats:
        g = f.geometry()
        line = _line_coords(g)
        if len(line) < 2:
            p = g.centroid().asPoint()
            up_end[f.id()] = (p.x(), p.y()); down_end[f.id()] = (p.x(), p.y()); az[f.id()] = 0.0
            continue
        a, b = line[0], line[-1]
        fa_a = _sample_raster_val(acc_path, a.x(), a.y(), 0.0)
        fa_b = _sample_raster_val(acc_path, b.x(), b.y(), 0.0)
        if fa_b > fa_a:
            up_end[f.id()] = (a.x(), a.y()); down_end[f.id()] = (b.x(), b.y())
            az[f.id()] = _azimuth_0_360_from(a, b)
        else:
            up_end[f.id()] = (b.x(), b.y()); down_end[f.id()] = (a.x(), a.y())
            az[f.id()] = _azimuth_0_360_from(b, a)

    sidx = QgsSpatialIndex()
    for fid, (ux, uy) in up_end.items():
        feat = QgsFeature(fid); feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(ux, uy)))
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
            if fa_cand + 1e-9 < fa_here: continue
            if _deflection(az[cur_fid], az[cid]) >= 90.0: continue
            if d < best_d: best_d = d; best = cid
        return best

    down_map = {fid: -1 for fid in fid_to_feat.keys()}
    indeg    = defaultdict(int)

    for fid in fid_to_feat.keys():
        x, y = down_end[fid]
        j = _best_neighbor(fid, x, y)
        if j is None:
            for step_cells in AHEAD_STEPS:
                x2, y2 = _step_along_az(x, y, az[fid], step_cells * max(px, py))
                j = _best_neighbor(fid, x2, y2)
                if j is not None: break
        if j is not None:
            down_map[fid] = j; indeg[j] += 1

    sinks = [fid for fid, dv in down_map.items() if dv == -1]
    if len(sinks) > 1:
        best_sink = max(sinks, key=lambda s: _sample_raster_val(acc_path, *down_end[s], 0.0))
        for s in sinks:
            if s != best_sink: down_map[s] = best_sink; indeg[best_sink] += 1

    upstream = defaultdict(list)
    for u, v in down_map.items():
        if v != -1: upstream[v].append(u)
    memo = {}
    def _shreve(u):
        if u in memo: return memo[u]
        ups = upstream.get(u, [])
        memo[u] = 1 if not ups else sum(_shreve(x) for x in ups)
        return memo[u]
    for k in set(list(upstream.keys()) + list(down_map.keys())): _shreve(k)

    def _toposort(feats_ids, down):
        nbrs = defaultdict(set)
        for u in feats_ids:
            v = down.get(u, -1)
            if v != -1: nbrs[u].add(v); nbrs[v].add(u)
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
            if len(od)!=len(comp): rest=[u for u in comp if u not in od]; od.extend(rest)
            order.extend(od)
        return order

    topo_fids = _toposort(list(fid_to_feat.keys()), down_map)
    new_id = {fid: i+1 for i, fid in enumerate(topo_fids)}
    n_outlets = sum(1 for fid in topo_fids if down_map.get(fid, -1) == -1)

    return fid_to_feat, down_map, indeg, memo, az, topo_fids, new_id, n_outlets

def _sample_mean_slope_tan_along_line(geom, slope_raster_path: str, n_samples=BASE_SLOPE_SAMPLES):
    if (geom is None) or geom.isEmpty(): return 0.01
    L = geom.length(); prov = _provider_for(slope_raster_path); vals=[]
    for i in range(1, n_samples+1):
        d=(i/(n_samples+1.0))*L
        p=geom.interpolate(d).asPoint(); val,ok = prov.sample(p,1)
        if ok and val is not None:
            try: vals.append(max(0.0, tan(radians(float(val)))))
            except: pass
    return sum(vals)/len(vals) if vals else 0.01

def _first_existing_field(vl: QgsVectorLayer, candidates):
    names = {n.lower(): n for n in vl.fields().names()}
    for cand in candidates:
        if cand.lower() in names: return names[cand.lower()]
    return None

def _write_stream_network_FA(vl: QgsVectorLayer, slope_raster_path: str, out_dir: Path, nin_mode=NIN_MODE):
    fid_to_feat, down_map, indeg_map, shreve_map, az, topo_fids, new_id, n_out = _build_directed_by_FA(vl, flow_acc, elev)
    class_field  = _first_existing_field(vl, ["chanclass","class_id","class"])
    length_field = _first_existing_field(vl, ["Shape_Leng","length","len"])

    out_path = out_dir / "stream.network.dat"
    with open(out_path, "w") as f:
        for fid in topo_fids:
            ft = fid_to_feat[fid]; geom = ft.geometry()
            length_m = _to_float(ft[length_field]) if length_field else (geom.length() if geom else 0.0)
            slope_tan = _sample_mean_slope_tan_along_line(geom, slope_raster_path)
            class_id = int(ft[class_field]) if class_field and _to_float(ft[class_field]) else 1
            nin_raw = int(shreve_map.get(fid, 1)) if nin_mode == "shreve" else int(indeg_map.get(fid, 0))
            nin = max(1, nin_raw)
            dwn_fid = down_map.get(fid,-1)
            dwn_new = -1 if dwn_fid==-1 else int(new_id.get(dwn_fid,-1))
            segid = new_id[fid]
            f.write(f"{segid:d} {nin:d} {slope_tan:0.5f} {length_m:0.5f} {class_id:d} {dwn_new:d}\n")
    return fid_to_feat, az, topo_fids, new_id

def _write_stream_map(vl: QgsVectorLayer, out_dir: Path, topo_fids, new_id, azimuth_map):
    rl_dem = _raster_layer_for(elev)
    px = abs(rl_dem.rasterUnitsPerPixelX()); py = abs(rl_dem.rasterUnitsPerPixelY())
    xmin = rl_dem.extent().xMinimum(); ymax = rl_dem.extent().yMaximum()
    ncols = rl_dem.width(); nrows = rl_dem.height()
    diag = (px * px + py * py) ** 0.5
    base_step = max(px, py) * 0.5

    len_field = _first_existing_field(vl, ["Shape_Leng", "length", "len"])
    w_field   = _first_existing_field(vl, ["Width", "effwidth", "hydwidth", "w", "width"])
    h_field   = _first_existing_field(vl, ["Height", "effdepth", "hyddepth", "d", "depth", "bankht", "bankheight"])

    out_path = out_dir / "stream.map.dat"
    with open(out_path, "w") as f:
        f.write("###### This file has been automatically generated #####\n")
        f.write("#  Col  Row  ID     Length  Height    Width    Aspect   SINK?\n")
        fid_to_ft = {ft.id(): ft for ft in vl.getFeatures()}

        for fid in sorted(topo_fids, key=lambda x: new_id[x]):
            ft = fid_to_ft.get(fid)
            geom = ft.geometry()
            if not geom or geom.isEmpty(): continue
            segid = new_id[fid]
            L_total = _to_float(ft[len_field]) if len_field else (geom.length() or 0.0)
            width = _to_float(ft[w_field]) if w_field else 0.50
            height = _to_float(ft[h_field]) if h_field else 0.10
            aspect = int(round(azimuth_map.get(fid, 0.0))) % 360

            n_steps = max(1, int(math.ceil(L_total / max(base_step, 1e-6))))
            step_len = L_total / n_steps
            cell_lengths = {}

            for k in range(n_steps):
                pt = geom.interpolate(((k + 0.5) / n_steps) * L_total).asPoint()
                col = int((pt.x() - xmin) / px) if px > 0 else 0
                row = int((ymax - pt.y()) / py) if py > 0 else 0
                if 0 <= col < ncols and 0 <= row < nrows:
                    cell_lengths[(row, col)] = cell_lengths.get((row, col), 0.0) + step_len

            sum_cells = sum(cell_lengths.values())
            scale = (L_total / sum_cells) if sum_cells > 0 else 1.0
            for (row, col), L_cell_raw in sorted(cell_lengths.items()):
                L_cell = min(L_cell_raw * scale, diag)
                f.write(f"{col:5d}{row:6d}{segid:6d}{L_cell:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n")
    return out_path

# ----------------------------- #
#  Run writers & Submodules     #
# ----------------------------- #
print("[step] Writing topologic stream files to DHSVM_input_streams/...")
fid2f, az_map, topo_fids, new_id = _write_stream_network_FA(vl_lines, stream_slope, DIR_STREAMS)
_ = _write_stream_map(vl_lines, DIR_STREAMS, topo_fids, new_id, az_map)

print("[step] Compute Soil Depth Binary (soildepthscript.py)...")
try:
    spec = importlib.util.spec_from_file_location("soildepthscript", str(SCRIPT_DIR/"soildepthscript.py"))
    mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    if hasattr(mod, "generate_soildepth"):
        bin_out = str(DIR_BINARIES / "soildepth.bin")
        tif_out = str(DIR_INTERMED / "soildepth.tif")
        mod.generate_soildepth(elev, stream_slope, flow_acc, bin_out, tif_out)
except Exception as e:
    print(f"[error] Soil depth generation failed: {e}")

print("[step] Generate DHSVM Base Maps (dem_to_dhsvm_bins.py)...")
try:
    spec_base = importlib.util.spec_from_file_location("dem_to_dhsvm_bins", str(SCRIPT_DIR/"dem_to_dhsvm_bins.py"))
    mod_base  = importlib.util.module_from_spec(spec_base)
    spec_base.loader.exec_module(mod_base) 
    # The script automatically runs generate_basemaps() upon execution.
except Exception as e:
    print(f"[error] Base map generation failed: {e}")

print("[step] Generate DHSVM Initial States (generate_dhsvm_states.py)...")
try:
    spec_states = importlib.util.spec_from_file_location("generate_dhsvm_states", str(SCRIPT_DIR/"generate_dhsvm_states.py"))
    mod_states  = importlib.util.module_from_spec(spec_states)
    spec_states.loader.exec_module(mod_states)
    # The script automatically runs generate_grid_states() and generate_channel_state() upon execution.
except Exception as e:
    print(f"[error] Initial state generation failed: {e}")

# ----------------------------- #
#              Done             #
# ----------------------------- #
print("\n=======================================================")
print("  PIPELINE FINISHED SUCCESSFULLY ")
print("=======================================================")
print(f"STREAM FILES:  {DIR_STREAMS.name}/")
print(f"BINARY GRIDS:  {DIR_BINARIES.name}/")
print(f"MODEL STATES:  {DIR_STATES.name}/")
print(f"INTERMEDIATE:  {DIR_INTERMED.name}/")