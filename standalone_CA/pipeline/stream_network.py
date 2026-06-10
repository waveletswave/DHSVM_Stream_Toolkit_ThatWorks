# -*- coding: utf-8 -*-
# =====================================================================
# stream_network.py  -  standalone port of prep's directed-network + writers
#
# Sub-step C of the #6 vector-IO rebuild. Reproduces:
#   _build_directed_by_FA   (prep line 434): FA-directed topology, spatial-index
#       nearest-neighbour linking, sink resolution, Shreve + propagated order +
#       topological sort + new_id renumbering.
#   _write_stream_network_FA (prep line 578) -> stream.network.dat
#   _write_stream_map        (prep line 612) -> stream.map.dat
#
# Reads streamfile_attr.shp (geometry + arcid/Shape_Leng/Row/Col/slope_deg/
# meanmsq + chanclass written back by channelclass_standalone.py). geopandas/
# shapely/numpy replace QgsVectorLayer/QgsGeometry/QgsSpatialIndex; raster
# sampling via InvGeoTransform+floor (verified == provider.sample in sub-step B).
#
# Verbatim constants from prep:
#   NIN_MODE="propagated", FA_TOL_MULT=1.10, AHEAD_STEPS=(1.0,2.0),
#   BASE_SLOPE_SAMPLES=12. Tier A1: outlet down-id is 0 (not -1).
#   Tier D2: Row/Col top-left origin. azimuth = north-azimuth (CW from North).
#
# Row order of streamfile_attr.shp aligns 1:1 with the reference (sub-step A),
# so fid == row index; new_id renumbering should match the reference.
#
# Run:  python3 stream_network.py
# =====================================================================

import math
from math import hypot, radians, cos, sin, tan, atan2, degrees
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
from osgeo import gdal

# ----------------------------- paths (fold into paths.py later) -------------
OUT_DIR     = Path("/work/ys451/dhsvm_ca/standalone_dev/outputs")
STREAM_IN   = OUT_DIR / "streamfile_attr.shp"          # attrs + chanclass written back
STREAMS_OUT = OUT_DIR / "DHSVM_input_streams"

ELEV_TIF    = OUT_DIR / "elev_clipped.tif"
SLOPE_TIF   = OUT_DIR / "slope_filled.tif"
FACC_TIF    = OUT_DIR / "flow_acc.tif"

# verbatim prep constants
NIN_MODE           = "propagated"
FA_TOL_MULT        = 1.10
AHEAD_STEPS        = (1.0, 2.0)
BASE_SLOPE_SAMPLES = 12


# ----------------------------- raster sampling ------------------------------
def _open_raster(path):
    ds = gdal.Open(str(path))
    arr = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    gt = ds.GetGeoTransform()
    inv = gdal.InvGeoTransform(gt)
    return arr, inv


def _sample(arr, inv, x, y, default=0.0):
    col = int(inv[0] + inv[1] * x + inv[2] * y)
    row = int(inv[3] + inv[4] * x + inv[5] * y)
    h, w = arr.shape
    if 0 <= row < h and 0 <= col < w:
        v = arr[row, col]
        if v is not None and np.isfinite(v):
            return float(v)
    return default


# ----------------------------- geometry helpers (verbatim) ------------------
def _line_coords(geom):
    """Mirror prep _line_coords: longest part of a multilinestring."""
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, MultiLineString):
        parts = list(geom.geoms)
        if not parts:
            return []
        geom = max(parts, key=lambda g: g.length)
    return list(geom.coords)


def _azimuth_0_360(p0, p1):
    """North azimuth (CW from North), verbatim from prep _azimuth_0_360_from."""
    theta = degrees(atan2(p1[1] - p0[1], p1[0] - p0[0]))
    return (90.0 - theta) % 360.0


def _deflection(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def _step_along_az(x, y, az_deg, step_m):
    theta = radians((90.0 - az_deg) % 360.0)
    return (x + step_m * cos(theta), y + step_m * sin(theta))


def _interp(geom_longest_coords, full_geom, d):
    """Interpolate distance d along the geometry (use the shapely line directly)."""
    return full_geom.interpolate(d)


# ----------------------------- directed network (verbatim logic) ------------
def build_directed_by_FA(gdf, facc_arr, facc_inv, px, py, diag):
    tol = FA_TOL_MULT * diag

    # fid == row index (sub-step A: order aligned with reference)
    fids = list(range(len(gdf)))
    geoms = list(gdf.geometry)

    up_end, down_end, az = {}, {}, {}
    for fid in fids:
        g = geoms[fid]
        line = _line_coords(g)
        if len(line) < 2:
            c = g.centroid
            up_end[fid] = (c.x, c.y); down_end[fid] = (c.x, c.y); az[fid] = 0.0
            continue
        a, b = line[0], line[-1]
        fa_a = _sample(facc_arr, facc_inv, a[0], a[1], 0.0)
        fa_b = _sample(facc_arr, facc_inv, b[0], b[1], 0.0)
        if fa_b > fa_a:
            up_end[fid] = (a[0], a[1]); down_end[fid] = (b[0], b[1])
            az[fid] = _azimuth_0_360(a, b)
        else:
            up_end[fid] = (b[0], b[1]); down_end[fid] = (a[0], a[1])
            az[fid] = _azimuth_0_360(b, a)

    # spatial index over UP endpoints. QGIS used nearestNeighbor(pt, 8); shapely
    # STRtree.query_nearest does not take k= in 2.1 (returns only equidistant
    # nearest), so we use an exact brute-force k-NN -- trivial for ~41 points and
    # dependency-free. Returns up to 8 nearest up-endpoint fids by distance, which
    # _best_neighbor then filters by tol/fa/deflection, matching prep exactly.
    up_xy = np.array([up_end[fid] for fid in fids], dtype=float)

    def _knn(x, y, k=8):
        d2 = (up_xy[:, 0] - x) ** 2 + (up_xy[:, 1] - y) ** 2
        order = np.argsort(d2, kind="stable")[:k]
        return [fids[int(i)] for i in order]

    def _best_neighbor(cur_fid, x, y):
        cand = _knn(x, y, 8)
        best = None; best_d = tol
        fa_here = _sample(facc_arr, facc_inv, x, y, 0.0)
        for cid in cand:
            if cid == cur_fid:
                continue
            ux, uy = up_end[cid]
            d = hypot(ux - x, uy - y)
            if d > tol:
                continue
            fa_cand = _sample(facc_arr, facc_inv, ux, uy, 0.0)
            if fa_cand + 1e-9 < fa_here:
                continue
            if _deflection(az[cur_fid], az[cid]) >= 90.0:
                continue
            # Strictly-closer wins. At a junction two upstream endpoints can
            # coincide with this downstream endpoint (d=0); the link is then
            # geometrically undefined and any tie-break is arbitrary. We do NOT
            # reverse-engineer QGIS's R-tree order (it has no ground truth here);
            # the k-NN is distance-sorted, so on a d=0 tie the first-seen
            # candidate wins, which happens to be the straighter (lower-
            # deflection) neighbour -- consistent with flow momentum. The
            # resulting network still satisfies DHSVM's requirements (single
            # outlet, all-reachable, dense order); see validation_log.
            if d < best_d:
                best_d = d; best = cid
        return best

    down_map = {fid: -1 for fid in fids}
    indeg = defaultdict(int)

    for fid in fids:
        x, y = down_end[fid]
        j = _best_neighbor(fid, x, y)
        if j is None:
            for step_cells in AHEAD_STEPS:
                x2, y2 = _step_along_az(x, y, az[fid], step_cells * max(px, py))
                j = _best_neighbor(fid, x2, y2)
                if j is not None:
                    break
        if j is not None:
            down_map[fid] = j; indeg[j] += 1

    sinks = [fid for fid, dv in down_map.items() if dv == -1]
    if len(sinks) > 1:
        best_sink = max(sinks, key=lambda s: _sample(facc_arr, facc_inv, *down_end[s], 0.0))
        for s in sinks:
            if s != best_sink:
                down_map[s] = best_sink; indeg[best_sink] += 1

    # Shreve (kept as-is)
    upstream = defaultdict(list)
    for u, v in down_map.items():
        if v != -1:
            upstream[v].append(u)
    memo = {}
    def _shreve(u):
        if u in memo:
            return memo[u]
        ups = upstream.get(u, [])
        memo[u] = 1 if not ups else sum(_shreve(x) for x in ups)
        return memo[u]
    for k in set(list(upstream.keys()) + list(down_map.keys())):
        _shreve(k)

    # toposort (verbatim)
    def _toposort(feats_ids, down):
        nbrs = defaultdict(set)
        for u in feats_ids:
            v = down.get(u, -1)
            if v != -1:
                nbrs[u].add(v); nbrs[v].add(u)
        seen = set(); comps = []
        for u in feats_ids:
            if u in seen:
                continue
            q = [u]; seen.add(u); comp = [u]
            for x in q:
                for y in nbrs[x]:
                    if y not in seen:
                        seen.add(y); q.append(y); comp.append(y)
            comps.append(comp)
        order = []
        for comp in comps:
            indeg_local = defaultdict(int)
            for u in comp:
                v = down.get(u, -1)
                if v in comp and v != -1:
                    indeg_local[v] += 1
            dq = deque([u for u in comp if indeg_local[u] == 0])
            od = []
            while dq:
                u = dq.popleft(); od.append(u)
                v = down.get(u, -1)
                if v in comp and v != -1:
                    indeg_local[v] -= 1
                    if indeg_local[v] == 0:
                        dq.append(v)
            if len(od) != len(comp):
                rest = [u for u in comp if u not in od]; od.extend(rest)
            order.extend(od)
        return order

    topo_fids = _toposort(fids, down_map)

    # propagated order (Tier A): headwater first
    prop_order = {}
    for fid in topo_fids:
        ups = upstream.get(fid, [])
        prop_order[fid] = 1 if not ups else max(prop_order[u] for u in ups) + 1

    new_id = {fid: i + 1 for i, fid in enumerate(topo_fids)}

    return (geoms, down_map, indeg, memo, prop_order, az, topo_fids, new_id, up_end, down_end)


# ----------------------------- slope tan sampler (verbatim, line 561) -------
def _sample_mean_slope_tan(full_geom, slope_arr, slope_inv, n_samples=BASE_SLOPE_SAMPLES):
    if full_geom is None or full_geom.is_empty:
        return 0.01
    L = full_geom.length
    vals = []
    for i in range(1, n_samples + 1):
        d = (i / (n_samples + 1.0)) * L
        p = full_geom.interpolate(d)
        v = _sample(slope_arr, slope_inv, p.x, p.y, default=None)
        if v is not None:
            try:
                vals.append(max(0.0, tan(radians(float(v)))))
            except Exception:
                pass
    return sum(vals) / len(vals) if vals else 0.01


def _first_existing(cols, candidates):
    names = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in names:
            return names[cand.lower()]
    return None


# ----------------------------- writers --------------------------------------
def write_stream_network(gdf, topo, dn, new_id, prop, shreve, indeg,
                         slope_arr, slope_inv, out_dir, nin_mode=NIN_MODE):
    class_field = _first_existing(gdf.columns, ["chanclass", "class_id", "class"])
    length_field = _first_existing(gdf.columns, ["Shape_Leng", "length", "len"])
    geoms = list(gdf.geometry)

    out_path = out_dir / "stream.network.dat"
    with open(out_path, "w") as f:
        for fid in topo:
            geom = geoms[fid]
            if length_field is not None:
                length_m = float(gdf.iloc[fid][length_field])
            else:
                length_m = geom.length if geom is not None else 0.0
            slope_tan = _sample_mean_slope_tan(geom, slope_arr, slope_inv)
            if class_field is not None:
                cv = gdf.iloc[fid][class_field]
                class_id = int(cv) if (cv is not None and float(cv)) else 1
            else:
                class_id = 1

            if nin_mode == "propagated":
                nin_raw = int(prop.get(fid, 1))
            elif nin_mode == "shreve":
                nin_raw = int(shreve.get(fid, 1))
            else:
                nin_raw = int(indeg.get(fid, 0))
            nin = max(1, nin_raw)

            dwn_fid = dn.get(fid, -1)
            dwn_new = 0 if dwn_fid == -1 else int(new_id.get(dwn_fid, 0))   # Tier A1: 0 for outlet

            segid = new_id[fid]
            f.write(f"{segid:d} {nin:d} {slope_tan:0.5f} {length_m:0.5f} {class_id:d} {dwn_new:d}\n")
    return out_path


def write_stream_map(gdf, topo, new_id, az_map, px, py, xmin, ymax, ncols, nrows, out_dir):
    diag = (px * px + py * py) ** 0.5
    base_step = max(px, py) * 0.5
    geoms = list(gdf.geometry)

    len_field = _first_existing(gdf.columns, ["Shape_Leng", "length", "len"])
    w_field = _first_existing(gdf.columns, ["Width", "effwidth", "hydwidth", "w", "width"])
    h_field = _first_existing(gdf.columns, ["Height", "effdepth", "hyddepth", "d", "depth", "bankht", "bankheight"])

    out_path = out_dir / "stream.map.dat"
    with open(out_path, "w") as f:
        f.write("###### This file has been automatically generated #####\n")
        f.write("#  Col  Row  ID     Length  Height    Width    Aspect   SINK?\n")
        for fid in sorted(topo, key=lambda x: new_id[x]):
            geom = geoms[fid]
            if geom is None or geom.is_empty:
                continue
            segid = new_id[fid]
            L_total = float(gdf.iloc[fid][len_field]) if len_field is not None else (geom.length or 0.0)
            width = float(gdf.iloc[fid][w_field]) if w_field is not None else 0.50
            height = float(gdf.iloc[fid][h_field]) if h_field is not None else 0.10
            aspect = int(round(az_map.get(fid, 0.0))) % 360

            n_steps = max(1, int(math.ceil(L_total / max(base_step, 1e-6))))
            cell_lengths = {}
            for k in range(n_steps):
                pt = geom.interpolate(((k + 0.5) / n_steps) * L_total)
                col = int((pt.x - xmin) / px) if px > 0 else 0
                row = int((ymax - pt.y) / py) if py > 0 else 0
                if 0 <= col < ncols and 0 <= row < nrows:
                    cell_lengths[(row, col)] = cell_lengths.get((row, col), 0.0) + (L_total / n_steps)

            sum_cells = sum(cell_lengths.values())
            scale = (L_total / sum_cells) if sum_cells > 0 else 1.0
            for (row, col), L_cell_raw in sorted(cell_lengths.items()):
                L_cell = min(L_cell_raw * scale, diag)
                f.write(f"{col:5d}{row:6d}{segid:6d}{L_cell:11.4f}{height:11.4f}{width:10.4f}{aspect:11d}\n")
    return out_path


# ----------------------------- driver ---------------------------------------
def main():
    for tag, p in [("streamfile_attr", STREAM_IN), ("elev", ELEV_TIF),
                   ("slope", SLOPE_TIF), ("flow_acc", FACC_TIF)]:
        if not p.exists():
            raise FileNotFoundError(f"[error] missing {tag}: {p}")
    STREAMS_OUT.mkdir(parents=True, exist_ok=True)

    gdf = gpd.read_file(STREAM_IN)
    print(f"[stream_network] read {len(gdf)} features")

    # DEM grid geometry (top-left)
    import rasterio
    dem = rasterio.open(ELEV_TIF)
    px = abs(dem.transform.a); py = abs(dem.transform.e)
    xmin = dem.bounds.left; ymax = dem.bounds.top
    ncols = dem.width; nrows = dem.height
    diag = (px * px + py * py) ** 0.5

    facc_arr, facc_inv = _open_raster(FACC_TIF)
    slope_arr, slope_inv = _open_raster(SLOPE_TIF)

    (geoms, dn, indeg, shreve, prop, az, topo, new_id, up_end, down_end) = \
        build_directed_by_FA(gdf, facc_arr, facc_inv, px, py, diag)

    n_outlets = sum(1 for fid in topo if dn.get(fid, -1) == -1)
    print(f"[stream_network] segments={len(topo)}  outlets={n_outlets}  "
          f"max_propagated_order={max(prop.values()) if prop else 0}")

    p_net = write_stream_network(gdf, topo, dn, new_id, prop, shreve, indeg,
                                 slope_arr, slope_inv, STREAMS_OUT)
    print(f"[stream_network] wrote {p_net.name}")

    p_map = write_stream_map(gdf, topo, new_id, az, px, py, xmin, ymax, ncols, nrows, STREAMS_OUT)
    print(f"[stream_network] wrote {p_map.name}")


if __name__ == "__main__":
    main()
