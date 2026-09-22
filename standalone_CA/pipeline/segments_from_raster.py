# -*- coding: utf-8 -*-
# =====================================================================
# segments_from_raster.py  -  channel segments from the stream raster
#
# Tier E replacement for vector_attrs.py. The stream network is built
# from the two rasters r.stream.extract wrote, stream_raster.tif (stream
# cells) and stream_dir.tif (the D8 direction it followed), instead of
# from r.to.vect lines re-linked by geometry. Orientation and topology
# are then intrinsic: every stream cell has one successor, a segment
# starts at a head (no stream inflow) or just below a confluence (two or
# more stream inflows) and ends at the next start or at an outlet cell.
# There is no tie-break, no tolerance and no sink merge (validation_log,
# Tier E).
#
# Conventions kept from the previous stage:
#   slope_deg  mean of the filled slope raster (degrees) over the cells
#   slope_tan  mean of tan(slope) over the same cells, default 0.01
#   meanmsq    contributing area in m2: max |flow_acc| on the segment
#              times the cell area (r.watershed writes negative values
#              for cells that may receive flow from outside the region;
#              only the magnitude is the accumulation)
#   Row, Col   the middle cell of the segment, top-left origin
# Per cell: length = the D8 step to the successor (cell size, or cell
# size times sqrt 2 on a diagonal); an outlet cell takes the cell size;
# aspect = azimuth of that step, degrees clockwise from north.
#
# Outputs (paths.py): STREAMFILE_ATTR (segment polylines through cell
# centres, upstream to downstream, with the attributes channelclass
# needs), SEGMENTS_CSV, STREAM_CELLS (one record per stream cell).
#
# Run:  python3 segments_from_raster.py
# =====================================================================

import csv
import math
from collections import Counter, deque

import numpy as np

# GRASS D8 codes: 1..8 counter-clockwise from north-east; a negative
# code is the same direction for a cell that drains out of the region.
D8 = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1),
      5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}
AZIMUTH = {1: 45.0, 2: 0.0, 3: 315.0, 4: 270.0,
           5: 225.0, 6: 180.0, 7: 135.0, 8: 90.0}
SLOPE_TAN_DEFAULT = 0.01


# ----------------------------- core (pure numpy) ----------------------
def successors(stream, fdir):
    """Successor of every stream cell.

    Returns (succ, code): succ is an int array of flat indices, -1 where
    the successor is off-grid, not a stream cell, or the code is not a
    D8 code (an outlet cell); code is the absolute D8 code (0 if none).
    """
    nr, nc = stream.shape
    succ = np.full(stream.shape, -1, dtype=np.int64)
    code = np.zeros(stream.shape, dtype=np.int32)
    for r, c in zip(*np.nonzero(stream)):
        d = fdir[r, c]
        if not np.isfinite(d):
            continue
        k = int(abs(int(d)))
        if k not in D8:
            continue
        code[r, c] = k
        dr, dc = D8[k]
        r2, c2 = r + dr, c + dc
        if 0 <= r2 < nr and 0 <= c2 < nc and stream[r2, c2]:
            succ[r, c] = r2 * nc + c2
    return succ, code


def build_segments(stream, fdir):
    """Split the stream cells into segments.

    Returns (segments, seg_of). segments is a list of dicts with keys
    cells (list of (row, col), upstream to downstream), down (index into
    segments, or -1 for an outlet segment) and succ (the (row, col) the
    last cell drains to, or None). seg_of maps every stream cell to its
    segment index (-1 elsewhere). Raises if a stream cell is left
    unassigned, which would mean the direction raster has a cycle.
    """
    nr, nc = stream.shape
    succ, _ = successors(stream, fdir)
    inflow = np.zeros(stream.shape, dtype=np.int32)
    for r, c in zip(*np.nonzero(stream)):
        s = succ[r, c]
        if s >= 0:
            inflow[s // nc, s % nc] += 1
    is_start = stream & ((inflow == 0) | (inflow >= 2))
    seg_of = np.full(stream.shape, -1, dtype=np.int64)
    segments = []
    for r0, c0 in zip(*np.nonzero(is_start)):      # row-major, so stable
        cells = []
        r, c = int(r0), int(c0)
        while True:
            cells.append((r, c))
            seg_of[r, c] = len(segments)
            s = succ[r, c]
            if s < 0:
                nxt = None
                break
            r2, c2 = int(s // nc), int(s % nc)
            if is_start[r2, c2]:
                nxt = (r2, c2)
                break
            r, c = r2, c2
        segments.append(dict(cells=cells, succ=nxt, down=-1))
    n_unassigned = int((stream & (seg_of < 0)).sum())
    if n_unassigned:
        raise RuntimeError(f"{n_unassigned} stream cells not reachable "
                           "from any segment start; the direction "
                           "raster is not a tree")
    for seg in segments:
        if seg["succ"] is not None:
            seg["down"] = int(seg_of[seg["succ"]])
    return segments, seg_of


def rank_and_number(segments):
    """Propagated routing rank (1 + max over upstream) and a dense
    topological numbering, upstream first. Raises on a cycle."""
    n = len(segments)
    ups = [[] for _ in range(n)]
    indeg = [0] * n
    for i, seg in enumerate(segments):
        d = seg["down"]
        if d >= 0:
            ups[d].append(i)
            indeg[d] += 1
    rank = [0] * n
    queue = deque(i for i in range(n) if indeg[i] == 0)
    for i in queue:
        rank[i] = 1
    order = []
    while queue:
        i = queue.popleft()
        order.append(i)
        d = segments[i]["down"]
        if d >= 0:
            rank[d] = max(rank[d], rank[i] + 1)
            indeg[d] -= 1
            if indeg[d] == 0:
                queue.append(d)
    if len(order) != n:
        raise RuntimeError("cycle in the segment graph")
    segid = [0] * n
    for k, i in enumerate(order):
        segid[i] = k + 1
    ranks = sorted(set(rank))
    if ranks != list(range(1, max(rank) + 1)):
        raise RuntimeError(f"routing ranks not dense: {ranks}")
    return rank, segid, ups


def step_geometry(code_k, px, py):
    """Length and azimuth of the D8 step with code k (0: no step)."""
    if code_k in D8:
        dr, dc = D8[code_k]
        return math.hypot(dc * px, dr * py), AZIMUTH[code_k]
    return px, 0.0


def segment_attributes(segments, stream, fdir, dem, slope_deg, acc,
                       px, py):
    """Per-cell and per-segment attributes (see the header)."""
    succ, code = successors(stream, fdir)
    cell_area = px * py
    rank, segid, ups = rank_and_number(segments)
    cells_out = []
    segs_out = []
    for i, seg in enumerate(segments):
        cells = seg["cells"]
        length = 0.0
        tans, degs, accs = [], [], []
        for k, (r, c) in enumerate(cells):
            code_k = int(code[r, c])
            if succ[r, c] >= 0:
                step, az = step_geometry(code_k, px, py)
            else:
                step, az = px, (AZIMUTH.get(code_k, 0.0))
            length += step
            s = float(slope_deg[r, c])
            cell_tan = float("nan")
            if np.isfinite(s) and s >= 0.0:
                cell_tan = math.tan(math.radians(s))
                degs.append(s)
                tans.append(cell_tan)
            a = float(acc[r, c])
            cell_acc = abs(a) if np.isfinite(a) else float("nan")
            if np.isfinite(a):
                accs.append(cell_acc)
            cells_out.append(dict(segid=segid[i], seq=k, row=r, col=c,
                                  length_m=step, azimuth_deg=az,
                                  slope_tan=cell_tan, abs_acc=cell_acc))
        slope_tan = (sum(tans) / len(tans)) if tans else SLOPE_TAN_DEFAULT
        if slope_tan <= 0.0:
            slope_tan = SLOPE_TAN_DEFAULT
        sdeg = (sum(degs) / len(degs)) if degs else 0.0
        meanmsq = (max(accs) * cell_area) if accs else cell_area
        mid = cells[len(cells) // 2]
        head, tail = cells[0], cells[-1]
        d = seg["down"]
        segs_out.append(dict(
            segid=segid[i], downid=(segid[d] if d >= 0 else 0),
            rank=rank[i], ncells=len(cells), length_m=length,
            slope_deg=sdeg, slope_tan=slope_tan, meanmsq_m2=meanmsq,
            row=mid[0], col=mid[1], is_outlet=int(d < 0),
            head_row=head[0], head_col=head[1],
            tail_row=tail[0], tail_col=tail[1],
            z_head=float(dem[head]), z_tail=float(dem[tail]),
            indeg=len(ups[i]),
            succ=seg["succ"], tail_code=int(code[tail])))
    return segs_out, cells_out


# ----------------------------- I/O ------------------------------------
def read_raster(path):
    import rasterio
    with rasterio.open(path) as ds:
        arr = ds.read(1).astype(np.float64)
        nd = ds.nodata
        valid = np.isfinite(arr)
        if nd is not None and np.isfinite(nd):
            valid &= arr != nd
        return arr, valid, ds.transform, ds.crs


def polyline(seg, transform, px, py):
    """Cell-centre polyline, upstream to downstream. A non-outlet
    segment ends at the successor cell centre (the junction); a
    one-cell outlet segment is extended half a cell along its D8 step
    so the line has two distinct vertices."""
    from shapely.geometry import LineString

    def centre(r, c):
        x, y = transform * (c + 0.5, r + 0.5)
        return (float(x), float(y))
    pts = [centre(*rc) for rc in seg["cells"]]
    if seg["succ"] is not None:
        pts.append(centre(*seg["succ"]))
    if len(pts) < 2:
        k = seg["tail_code"]
        dr, dc = D8.get(k, (1, 0))
        x, y = pts[0]
        pts.append((x + 0.5 * dc * px, y - 0.5 * dr * py))
    return LineString(pts)


def write_outputs(segs, cells, segments, transform, crs, px, py,
                  attr_path, segments_csv, cells_csv, epsg):
    import geopandas as gpd
    from rasterio.crs import CRS
    order = sorted(range(len(segs)), key=lambda i: segs[i]["segid"])
    rows = []
    geoms = []
    for i in order:
        s = segs[i]
        rows.append(dict(arcid=s["segid"], Shape_Leng=s["length_m"],
                         Row=s["row"], Col=s["col"],
                         slope_deg=s["slope_deg"], slope_tan=s["slope_tan"],
                         meanmsq=s["meanmsq_m2"], segid=s["segid"],
                         downid=s["downid"], rank=s["rank"],
                         ncells=s["ncells"], is_outlet=s["is_outlet"]))
        geoms.append(polyline(dict(segments[i], tail_code=s["tail_code"]),
                              transform, px, py))
    gdf = gpd.GeoDataFrame(rows, geometry=geoms,
                           crs=(crs if crs is not None
                                else CRS.from_epsg(epsg)))
    gdf.to_file(str(attr_path))
    keys = ["segid", "downid", "rank", "ncells", "length_m", "slope_deg",
            "slope_tan", "meanmsq_m2", "row", "col", "is_outlet",
            "head_row", "head_col", "tail_row", "tail_col", "z_head",
            "z_tail", "indeg"]
    with open(segments_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for i in order:
            w.writerow(segs[i])
    ckeys = ["segid", "seq", "row", "col", "length_m", "azimuth_deg",
             "slope_tan", "abs_acc"]
    with open(cells_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ckeys)
        w.writeheader()
        for rec in sorted(cells, key=lambda d: (d["segid"], d["seq"])):
            w.writerow(rec)


def summarize(segs, stream):
    n_cells = int(stream.sum())
    n_outlets = sum(s["is_outlet"] for s in segs)
    indeg = Counter(s["indeg"] for s in segs)
    ranks = Counter(s["rank"] for s in segs)
    total = sum(s["length_m"] for s in segs)
    print(f"[segments] stream cells {n_cells}, segments {len(segs)}, "
          f"outlet segments {n_outlets}, total length {total:.1f} m")
    print(f"[segments] rank histogram {dict(sorted(ranks.items()))}")
    print(f"[segments] in-degree histogram {dict(sorted(indeg.items()))}")
    if max(indeg) > 3:
        print(f"[segments] WARNING in-degree above 3 on "
              f"{sum(v for k, v in indeg.items() if k > 3)} segment(s)")
    for s in segs:
        if s["is_outlet"]:
            print(f"[segments] outlet segment {s['segid']}: tail cell "
                  f"(row {s['tail_row']}, col {s['tail_col']}) "
                  f"z {s['z_tail']:.2f} m, {s['ncells']} cells, "
                  f"rank {s['rank']}, meanmsq {s['meanmsq_m2']:.0f} m2")


def run():
    from paths import (ELEV_CLIPPED, SLOPE_FILLED, FLOW_ACC, STREAM_RASTER,
                       STREAM_DIR, STREAMFILE_ATTR, STREAM_CELLS,
                       SEGMENTS_CSV, EPSG)
    for tag, p in [("dem", ELEV_CLIPPED), ("slope", SLOPE_FILLED),
                   ("flow_acc", FLOW_ACC), ("stream_raster", STREAM_RASTER),
                   ("stream_dir", STREAM_DIR)]:
        if not p.exists():
            raise FileNotFoundError(f"[error] missing {tag}: {p}")
    dem, vdem, transform, crs = read_raster(ELEV_CLIPPED)
    slope, _, _, _ = read_raster(SLOPE_FILLED)
    acc, _, _, _ = read_raster(FLOW_ACC)
    srast, vsr, _, _ = read_raster(STREAM_RASTER)
    fdir, vfd, _, _ = read_raster(STREAM_DIR)
    for name, arr in (("slope", slope), ("flow_acc", acc),
                      ("stream_raster", srast), ("stream_dir", fdir)):
        if arr.shape != dem.shape:
            raise RuntimeError(f"{name} shape {arr.shape} != DEM "
                               f"{dem.shape}")
    stream = vsr & (srast != 0) & vdem
    fdir = np.where(vfd, fdir, 0.0)
    px, py = abs(transform.a), abs(transform.e)
    segments, _ = build_segments(stream, fdir)
    segs, cells = segment_attributes(segments, stream, fdir, dem, slope,
                                     acc, px, py)
    assert len(cells) == int(stream.sum()), "cell records != stream cells"
    write_outputs(segs, cells, segments, transform, crs, px, py,
                  STREAMFILE_ATTR, SEGMENTS_CSV, STREAM_CELLS, EPSG)
    summarize(segs, stream)
    print(f"[segments] wrote {STREAMFILE_ATTR.name}, {SEGMENTS_CSV.name}, "
          f"{STREAM_CELLS.name}")


if __name__ == "__main__":
    run()
