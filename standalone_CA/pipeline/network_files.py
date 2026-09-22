# -*- coding: utf-8 -*-
# =====================================================================
# network_files.py  -  stream.network.dat and stream.map.dat writers
#
# Tier E replacement for stream_network.py. Reads the segment table that
# segments_from_raster.py wrote and channelclass_standalone.py extended
# (chanclass, hydwidth, hyddepth written back), plus the per-cell table,
# and writes the two DHSVM routing files. Topology comes from the table
# (segid, downid, rank); nothing is re-derived from geometry.
#
# stream.network.dat  ID  order  slope  length  class  down  [SAVE "name"]
#   order  = propagated routing rank, dense (channel_route_network stops
#            at the first empty rank)
#   down   = 0 on an outlet segment (channel_read_network looks up any
#            non-zero id), followed by SAVE "outlet_<k>" so DHSVM writes
#            the segment's hydrograph to Streamflow.Only (the Lawler test
#            case form:  0  SAVE  "OUTLET")
# stream.map.dat      col  row  ID  length  height  width  aspect
#   one record per stream cell; height and width from the class table
#   through the write-back columns; aspect = azimuth of the cell's D8
#   step, degrees clockwise from north
#
# Checks before writing (fail loudly): every class id exists in
# stream.class.dat; slopes and lengths positive; ranks dense; every
# downid names a segment; the number of outlet rows equals the number
# of outlet cells in the direction raster; the lowest stream cell is the
# tail of an outlet segment (reported, not fatal); in-degree histogram.
#
# Run:  python3 network_files.py
# =====================================================================

import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def read_segments(attr_path):
    import geopandas as gpd
    gdf = gpd.read_file(str(attr_path))
    need = ["segid", "downid", "rank", "Shape_Leng", "slope_tan",
            "chanclass", "hydwidth", "hyddepth", "is_outlet"]
    missing = [c for c in need if c not in gdf.columns]
    if missing:
        raise RuntimeError(f"{Path(attr_path).name} lacks columns "
                           f"{missing}; run segments_from_raster.py and "
                           "channelclass_standalone.py first")
    segs = {}
    for _, row in gdf.iterrows():
        sid = int(row["segid"])
        segs[sid] = dict(down=int(row["downid"]), rank=int(row["rank"]),
                         length=float(row["Shape_Leng"]),
                         slope=float(row["slope_tan"]),
                         cls=int(row["chanclass"]),
                         width=float(row["hydwidth"]),
                         height=float(row["hyddepth"]),
                         is_outlet=int(row["is_outlet"]))
    return segs


def read_cells(cells_csv):
    cells = defaultdict(list)
    with open(cells_csv, newline="") as f:
        for rec in csv.DictReader(f):
            cells[int(rec["segid"])].append(dict(
                seq=int(rec["seq"]), row=int(rec["row"]),
                col=int(rec["col"]), length=float(rec["length_m"]),
                aspect=float(rec["azimuth_deg"])))
    for sid in cells:
        cells[sid].sort(key=lambda d: d["seq"])
    return cells


def read_class_ids(class_path):
    ids = set()
    with open(class_path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.add(int(s.split()[0]))
    return ids


def check(segs, cells, class_ids):
    ids = set(segs)
    for sid, s in segs.items():
        if s["down"] != 0 and s["down"] not in ids:
            raise RuntimeError(f"segment {sid} drains to unknown {s['down']}")
        if (s["down"] == 0) != bool(s["is_outlet"]):
            raise RuntimeError(f"segment {sid}: down/is_outlet disagree")
        if s["slope"] <= 0 or s["length"] <= 0:
            raise RuntimeError(f"segment {sid}: slope {s['slope']} "
                               f"length {s['length']} must be positive")
        if s["cls"] not in class_ids:
            raise RuntimeError(f"segment {sid}: class {s['cls']} not in "
                               "stream.class.dat")
        if sid not in cells:
            raise RuntimeError(f"segment {sid} has no cell records")
    ranks = sorted(set(s["rank"] for s in segs.values()))
    if ranks != list(range(1, max(ranks) + 1)):
        raise RuntimeError(f"ranks not dense: {ranks}")
    for sid, s in segs.items():
        if s["down"] and segs[s["down"]]["rank"] <= s["rank"]:
            raise RuntimeError(f"rank not increasing {sid} -> {s['down']}")
    seen = Counter((c["row"], c["col"]) for v in cells.values() for c in v)
    dup = [k for k, n in seen.items() if n > 1]
    if dup:
        raise RuntimeError(f"{len(dup)} cells listed under two segments")
    indeg = Counter(s["down"] for s in segs.values() if s["down"])
    return indeg


def write_network(segs, out_path):
    n_out = 0
    with open(out_path, "w") as f:
        for sid in sorted(segs):
            s = segs[sid]
            line = (f"{sid:d} {s['rank']:d} {s['slope']:0.5f} "
                    f"{s['length']:0.5f} {s['cls']:d} {s['down']:d}")
            if s["down"] == 0:
                n_out += 1
                line += f' SAVE "outlet_{n_out}"'
            f.write(line + "\n")
    return n_out


def write_map(segs, cells, out_path):
    n = 0
    with open(out_path, "w") as f:
        f.write("###### This file has been automatically generated #####\n")
        f.write("#  Col  Row  ID     Length  Height    Width    Aspect\n")
        for sid in sorted(segs):
            s = segs[sid]
            for c in cells[sid]:
                aspect = int(round(c["aspect"])) % 360
                f.write(f"{c['col']:5d}{c['row']:6d}{sid:6d}"
                        f"{c['length']:11.4f}{s['height']:11.4f}"
                        f"{s['width']:10.4f}{aspect:11d}\n")
                n += 1
    return n


def lowest_stream_cell(dem_path, stream_path):
    import rasterio
    with rasterio.open(dem_path) as d, rasterio.open(stream_path) as s:
        dem = d.read(1).astype(np.float64)
        sr = s.read(1).astype(np.float64)
        vdem = np.isfinite(dem)
        if d.nodata is not None and np.isfinite(d.nodata):
            vdem &= dem != d.nodata
        vsr = np.isfinite(sr) & (sr != 0)
        if s.nodata is not None and np.isfinite(s.nodata):
            vsr &= sr != s.nodata
    stream = vdem & vsr
    z = np.where(stream, dem, np.inf)
    r, c = np.unravel_index(int(np.argmin(z)), z.shape)
    return int(r), int(c), float(dem[r, c]), int(stream.sum())


def run():
    from paths import (STREAMFILE_ATTR, STREAM_CELLS, STREAMS_DIR,
                       ELEV_CLIPPED, STREAM_RASTER)
    STREAMS_DIR.mkdir(parents=True, exist_ok=True)
    class_path = STREAMS_DIR / "stream.class.dat"
    for tag, p in [("segments", STREAMFILE_ATTR), ("cells", STREAM_CELLS),
                   ("classes", class_path)]:
        if not p.exists():
            raise FileNotFoundError(f"[error] missing {tag}: {p}")
    segs = read_segments(STREAMFILE_ATTR)
    cells = read_cells(STREAM_CELLS)
    indeg = check(segs, cells, read_class_ids(class_path))

    net_path = STREAMS_DIR / "stream.network.dat"
    map_path = STREAMS_DIR / "stream.map.dat"
    n_out = write_network(segs, net_path)
    n_map = write_map(segs, cells, map_path)

    r, c, z, n_stream = lowest_stream_cell(ELEV_CLIPPED, STREAM_RASTER)
    tails = {}
    for sid, s in segs.items():
        last = cells[sid][-1]
        tails[(last["row"], last["col"])] = sid
    owner = tails.get((r, c))
    ok_low = owner is not None and segs[owner]["down"] == 0
    hist = Counter(indeg[sid] for sid in segs)
    print(f"[network] segments {len(segs)}, outlet rows {n_out}, "
          f"map records {n_map} (stream cells {n_stream}), "
          f"max rank {max(s['rank'] for s in segs.values())}")
    print(f"[network] in-degree histogram {dict(sorted(hist.items()))}")
    where = (f"outlet segment {owner}" if ok_low
             else "NO outlet segment <-- check the mouth")
    print(f"[network] lowest stream cell (row {r}, col {c}) z {z:.2f} m "
          f"is the tail of {where}")
    if n_map != n_stream:
        raise RuntimeError(f"map records {n_map} != stream cells "
                           f"{n_stream}")
    print(f"[network] wrote {net_path.name}, {map_path.name}")


if __name__ == "__main__":
    run()
