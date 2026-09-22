# -*- coding: utf-8 -*-
# =====================================================================
# test_segments_from_raster.py  -  Tier E network stage, synthetic tests
#
# Runs without GRASS. The core tests need numpy only; the end-to-end
# test needs rasterio and geopandas and writes to a temporary directory.
#
# Run:  python3 tests/test_segments_from_raster.py      (from standalone_CA)
#   or  python3 -m pytest tests/test_segments_from_raster.py
# =====================================================================

import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pipeline"))

import segments_from_raster as sfr  # noqa: E402

CODE = {(-1, 1): 1, (-1, 0): 2, (-1, -1): 3, (0, -1): 4,
        (1, -1): 5, (1, 0): 6, (1, 1): 7, (0, 1): 8}
N = 10
PX = 10.0


def chain(a, b):
    """Cells from a to b along a straight D8 path (inclusive)."""
    (r0, c0), (r1, c1) = a, b
    n = max(abs(r1 - r0), abs(c1 - c0))
    return [(r0 + round(k * (r1 - r0) / n), c0 + round(k * (c1 - c0) / n))
            for k in range(n + 1)]


def y_network():
    """Five segments, one outlet at the south edge.

    H1 (0,0)->(3,3), H2 (0,6)->(3,3), M1 (3,3)->(6,3), H3 (3,6)->(6,3),
    M2 (6,3)->(9,3) draining south off the grid.
    """
    stream = np.zeros((N, N), dtype=bool)
    fdir = np.zeros((N, N), dtype=np.int32)
    chains = [chain((0, 0), (3, 3)), chain((0, 6), (3, 3)),
              chain((3, 3), (6, 3)), chain((3, 6), (6, 3)),
              chain((6, 3), (9, 3))]
    for cells in chains:
        for rc in cells:
            stream[rc] = True
        for u, v in zip(cells[:-1], cells[1:]):
            fdir[u] = CODE[(v[0] - u[0], v[1] - u[1])]
    fdir[9, 3] = -6            # GRASS marks a cell draining out as negative
    return stream, fdir


def test_y_network_topology():
    stream, fdir = y_network()
    segments, seg_of = sfr.build_segments(stream, fdir)
    assert len(segments) == 5
    assert int((seg_of >= 0).sum()) == int(stream.sum())
    heads = {tuple(s["cells"][0]) for s in segments}
    assert heads == {(0, 0), (0, 6), (3, 3), (3, 6), (6, 3)}
    by_head = {tuple(s["cells"][0]): s for s in segments}
    assert by_head[(0, 0)]["succ"] == (3, 3)
    assert by_head[(0, 6)]["succ"] == (3, 3)
    assert by_head[(3, 3)]["succ"] == (6, 3)
    assert by_head[(3, 6)]["succ"] == (6, 3)
    assert by_head[(6, 3)]["succ"] is None
    assert by_head[(6, 3)]["down"] == -1
    assert segments[by_head[(0, 0)]["down"]]["cells"][0] == (3, 3)
    rank, segid, ups = sfr.rank_and_number(segments)
    ranks = {tuple(s["cells"][0]): rank[i] for i, s in enumerate(segments)}
    assert ranks == {(0, 0): 1, (0, 6): 1, (3, 3): 2, (3, 6): 1, (6, 3): 3}
    assert sorted(segid) == [1, 2, 3, 4, 5]
    assert segid[segments.index(by_head[(6, 3)])] == 5   # outlet is last


def test_y_network_attributes():
    stream, fdir = y_network()
    segments, _ = sfr.build_segments(stream, fdir)
    rows, cols = np.mgrid[0:N, 0:N]
    dem = 200.0 - 5.0 * rows + 0.5 * np.abs(cols - 3)
    slope = np.full((N, N), 10.0)
    slope[0, 0] = np.nan                       # a missing slope cell
    acc = -(1.0 + rows * 3.0)                  # negative, as r.watershed
    acc[:, 3] -= 20.0
    segs, cells = sfr.segment_attributes(segments, stream, fdir, dem,
                                         slope, acc, PX, PX)
    assert len(cells) == int(stream.sum())
    by_head = {(s["head_row"], s["head_col"]): s for s in segs}
    h1 = by_head[(0, 0)]
    assert abs(h1["length_m"] - 3 * math.sqrt(2) * PX) < 1e-9
    m1 = by_head[(3, 3)]
    assert abs(m1["length_m"] - 3 * PX) < 1e-9
    m2 = by_head[(6, 3)]
    assert abs(m2["length_m"] - 4 * PX) < 1e-9      # 3 steps + outlet cell
    assert m2["is_outlet"] == 1 and m2["downid"] == 0
    assert m2["rank"] == 3 and m2["indeg"] == 2
    # meanmsq: max |acc| on the segment times the cell area, sign ignored
    assert abs(m2["meanmsq_m2"] - (1 + 9 * 3 + 20) * PX * PX) < 1e-9
    # slope: mean tan over cells with a valid slope
    assert abs(h1["slope_tan"] - math.tan(math.radians(10.0))) < 1e-12
    assert abs(h1["slope_deg"] - 10.0) < 1e-12
    # azimuths: H1 steps south-east, M2 south
    az = {(c["row"], c["col"]): c["azimuth_deg"] for c in cells}
    assert az[(0, 0)] == 135.0 and az[(6, 3)] == 180.0
    assert az[(9, 3)] == 180.0                  # outlet cell keeps its code
    assert all(s["slope_tan"] > 0 for s in segs)


def test_two_outlets_and_negative_codes():
    stream, fdir = y_network()
    extra = chain((7, 7), (9, 9))               # a second stream, exits SE
    for rc in extra:
        stream[rc] = True
    for u, v in zip(extra[:-1], extra[1:]):
        fdir[u] = CODE[(v[0] - u[0], v[1] - u[1])]
    fdir[9, 9] = -7
    segments, _ = sfr.build_segments(stream, fdir)
    assert len(segments) == 6
    assert sum(s["down"] == -1 for s in segments) == 2
    rank, segid, _ = sfr.rank_and_number(segments)
    assert sorted(set(rank)) == [1, 2, 3]


def test_missing_direction_is_an_outlet():
    stream, fdir = y_network()
    fdir[9, 3] = 0                              # no direction at the end
    segments, _ = sfr.build_segments(stream, fdir)
    assert sum(s["down"] == -1 for s in segments) == 1


def test_cycle_is_rejected():
    stream = np.zeros((N, N), dtype=bool)
    fdir = np.zeros((N, N), dtype=np.int32)
    stream[4, 4] = stream[4, 5] = True
    fdir[4, 4] = 8                              # east
    fdir[4, 5] = 4                              # west: a two-cell loop
    try:
        sfr.build_segments(stream, fdir)
    except RuntimeError:
        return
    raise AssertionError("a cycle must raise")


def test_end_to_end_files():
    try:
        import rasterio
        from rasterio.transform import from_origin
        import geopandas  # noqa: F401
    except ImportError:
        print("  (rasterio/geopandas not available; end-to-end skipped)")
        return
    stream, fdir = y_network()
    rows, cols = np.mgrid[0:N, 0:N]
    dem = (200.0 - 5.0 * rows + 0.5 * np.abs(cols - 3)).astype("float32")
    slope = np.full((N, N), 10.0, dtype="float32")
    acc = -(1.0 + rows * 3.0)
    acc[:, 3] -= 20.0
    T = from_origin(500000.0, 4000000.0, PX, PX)
    tmp = Path(tempfile.mkdtemp(prefix="tierE_"))
    os.environ["DHSVM_OUT"] = str(tmp)
    os.environ["DHSVM_EPSG"] = "32617"

    def write(name, arr, dtype, nodata):
        with rasterio.open(tmp / name, "w", driver="GTiff", height=N,
                           width=N, count=1, dtype=dtype, crs="EPSG:32617",
                           transform=T, nodata=nodata) as d:
            d.write(arr.astype(dtype), 1)
    write("elev_clipped.tif", dem, "float32", -9999.0)
    write("slope_filled.tif", slope, "float32", None)
    write("flow_acc.tif", acc, "float64", None)
    write("stream_raster.tif", stream.astype("int32"), "int32", 0)
    write("stream_dir.tif", fdir, "int32", None)

    import importlib
    import paths
    importlib.reload(paths)
    assert Path(paths.OUT) == tmp
    sfr.run()
    import channelclass_standalone as cc
    importlib.reload(cc)
    cc.channelclassfun(paths.STREAMFILE_ATTR, paths.STREAMS_DIR)
    import network_files as nf
    nf.run()

    net = [ln.split() for ln in open(paths.STREAMS_DIR / "stream.network.dat")
           if ln.strip() and not ln.startswith("#")]
    assert len(net) == 5
    outlet_rows = [ln for ln in net if ln[5] == "0"]
    assert len(outlet_rows) == 1 and outlet_rows[0][6] == "SAVE"
    assert outlet_rows[0][7] == '"outlet_1"'
    assert all(len(ln) == 6 for ln in net if ln[5] != "0")
    ranks = sorted(int(ln[1]) for ln in net)
    assert ranks == [1, 1, 1, 2, 3]
    ids = {int(ln[0]) for ln in net}
    assert ids == {1, 2, 3, 4, 5}
    assert all(int(ln[5]) in ids | {0} for ln in net)
    recs = [ln.split() for ln in open(paths.STREAMS_DIR / "stream.map.dat")
            if ln.strip() and not ln.startswith("#")]
    assert len(recs) == int(stream.sum())
    cells = {(int(r[1]), int(r[0])) for r in recs}
    assert cells == set(zip(*np.nonzero(stream)))
    assert all(len(r) == 7 for r in recs)
    print(f"  end-to-end files written under {tmp}")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("OK")
