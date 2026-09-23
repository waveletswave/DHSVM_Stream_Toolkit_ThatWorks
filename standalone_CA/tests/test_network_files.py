# -*- coding: utf-8 -*-
# =====================================================================
# test_network_files.py  -  the invariant check and the writers of the
#                           DHSVM stream files, on synthetic networks
#
# network_files.check() is the gate between the segment tables and
# stream.network.dat / stream.map.dat: every failure mode it guards
# against is exercised here, and the two writers are checked for the
# format DHSVM reads (SAVE on every outlet row, six fields elsewhere,
# fixed-width map records, one record per stream cell). numpy only,
# except the lowest-stream-cell test, which writes two tiny rasters.
# =====================================================================

import copy
from collections import Counter

import numpy as np
import pytest

import network_files as nf


def seg(down, rank, is_outlet=0, length=100.0, slope=0.1, cls=13,
        width=1.0, height=0.1):
    return dict(down=down, rank=rank, length=length, slope=slope, cls=cls,
                width=width, height=height, is_outlet=is_outlet)


def cell(seq, row, col, length=10.0, aspect=180.0):
    return dict(seq=seq, row=row, col=col, length=length, aspect=aspect)


def y_network():
    """Two heads (1, 2) joining segment 3, which drains out."""
    segs = {1: seg(down=3, rank=1), 2: seg(down=3, rank=1),
            3: seg(down=0, rank=2, is_outlet=1, cls=14)}
    cells = {1: [cell(0, 0, 0, aspect=135.0), cell(1, 1, 1, aspect=135.0)],
             2: [cell(0, 0, 4, aspect=225.0), cell(1, 1, 3, aspect=225.0)],
             3: [cell(0, 2, 2), cell(1, 3, 2), cell(2, 4, 2, aspect=359.6)]}
    return segs, cells, {13, 14}


def test_check_passes_and_reports_in_degree():
    segs, cells, ids = y_network()
    indeg = nf.check(segs, cells, ids)
    assert indeg == Counter({3: 2})


@pytest.mark.parametrize("edit, message", [
    (lambda s, c: s[1].update(down=9), "unknown"),
    (lambda s, c: s[3].update(is_outlet=0), "down/is_outlet"),
    (lambda s, c: s[1].update(down=0), "down/is_outlet"),
    (lambda s, c: s[2].update(slope=0.0), "must be positive"),
    (lambda s, c: s[2].update(length=-1.0), "must be positive"),
    (lambda s, c: s[1].update(cls=99), "not in"),
    (lambda s, c: c.pop(2), "no cell records"),
    (lambda s, c: s[3].update(rank=3), "not dense"),
    (lambda s, c: s[1].update(rank=2), "increasing"),
    (lambda s, c: c[2].append(cell(2, 1, 1)), "two segments"),
])
def test_check_rejects(edit, message):
    segs, cells, ids = y_network()
    edit(segs, cells)
    with pytest.raises(RuntimeError, match=message):
        nf.check(segs, cells, ids)


def test_write_network_format(tmp_path):
    segs, cells, ids = y_network()
    out = tmp_path / "stream.network.dat"
    assert nf.write_network(segs, out) == 1
    rows = [ln.rstrip("\n") for ln in open(out)]
    assert rows == ["1 1 0.10000 100.00000 13 3",
                    "2 1 0.10000 100.00000 13 3",
                    '3 2 0.10000 100.00000 14 0 SAVE "outlet_1"']


def test_write_network_numbers_outlets_in_id_order(tmp_path):
    segs = {1: seg(down=0, rank=1, is_outlet=1),
            2: seg(down=4, rank=1), 4: seg(down=0, rank=2, is_outlet=1),
            3: seg(down=4, rank=1)}
    out = tmp_path / "stream.network.dat"
    assert nf.write_network(segs, out) == 2
    rows = [ln.split() for ln in open(out)]
    assert [r[0] for r in rows] == ["1", "2", "3", "4"]
    assert rows[0][6:] == ["SAVE", '"outlet_1"']
    assert rows[3][6:] == ["SAVE", '"outlet_2"']
    assert all(len(r) == 6 for r in rows[1:3])


def test_write_map_format(tmp_path):
    segs, cells, ids = y_network()
    out = tmp_path / "stream.map.dat"
    assert nf.write_map(segs, cells, out) == 7
    lines = open(out).read().splitlines()
    assert lines[0].startswith("######")
    assert lines[1].split() == ["#", "Col", "Row", "ID", "Length", "Height",
                                "Width", "Aspect"]
    recs = [ln for ln in lines[2:] if ln.strip()]
    assert len(recs) == 7
    # fixed-width columns: 5 + 6 + 6 + 11 + 11 + 10 + 11 characters
    assert all(len(ln) == 60 for ln in recs)
    first = recs[0]
    assert (int(first[0:5]), int(first[5:11]), int(first[11:17])) == (0, 0, 1)
    assert float(first[17:28]) == 10.0
    assert float(first[28:39]) == 0.1 and float(first[39:49]) == 1.0
    assert int(first[49:60]) == 135
    # the aspect is rounded to whole degrees modulo 360
    last = recs[-1]
    assert int(last[49:60]) == 0
    # one record per cell, segments in id order, cells in sequence
    order = [(int(ln[11:17]), int(ln[5:11]), int(ln[0:5])) for ln in recs]
    assert order == [(1, 0, 0), (1, 1, 1), (2, 0, 4), (2, 1, 3),
                     (3, 2, 2), (3, 3, 2), (3, 4, 2)]


def test_read_class_ids_and_cells(tmp_path):
    cls = tmp_path / "stream.class.dat"
    cls.write_text("#ID W  D   n    inf\n13   0.5 0.100 0.0450 0.0\n"
                   "14   1.0 0.100 0.0450 0.0\n")
    assert nf.read_class_ids(cls) == {13, 14}
    csvp = tmp_path / "stream_cells.csv"
    csvp.write_text("segid,seq,row,col,length_m,azimuth_deg,slope_tan,abs_acc\n"
                    "1,1,5,6,28.2,180.0,0.1,3.0\n"
                    "1,0,4,6,39.8,135.0,0.2,2.0\n"
                    "2,0,0,0,28.2,90.0,0.3,1.0\n")
    cells = nf.read_cells(csvp)
    assert [c["seq"] for c in cells[1]] == [0, 1]
    assert cells[1][0]["row"] == 4 and cells[1][1]["length"] == 28.2
    assert cells[2][0]["aspect"] == 90.0


def test_lowest_stream_cell(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin
    dem = np.array([[300.0, 290.0, 280.0],
                    [260.0, 250.0, 240.0],
                    [230.0, 210.0, -9999.0]], dtype="float32")
    stream = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]], dtype="uint8")
    kw = dict(driver="GTiff", height=3, width=3, count=1, crs="EPSG:32617",
              transform=from_origin(0.0, 30.0, 10.0, 10.0))
    with rasterio.open(tmp_path / "dem.tif", "w", dtype="float32",
                       nodata=-9999.0, **kw) as d:
        d.write(dem, 1)
    with rasterio.open(tmp_path / "stream.tif", "w", dtype="uint8",
                       nodata=0, **kw) as d:
        d.write(stream, 1)
    r, c, z, n = nf.lowest_stream_cell(tmp_path / "dem.tif",
                                       tmp_path / "stream.tif")
    assert (r, c, z, n) == (2, 1, 210.0, 3)


def test_check_does_not_modify_inputs():
    segs, cells, ids = y_network()
    before = copy.deepcopy((segs, cells))
    nf.check(segs, cells, ids)
    assert (segs, cells) == before
