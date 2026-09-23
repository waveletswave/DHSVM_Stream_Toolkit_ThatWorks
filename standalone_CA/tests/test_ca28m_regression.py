# -*- coding: utf-8 -*-
# =====================================================================
# test_ca28m_regression.py  -  the Tier E network stage on the real
#                              CA 28 m rasters, against the audited files
#
# Runs segments_from_raster, channelclass_standalone and network_files
# on the five GRASS-stage rasters in fixtures/CA_28m (elev_clipped,
# slope_filled, flow_acc, stream_raster, stream_dir; provenance in
# fixtures/CA_28m/README.md) and compares the results with the files the
# same stages wrote on DCC on 2026-09-22 and the DHSVM reruns used:
#
#   segments.csv, stream_cells.csv   value by value
#   stream.class.dat                 byte identical, sha256 26983e53..
#   stream.map.dat                   byte identical, sha256 a4f0f02e..
#   stream.network.dat               byte identical, sha256 6a01bd2f..
#
# plus the network facts recorded in the audit
# (docs/audit/tier_e_network_orientation_2026_09_22.md). No GRASS is
# needed; rasterio and geopandas are.
# =====================================================================

import csv
import hashlib
import importlib
import math
import os
import shutil
from collections import Counter
from pathlib import Path

import pytest

rasterio = pytest.importorskip("rasterio")
pytest.importorskip("geopandas")

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "fixtures" / "CA_28m"
EXPECTED = FIXTURE / "expected"
RASTERS = ["elev_clipped.tif", "slope_filled.tif", "flow_acc.tif",
           "stream_raster.tif", "stream_dir.tif"]
RASTER_SHA = {
    "elev_clipped.tif":
    "181539c95da679c91fb9aea2dbba7a04db5ca18b81d3ff632a3f8e4114fb7d4e",
    "slope_filled.tif":
    "c6c55348a735ead9a00f35b049e79c4191086c4cd20f217087417bb77a3a7321",
    "flow_acc.tif":
    "f0126f0d2e24a90c8bb4be24e010ca4fbc076c50b3ceff06f7cf30a1213a3c65",
    "stream_raster.tif":
    "6a6d112c238ebaa661174a5dd42f6d01108df940893bf0ed08083365f17ee9a0",
    "stream_dir.tif":
    "7e1c7acbcafb26219da3a03c6b0bbc787c3a1fcadfee99197da7810e8ec40584",
}
# the stream files of DEM_CA_tierE, the inputs of the audit's nA_ and
# new_ DHSVM reruns
STREAM_SHA = {
    "stream.class.dat":
    "26983e5342993a033cba391d720c81ce664c773cd3f0c554e16059e1999f4751",
    "stream.map.dat":
    "a4f0f02ea78d8ab5c6c28d4d38c152f420bc3b772aff442ab2ec55b6d98dd442",
    "stream.network.dat":
    "6a01bd2f6c1d64049ca85c920dcbaed09aa2a37ef380b9791b6af3b36a474595",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def same_value(a, b):
    """Integers exactly, floats to 1e-9 relative, other strings exactly."""
    try:
        return int(a) == int(b)
    except ValueError:
        pass
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-9)
    except ValueError:
        return a == b


def assert_same_table(got_path, exp_path):
    got, exp = read_csv(got_path), read_csv(exp_path)
    assert len(got) == len(exp), f"{got_path.name}: {len(got)} rows, " \
                                 f"expected {len(exp)}"
    assert got[0].keys() == exp[0].keys(), f"{got_path.name}: columns differ"
    for i, (g, e) in enumerate(zip(got, exp)):
        for k in e:
            assert same_value(g[k], e[k]), \
                f"{got_path.name} row {i} column {k}: {g[k]} != {e[k]}"


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """Run the three network stages on the fixture; return (paths, out)."""
    out = tmp_path_factory.mktemp("ca28m")
    for name in RASTERS:
        shutil.copy2(FIXTURE / name, out / name)
    os.environ["DHSVM_OUT"] = str(out)
    os.environ["DHSVM_EPSG"] = "32617"
    import paths
    importlib.reload(paths)
    assert Path(paths.OUT) == out
    import segments_from_raster as sfr
    sfr.run()
    import channelclass_standalone as cc
    importlib.reload(cc)
    cc.channelclassfun(paths.STREAMFILE_ATTR, paths.STREAMS_DIR)
    import network_files as nf
    nf.run()
    return paths, out


def test_fixture_rasters_are_the_audited_ones():
    for name, want in RASTER_SHA.items():
        assert sha256(FIXTURE / name) == want, f"{name} is not the fixture"
    with rasterio.open(FIXTURE / "elev_clipped.tif") as d:
        assert d.shape == (74, 82)
        assert d.crs.to_epsg() == 32617
        assert abs(d.res[0] - 28.15774) < 1e-4


def test_segments_table(run):
    paths, out = run
    assert_same_table(paths.SEGMENTS_CSV, EXPECTED / "segments.csv")


def test_stream_cells_table(run):
    paths, out = run
    assert_same_table(paths.STREAM_CELLS, EXPECTED / "stream_cells.csv")


def test_stream_files_byte_identical(run):
    paths, out = run
    for name, want in STREAM_SHA.items():
        got = paths.STREAMS_DIR / name
        assert got.read_bytes() == (EXPECTED / name).read_bytes(), \
            f"{name} differs from the audited file"
        assert sha256(got) == want


def test_network_facts_from_the_audit(run):
    paths, out = run
    rows = [ln.split() for ln in open(paths.STREAMS_DIR / "stream.network.dat")
            if ln.strip() and not ln.startswith("#")]
    assert len(rows) == 22
    ids = [int(r[0]) for r in rows]
    assert ids == list(range(1, 23))
    outlets = [int(r[0]) for r in rows if r[5] == "0"]
    assert outlets == [12, 22]
    assert [r[6:8] for r in rows if r[5] == "0"] == \
        [["SAVE", '"outlet_1"'], ["SAVE", '"outlet_2"']]
    assert all(len(r) == 6 for r in rows if r[5] != "0")
    ranks = Counter(int(r[1]) for r in rows)
    assert ranks == {1: 12, 2: 2, 3: 2, 4: 2, 5: 1, 6: 1, 7: 1, 8: 1}
    classes = Counter(int(r[4]) for r in rows)
    assert classes == {13: 17, 14: 5}
    indeg = Counter(int(r[5]) for r in rows if r[5] != "0")
    assert max(indeg.values()) == 2 and len(indeg) == 10
    assert all(float(r[2]) > 0 and float(r[3]) > 0 for r in rows)
    assert abs(sum(float(r[3]) for r in rows) - 7600.8) < 0.1

    recs = [ln.split() for ln in open(paths.STREAMS_DIR / "stream.map.dat")
            if ln.strip() and not ln.startswith("#")]
    assert len(recs) == 231
    cells = [(int(r[1]), int(r[0])) for r in recs]
    assert len(set(cells)) == 231, "a cell listed under two segments"

    segs = {int(r["segid"]): r for r in read_csv(paths.SEGMENTS_CSV)}
    assert (int(segs[12]["tail_row"]), int(segs[12]["tail_col"])) == (67, 68)
    assert abs(float(segs[12]["z_tail"]) - 829.61) < 0.01
    assert (int(segs[22]["tail_row"]), int(segs[22]["tail_col"])) == (65, 68)
    assert abs(float(segs[22]["z_tail"]) - 838.00) < 0.01
    assert int(segs[22]["rank"]) == 8 and int(segs[12]["rank"]) == 1
