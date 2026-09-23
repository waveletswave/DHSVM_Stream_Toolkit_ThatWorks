# -*- coding: utf-8 -*-
# =====================================================================
# test_drop_analysis.py  -  the constant-drop threshold sweep recomputes
#                           the Strahler order for every threshold
#
# pyflwdir's FlwdirRaster.stream_order (0.5.5 to 0.5.12) caches the
# Strahler map on the object and ignores the mask on later calls, which
# silently gave every threshold of the sweep the orders of the first one
# (docs/audit/drop_analysis_strahler_cache_2026_09_23.md). These tests
# pin the corrected behaviour on the CA 28 m fixture: the maximum order
# falls as the threshold rises, and the objective is 120 cells.
# Needs pyflwdir and scipy (the [dev] extras).
# =====================================================================

import importlib.util
from pathlib import Path

import numpy as np
import pytest

pyflwdir = pytest.importorskip("pyflwdir")
pytest.importorskip("scipy")

HERE = Path(__file__).resolve().parent
DEM = HERE / "fixtures" / "CA_28m" / "elev_clipped.tif"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "drop_analysis", HERE.parent / "diagnostics" / "drop_analysis.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def flow():
    da = load_module()
    elev, nodata, transform, cell_area = da.load_dem(DEM)
    flw = da.build_flow(elev, nodata, transform)
    acc = np.asarray(flw.upstream_area(unit="cell"))
    return da, flw, elev, transform, acc, cell_area


def test_strahler_order_follows_the_mask(flow):
    da, flw, elev, transform, acc, cell_area = flow
    dense = da.strahler_order(flw, acc > 10)
    sparse = da.strahler_order(flw, acc > 300)
    assert int(dense.max()) == 4
    assert int(sparse.max()) == 2
    # the order map is zero off the streams
    assert int(sparse[~(acc > 300)].max()) == 0
    # recomputing in the other order gives the same maps
    assert np.array_equal(da.strahler_order(flw, acc > 300), sparse)
    assert np.array_equal(da.strahler_order(flw, acc > 10), dense)


def test_max_order_falls_along_the_sweep(flow):
    da, flw, elev, transform, acc, cell_area = flow
    basin_cells = int(np.isfinite(elev).sum() - (elev == -9999.0).sum())
    orders = []
    for thr in (10, 60, 120, 300):
        r = da.evaluate_threshold(flw, elev, transform, acc, cell_area,
                                  thr, basin_cells)
        orders.append(r["max_order"])
    assert orders == [4, 3, 3, 2]
    assert orders == sorted(orders, reverse=True)


def test_ca28m_objective_is_120_cells(flow):
    da, flw, elev, transform, acc, cell_area = flow
    basin_cells = int(((elev != -9999.0) & np.isfinite(elev)).sum())
    assert basin_cells == 4334
    results = [da.evaluate_threshold(flw, elev, transform, acc, cell_area,
                                     thr, basin_cells)
               for thr in range(10, 301, 10)]
    results = [r for r in results if r is not None]
    obj, band = da.first_sustained_band(results, 3)
    assert obj["cells"] == 120
    assert (band[0]["cells"], band[-1]["cells"]) == (120, 300)
    assert not any(r["passes"] for r in results if r["cells"] < 120)
    by = {r["cells"]: r for r in results}
    assert by[60]["t_abs"] > 2.0 and by[120]["t_abs"] < 2.0
    assert by[120]["n_order1"] == 9 and by[120]["n_higher"] == 8
