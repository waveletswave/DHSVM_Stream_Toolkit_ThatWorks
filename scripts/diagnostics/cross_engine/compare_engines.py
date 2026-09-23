# -*- coding: utf-8 -*-
# =====================================================================
# compare_engines.py  -  the Toolkit's stream network against WW-DHSVM's,
#                        on the same grid, DEM and support area
#
# Two independent engines build the channel network of a basin:
#   Toolkit   GRASS r.watershed (MFD) + r.stream.extract, then the
#             raster-native segment stage (Tier E)
#   WW-DHSVM  pyflwdir priority-flood fill + D8, then its own raster-
#             native segment builder (Zhi Li, UConn-EFC/WW-DHSVM)
# This script runs WW-DHSVM's terrain and network stages on the Toolkit's
# grid (from the Toolkit's elev_clipped.tif), with the elevations read
# from the unclipped reprojected tile and the basin mask from the clip,
# at the Toolkit's support area A_c, and compares the result with the
# Toolkit's own outputs at the level DHSVM sees: channel cells, outlets,
# segments, lengths, ranks, classes, and the outlet invariants. It also
# writes WW-DHSVM's three stream files and a Channel.State for them, so
# the network can be swapped into a DHSVM run, and rasters in the
# Toolkit's conventions so check_network_orientation.py can be run on
# the WW-DHSVM network as well.
#
# Needs the WW-DHSVM fork with terrain.demFromRaster, streams.checkOutlets
# and channel_initiation (branch feat-channel-initiation) and its
# import-time dependencies (pyflwdir, xarray, rioxarray, cftime, ...).
#
# Run (paths are examples; every input is an argument):
#   python3 compare_engines.py --case CA \
#     --toolkit-rasters standalone_CA/tests/fixtures/CA_28m \
#     --toolkit-streams standalone_CA/tests/fixtures/CA_28m/expected \
#     --toolkit-segments \
#         standalone_CA/tests/fixtures/CA_28m/expected/segments.csv \
#     --tile /path/USGS_1_n36w084_20220725_UTM17.tif \
#     --ww-repo /path/WW-DHSVM --out /path/out/CA
# =====================================================================

import argparse
import csv
import hashlib
import json
import logging
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import rasterio

# GRASS D8 codes, counter-clockwise from north-east, as the Toolkit's
# rasters carry them: (drow, dcol) -> code
GRASS_CODE = {(-1, 1): 1, (-1, 0): 2, (-1, -1): 3, (0, -1): 4,
              (1, -1): 5, (1, 0): 6, (1, 1): 7, (0, 1): 8}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_raster(path):
    with rasterio.open(path) as d:
        arr = d.read(1).astype("float64")
        nd = d.nodata
        valid = np.isfinite(arr)
        if nd is not None and np.isfinite(nd):
            valid &= arr != nd
        return arr, valid, d.transform, d.crs


def write_raster(path, arr, transform, crs, dtype, nodata):
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0],
                       width=arr.shape[1], count=1, dtype=dtype, crs=crs,
                       transform=transform, nodata=nodata) as d:
        d.write(arr.astype(dtype), 1)


def read_network(path):
    rows = []
    for line in open(path):
        p = line.split()
        if not p or p[0].startswith("#"):
            continue
        rows.append(dict(id=int(p[0]), rank=int(p[1]), slope=float(p[2]),
                         length=float(p[3]), cls=int(p[4]), down=int(p[5]),
                         save="SAVE" in line))
    return rows


def read_map_cells(path):
    cells = {}
    for line in open(path):
        p = line.split()
        if not p or p[0].startswith("#"):
            continue
        cells.setdefault(int(p[2]), []).append((int(p[1]), int(p[0])))
    return cells


def read_classes(path):
    out = {}
    for line in open(path):
        p = line.split()
        if not p or p[0].startswith("#"):
            continue
        out[int(p[0])] = dict(width=float(p[1]), depth=float(p[2]),
                              n=float(p[3]))
    return out


def network_summary(rows, cells, classes, dem, valid, acc_abs, cell_area):
    """Facts DHSVM sees, from the three stream files and the rasters."""
    ids = {r["id"] for r in rows}
    down = {r["id"]: r["down"] for r in rows}
    outlets = sorted(r["id"] for r in rows if r["down"] == 0)
    indeg = Counter(d for d in down.values() if d)
    all_cells = [c for cs in cells.values() for c in cs]
    cellset = set(all_cells)
    # tails: the last map record of each outlet segment
    tails = {sid: cells[sid][-1] for sid in outlets if sid in cells}
    r, c = np.unravel_index(
        int(np.argmax(np.where(valid, acc_abs, -np.inf))), acc_abs.shape)
    seg_of = {cell: sid for sid, cs in cells.items() for cell in cs}
    max_seg = seg_of.get((int(r), int(c)), 0)
    lowest = min(cellset, key=lambda rc: dem[rc]) if cellset else None
    low_seg = seg_of.get(lowest, 0)
    length = sum(r["length"] for r in rows)
    return dict(
        segments=len(rows), stream_cells=len(cellset),
        map_records=len(all_cells),
        cells_under_two_segments=len(all_cells) - len(cellset),
        outlets=outlets,
        outlet_tails={sid: (int(t[0]), int(t[1]), round(float(dem[t]), 2))
                      for sid, t in tails.items()},
        save_rows=sum(1 for r in rows if r["save"]),
        total_length_m=round(length, 1),
        mean_segment_length_m=round(length / max(len(rows), 1), 1),
        rank_hist=dict(sorted(Counter(r["rank"] for r in rows).items())),
        max_rank=max(r["rank"] for r in rows),
        indegree_hist=dict(sorted(Counter(indeg[s] for s in ids).items())),
        max_indegree=max(indeg.values()) if indeg else 0,
        class_hist=dict(sorted(Counter(r["cls"] for r in rows).items())),
        class_table={k: v for k, v in classes.items()
                     if k in {r["cls"] for r in rows}},
        max_acc_cell=(int(r), int(c)),
        max_acc_cells=round(float(acc_abs[r, c]), 1),
        max_acc_in_outlet=(max_seg in outlets), max_acc_segment=max_seg,
        lowest_channel_cell=((int(lowest[0]), int(lowest[1]))
                             if lowest else None),
        lowest_channel_z=round(float(dem[lowest]), 2) if lowest else None,
        lowest_in_outlet=(low_seg in outlets), lowest_segment=low_seg,
        drainage_density_km_per_km2=round(
            length / (int(valid.sum()) * cell_area) * 1000.0, 3),
    )


def main():
    ap = argparse.ArgumentParser(
        description="Toolkit versus WW-DHSVM network comparison")
    ap.add_argument("--case", required=True)
    ap.add_argument("--toolkit-rasters", required=True,
                    help="folder with elev_clipped, slope_filled, flow_acc, "
                         "stream_raster, stream_dir .tif")
    ap.add_argument("--toolkit-streams", required=True,
                    help="folder with stream.class.dat, stream.network.dat, "
                         "stream.map.dat")
    ap.add_argument("--toolkit-segments", required=True,
                    help="segments.csv of the Toolkit run")
    ap.add_argument("--tile", required=True,
                    help="the reprojected DEM tile the clip was cut from")
    ap.add_argument("--ww-repo", required=True,
                    help="WW-DHSVM checkout (feat-channel-initiation)")
    ap.add_argument("--area-m2", type=float, default=47571.5)
    ap.add_argument("--initial-depth", type=float, default=0.05,
                    help="initial channel water depth for Channel.State, m "
                         "(Toolkit states.py)")
    ap.add_argument("--state-date", default="01.01.2016.00.00.00")
    ap.add_argument("--drop-tmax", type=int, default=600)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="  %(message)s")
    sys.path.insert(0, str(Path(args.ww_repo).resolve()))
    import ww_dhsvm.terrain as T
    import ww_dhsvm.streams as S
    import ww_dhsvm.channel_initiation as CI

    out = Path(args.out).resolve()
    (out / "ww_streams").mkdir(parents=True, exist_ok=True)
    (out / "ww_rasters").mkdir(parents=True, exist_ok=True)
    tk = Path(args.toolkit_rasters).resolve()
    tks = Path(args.toolkit_streams).resolve()

    print("=" * 72)
    print(f"cross-engine comparison, {args.case}, "
          f"A_c = {args.area_m2:.1f} m2")
    print("=" * 72)
    print("inputs:")
    inputs = {}
    for name in ["elev_clipped.tif", "slope_filled.tif", "flow_acc.tif",
                 "stream_raster.tif", "stream_dir.tif"]:
        inputs[name] = sha256(tk / name)
    for name in ["stream.class.dat", "stream.network.dat", "stream.map.dat"]:
        inputs[name] = sha256(tks / name)
    inputs["segments.csv"] = sha256(args.toolkit_segments)
    inputs["tile"] = sha256(args.tile)
    for k, v in inputs.items():
        print(f"  {k:22s} {v[:16]}")

    # ---------------------------------------------------------- WW-DHSVM
    print("WW-DHSVM on the Toolkit's grid:")
    grid, _ = T.demFromRaster(str(tk / "elev_clipped.tif"))
    _, dem = T.demFromRaster(args.tile, grid=grid, mask_from_nodata=False)
    terrain = T.conditionDEM(dem.astype("float32"), grid)
    net = S.extractNetwork(terrain, grid,
                           channel_threshold_km2=args.area_m2 / 1e6)
    chk = S.checkTopology(net["segments"], net, terrain, grid)
    drop = CI.dropAnalysis(terrain, grid, tmin_cells=10,
                           tmax_cells=args.drop_tmax, step_cells=10)
    segs = net["segments"]
    files = S.writeAll(str(out / "ww_streams"), net)
    # Channel.State as the Toolkit's states.py writes it:
    # id, width x length x depth
    cls = net["class_table"].set_index("ID")
    state = out / "ww_streams" / f"Channel.State.{args.state_date}"
    with open(state, "w") as f:
        for _, s in segs.sort_values("ID").iterrows():
            width = float(cls.loc[int(s["class"]), "width"])
            vol = width * float(s["length"]) * args.initial_depth
            f.write(f"{int(s['ID'])} {vol:.6f}\n")
    # rasters in the Toolkit's conventions for check_network_orientation.py
    tr, crs = grid.transform, grid.crs.to_wkt()
    mask = grid.mask != 0
    nrows, ncols = grid.shape
    idxs = np.asarray(terrain["flowdir"].idxs_ds).reshape(grid.shape)
    fdir = np.zeros(grid.shape, dtype="int16")
    for rr in range(nrows):
        for cc in range(ncols):
            if not mask[rr, cc]:
                continue
            j = int(idxs[rr, cc])
            if j < 0 or j == rr * ncols + cc:
                continue
            fdir[rr, cc] = GRASS_CODE.get((j // ncols - rr, j % ncols - cc),
                                          0)
    wr = out / "ww_rasters"
    write_raster(wr / "elev_clipped.tif",
                 np.where(mask, terrain["raw_dem"], -9999.0), tr, crs,
                 "float32", -9999.0)
    write_raster(wr / "flow_acc.tif",
                 np.where(mask, terrain["flowacc"], -9999.0), tr, crs,
                 "float32", -9999.0)
    write_raster(wr / "stream_raster.tif",
                 net["channel_mask"].astype("uint8"), tr, crs, "uint8", 0)
    write_raster(wr / "flow_dir.tif", fdir, tr, crs, "int16", -32768)
    write_raster(wr / "segment_id.tif",
                 net["segment_id_grid"].astype("int32"), tr, crs, "int32", 0)
    write_raster(wr / "slope_filled.tif",
                 np.where(mask, np.degrees(terrain["slope"]), np.nan), tr,
                 crs, "float32", np.nan)
    drop["table"].to_csv(out / f"drop_sweep_ww_{args.case}.csv", index=False)

    # ---------------------------------------------------------- summaries
    dem_tk, valid_tk, _, _ = read_raster(tk / "elev_clipped.tif")
    acc_tk, _, _, _ = read_raster(tk / "flow_acc.tif")
    acc_tk = np.abs(np.nan_to_num(acc_tk, nan=0.0))
    sr_tk, vsr, _, _ = read_raster(tk / "stream_raster.tif")
    cells_tk = set(map(tuple, np.argwhere(vsr & (sr_tk != 0))))
    tk_sum = network_summary(read_network(tks / "stream.network.dat"),
                             read_map_cells(tks / "stream.map.dat"),
                             read_classes(tks / "stream.class.dat"),
                             dem_tk, valid_tk, acc_tk, grid.cell_area)
    ws = out / "ww_streams"
    ww_sum = network_summary(read_network(ws / "stream.network.dat"),
                             read_map_cells(ws / "stream.map.dat"),
                             read_classes(ws / "stream.class.dat"),
                             np.where(mask, terrain["raw_dem"], np.nan), mask,
                             np.where(mask, terrain["flowacc"], 0.0),
                             grid.cell_area)
    cells_ww = set(map(tuple, np.argwhere(net["channel_mask"])))
    common = cells_tk & cells_ww
    overlap = dict(
        toolkit_cells=len(cells_tk), ww_cells=len(cells_ww),
        common=len(common),
        toolkit_only=len(cells_tk - cells_ww),
        ww_only=len(cells_ww - cells_tk),
        jaccard=round(len(common) / len(cells_tk | cells_ww), 3),
        toolkit_max_acc_cell_is_ww_channel=(
            tk_sum["max_acc_cell"] in cells_ww),
        ww_max_acc_cell_is_toolkit_channel=(
            ww_sum["max_acc_cell"] in cells_tk),
        toolkit_outlet_tails_in_ww=[t[:2] in cells_ww
                                    for t in tk_sum["outlet_tails"].values()],
        # agreement of the channel sets by upstream-area band (Toolkit MFD
        # |acc|, in cells)
        common_by_acc_band={},
    )
    bands = [(0, 60), (60, 120), (120, 300), (300, 1000), (1000, 1e9)]
    for lo, hi in bands:
        tkb = {c for c in cells_tk if lo <= acc_tk[c] < hi}
        wwb = {c for c in cells_ww if lo <= acc_tk[c] < hi}
        if tkb or wwb:
            key = f"{lo}-{int(hi) if hi < 1e9 else 'max'}"
            overlap["common_by_acc_band"][key] = dict(
                toolkit=len(tkb), ww=len(wwb), common=len(tkb & wwb))
    outlet_keys = ("ok", "max_area_cell", "max_area_segment", "lowest_cell",
                   "lowest_segment", "outlet_segments", "outlet_tails",
                   "n_edge_outlets")
    ww_checks = dict(topology_ok=chk["ok"], errors=chk["errors"],
                     warnings=chk["warnings"],
                     outlets=chk["outlets"] and {k: chk["outlets"][k]
                                                 for k in outlet_keys},
                     drop_objective_cells=drop["objective_cells"],
                     drop_objective_km2=drop["objective_km2"],
                     drop_band_cells=drop["band_cells"],
                     fill_cells_raised=terrain["n_pits_filled"],
                     stream_cells_at_objective=None)
    if drop["objective_km2"] is not None:
        net_obj = S.extractNetwork(
            terrain, grid, channel_threshold_km2=drop["objective_km2"])
        ww_checks["stream_cells_at_objective"] = int(
            net_obj["channel_mask"].sum())
        ww_checks["segments_at_objective"] = len(net_obj["segments"])

    result = dict(case=args.case, area_m2=args.area_m2, grid=dict(
        nrows=grid.nrows, ncols=grid.ncols, cellsize=grid.cellsize,
        basin_cells=int(mask.sum()), basin_km2=round(grid.basin_area_km2, 4)),
        inputs=inputs, toolkit=tk_sum, ww_dhsvm=ww_sum, overlap=overlap,
        ww_checks=ww_checks,
        ww_files={k: sha256(v) for k, v in files.items()})
    (out / f"comparison_{args.case}.json").write_text(
        json.dumps(result, indent=1, default=str))

    # ---------------------------------------------------------- report
    L = []
    L.append(f"# {args.case}: Toolkit versus WW-DHSVM, "
             f"A_c = {args.area_m2:.1f} m2 "
             f"({args.area_m2 / grid.cell_area:.0f} cells)\n")
    L.append(f"Grid {grid.nrows} x {grid.ncols} at {grid.cellsize:.3f} m, "
             f"{int(mask.sum())} basin cells ({grid.basin_area_km2:.3f} "
             f"km2). Same grid, mask, DEM and support area for both "
             f"engines; no stream burning.\n")
    L.append("| | Toolkit (GRASS MFD, r.stream.extract) | "
             "WW-DHSVM (pyflwdir D8) |")
    L.append("|---|---|---|")
    labels = [("stream_cells", "channel cells"), ("segments", "segments"),
              ("map_records", "map records"),
              ("cells_under_two_segments", "cells under two segments"),
              ("total_length_m", "total length (m)"),
              ("mean_segment_length_m", "mean segment length (m)"),
              ("drainage_density_km_per_km2", "drainage density (km/km2)"),
              ("outlets", "outlet segments"),
              ("outlet_tails", "outlet tail cells (row, col, z)"),
              ("save_rows", "SAVE rows"), ("max_rank", "max routing rank"),
              ("rank_hist", "rank histogram"),
              ("max_indegree", "max in-degree"),
              ("indegree_hist", "in-degree histogram"),
              ("class_hist", "classes used"),
              ("max_acc_cell", "max-|acc| cell"),
              ("max_acc_cells", "its |acc| (cells)"),
              ("max_acc_in_outlet", "max-|acc| cell in an outlet segment"),
              ("lowest_channel_cell", "lowest channel cell"),
              ("lowest_channel_z", "its elevation (m)"),
              ("lowest_in_outlet",
               "lowest channel cell in an outlet segment")]
    for key, label in labels:
        L.append(f"| {label} | {tk_sum[key]} | {ww_sum[key]} |")
    L.append("")
    bands_txt = "; ".join(f"{k}: {v['common']} of {v['toolkit']} / {v['ww']}"
                          for k, v in overlap["common_by_acc_band"].items())
    L.append(f"Channel-cell overlap: {overlap['common']} common, "
             f"{overlap['toolkit_only']} Toolkit only, {overlap['ww_only']} "
             f"WW-DHSVM only, Jaccard {overlap['jaccard']}. By Toolkit |acc| "
             f"band (cells): {bands_txt}.")
    L.append(f"The Toolkit's max-|acc| cell {tk_sum['max_acc_cell']} is a "
             f"WW-DHSVM channel cell: "
             f"{overlap['toolkit_max_acc_cell_is_ww_channel']}; WW-DHSVM's "
             f"max-area cell {ww_sum['max_acc_cell']} is a Toolkit channel "
             f"cell: {overlap['ww_max_acc_cell_is_toolkit_channel']}.")
    L.append("")
    L.append("Class tables (id: width, depth, Manning n):")
    for name, s in (("Toolkit", tk_sum), ("WW-DHSVM", ww_sum)):
        L.append(f"- {name}: " + "; ".join(
            f"{k}: {v['width']:.2f}, {v['depth']:.3f}, {v['n']:.3f}"
            for k, v in sorted(s["class_table"].items())))
    L.append("")
    at_obj = ""
    if drop["objective_km2"]:
        at_obj = (f", giving {ww_checks['stream_cells_at_objective']} "
                  f"channel cells in {ww_checks['segments_at_objective']} "
                  f"segments")
    L.append(f"WW-DHSVM checks: topology ok {chk['ok']}, outlet invariants "
             f"ok {chk['outlets']['ok']}, errors {chk['errors']}, warnings "
             f"{chk['warnings']}; depression fill raised "
             f"{terrain['n_pits_filled']} cells; drop analysis objective "
             f"{drop['objective_cells']} cells ({drop['objective_km2']} km2), "
             f"band {drop['band_cells']}{at_obj}.")
    L.append("")
    shas = ", ".join(f"{k} {v[:16]}" for k, v in result["ww_files"].items())
    L.append(f"WW-DHSVM stream files (sha256): {shas}; Channel.State written "
             f"for {len(segs)} segments at {args.initial_depth} m initial "
             f"depth.")
    (out / f"comparison_{args.case}.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwritten: {out / f'comparison_{args.case}.md'}, .json, "
          f"ww_streams/, ww_rasters/")


if __name__ == "__main__":
    main()
