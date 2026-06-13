# -*- coding: utf-8 -*-
# iso_check_streamfile.py
# Sub-step A of the #6 vector-IO rebuild: characterize the standalone
# streamfile.shp (GRASS r.to.vect) against the qgis_CA reference streamfile.shp
# (QGIS r.to.vect), since the whole directed-network + text-file chain is built
# on this layer's geometry and attributes.
#
# Four questions:
#   1. feature count (known 41 vs 41, reconfirm)
#   2. dbf field list on each side (what attributes r.to.vect emits)
#   3. geometry: match lines across the two layers and compare vertex
#      coordinates -- are the lines bit/near identical, and is vertex ORDER
#      the same (endpoint up/down assignment depends on it)
#   4. feature ordering / id: does record i on one side correspond to record i
#      on the other (new_id renumbering depends on iteration order)
#
# Verdict this informs:
#   - geometry identical + same order  -> downstream text files can be byte-identical
#   - geometry equivalent, order/ids differ -> only topological equivalence;
#     stream.network/map/class will differ in segment IDs though the network matches
#
# Run:  python3 iso_check_streamfile.py

import sys
from pathlib import Path
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import REF_INTERMED, STREAMFILE

# Reference (QGIS) and standalone stream vector, from paths.py.
REF = str(REF_INTERMED / "streamfile.shp")
STD = str(STREAMFILE)

ref = gpd.read_file(REF)
std = gpd.read_file(STD)

print("=== 1) feature count ===")
print(f"  ref={len(ref)}  std={len(std)}")

print("\n=== 2) dbf fields ===")
print(f"  ref cols: {list(ref.columns)}")
print(f"  std cols: {list(std.columns)}")

print("\n=== 3+4) geometry + ordering ===")
# Build endpoint signatures for each line: (first_xy, last_xy), rounded.
def endpoints(geom):
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "MultiLineString":
        # take the longest part, mirroring _line_coords
        parts = list(geom.geoms)
        geom = max(parts, key=lambda g: g.length) if parts else None
        if geom is None:
            return None
    cs = list(geom.coords)
    return cs[0], cs[-1], len(cs)

def round_xy(xy, nd=3):
    return (round(xy[0], nd), round(xy[1], nd))

# match each ref line to the nearest std line by midpoint, then compare vertices
def midpoint(geom):
    return geom.interpolate(0.5, normalized=True)

ref_mid = [midpoint(g) for g in ref.geometry]
std_mid = [midpoint(g) for g in std.geometry]

# nearest std line for each ref line (by midpoint distance)
from shapely.strtree import STRtree
std_tree = STRtree(list(std.geometry))

matched = 0
order_same = 0
vtx_count_same = 0
max_vertex_diff = 0.0
endpoint_order_same = 0
unmatched = []

for i, rg in enumerate(ref.geometry):
    rmid = ref_mid[i]
    # query nearest std geometry
    nearest_idx = std_tree.nearest(rmid)
    sg = std.geometry.iloc[nearest_idx]

    rep = endpoints(rg); sep = endpoints(sg)
    if rep is None or sep is None:
        unmatched.append(i); continue

    (r0, r1, rn) = rep; (s0, s1, sn) = sep
    # vertex count
    if rn == sn:
        vtx_count_same += 1
    # endpoint order: do first->first and last->last align (same direction)?
    same_dir = (round_xy(r0) == round_xy(s0)) and (round_xy(r1) == round_xy(s1))
    rev_dir  = (round_xy(r0) == round_xy(s1)) and (round_xy(r1) == round_xy(s0))
    if same_dir:
        endpoint_order_same += 1
    if same_dir or rev_dir:
        matched += 1
        # max coordinate diff over the shared endpoints (in matched orientation)
        if same_dir:
            d = max(abs(r0[0]-s0[0]), abs(r0[1]-s0[1]), abs(r1[0]-s1[0]), abs(r1[1]-s1[1]))
        else:
            d = max(abs(r0[0]-s1[0]), abs(r0[1]-s1[1]), abs(r1[0]-s0[0]), abs(r1[1]-s0[1]))
        max_vertex_diff = max(max_vertex_diff, d)
    else:
        unmatched.append(i)
    # ordering: is the nearest std index the same as ref index?
    if nearest_idx == i:
        order_same += 1

print(f"  lines matched by endpoints (either direction): {matched} / {len(ref)}")
print(f"  of matched: same vertex direction (first->first): {endpoint_order_same}")
print(f"  vertex COUNT identical: {vtx_count_same} / {len(ref)}")
print(f"  max endpoint coord diff (matched lines): {max_vertex_diff:.3e}")
print(f"  record-order aligned (ref[i] nearest == std[i]): {order_same} / {len(ref)}")
if unmatched:
    print(f"  UNMATCHED ref line indices: {unmatched[:20]}")

print("\n=== interpretation hint ===")
if matched == len(ref) and max_vertex_diff < 1e-3 and order_same == len(ref) and endpoint_order_same == len(ref):
    print("  geometry identical AND same order/direction -> text files may be byte-identical")
elif matched == len(ref) and max_vertex_diff < 1.0:
    print("  geometry equivalent; check order/direction columns above -> likely topological equivalence only")
else:
    print("  geometry differs non-trivially -> investigate before porting the chain")
