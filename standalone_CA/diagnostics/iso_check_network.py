# -*- coding: utf-8 -*-
# iso_check_network.py
# Prove the standalone stream.network.dat / stream.map.dat are the SAME directed
# network as the reference under a segment-ID relabeling. The toposort tie-break
# differs, so new_id numbering is permuted; topology should be identical.
#
# Method:
#   1. parse both network.dat. Each row: segid nin slope_tan length class down.
#   2. build a bijection ref_id <-> std_id by matching on (round(length,5),
#      round(slope_tan,5), class). If unique, we have the relabeling.
#   3. verify: for every ref segment A with down=B, the mapped std segment
#      map[A] has down == map[B]. If all hold -> isomorphic directed network.
#   4. also check nin (order) agrees under the mapping.
#   5. for map.dat: group rows by segid into a set of (col,row,len,aspect),
#      and check ref segment A's cell-set equals std map[A]'s cell-set.
#
# Run:  python3 iso_check_network.py

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
from paths import REF_STREAMS, STREAMS_DIR

# Reference (QGIS) and standalone stream text-file directories, from paths.py.
REF = str(REF_STREAMS)
STD = str(STREAMS_DIR)


def parse_network(path):
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            segid = int(p[0]); nin = int(p[1])
            slope = round(float(p[2]), 5); length = round(float(p[3]), 5)
            cls = int(p[4]); down = int(p[5])
            rows[segid] = dict(nin=nin, slope=slope, length=length, cls=cls, down=down)
    return rows


def main():
    ref = parse_network(f"{REF}/stream.network.dat")
    std = parse_network(f"{STD}/stream.network.dat")
    print(f"segments: ref={len(ref)} std={len(std)}")

    # signature -> id, on each side
    def sig(r):
        return (r["length"], r["slope"], r["cls"])

    ref_by_sig = {}
    for sid, r in ref.items():
        ref_by_sig.setdefault(sig(r), []).append(sid)
    std_by_sig = {}
    for sid, r in std.items():
        std_by_sig.setdefault(sig(r), []).append(sid)

    # build bijection where signature is unique on both sides
    ref2std = {}
    ambiguous = []
    for s, ids in ref_by_sig.items():
        sids = std_by_sig.get(s, [])
        if len(ids) == 1 and len(sids) == 1:
            ref2std[ids[0]] = sids[0]
        else:
            ambiguous.append((s, ids, sids))

    print(f"unique-signature segments mapped: {len(ref2std)} / {len(ref)}")
    if ambiguous:
        print(f"ambiguous (non-unique signature) groups: {len(ambiguous)}")
        for s, ids, sids in ambiguous[:10]:
            print(f"   sig={s}  ref_ids={ids}  std_ids={sids}")

    # verify down-pointer consistency under the mapping (outlet down=0 maps to 0)
    ref2std_full = dict(ref2std)
    # try to resolve ambiguous by down-consistency later; first report on unique part
    ok_down = 0
    bad_down = []
    checkable = 0
    for rid, r in ref.items():
        if rid not in ref2std_full:
            continue
        rdown = r["down"]
        if rdown == 0:
            mapped_down_expected = 0
        elif rdown in ref2std_full:
            mapped_down_expected = ref2std_full[rdown]
        else:
            continue  # down target ambiguous, skip
        checkable += 1
        actual_down = std[ref2std_full[rid]]["down"]
        if actual_down == mapped_down_expected:
            ok_down += 1
        else:
            bad_down.append((rid, ref2std_full[rid], rdown, mapped_down_expected, actual_down))

    print(f"\ndown-pointer consistency (unique-mapped, resolvable targets):")
    print(f"  checked={checkable}  consistent={ok_down}  inconsistent={len(bad_down)}")
    if bad_down:
        for b in bad_down[:10]:
            print(f"   ref{b[0]}->std{b[1]}: ref.down={b[2]} expect std.down={b[3]} got {b[4]}")

    # order (nin) agreement under mapping
    ok_nin = 0; bad_nin = 0
    for rid, sid in ref2std_full.items():
        if ref[rid]["nin"] == std[sid]["nin"]:
            ok_nin += 1
        else:
            bad_nin += 1
    print(f"\nnin/order agreement (unique-mapped): ok={ok_nin} bad={bad_nin}")

    # ---- map.dat cell-set check under the same mapping ----
    def parse_map(path):
        seg_cells = {}
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                col = int(p[0]); row = int(p[1]); segid = int(p[2])
                length = round(float(p[3]), 3); aspect = int(p[6])
                seg_cells.setdefault(segid, set()).add((col, row, length, aspect))
        return seg_cells

    rmap = parse_map(f"{REF}/stream.map.dat")
    smap = parse_map(f"{STD}/stream.map.dat")
    print(f"\nmap.dat: ref segs={len(rmap)} std segs={len(smap)}")
    same_cells = 0; diff_cells = 0; skipped = 0
    for rid, sid in ref2std_full.items():
        if rid in rmap and sid in smap:
            if rmap[rid] == smap[sid]:
                same_cells += 1
            else:
                diff_cells += 1
        else:
            skipped += 1
    print(f"  per-segment cell-set identical under mapping: {same_cells}")
    print(f"  per-segment cell-set DIFFERS: {diff_cells}   skipped: {skipped}")

    print("\n=== verdict ===")
    if len(ref2std) == len(ref) and len(bad_down) == 0 and bad_nin == 0 and diff_cells == 0:
        print("  ISOMORPHIC: same directed network + same order + same cells, "
              "only segment-ID numbering permuted (toposort tie-break).")
    elif len(bad_down) == 0 and len(ambiguous) > 0:
        print("  unique part consistent; ambiguous groups need manual look (see above).")
    else:
        print("  NOT clean -- inspect down/nin/cell mismatches above.")


if __name__ == "__main__":
    main()
