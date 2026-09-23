# -*- coding: utf-8 -*-
# =====================================================================
# tierE_rerun_AR.py  -  rerun the manuscript AR configuration (AR S4h)
#                       with the April inputs, with only the network
#                       swapped, and with the current pipeline outputs
#
# The template AR_0416_S4h_UA.dhs points at DEM_AR_0406, the April 2026
# input set the manuscript AR run used (unchanged since; the June rewrite
# touched CA only). Three configurations are written from it:
#
#   control  input/oA_S4h.dhs   one line changed, Output Directory ->
#            TestCase/AR/output/oA_S4h; must reproduce the manuscript
#            output S4h* byte for byte
#   Tier E   input/nA_S4h.dhs   five lines changed: the three stream files
#            and the Initial State Directory -> DEM_AR_tierE (the pipeline
#            rerun on the April grid; its grid states are byte-identical
#            to the April ones, its Channel.State belongs to the new
#            network), Output Directory -> output/nA_S4h
#   new      input/new_S4h.dhs  ten lines changed: every input file and
#            the state directory -> DEM_AR_tierE (the current pipeline
#            outputs: conditioned soil depth and Tier E network),
#            Output Directory -> output/new_S4h
#
# DHSVM appends the file names to the Output Directory value; the value
# must stay at or under 78 characters (100-byte file name buffers in
# InitDump.c and RouteSubSurface.c), so the prefixes are short. --smoke
# runs nine days under oAsm_S4h, nAsm_S4h, newsm_S4h.
#
# Usage:
#   python3 tierE_rerun_AR.py --smoke
#   python3 tierE_rerun_AR.py
#   python3 tierE_rerun_AR.py --kinds tierE new
# =====================================================================

import argparse
import difflib
import hashlib
import os
import re
import subprocess
import sys
import time
from pathlib import Path

DHSVM_EXE = Path(os.environ.get(
    "TIERE_DHSVM_EXE",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/DHSVM/sourcecode/"
    "DHSVM3.2"))
CASE = Path(os.environ.get(
    "TIERE_CASE",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/AR"))
OLD_ROOT = Path(os.environ.get(
    "TIERE_OLD_ROOT",
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_AR_0406"))
NEW_ROOT = Path(os.environ.get(
    "TIERE_NEW_ROOT",
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_AR_tierE"))
# the WW-DHSVM network built on the same grid, DEM and A_c (cross-engine
# comparison, scripts/diagnostics/cross_engine/compare_engines.py,
# 2026-09-23): three stream files, its Channel.State, and copies of the
# grid states
WW_ROOT = Path(os.environ.get(
    "TIERE_WW_ROOT",
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/AR/DEM_AR_ww"))

RUNS = {"AR_S4h": "AR_0416_S4h_UA.dhs"}
TAGS = {"AR_S4h": "S4h"}
KINDS = ["ctrl", "tierE", "new", "ww"]
WORDS = {"ctrl": "oA", "tierE": "nA", "new": "new", "ww": "wA"}
MAX_OUTPUT_PATH = 78        # see the header


def prefix_for(kind, run, smoke):
    """Output prefix (the last path element of Output Directory), at
    most 9 characters: oA_S4h, nA_S4h, new_S4h, oAsm_S4h, ..."""
    return f"{WORDS[kind]}{'sm' if smoke else ''}_{TAGS[run]}"


STREAM_KEYS = (("Stream Map File", "DHSVM_input_streams/stream.map.dat"),
               ("Stream Network File",
                "DHSVM_input_streams/stream.network.dat"),
               ("Stream Class File",
                "DHSVM_input_streams/stream.class.dat"))
BINARY_KEYS = ("DEM File", "Basin Mask File", "Soil Map File",
               "Soil Depth File", "Vegetation Map File")

# sha256 of the April AR input set (DEM_AR_0406, files dated 2026-04-06;
# the manuscript AR run points at it)
EXPECTED_SHA_OLD = {
    "DHSVM_input_streams/stream.class.dat":
        "fd4ad1df567e612541c7f8e392e4501c83a6800bd9103c49bbe18822a550bca3",
    "DHSVM_input_streams/stream.map.dat":
        "84ba7709014b154db700c5be97d697efaf95d35766258a9bef7607352a0e5d78",
    "DHSVM_input_streams/stream.network.dat":
        "f5626c43379ac01cf0f43c51072d4d66a1b4ddfa8dd244a5f98f0770685e183c",
    "modelstate/Channel.State.01.01.2016.00.00.00":
        "583bbc8dc913febf4f9d763d2033ddd5e5bea06f12c6bd00d0b8eb3d13a3674a",
    "modelstate/Interception.State.01.01.2016.00.00.00.bin":
        "433daa6185c417d3bac1ff5b40f82c3bc325bf33cbfa0d0105abf3a0db6a625b",
    "modelstate/Snow.State.01.01.2016.00.00.00.bin":
        "ff8c02764582ef48a4ddbbc80332d6b57b666ec6ca9bffc58fc7fa182daf6708",
    "modelstate/Soil.State.01.01.2016.00.00.00.bin":
        "55b0bc94441fe8244f30d65f7bf337eecef831a7a94208b0f8251b1dc6128345",
    "DHSVM_input_binaries/dem.bin":
        "25dfbd411e34e16fb78ea91a179c6b083a2784a62e3c41a8db56c5e8d34656b5",
    "DHSVM_input_binaries/mask.bin":
        "6c4165528f417979c8d39c1ff77946f1300995a31e1294f078a9b9dbdbb0a5ba",
    "DHSVM_input_binaries/soil.bin":
        "6c4165528f417979c8d39c1ff77946f1300995a31e1294f078a9b9dbdbb0a5ba",
    "DHSVM_input_binaries/veg.bin":
        "6c4165528f417979c8d39c1ff77946f1300995a31e1294f078a9b9dbdbb0a5ba",
    "DHSVM_input_binaries/soildepth.bin":
        "aa84ca0a0d7bd8c67c4eedaacac1102a95c5d26f97baef28f85f495cf8918a30",
}
# sha256 of the pipeline rerun on the April AR grid (DCC
# /work/ys451/dhsvm_ca/tierE/fixed_AR_28m, main at eecd961, 2026-09-23):
# dem, mask, soil, veg and the grid states equal the April files;
# soildepth.bin, the three stream files and Channel.State are new
EXPECTED_SHA_NEW = dict(EXPECTED_SHA_OLD)
EXPECTED_SHA_NEW.update({
    "DHSVM_input_streams/stream.class.dat":
        "26983e5342993a033cba391d720c81ce664c773cd3f0c554e16059e1999f4751",
    "DHSVM_input_streams/stream.map.dat":
        "94cf497545aa1fad81c4d4c6d6986f39da942b9983ee78163715985b04445e78",
    "DHSVM_input_streams/stream.network.dat":
        "9d5fdca8f1222bbf976e644b9581bd69621f5a2fe21dc26df1a3cb68fc9e0dc0",
    "modelstate/Channel.State.01.01.2016.00.00.00":
        "b7276040a2d2b0dd2af32524702cabc226ff10e35fc274cd4648b407361aac99",
    "DHSVM_input_binaries/soildepth.bin":
        "262ca509decd739aef37f46a2ed07297d3dbdf9afa172ee9618810c48c922a34",
})
# sha256 of the WW-DHSVM network on the AR grid at A_c 47571.5 m2
# (compare_engines.py, WW-DHSVM fork feat-channel-initiation, 2026-09-23)
EXPECTED_SHA_WW = {
    "DHSVM_input_streams/stream.class.dat":
        "696ebf8beab15040fc16d8d53fb50c57c6ac2ef8a1ddc99a187ba3c7fa5a26aa",
    "DHSVM_input_streams/stream.map.dat":
        "0413270a3e800c4a35983a6a01117efbebe8ed6011136afecea3a8afd761dfc7",
    "DHSVM_input_streams/stream.network.dat":
        "3ce3cac9365a968225f3e2bcc6602a0fce3d66a8c21b17c396de96fa00d8c82c",
    "modelstate/Channel.State.01.01.2016.00.00.00":
        "b6ec8f314c5fd517c0c6221ce29dfe01035944ad106c70e504871a59d3e6afb3",
}
GRID_STATES = ["Interception.State.01.01.2016.00.00.00.bin",
               "Snow.State.01.01.2016.00.00.00.bin",
               "Soil.State.01.01.2016.00.00.00.bin"]
SMOKE_END = "01/10/2016-00"
OUTPUT_FILES = ["Aggregated.Values", "Mass.Balance", "Mass.Final.Balance",
                "Streamflow.Only", "Stream.Flow"]
FINAL_KEYS = ("Initial Storage", "Precip/Inflow", "ET ", "ChannelInt",
              "Final Storage", "Mass Error")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_tree(root, expected, label, sha):
    print(f"[inputs] {label} under {root}")
    bad = []
    for rel, want in expected.items():
        p = root / rel
        assert p.exists(), f"missing {p}"
        got = sha256(p)
        ok = got == want
        if not ok:
            bad.append(rel)
        print(f"  {'ok ' if ok else 'DIFFERS'} {rel}  {got[:12]}"
              f"{'' if ok else '  expected ' + want[:12]}")
    if sha:
        assert not bad, f"{label}: sha256 differs for {bad}"


def check_inputs(kinds, sha=True):
    check_tree(OLD_ROOT, EXPECTED_SHA_OLD,
               "DEM_AR_0406 (April inputs, the manuscript AR run)", sha)
    if "tierE" in kinds or "new" in kinds:
        check_tree(NEW_ROOT, EXPECTED_SHA_NEW,
                   "DEM_AR_tierE (pipeline rerun, DCC fixed_AR_28m)", sha)
        for name in GRID_STATES:
            assert sha256(OLD_ROOT / "modelstate" / name) == \
                sha256(NEW_ROOT / "modelstate" / name), \
                f"{name} differs between DEM_AR_0406 and DEM_AR_tierE"
            print(f"  ok  modelstate/{name} identical in both trees")
    if "ww" in kinds:
        check_tree(WW_ROOT, EXPECTED_SHA_WW,
                   "DEM_AR_ww (WW-DHSVM network, compare_engines.py)", sha)
        for name in GRID_STATES:
            new = WW_ROOT / "modelstate" / name
            assert new.exists(), \
                f"missing {new}: copy it from {OLD_ROOT / 'modelstate'}"
            assert sha256(OLD_ROOT / "modelstate" / name) == sha256(new), \
                f"{name} differs between DEM_AR_0406 and DEM_AR_ww"
            print(f"  ok  DEM_AR_ww/modelstate/{name} identical to "
                  "DEM_AR_0406")
    assert DHSVM_EXE.exists(), f"DHSVM executable not found: {DHSVM_EXE}"


def set_line(text, key, value, tag):
    """Replace the value of `key = ...` on its line; exactly one match."""
    pat = re.compile(rf"^({re.escape(key)}\s*=\s*)(\S+)(.*)$",
                     re.MULTILINE)
    hits = pat.findall(text)
    assert len(hits) == 1, \
        f"{tag}: expected one line for {key!r}, got {len(hits)}"
    new_text = pat.sub(lambda m: m.group(1) + value + m.group(3), text)
    return new_text, hits[0][1]


def get_value(text, key, tag):
    m = re.search(rf"^{re.escape(key)}\s*=\s*(\S+)", text, re.MULTILINE)
    assert m, f"{tag}: no line for {key!r}"
    return m.group(1)


def make_config(run, template, kind, smoke):
    text = template.read_text()
    tag = template.name
    prefix = prefix_for(kind, run, smoke)
    changed = {}
    network_keys = [k for k, _ in STREAM_KEYS] + ["Initial State Directory"]
    for key in network_keys + list(BINARY_KEYS):
        assert get_value(text, key, tag).startswith(str(OLD_ROOT)), \
            f"{tag}: {key} does not point at {OLD_ROOT}"
    n_expected = 1
    if kind == "new":
        for key in BINARY_KEYS:
            rel = get_value(text, key, tag)[len(str(OLD_ROOT)) + 1:]
            text, changed[key] = set_line(text, key, str(NEW_ROOT / rel),
                                          tag)
        n_expected += 5
    if kind in ("tierE", "new", "ww"):
        net_root = WW_ROOT if kind == "ww" else NEW_ROOT
        for key, rel in STREAM_KEYS:
            text, changed[key] = set_line(text, key, str(net_root / rel),
                                          tag)
        text, changed["Initial State Directory"] = set_line(
            text, "Initial State Directory",
            str(net_root / "modelstate") + "/", tag)
        n_expected += 4
    out_prefix = CASE / "output" / prefix
    assert len(str(out_prefix)) <= MAX_OUTPUT_PATH, (
        f"Output Directory {out_prefix} is {len(str(out_prefix))} "
        f"characters; this DHSVM build accepts at most {MAX_OUTPUT_PATH} "
        f"(100-byte file name buffers, see the header)")
    text, changed["Output Directory"] = set_line(
        text, "Output Directory", str(out_prefix), tag)
    if smoke:
        text, changed["Model End"] = set_line(text, "Model End", SMOKE_END,
                                              tag)
        n_expected += 1
    assert get_value(text, "Model Start", tag) == "01/01/2016-00", \
        f"{tag}: Model Start is not 01/01/2016-00"
    diff = [ln for ln in difflib.unified_diff(
        template.read_text().splitlines(), text.splitlines(), lineterm="")
        if ln.startswith("-") and not ln.startswith("---")]
    assert len(diff) == n_expected, (
        f"{tag}: {len(diff)} lines changed, expected {n_expected}:\n"
        + "\n".join(diff))
    cfg = CASE / "input" / f"{prefix}.dhs"
    cfg.write_text(text)
    print(f"[config] {cfg.name}: {n_expected} line(s) changed from {tag}")
    for key, old in changed.items():
        print(f"         {key}: {old}")
        print(f"           -> {get_value(text, key, tag)}")
    return cfg, out_prefix


def report_final(out_prefix):
    p = Path(str(out_prefix) + "Mass.Final.Balance")
    if not p.exists():
        return
    for ln in p.read_text().splitlines():
        if any(k in ln for k in FINAL_KEYS):
            print("      " + ln.strip())


def run_model(cfg, out_prefix):
    out_prefix.mkdir(parents=True, exist_ok=True)
    log = out_prefix / "terminal_log.txt"
    t0 = time.time()
    with open(log, "w") as f:
        proc = subprocess.run([str(DHSVM_EXE), str(cfg)], stdout=f,
                              stderr=subprocess.STDOUT, text=True)
    dt = time.time() - t0
    status = ("SUCCESS" if proc.returncode == 0
              else f"FAILED (rc {proc.returncode})")
    print(f"[run] {cfg.stem}: {status} in {dt / 60:.1f} min; log {log}")
    if proc.returncode != 0:
        tail = log.read_text().splitlines()[-25:]
        print("\n".join("      " + ln for ln in tail))
        return False
    for name in OUTPUT_FILES:
        p = Path(str(out_prefix) + name)
        size = f"  {p.stat().st_size} bytes" if p.exists() else ""
        print(f"      {'ok ' if p.exists() else 'MISSING'} {p.name}{size}")
    so = Path(str(out_prefix) + "Streamflow.Only")
    if so.exists():
        lines = so.read_text().splitlines()
        print(f"      Streamflow.Only header: {lines[0].strip()!r}  "
              f"rows {len(lines) - 1}")
    report_final(out_prefix)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=list(RUNS),
                    choices=list(RUNS))
    ap.add_argument("--kinds", nargs="+", default=KINDS[:3], choices=KINDS,
                    help="ctrl (April inputs), tierE (April inputs, Tier E "
                         "network), new (current pipeline outputs)")
    ap.add_argument("--smoke", action="store_true",
                    help="nine-day runs under the smoke prefixes")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-sha-check", action="store_true",
                    help="report the sha256 checks without stopping (tests)")
    args = ap.parse_args()
    check_inputs(args.kinds, sha=not args.no_sha_check)
    runs = args.runs
    ok = True
    for run in runs:
        template = CASE / "input" / RUNS[run]
        assert template.exists(), f"template missing: {template}"
        for kind in args.kinds:
            cfg, out_prefix = make_config(run, template, kind, args.smoke)
            if args.dry_run:
                continue
            ok = run_model(cfg, out_prefix) and ok
    print("done" if ok else "some runs failed")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
