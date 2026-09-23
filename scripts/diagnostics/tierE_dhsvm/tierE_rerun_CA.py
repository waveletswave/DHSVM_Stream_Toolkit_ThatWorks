# -*- coding: utf-8 -*-
# =====================================================================
# tierE_rerun_CA.py  -  rerun the CA manuscript configurations today,
#                       once with the old network and once with the
#                       Tier E network, everything else unchanged
#
# For each run (CA_S4h from CA_prefire_S4h.dhs, CA_LAI70, CA_LAI20) two
# configurations are written from the template:
#
#   control  input/old_<tag>.dhs    the template with one line changed,
#            Output Directory -> TestCase/CA/output/old_<tag>
#            (the old DEM_CA_0406 network runs again with today's files;
#            comparing it with the April output shows whether the
#            manuscript runs are reproducible from today's inputs)
#
#   Tier E   input/new_<tag>.dhs    the template with five lines changed:
#            Stream Map / Network / Class File -> DEM_CA_tierE/...
#            Initial State Directory -> DEM_CA_tierE/modelstate/ (the
#            grid states are byte copies of the old ones, Channel.State
#            is the one rebuilt for the new network)
#            Output Directory -> TestCase/CA/output/new_<tag>
#
# <tag> is S4h, LAI70 or LAI20. DHSVM appends the file names to the
# Output Directory value, as in the existing runs
# (output/CA_LAI70Aggregated.Values); the terminal log goes to
# output/<prefix>/terminal_log.txt as before. The DEM, mask, soil, soil
# depth and vegetation binaries stay the DEM_CA_0406 files in both
# configurations, so control versus Tier E isolates the network.
#
# The Output Directory value must be at most 78 characters for this
# DHSVM3.2 build: the value is written into 100-byte buffers with a
# file name appended (InitDump.c `char sumoutfile[100]` with
# "failure_summary.txt", and the fork's "saturation_extent.txt", 21
# characters), and a longer value trips the fortified sprintf
# (__chk_fail_overflow, SIGTRAP). Observed: every prefix of 78
# characters ran (CA_LAIo90 and the like), 79 and 80 crashed after the
# five output files were opened, 86 crashed before. Under
# TestCase/CA/output/ (69 characters) that leaves 9 for the prefix.
#
# --base apr uses the April input set instead (DEM_CA_apr, a copy of
# "DEM_CA_0406 copy": the manuscript's soildepth.bin 7443b833.. and the
# pre-Tier-A stream.network.dat d0973965..; everything else is byte
# identical to DEM_CA_0406). Its control, oA_<tag>, points all nine
# input files and the state directory at DEM_CA_apr (10 lines changed)
# and must reproduce the April outputs byte for byte; its Tier E run,
# nA_<tag>, takes the binaries from DEM_CA_apr and the network and
# states from DEM_CA_tierE (10 lines). That is the strict R11: the
# manuscript inputs with only the network swapped.
#
# --smoke runs CA_S4h only, nine days (Model End 01/10/2016-00), under
# the prefixes oldsm_S4h and newsm_S4h (oAsm_S4h, nAsm_S4h for --base
# apr): it confirms that DHSVM reads the two SAVE outlet rows and shows
# the Initial Storage of the state files (606.150 mm in the April runs).
# --dry-run writes the configurations and prints what changed, no run.
#
# Usage (from anywhere):
#   python3 tierE_rerun_CA.py --smoke
#   python3 tierE_rerun_CA.py
#   python3 tierE_rerun_CA.py --base apr
#   python3 tierE_rerun_CA.py --runs CA_LAI70 CA_LAI20 --kinds tierE
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
    "/Users/benthosyy/Desktop/CodeBits/DHSVM-PNNL-2025/TestCase/CA"))
OLD_ROOT = Path(os.environ.get(
    "TIERE_OLD_ROOT",
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_0406"))
NEW_ROOT = Path(os.environ.get(
    "TIERE_NEW_ROOT",
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_tierE"))
APR_ROOT = Path(os.environ.get(
    "TIERE_APR_ROOT",
    "/Users/benthosyy/Desktop/CreateStreamNetwork_PythonV/DEM_CA_apr"))

RUNS = {"CA_S4h": "CA_prefire_S4h.dhs",
        "CA_LAI70": "CA_LAI70.dhs",
        "CA_LAI20": "CA_LAI20.dhs"}
TAGS = {"CA_S4h": "S4h", "CA_LAI70": "LAI70", "CA_LAI20": "LAI20"}
KINDS = ["ctrl", "tierE"]
# prefix words per base: old/new network on today's inputs (jun), oA/nA
# on the April inputs (apr)
WORDS = {"jun": {"ctrl": "old", "tierE": "new"},
         "apr": {"ctrl": "oA", "tierE": "nA"}}
MAX_OUTPUT_PATH = 78        # see the header


def prefix_for(kind, run, smoke, base="jun"):
    """Output prefix (the last path element of Output Directory), at
    most 9 characters: old_S4h, new_LAI70, oldsm_S4h, oA_S4h, nA_LAI20"""
    return f"{WORDS[base][kind]}{'sm' if smoke else ''}_{TAGS[run]}"


STREAM_KEYS = (("Stream Map File", "DHSVM_input_streams/stream.map.dat"),
               ("Stream Network File",
                "DHSVM_input_streams/stream.network.dat"),
               ("Stream Class File",
                "DHSVM_input_streams/stream.class.dat"))
BINARY_KEYS = ("DEM File", "Basin Mask File", "Soil Map File",
               "Soil Depth File", "Vegetation Map File")

# sha256 of the Tier E files as produced on DCC (fixed_CA_28m, 2026-09-22)
EXPECTED_SHA_NEW = {
    "DHSVM_input_streams/stream.class.dat":
        "26983e5342993a033cba391d720c81ce664c773cd3f0c554e16059e1999f4751",
    "DHSVM_input_streams/stream.map.dat":
        "a4f0f02ea78d8ab5c6c28d4d38c152f420bc3b772aff442ab2ec55b6d98dd442",
    "DHSVM_input_streams/stream.network.dat":
        "6a01bd2f6c1d64049ca85c920dcbaed09aa2a37ef380b9791b6af3b36a474595",
    "modelstate/Channel.State.01.01.2016.00.00.00":
        "819fa1f74422bb8e3404f47cdd02221a8957fcf3cb05acb89ca21cd41a24cb1f",
}
# sha256 of the DEM_CA_0406 files the manuscript runs point at, equal to
# the QGIS reference tree on DCC (qgis_CA_ref), checked 2026-09-22
EXPECTED_SHA_OLD = {
    "DHSVM_input_streams/stream.class.dat":
        "fd4ad1df567e612541c7f8e392e4501c83a6800bd9103c49bbe18822a550bca3",
    "DHSVM_input_streams/stream.map.dat":
        "a5f4bd32073bec4311cd548403fdade29af7801c21c8ac3d69589904a0079472",
    "DHSVM_input_streams/stream.network.dat":
        "e6ab190f0a683ff3f780124405885e137769c9bf1390a97db824819449ed8f4e",
    "modelstate/Channel.State.01.01.2016.00.00.00":
        "dc89adfb1acccf123f9b9466155a217b7d2807aa10bfaeda371feb07ec166826",
    "modelstate/Interception.State.01.01.2016.00.00.00.bin":
        "ba1b6c77b5f01352662febb861ab21d3603fa0a36d6d505476367e9b2f1ca37e",
    "modelstate/Snow.State.01.01.2016.00.00.00.bin":
        "c763d1f5826b63b85564c4402a69701eff59e95d06048b92682a2277a7f89c75",
    "modelstate/Soil.State.01.01.2016.00.00.00.bin":
        "3e71f33df8aa94a61ecccb2f0118aec6cf5c8a995d33394a640483fd82d05cf5",
    "DHSVM_input_binaries/dem.bin":
        "54539d716e807c8246cd139c1c687894f3c1730a36cb5ed0235c91ee8c63a454",
    "DHSVM_input_binaries/mask.bin":
        "3ef38adf096e43709d07afaa9f684906044f951c1fa5d3ee6dcccf80c00be1f6",
    "DHSVM_input_binaries/soil.bin":
        "3ef38adf096e43709d07afaa9f684906044f951c1fa5d3ee6dcccf80c00be1f6",
    "DHSVM_input_binaries/veg.bin":
        "3ef38adf096e43709d07afaa9f684906044f951c1fa5d3ee6dcccf80c00be1f6",
    "DHSVM_input_binaries/soildepth.bin":
        "ff76ec11d604edc726a2184e0d20354eef1b8658b3bd78988ccb2a674fc48410",
}
# sha256 of DEM_CA_apr, the April input set ("DEM_CA_0406 copy", files
# dated 2026-04-06 16:17): differs from DEM_CA_0406 in soildepth.bin and
# stream.network.dat only
EXPECTED_SHA_APR = dict(EXPECTED_SHA_OLD)
EXPECTED_SHA_APR.update({
    "DHSVM_input_binaries/soildepth.bin":
        "7443b83360e2251a6a865212b1380749080c73f61553ed33d717e775a4849ed8",
    "DHSVM_input_streams/stream.network.dat":
        "d0973965f82d344373a4c13f22610d7d9e2833e593deb019df54e97d7ac491d1",
})
BASES = {"jun": (OLD_ROOT, EXPECTED_SHA_OLD,
                 "DEM_CA_0406 (today's inputs, DCC qgis_CA_ref)"),
         "apr": (APR_ROOT, EXPECTED_SHA_APR,
                 "DEM_CA_apr (April inputs, the manuscript runs)")}
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


def check_inputs(kinds, base, sha=True):
    root, expected, label = BASES[base]
    check_tree(root, expected, label, sha)
    if "tierE" in kinds:
        check_tree(NEW_ROOT, EXPECTED_SHA_NEW,
                   "Tier E network (DCC fixed_CA_28m)", sha)
        for name in GRID_STATES:
            old = root / "modelstate" / name
            new = NEW_ROOT / "modelstate" / name
            assert new.exists(), \
                f"missing {new}: copy it from {old.parent}"
            assert sha256(old) == sha256(new), \
                f"{name} differs from the {root.name} one"
            print(f"  ok  modelstate/{name} identical to {root.name}")
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


def make_config(run, template, kind, smoke, base="jun"):
    text = template.read_text()
    tag = template.name
    prefix = prefix_for(kind, run, smoke, base)
    changed = {}
    network_keys = [k for k, _ in STREAM_KEYS] + ["Initial State Directory"]
    for key in network_keys + list(BINARY_KEYS):
        assert get_value(text, key, tag).startswith(str(OLD_ROOT)), \
            f"{tag}: {key} does not point at {OLD_ROOT}"
    n_expected = 1
    if base == "apr":
        # the manuscript binaries; the control also takes the April
        # network and states
        for key in BINARY_KEYS:
            rel = get_value(text, key, tag)[len(str(OLD_ROOT)) + 1:]
            text, changed[key] = set_line(text, key, str(APR_ROOT / rel),
                                          tag)
        n_expected += 5
        if kind == "ctrl":
            for key, rel in STREAM_KEYS:
                text, changed[key] = set_line(text, key,
                                              str(APR_ROOT / rel), tag)
            text, changed["Initial State Directory"] = set_line(
                text, "Initial State Directory",
                str(APR_ROOT / "modelstate") + "/", tag)
            n_expected += 4
    if kind == "tierE":
        for key, rel in STREAM_KEYS:
            text, changed[key] = set_line(text, key, str(NEW_ROOT / rel),
                                          tag)
        text, changed["Initial State Directory"] = set_line(
            text, "Initial State Directory",
            str(NEW_ROOT / "modelstate") + "/", tag)
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
    ap.add_argument("--kinds", nargs="+", default=KINDS, choices=KINDS,
                    help="ctrl (old network) and/or tierE")
    ap.add_argument("--base", default="jun", choices=list(BASES),
                    help="jun: today's DEM_CA_0406; apr: DEM_CA_apr")
    ap.add_argument("--smoke", action="store_true",
                    help="nine-day runs of CA_S4h under the smoke prefixes")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-sha-check", action="store_true",
                    help="report the sha256 checks without stopping (tests)")
    args = ap.parse_args()
    check_inputs(args.kinds, args.base, sha=not args.no_sha_check)
    runs = ["CA_S4h"] if args.smoke else args.runs
    ok = True
    for run in runs:
        template = CASE / "input" / RUNS[run]
        assert template.exists(), f"template missing: {template}"
        for kind in args.kinds:
            cfg, out_prefix = make_config(run, template, kind, args.smoke,
                                          args.base)
            if args.dry_run:
                continue
            ok = run_model(cfg, out_prefix) and ok
    print("done" if ok else "some runs failed")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
