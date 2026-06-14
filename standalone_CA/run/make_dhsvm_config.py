"""Render a DHSVM input file (.dhs) from CA.dhs.template using paths.py.

CA.dhs.template is the calibrated baseline CA_LAI70.dhs with every input/output
path replaced by an @TOKEN@, and the [AREA] block replaced by @AREA_BLOCK@. The
paths and [AREA] are filled in here; all other model parameters and options stay
identical to the baseline. [AREA] is computed from the run DEM (area.py), so the
same template fits any resolution or watershed; on the 28 m CA DEM it reproduces
the calibrated block, so the standalone run still differs from the baseline only
in the pipeline-produced inputs, which is the point of the comparison.

The same template + paths.py generate a config for any input set: leave
DHSVM_OUT at the standalone outputs (default) or point it at qgis_CA_ref (which
has the same DHSVM_input_binaries/ DHSVM_input_streams/ modelstate/ layout), set
DHSVM_RUN_PREFIX so the outputs do not collide, and re-run.

Usage:
    python make_dhsvm_config.py
        -> writes {RUN_OUT}/{CASE}.dhs ; DHSVM output prefix {RUN_OUT}/{PREFIX}

The two stale Chiwawa shading/skyview paths are left as-is in the template:
Shading = FALSE, so DHSVM never opens them.
"""
import os
import re
import sys
from pathlib import Path

# paths.py lives in the sibling pipeline/ dir (single source of truth); area.py
# sits next to this file. Make both importable regardless of the working dir.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pipeline"))
sys.path.insert(0, str(HERE))
import paths  # noqa: E402
import area   # noqa: E402

TEMPLATE = HERE / "CA.dhs.template"

# DHSVM output prefix: the "Output Directory" value is a directory + filename
# prefix, so outputs become <prefix>Stream.Flow etc. Default <CASE>_sa keeps the
# standalone run apart from the baseline (CA_LAI70*) and a future qgis-ref run.
RUN_PREFIX = os.environ.get("DHSVM_RUN_PREFIX", f"{paths.CASE}_sa")
CONFIG_OUT = paths.RUN_OUT / f"{paths.CASE}.dhs"

# @TOKEN@ -> absolute path. STATE_DIR keeps the baseline trailing slash (DHSVM
# concatenates dir + State filename); OUTPUT_PREFIX has none (it is a prefix).
SUBS = {
    "@DEM_BIN@":        str(paths.BIN_DIR / "dem.bin"),
    "@MASK_BIN@":       str(paths.BIN_DIR / "mask.bin"),
    "@SOIL_BIN@":       str(paths.BIN_DIR / "soil.bin"),
    "@SOILDEPTH_BIN@":  str(paths.BIN_DIR / "soildepth.bin"),
    "@VEG_BIN@":        str(paths.BIN_DIR / "veg.bin"),
    "@STREAM_MAP@":     str(paths.STREAMS_DIR / "stream.map.dat"),
    "@STREAM_NETWORK@": str(paths.STREAMS_DIR / "stream.network.dat"),
    "@STREAM_CLASS@":   str(paths.STREAMS_DIR / "stream.class.dat"),
    "@MET_FILE@":       str(paths.MET_FILE),
    "@STATE_DIR@":      str(paths.STATE_DIR) + "/",
    "@OUTPUT_PREFIX@":  str(paths.RUN_OUT / RUN_PREFIX),
}


def render():
    if not TEMPLATE.exists():
        raise FileNotFoundError(f"template not found: {TEMPLATE}")
    text = TEMPLATE.read_text()

    # [AREA] is computed from the run DEM (paths.ELEV_CLIPPED), not frozen in the
    # template, so the same template fits any resolution / watershed. The time
    # zone meridian comes from paths.TIME_ZONE_MERIDIAN (declared, not derived
    # from longitude), so the config needs no hand edit.
    subs = dict(SUBS)
    subs["@AREA_BLOCK@"] = area.format_area_block(area.area_fields(paths.ELEV_CLIPPED)).rstrip("\n")

    for token, value in subs.items():
        text = text.replace(token, value)

    # fail loud if any @TOKEN@ slipped through (typo / unfilled placeholder)
    stray = sorted(set(re.findall(r"@[A-Z_]+@", text)))
    if stray:
        raise SystemExit(f"[ERROR] unsubstituted tokens remain: {stray}")

    paths.RUN_OUT.mkdir(parents=True, exist_ok=True)
    CONFIG_OUT.write_text(text)

    print(f"[make_dhsvm_config] wrote {CONFIG_OUT}")
    print(f"  case          : {paths.CASE}")
    print(f"  area DEM       : {paths.ELEV_CLIPPED}")
    print(f"  tz meridian    : {paths.TIME_ZONE_MERIDIAN}")
    print(f"  grid binaries : {paths.BIN_DIR}")
    print(f"  stream files  : {paths.STREAMS_DIR}")
    print(f"  initial state : {paths.STATE_DIR}/")
    print(f"  met file      : {paths.MET_FILE}")
    print(f"  output prefix : {paths.RUN_OUT / RUN_PREFIX}")
    print(f"  run binary    : {paths.DHSVM_BIN}")


if __name__ == "__main__":
    render()
