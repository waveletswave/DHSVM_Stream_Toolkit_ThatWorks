"""Path constants for the standalone DHSVM preprocessing pipeline (CA case, DCC)."""
# NOTE: paths are hardcoded to the CA case and DCC absolute locations.
# Parameterization (env vars / config) is deferred until all trivial
# stages pass; see standalone_CA/docs/validation_log.md "Outstanding".
from pathlib import Path

# Inputs (persistent group space)
INPUTS = Path("/hpc/group/abmurraylab/ys451/dhsvm_ca/inputs")
SRC_DEM = INPUTS / "USGS_1_n36w084_20220725_UTM17.tif"
WATERSHED = INPUTS / "cabr_watershed_UTM17.shp"

# qgis_CA reference outputs, for byte-level validation
REF = Path("/hpc/group/abmurraylab/ys451/dhsvm_ca/qgis_CA_ref")
REF_ELEV_CLIPPED = REF / "elev_clipped.tif"

# Standalone outputs (scratch during development)
OUT = Path("/work/ys451/dhsvm_ca/standalone_dev/outputs")
OUT.mkdir(parents=True, exist_ok=True)
ELEV_CLIPPED = OUT / "elev_clipped.tif"

# Target CRS for this case (UTM Zone 17N). Used to re-stamp CRS after GRASS export.
EPSG = 32617

# Slope stage
REF_SLOPE = REF / "stream_slope_filled.tif"   # qgis_CA reference (already conditioned)
SLOPE_RAW = OUT / "slope_raw.tif"             # r.slope.aspect output before fill
SLOPE_FILLED = OUT / "slope_filled.tif"       # after slope_fill conditioning
