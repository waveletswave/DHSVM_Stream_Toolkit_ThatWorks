# DCC setup and project handoff (2026-06-08, updated 2026-06-09)

## Purpose

This document does two things: it snapshots where the DHSVM preprocessing toolkit stands on both tracks (the canonical QGIS pipeline and the planned DCC/cluster standalone), and it records a reproducible DCC environment runbook so another user, or a future session, can pick up without rediscovering the environment problems solved here. The standalone code-port plan lives in `docs/audit/standalone_rebuild_inventory_2026_06_08.md`; this document covers status and DCC operations, and points there for the port itself.

## 1. Where things stand

### Canonical pipeline

`qgis_CA/` is the canonical pipeline. It is the Tier A/B audited version and is what all new work builds from. It runs in QGIS, using QGIS processing for the GDAL clip and four GRASS algorithms (`r.watershed`, `r.stream.extract`, `r.to.vect`, `r.slope.aspect`), with the rest in NumPy/GDAL/rasterio plus pure-Python topology.

### Repo layout (as of 2026-06-09)

- `qgis_CA/` — canonical QGIS + GRASS pipeline.
- `standalone_CA/` — DCC standalone (no-QGIS) pipeline. Scaffold plus the clip/slope/bins trivial stages, validated against `qgis_CA` on the CA case (PR #5).
- `scripts/diagnostics/` — read-only audit and analysis tools (see `diagnostics_README.md` there). None are pipeline steps; none run DHSVM.
- `scripts/legacy/qgis/` and `scripts/legacy/standalone/` — retired pre-Tier copies, kept for reference (including the row/col-origin reference scripts). Not maintained.
- `docs/audit/` — audit and design docs, including the slope-conditioning and B-Q2 records and the standalone rebuild inventory.

### Tier status

- Tier A (channel network propagated/dense order): done and validated.
- Tier B unit fixes (slope `format=1` percent to degrees; `sample()` unpacking) plus the slope NoData conditioning: done. Tier B slope-correctness is closed. On CA, 318 boundary/void-ring cells now get interpolated slope instead of a spurious 0, removing the soil-depth boundary artifact (318 cells changed, 4016 byte-identical; channel files unaffected).
- Tier D: done and validated.
- B-Q2 (uniform vs variable soil depth): decision deferred. Evidence leans toward uniform, but undecided. Independent of the slope fix.
- B-Q1 (channel geometry): low priority for CA.

### Changes made this session

Four merged PRs: the slope conditioning (`slope_fill.py` plus the prep integration) with its audit doc and the B-Q2 report; the diagnostics plus the rebuild inventory; the legacy reorg; and the diagnostics README plus workflow-PDF rename.

## 2. Two tracks to finish

### QGIS version (multi-user)

Mostly working; it is the canonical pipeline and runs today. To make it cleanly usable by other people, the remaining work is: parameterize the hardcoded CA paths (the veg-class scripts `make_2class_veg_map.py`, `make_4class_veg_map.py`, `rdnbr_to_veg_bs.py`, and `dhsvm_area_config_from_dem.py` carry machine/case-specific paths), and document the QGIS version used and the run order for a new user. No algorithmic work is outstanding here.

### DCC standalone version

A true standalone (no QGIS) pipeline for the cluster. The decision is to keep GRASS for the hydrology (path a: headless GRASS, for parity with `qgis_CA`) rather than swap in a pure-Python hydrology stack. The per-script port plan, the qgis.core to portable mapping, the rebuild order, and the validation gate are in `docs/audit/standalone_rebuild_inventory_2026_06_08.md`. What that inventory does not cover, and what this document adds, is how to actually stand GRASS up on DCC (section 3).

## 3. DCC environment runbook

### Cluster basics

DCC uses Slurm for jobs and classic Tcl environment modules (not Lmod, so `module spider` is unavailable; use `module avail`). Modules are per-shell and per-node: after reconnecting, or on a newly allocated node, reload them. Relevant storage quotas for user ys451:

- `/hpc/home/ys451` — 25 GiB. Too small for conda environments.
- `/hpc/group/abmurraylab` — 1 TiB, ~616 GiB free. Persistent lab space. Put environments and data here.
- `/work/ys451` — 20 TiB. Scratch, purgeable. Not for persistent environments.

### Resume quick reference

To get a working shell with both the geo stack and GRASS after a fresh login, one command:

```bash
source /hpc/group/abmurraylab/ys451/env.sh
```

which runs:

```bash
module load Anaconda3/2024.02
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /hpc/group/abmurraylab/ys451/envs/dhsvm_geo
module load GRASS/7.6.0
export LC_ALL=C
```

Use `source`, not `bash` (conda activate and module load must affect the current shell).

All GRASS calls go through the Python 3 shim, not `grass76` directly:

```bash
python /hpc/group/abmurraylab/ys451/bin/grass76_py3.py -c <crs_raster.tif> <loc> --exec <module> ...
```

The shim patches the running launcher (lgettext, `file()`, collections ABCs) then runs it as `__main__`. See Step B for what it fixes. Login-node `--help` / `g.version` calls are light and fine without an allocation; a real `r.watershed` on a full DEM should use tmux + `srun --pty` to survive an SSH drop.

### Step A — Python geo stack (conda)  [VERIFIED]

The conda solve must run on a compute node. On the login node the solver is killed by the per-user memory limit (`Solving environment: Killed`). Allocate an interactive node first; a tmux session plus `srun --pty` survives an SSH drop until the walltime.

```bash
srun -p common --mem=16G -c 4 -t 2:00:00 --pty bash -i   # partition/limits per your allocation

module load Anaconda3/2024.02
source "$(conda info --base)/etc/profile.d/conda.sh"

ENVDIR=/hpc/group/abmurraylab/ys451/envs/dhsvm_geo
conda create -y --prefix $ENVDIR -c conda-forge --strict-channel-priority --solver=libmamba \
    python=3.10 rasterio geopandas pyproj shapely networkx numpy gdal

$ENVDIR/bin/python -c "import rasterio, geopandas, pyproj, shapely, networkx, numpy; print('geo stack ok')"
```

This is done. The environment exists at `/hpc/group/abmurraylab/ys451/envs/dhsvm_geo` and the import check prints `geo stack ok`. The env lives in shared group space, so the import check can be run from a login node without an allocation. Python is 3.10, which matters for GRASS (next step).

### Step B — GRASS 7.6  [VERIFIED]

Facts established about this GRASS:

- The module is `GRASS/7.6.0`; the launcher is `grass76`, not `grass`.
- The module only prepends PATH (`/opt/apps/rhel7/grass-7.6.0/bin`) and LD_LIBRARY_PATH (`/opt/apps/rhel7/compatlib`). It does not export GISBASE.
- GISBASE is `/opt/apps/rhel7/grass-7.6.0/grass-7.6.0` (from `grass76 --config path`).
- The four needed modules are present in the install: `r.slope.aspect` and `r.to.vect` were confirmed on disk, and `r.watershed` and `r.stream.extract` are core modules in the same install.
- GRASS 7.6 (2018) does not start under Python 3.12+ (the Anaconda base): it aborts at startup with `AttributeError: 'NullTranslations' object has no attribute 'lgettext'`, because `gettext.lgettext` was removed in 3.12. Running it under Python 3.10 (the `dhsvm_geo` env, 3.10.20) clears the startup abort, but is not sufficient on its own; two more Python 2/3 issues surface under `--exec` and at cleanup, fixed by the shim below.
- `export LC_ALL=C` is required, otherwise startup complains that default locale settings are missing.

#### Resolved: location creation and the GRASS launcher (shim)

**Location creation.** Building a location from an EPSG code fails on this install: GRASS 7.6 wants the old GDAL 2.x CSV EPSG files (`gcs.csv`) to translate `EPSG:32617`, but the conda GDAL ships modern PROJ (`proj.db`) with no `gcs.csv`, so `--tmp-location EPSG:32617` errors with `Unable to translate EPSG code`. Resolved with option (b): build the location from a CRS-bearing GeoTIFF, which bypasses the EPSG lookup entirely. Option (a) (point GDAL_DATA at a bundled `gcs.csv`) was not needed, and is closer to real usage anyway since the pipeline starts from a projected DEM.

```bash
python /hpc/group/abmurraylab/ys451/bin/grass76_py3.py -c /path/to/projected_dem.tif /tmp/gloc --exec g.version
```

**Launcher shim.** Python 3.10 clears the startup `lgettext` abort, but two further Python 2/3 incompatibilities surface and corrupt `--exec` runs. Both are patched by a shim at `/hpc/group/abmurraylab/ys451/bin/grass76_py3.py`, which patches the running process then runs the real launcher as `__main__`:

1. `lgettext()` returns bytes in Python 3, so the launcher's `_("Executing <%s> ...") % batch_job_string` is `bytes % str` and raises TypeError on every `--exec`. The shim forces the lgettext family to delegate to `gettext()` (str).
2. Python 2's `file()` builtin is gone; `grass/script/db.py` calls `file(os.devnull, 'w')` during post-exec cleanup (`clean_all` -> `clean_default_db`), which crashes after the job finishes and corrupts the exit code. The shim maps `file` -> `open`. It also aliases the collections ABCs (Iterable, Mapping, etc.) back onto `collections`, which Python 3.10 moved to `collections.abc`.

All pipeline GRASS calls go through the shim: `python /hpc/group/abmurraylab/ys451/bin/grass76_py3.py ... --exec ...`.

**Two harmless `r.out.gdal` quirks (raster data unaffected).** Every export on this install prints two errors that do not change the pixel values:

- `ERROR 4: Unable to open EPSG support file gcs.csv` — the output GeoTIFF ends up with `crs=None`. Fix: stamp the CRS back with rasterio in `r+` mode after export (`s.crs = CRS.from_epsg(32617)`), which rewrites only the header.
- `ERROR 6: SetColorTable() only supported for Byte or UInt16` — GRASS tries to attach a color table TIFF rejects for Float32. Harmless.

**Verified end to end through the shim, exit 0:** location creation, `g.version`, and the four hydrology modules (`r.watershed`, `r.stream.extract`, `r.to.vect`, `r.slope.aspect`) at `--help`. `r.slope.aspect` has also run for real on the CA DEM (the slope stage; output matches `qgis_CA` to the pixel). The heavy `r.watershed` / `r.stream.extract` / `r.to.vect` runs on the full DEM come with the hydrology core.

### Key insight: one env serves both

The same Python 3.10 conda env both launches `grass76` (clearing the lgettext crash) and runs rasterio/geopandas. There is no need for a separate GRASS Python or for grass-session, and no two-interpreter conflict to manage. The standalone architecture stays as file passing: run the GRASS stages via the shim with `--exec ... r.out.gdal` to write GeoTIFFs, then read those with rasterio in the same env. Rebuild scripts should activate `dhsvm_geo`, `module load GRASS/7.6.0`, and `export LC_ALL=C` (or just `source /hpc/group/abmurraylab/ys451/env.sh`), invoke GRASS through the `grass76_py3` shim, and stamp the CRS back onto exports with rasterio in `r+` mode (see Step B).

## 4. Standalone rebuild — next steps

Status (2026-06-09): the package is scaffolded and the clip, slope, and bins trivial stages are ported and validated against `qgis_CA` on the CA case, all with zero data difference (PR #5 merged). Next: the states stage (`generate_dhsvm_states.py`, NumPy+GDAL), then the hydrology core via the shim.

Follow section 7 of `docs/audit/standalone_rebuild_inventory_2026_06_08.md`: scaffold the package, port the trivial stages first (clip via `rasterio.mask`; slope via GRASS `r.slope.aspect` plus `slope_fill`; the already-portable compute modules), then the hydrology core via `grass76 --exec` (through the shim), then the geopandas vector I/O layer. Each stage is validated against `qgis_CA` on the CA case before moving on (the validation gate). Reuse the existing diagnostics for the soil-depth and raster comparisons.

## 5. Open items

- DCC: real `r.watershed` / `r.stream.extract` / `r.to.vect` runs on the full CA DEM still pending. Location creation and the four modules are confirmed exit 0 through the shim (`g.version`, `--help`); `r.slope.aspect` has run for real (slope stage). The heavy hydrology runs come with the hydrology core.
- Standalone: scaffold plus clip/slope/bins done and validated (PR #5). Next: port states, then the hydrology core.
- Standalone: hydrology core via the shim; validate flow accumulation and the stream network vs `qgis_CA`.
- Standalone: geopandas vector I/O layer; validate `stream.class.dat`, `stream.network.dat`, `stream.map.dat` vs `qgis_CA`.
- Both: parameterize the hardcoded CA paths (veg-class scripts and area config), and fold the per-stage local path blocks into `paths.py` in one pass.
- QGIS: document the QGIS version and the run order for a new user.
