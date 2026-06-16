# Running the pipeline on a new watershed

This toolkit turns a watershed polygon into a full DHSVM preprocessing input set: the gridded DEM, mask, slope, soil depth, stream network, and initial model state. You supply a watershed polygon and a target resolution. Everything below the DEM is derived from the DEM, so the per-basin inputs are small.

The pipeline is driven entirely by environment variables read in `standalone_CA/pipeline/paths.py`. Switching to a new basin is a matter of setting those variables; no code is edited.

## What you provide

- A watershed boundary polygon as a shapefile, in a projected CRS (a UTM zone). The polygon defines the basin; the DEM is fetched and masked to it.
- A target grid resolution in metres. 3DEP serves staged data at 10, 30, and 60 m; pick one of these per run.
- A target EPSG (the projected UTM zone of your basin).

## What you get

In the output directory: `DHSVM_input_binaries/` (dem, mask, soil, veg, soildepth, plus uniform soildepth variants), `DHSVM_input_streams/` (stream.class.dat, stream.network.dat, stream.map.dat), and `modelstate/` (Interception, Snow, Soil, and Channel state files), plus the intermediate GeoTIFFs the quick-look reads.

Note on scope. This toolkit produces the terrain, stream, and state inputs. The soil and vegetation binaries that `bins.py` writes are the toolkit defaults (uniform for the development basin). For a real simulation, replace the soil and vegetation/LAI inputs with maps for your basin, and supply your own meteorological forcing and run configuration. Those are separate from this preprocessing step.

## Prerequisites: two environments, step by step

The pipeline needs two separate environments. Keep them separate, and do not activate one inside the other (stacking a venv and a conda env in the same shell leads to the wrong Python being used). The only thing that crosses between them is a GeoTIFF, so their package versions do not have to agree.

You need conda first. If you do not have it, install Miniconda or Miniforge and confirm that `conda` works in your shell before continuing.

### Environment 1: the geometry environment

This runs the pipeline and the quick-look. Create it once:

```bash
conda create -n dhsvm_geo -c conda-forge python=3.10 \
    rasterio geopandas pyproj shapely numpy networkx pyflwdir scipy matplotlib gdal
```

cartopy is optional; the quick-look uses it for a basemap on the location panel and falls back to a plain graticule without it. Add `cartopy` to the create line if you want the basemap.

You also need GRASS GIS 7.x, which the slope and hydrology stages call. On a cluster, load it from the module system (for example `module load GRASS/7.6`); ask your admins for the exact module name. Off a cluster, install it from conda-forge or your system package manager.

One GRASS wrinkle. GRASS 7.6's launcher script is Python 2 on some builds. If `grass76 --exec ...` fails with a Python 2 or "bytes" error, the launcher needs a small Python 3 shim; the one used here is described in `standalone_CA/docs/DCC_SETUP.md`. Newer GRASS launchers (7.8 and 8.x) are Python 3 and need no shim. The pipeline was byte-validated against GRASS 7.6; other 7.x and 8.x versions should produce an equivalent network but have not been byte-checked.

### Environment 2: the fetch environment

This runs only `fetch_dem.py`, which pulls the DEM from 3DEP. A plain venv is enough:

```bash
python3 -m venv ~/fetch3dep
source ~/fetch3dep/bin/activate
pip install py3dep
deactivate
```

The fetch step needs outbound internet, so run it on a node that has it (on a cluster, a login node, not a compute node). The rest of the pipeline needs no network and can run on a compute node.

### A one-line session setup

So that each working session is a single command, put the GRASS load and the geometry-env activation in a small script and source it. On the development cluster this is `env.sh`, roughly:

```bash
# env.sh  (adjust the module name and the conda env name to your system)
module load GRASS/7.6
conda activate dhsvm_geo
```

Then every session in the geometry environment starts with `source /path/to/env.sh`.

### Storage

If you run across more than one node, put the fetched DEM and the output directory on shared storage, not node-local `/tmp`, so every node can read them.

## The flow, for one resolution

Set the shared environment first. Replace the placeholders with your values.

```bash
export DHSVM_EPSG=32617                                 # your basin's UTM zone EPSG
export DHSVM_WATERSHED=/shared/path/your_watershed.shp  # absolute path overrides the default
export DHSVM_CASE=YOURBASIN                              # a short tag
```

### Step 1: fetch the DEM at the target resolution (fetch env, networked node)

```bash
source ~/fetch3dep/bin/activate
cd /path/to/DHSVM_Stream_Toolkit_ThatWorks
python3 standalone_CA/prep/fetch_dem.py /shared/path/src_3dep_<RES>m.tif \
        --res <RES> --watershed "$DHSVM_WATERSHED"
deactivate
```

Fetch at the resolution you intend to model: 10, 30, or 60 m. For a second resolution, fetch again at that resolution into its own file. Do not fetch one resolution and resample it to another. The staged 3DEP product at each resolution is aggregated server-side from all available source data, which is the correct DEM for that grid; client-side resampling of a different product is not. The fetch output is in EPSG:5070; `prep_dem` reprojects it to your EPSG.

### Step 2: build the clipped DEM and set the stream support area A_c (geometry env)

Open a clean shell in the geometry env (do not keep the fetch venv active).

```bash
source /path/to/env.sh                    # loads GRASS and activates the geo env
cd /path/to/DHSVM_Stream_Toolkit_ThatWorks
export DHSVM_OUT=/shared/path/outputs_yourbasin_<RES>m
export DHSVM_DEM_SOURCE=/shared/path/src_3dep_<RES>m.tif
export DHSVM_DEM_RES=<RES>

# reproject and clip to the polygon so drop_analysis has a DEM to read
python3 standalone_CA/prep/prep_dem.py "$DHSVM_DEM_SOURCE"

# look at the constant-drop table and confirm there is a clean passing band
python3 standalone_CA/diagnostics/drop_analysis.py "$DHSVM_OUT/elev_clipped.tif"

# pin A_c from the objective and capture it into the environment
export DHSVM_STREAM_SOURCE_AREA_M2=$(python3 standalone_CA/diagnostics/drop_analysis.py \
       "$DHSVM_OUT/elev_clipped.tif" --emit-area)
echo "A_c = $DHSVM_STREAM_SOURCE_AREA_M2 m2"
```

Two things about this step.

First, the clip. `prep_dem` clips the grid tight to the basin by default, with a whole-cell margin of zero, so the grid hugs the basin and there is no wide nodata border. If you want a buffer of cells around the basin, pass `--pad-cells N` (or set `DHSVM_DEM_PAD_CELLS`).

Second, A_c. This is the one per-basin, per-resolution number you set once. `drop_analysis` sweeps the channel-initiation threshold and reports the objective as the start of the first sustained band of thresholds that satisfy the constant stream drop law, meaning at least three consecutive passes, which avoids picking an isolated noisy threshold. Run it without arguments first to read the table; `--emit-area` prints just the objective area so a script can capture it. The sweep range is in cells, so the useful range scales with resolution: a few tens of cells at 30 m is a few hundred at 10 m. If the band forms at the very edge of the swept range, widen it with `--tmin/--tmax` and rerun. On a new basin `drop_analysis` does not compare against any preset value; it reports the objective and a one-line hint to set `DHSVM_STREAM_SOURCE_AREA_M2` from it.

### Step 3: run the pipeline (geometry env)

```bash
bash standalone_CA/pipeline/run_pipeline.sh
```

With `DHSVM_DEM_SOURCE` set, stage 1 runs `prep_dem` (reproject to your EPSG at your resolution, mask to the polygon) and writes `elev_clipped.tif`. Stages 2 through 7 then build slope, base-map binaries, the hydrology and stream network, soil depth, and the initial states. The run is fail-fast and prints `PIPELINE COMPLETE` with the full input set at the end.

### Step 4: quick-look (geometry env)

Render a six-panel figure to eyeball the geometry: DEM with stream network and basin boundary, soil depth, slope, flow accumulation, flow direction, and basin location. Pass the watershed polygon with `--boundary` so the basin outline is drawn as a smooth line.

```bash
python3 standalone_CA/diagnostics/quicklook.py \
        --run-dir "$DHSVM_OUT" \
        --boundary "$DHSVM_WATERSHED" \
        --out "$DHSVM_OUT/quicklook/${DHSVM_CASE}_${DHSVM_DEM_RES}m_quicklook.png" \
        --title "$DHSVM_CASE  ${DHSVM_DEM_RES} m  pipeline quick-look"
```

This is a diagnostic for catching gross geometric errors, not a validation gate. Check that the stream network sits in the valleys, the boundary hugs the basin, flow accumulation concentrates toward the outlet, the basin fills most of the frame (the DEM panel title reports the basin's data fraction), and the location panel puts the basin where you expect. Run it on a node with internet if you want the cartopy basemap on the location panel; without internet or cartopy it falls back to a plain lon/lat graticule.

To run a second resolution, repeat all four steps with a new `DHSVM_OUT`, `DHSVM_DEM_RES`, and its own fetched DEM, and re-derive A_c for that resolution.

## Worked example: the AR (Arrowwood) watershed at 10 m and 30 m

Concrete commands, fetching each resolution natively from the AR polygon, with the values this basin produced.

```bash
# get the AR polygon onto shared storage (all shapefile parts)
scp /path/to/arwd_watershed_UTM17.* user@cluster:/work/ys451/dhsvm_ca/standalone_dev/

# shared environment
export DHSVM_EPSG=32617
export DHSVM_WATERSHED=/work/ys451/dhsvm_ca/standalone_dev/arwd_watershed_UTM17.shp
export DHSVM_CASE=AR

# ---- fetch both resolutions natively (fetch venv, login node) ----
source ~/fetch3dep/bin/activate
cd /hpc/group/abmurraylab/ys451/dhsvm_ca/DHSVM_Stream_Toolkit_ThatWorks
python3 standalone_CA/prep/fetch_dem.py /work/ys451/dhsvm_ca/standalone_dev/AR_src_3dep_10m.tif \
        --res 10 --watershed "$DHSVM_WATERSHED"
python3 standalone_CA/prep/fetch_dem.py /work/ys451/dhsvm_ca/standalone_dev/AR_src_3dep_30m.tif \
        --res 30 --watershed "$DHSVM_WATERSHED"
deactivate

# ---- process both (geometry env) ----
source /hpc/group/abmurraylab/ys451/env.sh
cd /hpc/group/abmurraylab/ys451/dhsvm_ca/DHSVM_Stream_Toolkit_ThatWorks

# AR at 10 m
export DHSVM_OUT=/work/ys451/dhsvm_ca/standalone_dev/outputs_AR_10m
export DHSVM_DEM_SOURCE=/work/ys451/dhsvm_ca/standalone_dev/AR_src_3dep_10m.tif
export DHSVM_DEM_RES=10
python3 standalone_CA/prep/prep_dem.py "$DHSVM_DEM_SOURCE"
python3 standalone_CA/diagnostics/drop_analysis.py "$DHSVM_OUT/elev_clipped.tif" --tmin 50 --tmax 800
export DHSVM_STREAM_SOURCE_AREA_M2=$(python3 standalone_CA/diagnostics/drop_analysis.py \
       "$DHSVM_OUT/elev_clipped.tif" --emit-area --tmin 50 --tmax 800)
bash standalone_CA/pipeline/run_pipeline.sh
python3 standalone_CA/diagnostics/quicklook.py --run-dir "$DHSVM_OUT" \
        --boundary "$DHSVM_WATERSHED" \
        --out "$DHSVM_OUT/quicklook/AR_10m_quicklook.png" \
        --title "AR (Arrowwood)  10 m  pipeline quick-look"

# AR at 30 m
export DHSVM_OUT=/work/ys451/dhsvm_ca/standalone_dev/outputs_AR_30m
export DHSVM_DEM_SOURCE=/work/ys451/dhsvm_ca/standalone_dev/AR_src_3dep_30m.tif
export DHSVM_DEM_RES=30
python3 standalone_CA/prep/prep_dem.py "$DHSVM_DEM_SOURCE"
python3 standalone_CA/diagnostics/drop_analysis.py "$DHSVM_OUT/elev_clipped.tif"
export DHSVM_STREAM_SOURCE_AREA_M2=$(python3 standalone_CA/diagnostics/drop_analysis.py \
       "$DHSVM_OUT/elev_clipped.tif" --emit-area)
bash standalone_CA/pipeline/run_pipeline.sh
python3 standalone_CA/diagnostics/quicklook.py --run-dir "$DHSVM_OUT" \
        --boundary "$DHSVM_WATERSHED" \
        --out "$DHSVM_OUT/quicklook/AR_30m_quicklook.png" \
        --title "AR (Arrowwood)  30 m  pipeline quick-look"
```

What this basin produced, for reference. At 10 m the tight grid is 159 x 209 (22849 valid cells, 69 percent of the grid), the drop objective is 220 cells (22000 m2), and the network is 67 segments with a single outlet. At 30 m the grid is 53 x 70 (2536 valid cells, 68 percent), the objective is 100 cells (90000 m2), and the network is 19 segments with a single outlet. Both reach PIPELINE COMPLETE with the full input set. Each resolution takes its own A_c from its own drop table; you set A_c as a physical area and the hydrology stage converts it to the cell count for that grid. The objective area is larger at 30 m than at 10 m, consistent with a coarser grid resolving fewer fine channels.

## Notes and rough edges

- EPSG and time zone. Set `DHSVM_EPSG` to your basin's UTM zone. For config generation later, set `DHSVM_TZ_MERIDIAN` to your standard meridian (default -75, US Eastern); it is declared, not derived from longitude.

- GRASS and OpenSSL. GRASS 7.6 links OpenSSL 1.1. If your system's OpenSSL has been upgraded to 3.x and the 1.1 compatibility library removed, GRASS fails to start with a `libssl.so.1.1: cannot open shared object file` error. Install an OpenSSL 1.1 build into its own location and add its lib directory to `LD_LIBRARY_PATH`; 1.1 and 3.x are different sonames, so they coexist. One way: `conda create -p /path/openssl11 -c conda-forge openssl=1.1.1`, then add `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/openssl11/lib` to your `env.sh`.

- The slope reference comparison. `slope.py` compares its output against a QGIS reference at `DHSVM_REF`. That reference is specific to the development basin. On a new basin the comparison reports no match, which is diagnostic only and does not affect the output or stop the run. If `DHSVM_REF` points at a path that does not exist on your system, the comparison may error; point it at any valid reference you have, or leave it at a path that exists.

- The GRASS messages. Lines about `gcs.csv` and `SetColorTable` from the GRASS stages are harmless on this GRASS 7.6 build; the pipeline stamps the CRS back with rasterio afterward.

- A_c sweep range. `drop_analysis` sweeps thresholds in cells, so the range scales with resolution. The objective is the start of the first sustained passing band, not the smallest single passing threshold, so an isolated noisy pass below the band is reported and skipped. If no band forms in the range, widen `--tmin/--tmax`.

- The tight clip. `prep_dem` clips to the basin with a zero-cell margin by default. `clip.py` is a separate entry that reproduces the development basin's 28 m grid exactly for the byte-level regression; use `prep_dem` for any new basin or resolution.

- The quick-look is diagnostic, not a gate. Use it to catch gross geometric errors. The byte-level regression on the development basin remains the source of truth for the pipeline itself.
