# standalone_CA validation log

Validation of the standalone (no-QGIS, DCC) DHSVM preprocessing pipeline against the Tier A/B audited `qgis_CA` reference, on the CA (Camp Branch) test case. Each stage is checked before the next is ported (the inventory's validation gate). Reference outputs live in `/hpc/group/abmurraylab/ys451/dhsvm_ca/qgis_CA_ref/`; standalone development outputs in `/work/ys451/dhsvm_ca/standalone_dev/outputs/`.

CRS: EPSG:32617 (UTM 17N). Grid: 74 rows x 82 cols = 6068 cells, of which 4334 are inside the watershed (data) and 1734 are nodata.

## Acceptance criteria

Two levels, by output type.

Rasters (GeoTIFF) are checked for data equivalence, not byte identity. GRASS and rasterio write different TIFF headers and GRASS emits a color table that does not round-trip, so byte-identity is not expected; what matters is that DHSVM reads pixel values, not the container. A raster passes when shape, transform, dtype, and the nodata layout all match the reference and the maximum absolute per-pixel difference over the shared data pixels is within tolerance.

Binaries (.bin) are checked for byte identity. These are raw `numpy.tofile()` dumps with no header, so sha256 equality is exact and is the standard applied.

## Methodology note: lock the grid to the reference

The recurring issue across stages is not the algorithm but how the output grid aligns to the source grid. `rasterio.mask(crop=True)` rounds the crop window outward (includes any touched edge pixel); GDAL warp / GRASS align differently. The first clip attempt came out 76x84 against the reference 74x82, a one-pixel ring on every side, with the transform origin offset by exactly one pixel.

The fix that generalizes: do not rely on either tool's rounding. Lock the output window explicitly to the reference footprint on the shared source grid, then mask. For GRASS stages the equivalent is importing the clipped DEM and running `g.region raster=<clipped>` before any raster op, so the computation lands on the same 74x82 grid. This is what made both clip and slope align to the pixel. Expect to apply the same lock to the hydrology stages.

## Stage results

### Stage 1 — clip (elev_clipped.tif)

Replaces `gdal:cliprasterbymasklayer`. rasterio window locked to the reference bounds on the source grid, then `rasterio.features.geometry_mask` sets pixels outside the watershed polygon to nodata=-9999, output dtype float32.

Result: PASS (data-equivalent).
- shape 74x82 = ref; transform match: True; dtype float32; nodata -9999.
- nodata layout identical: True (1734 / 1734).
- data pixels compared: 4334; max abs diff: 0.0; pixels diff > 0: 0.
- byte-identical: False (TIFF header only; data identical).

This stage fixes the valid-pixel extent that every downstream .bin inherits.

### Stage 2 — slope (slope_filled.tif)

Replaces `r.slope.aspect` + slope conditioning. Decision: keep GRASS `r.slope.aspect` (format=degrees) rather than a NumPy Horn's reimplementation, so the result is the same algorithm as the reference and aligns numerically. Region locked to the clipped DEM (g.region reported rows 74, cols 82, cells 6068). CRS stamped back after `r.out.gdal` (the gcs.csv export issue). Conditioned with the shared `slope_fill.fill_slope_nodata` using elev_clipped.tif as the valid-extent reference.

Result: PASS (data-equivalent, exact).
- shape 74x82 = ref; transform match: True; dtype float32; nodata NaN = ref.
- nodata layout identical: True (1734 / 1734).
- data pixels compared: 4334; max abs diff: 0.0; pixels diff > 0: 0.
- slope_fill: filled 318 of 318 missing-slope cells, 0 unfilled (matches the known boundary/void ring count); filled extent equals the DEM valid extent.
- byte-identical: False (TIFF header only; data identical).

Units confirmed as degrees on the reference before porting (max 44.9, p99 40.1; mountainous CA terrain), consistent with the Tier B fix of format=1 (percent) to format=0 (degrees).

### Stage 3 — base-map binaries (dem/mask/soil/veg + uniform soil depths)

Port of `dem_to_dhsvm_bins.py`. Reads elev_clipped.tif; valid mask from DEM nodata; writes dem.bin (Float32), mask/soil/veg.bin (Int8, 1 inside / 0 outside), and uniform soil-depth baselines at 2.0/2.5/3.0/3.5/4.0 m (Float32). Raw `tofile()`, C-order, little-endian, no header. soil.bin/veg.bin here are the uniform placeholders; per-class veg_bs binaries are a later stage.

Result: PASS (byte-identical, all nine files).
- valid pixels 4334 / 6068, inherited from clip unchanged.
- dem.bin and uniform soil depths: 24272 bytes = 6068 x 4 (Float32), sha256 identical.
- mask/soil/veg.bin: 6068 bytes = 6068 x 1 (Int8), sha256 identical.
- byte order, C-order, dtype all confirmed correct; .bin format locked for downstream soildepth/veg binaries.

## Known environment notes (carried from DCC_SETUP)

- GRASS 7.6 launcher needs the py3 shim (`grass76_py3.py`): lgettext returns bytes (bytes % str crash on every --exec) and Python 2 file() is gone (post-exec cleanup crash). Shim patches both, plus the collections ABC aliases.
- `r.out.gdal` cannot write the CRS into the output GeoTIFF on this install (ERROR 4, gcs.csv missing) and tries to set a color table TIFF rejects (ERROR 6). Both are harmless to the data; CRS is stamped back with rasterio in r+ mode after every export.

## Outstanding

- soildepth depends on flow_acc (hydrology output), so it cannot complete before the hydrology core; it is not a no-dependency trivial stage despite being NumPy+GDAL only.
- generate_dhsvm_states (states): no-dependency trivial, not yet ported.
- Hydrology core (r.watershed / r.stream.extract / r.to.vect): highest-risk stage, not started; needs the full processing.run parameter dicts from prep_dhsvm_inputs.py.
- Vector I/O layer (channelclass, network/order, stream.*.dat): after hydrology.
- paths.py is hardcoded to CA/DCC absolute paths; parameterization deferred until all trivial stages pass, then done in one pass before merge.

## states (modelstate: Interception / Snow / Soil / Channel)

Port of generate_dhsvm_states.py. Byte-parity port: only path resolution and the
run-on-import to run() guard changed; grid-state writes and channel-state
parsing/format copied verbatim. CA case, NY=74 NX=82 (6068 cells).

Grid binaries byte-identical to the qgis_CA reference:
- Interception.State.01.01.2016.00.00.00.bin  121360 B (5 x 6068 x 4)  identical
- Snow.State.01.01.2016.00.00.00.bin          194176 B (8 x 6068 x 4)  identical
- Soil.State.01.01.2016.00.00.00.bin          242720 B (10 x 6068 x 4) identical

Channel.State.01.01.2016.00.00.00 (text) identical via diff. Channel state reads
the qgis_CA reference stream files (stream.class.dat / stream.network.dat); this
validates the channel-state logic. It will track the standalone hydrology once
those stream files are ported and pass the byte gate.

Every grid-state array is spatially uniform, so output bytes are fixed by
(array count, write order, dtype, constants) alone; identical by construction
given matching grid dims. Reference: qgis_CA modelstate, regenerated 2026-06-08.
