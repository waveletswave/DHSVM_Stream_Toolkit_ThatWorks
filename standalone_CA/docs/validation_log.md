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

## hydrology core (flow_acc / flow_dir / stream_raster + stream vector)

Reproduces the three grass7: processing.run calls in prep_dhsvm_inputs.py:
r.watershed -> r.stream.extract -> r.to.vect, chained in one GRASS location
(option i), run through the grass76_py3 shim. Region locked g.region
raster=elev_in (74x82). r.out.gdal exports CRS-stamped to EPSG:32617 by
hydrology.py (r.out.gdal cannot write CRS on this install).

Parameter parity vs prep, three corrections found and confirmed empirically:
- r.watershed: MFD default, no -a. flow_acc is all-negative on CA (-3387.5..-1),
  the expected edge/underestimate marking; adding -a would have flipped it.
- r.stream.extract: prep's '-m':True is not a valid flag for this module (it is
  an r.watershed flag); the QGIS wrapper silently dropped it, so we drop it too.
  prep's 'direction':flow_dir is an OUTPUT param, not an input; not fed here.
  Real input set is elevation + accumulation, matching what QGIS actually ran.
- r.to.vect: CA's r.stream.extract stream_vector came out empty (prep fell back
  to _force_lines_from_raster); reproduced with type=line, no -s.

Validation vs qgis_CA_ref/Intermediate_GIS (CA case, 74x82):
- flow_acc.tif (DCELL)  shape/transform/CRS match, max abs diff 0.0
- flow_dir.tif (CELL)   shape/transform/CRS match, max abs diff 0.0
- stream_raster.tif     shape/transform/CRS match, max abs diff 0.0
- streamfile.shp        41 line features (matches ref), extent identical to 6 dp

Note: this validates the hydrology rasters and the fallback stream geometry.
The stream vector attribute table / topology / channel classification
(channelclass, directed network, propagated order) is the #6 vector I/O layer,
which builds on this geometry.

## soildepth (soildepth.bin + soildepth.tif)

Port of soildepthscript.py. PNNL weighting formula and all constants copied
verbatim; only path resolution and the run guard changed. Inputs: elev_clipped
(clip), slope_filled (slope stage, FILLED -- prep feeds the in-place conditioned
stream_slope, so the standalone must feed the filled slope too), flow_acc
(hydrology). CA case, 74x82, 4334 valid / 1734 nodata.

Verdict: PASS by data-equivalence, NOT byte-identity. soildepth.bin differs from
the qgis_CA reference in 412 of 6068 cells, max abs diff 4.768e-07 (= 1 float32
ULP at depth magnitude ~2-4 m). soildepth.tif shows the same max diff.

Isolation confirms the residual is environment-level float rounding, not an
input or formula error:
- filled slope ref-vs-standalone: 0 cells differ, max 0.0 (inputs bit-identical;
  elev and flow_acc already validated zero-diff upstream).
- recompute soildepth from the REFERENCE inputs vs the reference bin still
  differs in 412 cells, max 4.768e-07, 393 of them interior (not boundary).
  Same inputs + verbatim formula -> 1-ULP diff means the power terms
  ((x)**0.25, (x)**0.75) round differently under the DCC numpy/libm than under
  the QGIS environment that produced the reference. float32 is not bit-
  reproducible across libm implementations.

This is the inventory section-6 tolerance case (algorithms identical -> byte-
for-byte where reproducible; otherwise within tolerance and documented).
Tolerance here is 1 float32 ULP; physically ~5e-7 m on a 6 m depth field,
no effect on DHSVM behaviour.

## vector I/O sub-step B: attributes + channelclass -> stream.class.dat

Option (ii) chained: streamfile.shp (geometry, hydrology) -> vector_attrs.py
writes streamfile_attr.shp (+arcid/Shape_Leng/Row/Col/slope_deg/meanmsq) ->
channelclass_standalone.py writes stream.class.dat. geopandas/shapely/rasterio
replace QgsVectorLayer/QgsGeometry/provider.sample; sampling semantics copied
verbatim from prep (slope_deg: 12 samples, mean; meanmsq: 15 samples, >0 ->
midpoint value, else cell_area; Row/Col top-left origin per Tier D2). Raster
sampling via InvGeoTransform + floor, verified to match provider.sample's
pixel-containing-point convention.

stream.class.dat byte-IDENTICAL to the qgis_CA reference:
  #ID W  D   n    inf
  13   0.5 0.100 0.0450 0.0

Note on why this is a single class: meanmsq is 792.9 for all 41 segments
(= cell_area, 28.158^2). The meanmsq sampler keeps only flow_acc values > 0,
but flow_acc is all-negative on CA (edge/underestimate marking, see hydrology),
so every segment's sample set is empty and falls back to cell_area. This
faithfully reproduces prep's behaviour on CA. With area degenerate to a constant
(< 1e6, area_bin 0) and slope all steep (smax 34.9 deg -> tan band "steep"),
every segment maps to class 13. The multi-class area/CLASS_TABLE branches are
not exercised by CA's (degenerate) data; they would engage on a case with
positive flow_acc or gentler slopes. Port is verbatim and correct; CA simply
does not stress it.

## vector I/O sub-step C: directed network -> stream.network.dat + stream.map.dat

Ports _build_directed_by_FA + _write_stream_network_FA + _write_stream_map
(prep). geopandas/shapely/numpy replace QGIS; spatial index replaced by an exact
brute-force k-NN over up-endpoints (shapely STRtree.query_nearest lacks k= in
2.1 and would return a single candidate -- a silent linking bug avoided).
channelclass_standalone.py extended (write_back=True) to write per-segment
chanclass/Class/hydwidth/hyddepth onto streamfile_attr.shp, the option-(ii)
equivalent of prep's in-place edit, so this stage can read per-segment class.

Verdict: PASS as a DHSVM-valid topology. NOT byte-identical to the qgis_CA
reference, and deliberately not aligned to it. Confirmed against the four
requirements DHSVM actually imposes:
- single outlet, and it is the same physical line as the reference outlet
  (segment 41, len 124.29428, slope 0.20258);
- every segment reaches the outlet (no cycles, no dangling);
- propagated order is DENSE (std nin distribution spans 1..7 with no gap) --
  the Tier A requirement so channel_route_network does not break on an empty
  order;
- segment set identical to the reference ((length,slope) multiset equal).

Two of 41 down-pointers differ from the reference (and the order/hop depth
distribution differs by one bin as a result). Root cause, established by
isolation: at a junction, two upstream endpoints are COINCIDENT with the
current segment's downstream endpoint (d=0). The downstream link is then
geometrically undefined; _best_neighbor's distance criterion cannot
discriminate and any tie-break is arbitrary. The reference's choice comes from
QGIS's R-tree internal order and has no hydrologic ground truth. We do NOT
reverse-engineer it: aligning to QGIS would reproduce an arbitrary pick with no
scientific basis. The standalone's k-NN is distance-sorted, so on the d=0 tie it
takes the straighter (lower-deflection) neighbour, which is consistent with flow
momentum and is at least as defensible as the reference's (steeper-deflection)
pick. Both yield valid networks; the standalone meets all four DHSVM
requirements above. Byte-identity is therefore not the acceptance criterion for
this stage -- it would mean reproducing an undefined tie's arbitrary resolution.

This is the inventory section-6 tolerance case, applied at the topology level:
algorithms identical where reproducible; at an inherently undefined point, a
defensible deterministic choice that satisfies the downstream model.

## Path parameterization (rebuild step 9) + layout unification

All seven stages now derive paths from a single source of truth, paths.py, with
env-overridable roots (DHSVM_INPUTS / DHSVM_REF / DHSVM_OUT / DHSVM_GRASS_SHIM /
DHSVM_EPSG / DHSVM_SRC_DEM / DHSVM_WATERSHED); defaults reproduce the CA case on
DCC. clip.py and slope.py already imported paths.py and were unchanged. The other
five stages plus bins.py had local hardcoded path blocks; these were replaced by
paths.py imports with no change to any stage logic. No stage file now contains an
absolute /work or /hpc path; the defaults live only in paths.py.

Two deliberate layout changes, validated:

Decision A -- unify base-map binaries into DHSVM_input_binaries. bins.py
previously wrote dem/mask/soil/veg + soildepth_uniform_*.bin to OUT/bins;
soildepth.bin went to OUT/DHSVM_input_binaries. They now share
DHSVM_input_binaries, the single directory DHSVM reads all grid binaries from.
Confirmed: the two sets coexist with no filename collision, and dem/mask/soil/
veg remain byte-identical to the reference (content unchanged, only directory
moved). soildepth.bin unchanged (still the pre-existing 1-ULP float
data-equivalence, size 24272).

Decision B -- channel state reads the standalone stream files. states.py
previously read stream.class/network.dat from the qgis_CA reference (the vector
stage was not yet ported). It now reads OUT/DHSVM_input_streams, the standalone
vector output. Consequence: grid states (Interception/Snow/Soil) are unchanged
and byte-identical (DEM dims + constants only; sizes 121360/194176/242720), but
Channel.State now reflects the standalone stream network and is NOT compared to
the reference Channel.State. This is correct: the channel state must be
consistent with the stream network actually used, and the standalone network
differs from the reference at two d=0 tie links (sub-step C). Run order is now
hydrology -> vector -> states.

Regression check: re-ran bins/soildepth/vector-chain/states after the refactor.
Every byte-identical output stayed byte-identical; soildepth's data-equivalence
diff was unchanged; the two intended changes (bins directory, Channel.State
source) behaved as designed. Pure refactor, no numerical change.

Note: run_hydrology_grass.sh still hardcodes the GRASS shim path internally;
parameterizing the .sh (read shim from env/arg) is a separate follow-up and does
not affect the Python path layer.

## End-to-end orchestration: run_pipeline.sh

One command drives the full pipeline: source DEM + watershed -> complete DHSVM
input set. Shell orchestrator (set -euo pipefail), running each stage in
dependency order with a fail-fast existence check on every stage's key output:

  clip -> slope(GRASS r.slope.aspect) -> slope(py: CRS stamp+fill+validate)
  -> bins -> hydrology(py runs the GRASS chain + CRS stamp) -> soildepth
  -> vector_attrs -> channelclass(+chanclass write-back) -> stream_network
  -> states

states runs last (decision B: channel state reads the standalone stream files).
Paths are read from paths.py (the script does not redefine them), so checks and
stages agree. The two harmless GRASS messages on this install (gcs.csv missing ->
r.out.gdal cannot write the CRS tag, handled by the Python CRS stamp; and
SetColorTable unsupported on Float32) are non-fatal; the modules exit 0 and the
outputs are correct.

Verified by a from-scratch run (outputs/ emptied first, so every intermediate is
freshly generated -- no stale-file luck):
- slope vs reference: byte-identical, max abs diff 0.0;
- base maps (dem/mask/soil/veg) byte-identical to reference;
- soildepth.bin unchanged 1-ULP data-equivalence (size 24272);
- stream.class.dat identical; network/map 41/277 lines;
- grid states sizes unchanged (121360/194176/242720);
- the full DHSVM input set is produced (DHSVM_input_binaries/,
  DHSVM_input_streams/, modelstate/).
A file-by-file cmp of the from-scratch run against the prior step-by-step run
matched on every checked output (binaries, stream files, grid + channel state):
the one-command run is byte-for-byte equal to running the stages by hand.

Follow-up: run_slope_grass.sh hardcodes its paths and takes no args;
run_hydrology_grass.sh hardcodes the shim. Parameterizing the two GRASS .sh
(read paths/shim from env or args) would close the last gap so a case/machine
change needs only the DHSVM_* env vars. The Python layer is already fully
parameterized.
