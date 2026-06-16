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

Port of generate_dhsvm_states.py. Byte-parity port: only path resolution and the run-on-import to run() guard changed; grid-state writes and channel-state parsing/format copied verbatim. CA case, NY=74 NX=82 (6068 cells).

Grid binaries byte-identical to the qgis_CA reference:
- Interception.State.01.01.2016.00.00.00.bin  121360 B (5 x 6068 x 4)  identical
- Snow.State.01.01.2016.00.00.00.bin          194176 B (8 x 6068 x 4)  identical
- Soil.State.01.01.2016.00.00.00.bin          242720 B (10 x 6068 x 4) identical

Channel.State.01.01.2016.00.00.00 (text) identical via diff. Channel state reads the qgis_CA reference stream files (stream.class.dat / stream.network.dat); this validates the channel-state logic. It will track the standalone hydrology once those stream files are ported and pass the byte gate.

Every grid-state array is spatially uniform, so output bytes are fixed by (array count, write order, dtype, constants) alone; identical by construction given matching grid dims. Reference: qgis_CA modelstate, regenerated 2026-06-08.

## hydrology core (flow_acc / flow_dir / stream_raster + stream vector)

Reproduces the three grass7: processing.run calls in prep_dhsvm_inputs.py: r.watershed -> r.stream.extract -> r.to.vect, chained in one GRASS location (option i), run through the grass76_py3 shim. Region locked g.region raster=elev_in (74x82). r.out.gdal exports CRS-stamped to EPSG:32617 by hydrology.py (r.out.gdal cannot write CRS on this install).

Parameter parity vs prep, three corrections found and confirmed empirically:
- r.watershed: MFD default, no -a. flow_acc is all-negative on CA (-3387.5..-1), the expected edge/underestimate marking; adding -a would have flipped it.
- r.stream.extract: prep's '-m':True is not a valid flag for this module (it is an r.watershed flag); the QGIS wrapper silently dropped it, so we drop it too. prep's 'direction':flow_dir is an OUTPUT param, not an input; not fed here. Real input set is elevation + accumulation, matching what QGIS actually ran.
- r.to.vect: CA's r.stream.extract stream_vector came out empty (prep fell back to _force_lines_from_raster); reproduced with type=line, no -s.

Validation vs qgis_CA_ref/Intermediate_GIS (CA case, 74x82):
- flow_acc.tif (DCELL)  shape/transform/CRS match, max abs diff 0.0
- flow_dir.tif (CELL)   shape/transform/CRS match, max abs diff 0.0
- stream_raster.tif     shape/transform/CRS match, max abs diff 0.0
- streamfile.shp        41 line features (matches ref), extent identical to 6 dp

Note: this validates the hydrology rasters and the fallback stream geometry. The stream vector attribute table / topology / channel classification (channelclass, directed network, propagated order) is the #6 vector I/O layer, which builds on this geometry.

## soildepth (soildepth.bin + soildepth.tif)

Port of soildepthscript.py. PNNL weighting formula and all constants copied verbatim; only path resolution and the run guard changed. Inputs: elev_clipped (clip), slope_filled (slope stage, FILLED -- prep feeds the in-place conditioned stream_slope, so the standalone must feed the filled slope too), flow_acc (hydrology). CA case, 74x82, 4334 valid / 1734 nodata.

Verdict: PASS by data-equivalence, NOT byte-identity. soildepth.bin differs from the qgis_CA reference in 412 of 6068 cells, max abs diff 4.768e-07 (= 1 float32 ULP at depth magnitude ~2-4 m). soildepth.tif shows the same max diff.

Isolation confirms the residual is environment-level float rounding, not an input or formula error:
- filled slope ref-vs-standalone: 0 cells differ, max 0.0 (inputs bit-identical; elev and flow_acc already validated zero-diff upstream).
- recompute soildepth from the REFERENCE inputs vs the reference bin still differs in 412 cells, max 4.768e-07, 393 of them interior (not boundary). Same inputs + verbatim formula -> 1-ULP diff means the power terms ((x)**0.25, (x)**0.75) round differently under the DCC numpy/libm than under the QGIS environment that produced the reference. float32 is not bit-reproducible across libm implementations.

This is the inventory section-6 tolerance case (algorithms identical -> byte-for-byte where reproducible; otherwise within tolerance and documented). Tolerance here is 1 float32 ULP; physically ~5e-7 m on a 6 m depth field, no effect on DHSVM behaviour.

## vector I/O sub-step B: attributes + channelclass -> stream.class.dat

Option (ii) chained: streamfile.shp (geometry, hydrology) -> vector_attrs.py writes streamfile_attr.shp (+arcid/Shape_Leng/Row/Col/slope_deg/meanmsq) -> channelclass_standalone.py writes stream.class.dat. geopandas/shapely/rasterio replace QgsVectorLayer/QgsGeometry/provider.sample; sampling semantics copied verbatim from prep (slope_deg: 12 samples, mean; meanmsq: 15 samples, >0 -> midpoint value, else cell_area; Row/Col top-left origin per Tier D2). Raster sampling via InvGeoTransform + floor, verified to match provider.sample's pixel-containing-point convention.

stream.class.dat byte-IDENTICAL to the qgis_CA reference:
  #ID W  D   n    inf
  13   0.5 0.100 0.0450 0.0

Note on why this is a single class: meanmsq is 792.9 for all 41 segments (= cell_area, 28.158^2). The meanmsq sampler keeps only flow_acc values > 0, but flow_acc is all-negative on CA (edge/underestimate marking, see hydrology), so every segment's sample set is empty and falls back to cell_area. This faithfully reproduces prep's behaviour on CA. With area degenerate to a constant (< 1e6, area_bin 0) and slope all steep (smax 34.9 deg -> tan band "steep"), every segment maps to class 13. The multi-class area/CLASS_TABLE branches are not exercised by CA's (degenerate) data; they would engage on a case with positive flow_acc or gentler slopes. Port is verbatim and correct; CA simply does not stress it.

## vector I/O sub-step C: directed network -> stream.network.dat + stream.map.dat

Ports _build_directed_by_FA + _write_stream_network_FA + _write_stream_map (prep). geopandas/shapely/numpy replace QGIS; spatial index replaced by an exact brute-force k-NN over up-endpoints (shapely STRtree.query_nearest lacks k= in 2.1 and would return a single candidate -- a silent linking bug avoided). channelclass_standalone.py extended (write_back=True) to write per-segment chanclass/Class/hydwidth/hyddepth onto streamfile_attr.shp, the option-(ii) equivalent of prep's in-place edit, so this stage can read per-segment class.

Verdict: PASS as a DHSVM-valid topology. NOT byte-identical to the qgis_CA reference, and deliberately not aligned to it. Confirmed against the four requirements DHSVM actually imposes:
- single outlet, and it is the same physical line as the reference outlet (segment 41, len 124.29428, slope 0.20258);
- every segment reaches the outlet (no cycles, no dangling);
- propagated order is DENSE (std nin distribution spans 1..7 with no gap) -- the Tier A requirement so channel_route_network does not break on an empty order;
- segment set identical to the reference ((length,slope) multiset equal).

Two of 41 down-pointers differ from the reference (and the order/hop depth distribution differs by one bin as a result). Root cause, established by isolation: at a junction, two upstream endpoints are COINCIDENT with the current segment's downstream endpoint (d=0). The downstream link is then geometrically undefined; _best_neighbor's distance criterion cannot discriminate and any tie-break is arbitrary. The reference's choice comes from QGIS's R-tree internal order and has no hydrologic ground truth. We do NOT reverse-engineer it: aligning to QGIS would reproduce an arbitrary pick with no scientific basis. The standalone's k-NN is distance-sorted, so on the d=0 tie it takes the straighter (lower-deflection) neighbour, which is consistent with flow momentum and is at least as defensible as the reference's (steeper-deflection) pick. Both yield valid networks; the standalone meets all four DHSVM requirements above. Byte-identity is therefore not the acceptance criterion for this stage -- it would mean reproducing an undefined tie's arbitrary resolution.

This is the inventory section-6 tolerance case, applied at the topology level: algorithms identical where reproducible; at an inherently undefined point, a defensible deterministic choice that satisfies the downstream model.

## Path parameterization (rebuild step 9) + layout unification

All seven stages now derive paths from a single source of truth, paths.py, with env-overridable roots (DHSVM_INPUTS / DHSVM_REF / DHSVM_OUT / DHSVM_GRASS_SHIM / DHSVM_EPSG / DHSVM_SRC_DEM / DHSVM_WATERSHED); defaults reproduce the CA case on DCC. clip.py and slope.py already imported paths.py and were unchanged. The other five stages plus bins.py had local hardcoded path blocks; these were replaced by paths.py imports with no change to any stage logic. No stage file now contains an absolute /work or /hpc path; the defaults live only in paths.py.

Two deliberate layout changes, validated:

Decision A -- unify base-map binaries into DHSVM_input_binaries. bins.py previously wrote dem/mask/soil/veg + soildepth_uniform_*.bin to OUT/bins; soildepth.bin went to OUT/DHSVM_input_binaries. They now share DHSVM_input_binaries, the single directory DHSVM reads all grid binaries from. Confirmed: the two sets coexist with no filename collision, and dem/mask/soil/ veg remain byte-identical to the reference (content unchanged, only directory moved). soildepth.bin unchanged (still the pre-existing 1-ULP float data-equivalence, size 24272).

Decision B -- channel state reads the standalone stream files. states.py previously read stream.class/network.dat from the qgis_CA reference (the vector stage was not yet ported). It now reads OUT/DHSVM_input_streams, the standalone vector output. Consequence: grid states (Interception/Snow/Soil) are unchanged and byte-identical (DEM dims + constants only; sizes 121360/194176/242720), but Channel.State now reflects the standalone stream network and is NOT compared to the reference Channel.State. This is correct: the channel state must be consistent with the stream network actually used, and the standalone network differs from the reference at two d=0 tie links (sub-step C). Run order is now hydrology -> vector -> states.

Regression check: re-ran bins/soildepth/vector-chain/states after the refactor. Every byte-identical output stayed byte-identical; soildepth's data-equivalence diff was unchanged; the two intended changes (bins directory, Channel.State source) behaved as designed. Pure refactor, no numerical change.

Note: run_hydrology_grass.sh still hardcodes the GRASS shim path internally; parameterizing the .sh (read shim from env/arg) is a separate follow-up and does not affect the Python path layer.

## End-to-end orchestration: run_pipeline.sh

One command drives the full pipeline: source DEM + watershed -> complete DHSVM input set. Shell orchestrator (set -euo pipefail), running each stage in dependency order with a fail-fast existence check on every stage's key output:

  clip -> slope(GRASS r.slope.aspect) -> slope(py: CRS stamp+fill+validate)
  -> bins -> hydrology(py runs the GRASS chain + CRS stamp) -> soildepth
  -> vector_attrs -> channelclass(+chanclass write-back) -> stream_network
  -> states

states runs last (decision B: channel state reads the standalone stream files). Paths are read from paths.py (the script does not redefine them), so checks and stages agree. The two harmless GRASS messages on this install (gcs.csv missing -> r.out.gdal cannot write the CRS tag, handled by the Python CRS stamp; and SetColorTable unsupported on Float32) are non-fatal; the modules exit 0 and the outputs are correct.

Verified by a from-scratch run (outputs/ emptied first, so every intermediate is freshly generated -- no stale-file luck):
- slope vs reference: byte-identical, max abs diff 0.0;
- base maps (dem/mask/soil/veg) byte-identical to reference;
- soildepth.bin unchanged 1-ULP data-equivalence (size 24272);
- stream.class.dat identical; network/map 41/277 lines;
- grid states sizes unchanged (121360/194176/242720);
- the full DHSVM input set is produced (DHSVM_input_binaries/, DHSVM_input_streams/, modelstate/).
A file-by-file cmp of the from-scratch run against the prior step-by-step run matched on every checked output (binaries, stream files, grid + channel state): the one-command run is byte-for-byte equal to running the stages by hand.

Follow-up: run_slope_grass.sh hardcodes its paths and takes no args; run_hydrology_grass.sh hardcodes the shim. Parameterizing the two GRASS .sh (read paths/shim from env or args) would close the last gap so a case/machine change needs only the DHSVM_* env vars. The Python layer is already fully parameterized.

## GRASS .sh parameterization (closes the follow-up above)

Both GRASS shell scripts now take their paths and the shim as arguments, so nothing case- or machine-specific is hardcoded in them.

run_slope_grass.sh takes three positional args: clipped DEM, slope_raw out, and the shim. It builds its GRASS location in a throwaway /tmp dir (PID-tagged), matching run_hydrology_grass.sh, instead of a fixed /work path. run_hydrology_grass.sh now reads the shim from its third arg; it already took the DEM and out_dir. run_pipeline.sh resolves the three slope args from paths.py (the existing python3 -c idiom) and passes them in; hydrology.py imports SHIM from paths.py and passes it as the third arg.

Result: a case or machine change needs only the DHSVM_* env vars end to end. Both front-end layers, Python and GRASS, are now parameterized.

Regression: pure refactor. Under the CA defaults each arg equals the former hardcoded value (slope DEM = paths.ELEV_CLIPPED, slope out = paths.SLOPE_RAW, shim = paths.SHIM), so the GRASS-stage outputs are byte-identical by construction. The slope LOC change is a scratch GRASS location only and does not touch the exported slope_raw.tif. Re-ran run_pipeline.sh from a clean outputs/ and cmp'd the GRASS-stage outputs against the pre-change run: slope_raw.tif, flow_acc/flow_dir/stream_raster.tif, and streamfile.shp/.shx/.prj are byte-for-byte equal. streamfile.dbf differs only in the dBASE last-update date byte (the two runs were a day apart); no data byte differs.

The standalone preprocessing pipeline was validated against the QGIS reference pipeline at three levels. Hourly streamflow over the full 2016-2018 simulation is numerically equivalent (NSE 0.99999510, PBIAS -0.000059%). Evaluated against observed Camp Branch streamflow over 2017-2018, the two pipelines yield identical performance: NSE 0.641, PBIAS -15.4%, r 0.822 for both, with absolute differences of zero at reporting precision. Annual ET partitioning is identical across runs (transpiration 576 mm, interception 172 mm, total 748 mm). The two input sets differ only in two documented numerical respects: a 1-ULP float difference in soil depth on 412 cells, and a d=0 tie-break in stream network extraction. These produce a maximum single-step streamflow deviation of 0.7% at the largest flow event and do not affect any aggregate metric. The CA test case is validated end-to-end.

## Stream network d=0 tie-break: cell-set diagnosis

`iso_check_network.py` reports the standalone `stream.network.dat` as NOT clean against the reference: all 41 segments map, but 2 down-pointers, 7 nin values, and 2 map.dat cell-sets disagree under its segment-ID relabeling. That relabeling is built by matching segments on rounded (length, slope, class), so the open question was whether the 2 differences are real or an artifact of the signature match pairing the wrong segments.

The ref-to-std correspondence was rebuilt from the map.dat cell sets instead of the signature. The streamfile geometry is byte-identical and in the same record order (vector sub-step A), so corresponding segments occupy the same (col, row) cells. Matching by cell set with Jaccard overlap (order-independent, and robust to the ~40 shared junction cells that a last-writer cell-to-segment lookup mis-assigns) gives a clean bijection: 39 of 41 segments match exactly and the other 2 differ by a single cell. This cell-based correspondence agrees with the signature relabeling on all 41 segments, so the signature match was not pairing the wrong segments.

Under this true correspondence the 2 down-pointer differences remain. Reference segment 11 flows to ref 24, but the standalone segment on the same cells (std 15) flows to std 27; ref 36 flows to ref 38, but std 38 flows to std 40. Each is one headwater connecting to a different downstream segment at a confluence. This is the d=0 tie-break the extraction code already anticipates: where two downstream arcs are equidistant, the standalone takes the first-seen under a stable argsort and the QGIS reference takes the R-tree order. One cell of the rasterized footprint moves with the tie, (39, 37) in the reference versus (36, 36) in the standalone.

The 7 nin differences are not 7 independent topology differences. NIN_MODE is "propagated", so column 2 is a cumulative upstream count; the 2 root tie-break differences propagate down each path to the outlet and surface as 7 changed values. The in-degree recomputed directly from the down-pointers differs at only 4 segments, exactly the two that gain and the two that lose an inflow from the 2 reroutes.

Both networks are valid DHSVM drainage networks: a single outlet (segment 41, down 0), every segment reaches it by following down, and no cycles. The two differences sit at arbitrary equidistant ties with no hydrologic ground truth, which the acceptance criteria treat as acceptable, since the standard is satisfying DHSVM's topology rather than reproducing the QGIS tie-break. No code change follows from this. This is the network-side source of the 0.7% maximum single-step streamflow deviation in the end-to-end validation, which leaves every aggregate metric unchanged.

## Resolution test: A_c at 10 m

The physical-area stream threshold (line A) predicts that holding the support area A_c fixed, rather than the cell count, keeps the extracted network in the same geomorphic regime as grid resolution changes. This was tested by repeating the channel-initiation analysis on a 10 m DEM.

The 10 m DEM is the USGS 3DEP 1/3 arc-second tile n36w084 (NAD83 geographic, ~10 m native), reprojected to EPSG:32617 on a clean 10 m grid by bilinear resampling and masked to the watershed (prep/prep_dem_10m.py). Masking before routing matters: drop_analysis routes the whole raster, so an unmasked tile would draw in neighbouring drainage and the comparison would mix domain extent with resolution. The valid cell count scales as expected, 4334 at 28 m to 34347 at 10 m, a factor of 7.92 that matches (28.16/10)^2, confirming the same watershed at a finer grid.

drop_analysis on the 10 m DEM gives an objective threshold of 500 cells, 0.0500 km^2, with a passing band (|t| < 2) of roughly 0.050 to 0.10 km^2. At 28 m the objective was 50 cells, 0.0396 km^2, band 0.040 to 0.079 km^2. The two passing bands overlap over 0.050 to 0.079 km^2, so the support-area range the constant-drop law accepts is consistent across the 2.8x resolution change. The 10 m crossing is softer than the 28 m one: |t| hovers near 2 from about 0.038 to 0.048 km^2 before dropping, so the 10 m objective reads as a band rather than a sharp point.

The chosen value carries across cleanly. 60 cells at 28 m is 0.04757 km^2; the same physical area is 476 cells at 10 m. That value sits inside the 28 m band and at the lower edge of the 10 m band (|t| ~ 2.0 at 476 cells). Had the cell count been held fixed instead, 60 cells at 10 m would be 0.006 km^2, far denser than any threshold the law accepts, so it is the physical-area specification, not a cell count, that keeps the network valid across resolution. This is the resolution-portability the parameterization was built for.

There is a modest, expected drift: the objective A_c rises about a quarter from 28 m to 10 m, and drainage density at the chosen area is slightly higher at 10 m. Both are consistent with the known mild resolution-dependence of DEM channel initiation, where a finer grid resolves more fine-scale convergence. The drift stays within the overlap of the passing bands and does not unseat the chosen value for the validated 28 m CA case.

Pipeline-side confirmation: run on the 10 m DEM, the hydrology stage derives the threshold from STREAM_SOURCE_AREA_M2 and the 10 m cell size as 476 cells with no edit. It extracts a 55-feature network, 8.9 km total, against 41 features at 28 m. The extra segments and slightly greater length are the finer grid resolving more channel detail at the same support area, consistent with the drainage-density drift noted above. The network is dendritic and tracks the valleys, with no cross-ridge artifacts.

## DHSVM at 10 m: streamflow resolution sensitivity

This is the model-side companion to the earlier A_c resolution test. That test showed the channel-initiation threshold and the stream network are resolution-portable. This one runs DHSVM on the full 10 m input set and compares the simulated hydrograph at the Camp Branch outlet to the validated 28 m run, against observed CB over 2017 to 2018.

The 10 m input set was built by the standalone pipeline at 251 x 275 (34347 valid cells, A_c 47571.5 m2 giving 476 cells, a 55-segment single-outlet network) and run with DHSVM3.2 on the Mac. ET, soil, and vegetation are unchanged from the 28 m setup: a single soil type, a single vegetation type, the same met forcing, and the same calibrated parameters. The only thing that changes between the two runs is the DEM resolution, so the comparison isolates the grid.

Integrated yield is robust to resolution. Total simulated Q is 1125 mm at 28 m and 1134 mm at 10 m, within about one percent. PBIAS moves from -15.4 percent at 28 m to -14.8 percent at 10 m, the 10 m run slightly closer to observed. Annual ET partitioning is identical across resolutions (transpiration 576, canopy interception 172, soil evaporation 0, total 748 mm), as expected since vegetation is uniform and the met forcing does not depend on the grid. The fire-induced cumulative excess, observed minus the pre-fire-vegetation counterfactual, is about +400 mm at both resolutions, so the scientific signal does not depend on resolution.

Daily fit degrades modestly at 10 m. Overall NSE moves from 0.64 to 0.44 and r from 0.82 to 0.71. The drop is concentrated in 2018, where NSE moves from 0.69 to 0.47; 2017 is slightly better at 10 m, NSE -0.12 to -0.01. The 10 m hydrograph carries more sustained recession baseflow and larger errors at a few peak events.

The reading is conservative. The water balance, the ET partitioning, and the first-year deficit recovering to a second-year match are all resolution-robust. The lower daily NSE and r at 10 m are consistent with running the 28 m-calibrated parameters, in particular the soil depth scaling and the lateral conductivity, on a finer grid without recalibration, and with the finer grid resolving a flashier saturation and routing response. This is a reassuring but not surprising supplementary result. It does not bear on the soil-versus-vegetation question or the alpha-to-beta_u mapping, both of which are yield and ET arguments.

The 28 m pipeline reproducibility is unaffected. Standalone versus QGIS at 28 m remains identical at the reporting precision, NSE difference 0.000, PBIAS difference 0.0 percent, r difference 0.000. The config-generation rework that produced the 10 m configuration (DEM-derived [AREA], declared meridian) reproduces the 28 m [AREA] block byte for byte, so the 28 m baseline is intact.

One minor accounting note. The evaluation divides the outlet volume by a fixed basin area of 3.436e6 m2 for all runs. The 10 m valid area is 3.4347e6 m2, about 0.04 percent smaller, which shifts the 10 m Q by well under a tenth of a percent and does not change any conclusion.

## Portability: AR (Arrowwood), tight-clip, 10 m and 30 m

The pipeline was run end to end on a second watershed, AR (Arrowwood), a Coweeta basin adjacent to Camp Branch, at 10 m and 30 m. The aim is a portability check on the general entry point, fetch_dem and prep_dem, rather than a reproduction of the QGIS reference, so the acceptance criterion is a valid DHSVM input set on an arbitrary basin, not byte-identity against CA. The reference-comparison stages still run and still report no match against the CA reference, which is diagnostic only on a different basin and a different resolution.

Source DEMs were fetched natively per resolution from 3DEP, with no resampling between resolutions: AR_src_3dep_10m.tif (EPSG:5070, 351 x 309, 9.12 m) and AR_src_3dep_30m.tif (EPSG:5070, 117 x 104, 27.43 m). prep_dem reprojects each to clean UTM 17N (EPSG:32617) at the target resolution and masks to arwd_watershed_UTM17.shp. Outputs are in /work/ys451/dhsvm_ca/standalone_dev/outputs_AR_10m and outputs_AR_30m. The basin sits at 35.181 N, 83.516 W.

### Tight DEM clip (prep_dem PAD_CELLS = 0)

prep_dem previously padded the target grid by a fixed 200 m, which left a wide nodata margin around the basin: about 20 cells each side at 10 m and 7 at 30 m, so the cell margin even differed by resolution. CA, built by clip.py against the QGIS reference window, was tight; AR carried the margin. The model only computes masked cells, so the margin never changed any result, but it bloated the grid, the binaries, and the quick-look figures. prep_dem now pads by a whole-cell count, default 0, so the grid hugs the basin bounding box snapped to the grid, matching CA. The margin is the same ring of cells at any resolution, and a --pad-cells argument (env DHSVM_DEM_PAD_CELLS) restores a buffer if one is wanted. clip.py and the CA 28 m byte-match are untouched; this is the general prep_dem path only.

The tight clip preserves the basin DEM exactly. AR re-clipped at both resolutions keeps the same valid cell count as the padded version, 22849 at 10 m and 2536 at 30 m, with only the nodata margin removed: 10 m goes from 199 x 249 to 159 x 209 (of 33231 total), 30 m from 67 x 84 to 53 x 70 (of 3710 total). The tight grid is a subset of the padded grid on the same snapped lattice, so the basin cells fall at the same coordinates and reproject to the same values. The basin DEM is byte-preserved and only the surrounding nodata is dropped. Basin fraction rises from about 45 percent to 69 percent at 10 m and 68 percent at 30 m, close to CA's 71 percent.

### A_c is extent-independent

drop_analysis on the tight 10 m DEM gives the same objective as the padded version, 220 cells = 22000 m2, band 220 to 600, with the isolated sub-band passes at 140 and 190 rejected. Trimming the nodata margin does not change the basin's internal flow accumulation, so the support-area threshold is unchanged. The sweep matches the padded run except for one- or two-count differences in a few of the densest thresholds, where the basin touching the grid edge shifts the edge-cell D8 slightly; the band start is unaffected. At 30 m the objective is 100 cells = 90000 m2, band 100 to 140. The threshold is read from STREAM_SOURCE_AREA_M2, so the hydrology stage derives 220 cells at 10 m and 100 cells at 30 m with no edit.

### Both resolutions produce a valid input set

Both runs reach PIPELINE COMPLETE with the full input set (DHSVM_input_binaries, DHSVM_input_streams, modelstate). The stream network is a valid DHSVM topology at both resolutions: a single outlet, propagated order dense to 8, and a two-class channel network (10 m class counts 7:4 and 13:63; 30 m 7:2 and 13:17). The quick-look figure for each was reviewed. The stream network sits in the valleys, the boundary hugs the basin, flow accumulation concentrates toward the southern outlet, and the location panel places the basin at Coweeta.

Segment counts differ from the prior padded AR run, and the difference is the A_c change, not the clip. The padded run used the old drop objective, the smallest passing threshold: 140 cells at 10 m giving 94 segments at order 9, and 20 cells at 30 m giving 52 segments at order 10. The tight run uses the new sustained-band objective: 220 cells at 10 m giving 67 segments at order 8, and 100 cells at 30 m giving 19 segments at order 8. The larger support areas initiate channels later, so the networks are sparser and the order is one to two lower. This is the expected consequence of the drop-band change made this session, and it is orthogonal to the tight clip, which leaves the basin DEM and A_c unchanged.

The reading. A polygon plus a resolution yields a full DHSVM terrain, stream, and state input set on an arbitrary basin, with the DEM fetched server-side and the grid clipped tight. The tight AR input set is not byte-identical to the padded AR input set, by design: the grid dimensions and the stream.map.dat row and column indices shift with the grid origin, and A_c differs because the drop objective changed. The equivalence is at the level the model uses, a valid clipped DEM and a valid drainage network, with the basin DEM itself byte-preserved across the clip.
