# qgis_CA dependency inventory — for the standalone (DCC/cluster) rebuild

Date: 2026-06-08. Source of truth: `qgis_CA/` (Tier A/B audited). Goal: rebuild a true standalone (no QGIS) pipeline from this, for the Duke Compute Cluster, reproducing the qgis_CA outputs.

This is a planning inventory, not a port. It classifies every qgis_CA script and every stage of the orchestrator by external dependency, and names the portable replacement, so the GRASS-vs-pure-Python decision can be made before any code is written.

## 1. Per-script classification

| script | external dependency | rebuild verdict |
|---|---|---|
| `soildepthscript.py` | NumPy + GDAL only | port as-is (no change) |
| `dem_to_dhsvm_bins.py` | NumPy + GDAL only | port as-is |
| `generate_dhsvm_states.py` | NumPy + GDAL only | port as-is |
| `rdnbr_to_veg_bs.py` | NumPy + GDAL only | port as-is (veg phase) |
| `make_2class_veg_map.py` | NumPy + rasterio | port as-is, parameterize hardcoded paths (veg phase) |
| `make_4class_veg_map.py` | NumPy + rasterio (+ rasterio.warp) | port as-is, parameterize paths (veg phase) |
| `slope_fill.py` (new) | NumPy + GDAL only | already portable; used by prep |
| `channelclass.py` | qgis.core (QgsVectorLayer field editing) | rewrite I/O on geopandas; classification logic is pure Python |
| `dhsvm_area_config_from_dem.py` | qgis.core (QgsRasterLayer + QgsCoordinateTransform) | rewrite on rasterio + pyproj (small) |
| `plot_mainstem.py` | qgis.core + networkx + matplotlib | diagnostic only; port later or keep QGIS-only |
| `prep_dhsvm_inputs.py` | qgis.core + GRASS (5 processing.run) | the main rebuild target; see stage map below |

Bottom line: 7 of 11 scripts are already portable or trivially so. The work concentrates in `prep_dhsvm_inputs.py`, plus a geopandas rewrite of `channelclass.py` and `dhsvm_area_config_from_dem.py`.

## 2. prep_dhsvm_inputs.py — stage-by-stage

| stage | what it does | dependency | standalone replacement |
|---|---|---|---|
| config / paths | params, WS resolution | portable (drop the qgis imports) | unchanged |
| clip DEM | `gdal:cliprasterbymasklayer` | GDAL (via QGIS) | `rasterio.mask.mask` or `gdal.Warp` with cutline (deterministic, easy) |
| flow dir + accumulation | `grass7:r.watershed` | **GRASS hydrology** | decision a/b (see section 4) |
| stream extraction | `grass7:r.stream.extract` | **GRASS hydrology** | decision a/b |
| stream raster -> lines | `grass7:r.to.vect` | **GRASS** | decision a/b (raster-to-line is the least trivial in pure Python) |
| slope | `grass7:r.slope.aspect` (degrees) | GRASS | NumPy Horn's method (3x3), or richdem; trivial |
| slope conditioning | `slope_fill.fill_slope_nodata` | portable already | unchanged |
| row/col + length fields | QgsVectorLayer field editing, raster sample | qgis.core | geopandas GeoDataFrame + `rasterio.sample` |
| channel slope sampling | sample slope along lines (QgsPointXY) | qgis.core | `rasterio.sample.sample_gen` along shapely line points |
| channel classification | `channelclass.py` (importlib) | qgis.core | geopandas rewrite of channelclass |
| directed network + order | `_build_directed_by_FA`, `_toposort`, propagated order | **algorithm is pure Python**; geometry via QgsGeometry/QgsPointXY; nearest via QgsSpatialIndex | keep the graph logic; swap geometry to shapely, nearest to geopandas `.sindex` (rtree) or `shapely.STRtree` |
| write stream.network.dat / stream.map.dat | sample slope, write text | qgis.core sampling + plain text | `rasterio.sample` + text write (text part already portable) |
| soil depth | `soildepthscript.py` (importlib) | portable | unchanged |
| base maps | `dem_to_dhsvm_bins.py` (importlib) | portable | unchanged |
| initial states | `generate_dhsvm_states.py` (importlib) | portable | unchanged |

The directed-network and ordering logic (including the Tier A propagated/dense order fix) is pure-Python graph work and ports directly. The only QGIS surface there is geometry objects, raster sampling, and the spatial index.

## 3. qgis.core API -> portable library

| qgis.core | used for | portable replacement |
|---|---|---|
| QgsVectorLayer, QgsFeature, QgsField, QVariant | read/edit stream vector + attributes | geopandas GeoDataFrame |
| QgsGeometry, QgsPointXY | line geometry, endpoints, length | shapely LineString / Point |
| QgsSpatialIndex, nearestNeighbor | snap segment endpoints for topology | geopandas `.sindex` (rtree) or `shapely.STRtree` |
| QgsRasterLayer + provider.sample / identify | sample slope / elevation at points | `rasterio.open(...).sample([(x,y)])` |
| QgsCoordinateReferenceSystem, QgsCoordinateTransform | center lat/lon for [AREA] | pyproj.Transformer |
| QgsRasterLayer extent / rasterUnitsPerPixel | grid extent, cell size | rasterio `dataset.bounds`, `dataset.res` |
| processing `gdal:cliprasterbymasklayer` | clip DEM to watershed | rasterio.mask or gdal.Warp |

## 4. The central decision: GRASS hydrology replacement

Three stages are genuine GRASS hydrology: `r.watershed` (flow dir + accumulation), `r.stream.extract` (threshold + network + IDs), `r.to.vect` (stream raster to line topology). The other GRASS/GDAL calls (clip, slope) port trivially. Two paths:

a. Headless GRASS (grass-session, or conda `grass` in batch). Keeps r.watershed / r.stream.extract / r.slope.aspect / r.to.vect identical, so the standalone reproduces qgis_CA byte-for-byte. Lowest correctness risk. Cost: GRASS must be installable on DCC (it is, via conda, but it is a heavier dependency).

b. Pure-Python hydrology (pysheds / richdem / WhiteboxTools). No GRASS dependency, lightest install. Cost: different algorithms (depression handling, flow routing, stream-to-vector), so results will not match qgis_CA and the whole chain needs re-validation. `r.to.vect`'s raster-line-to-vector with topology is the least trivial to reproduce here (candidates: WhiteboxTools `RasterStreamsToVector`, or a skeleton-to-line builder).

Recommendation: confirm whether GRASS can run on DCC. If yes, prefer (a) for parity; the two interfaces (QGIS, cluster) then produce identical inputs, which matters for "same pipeline, two front ends" and for other users. Choose (b) only to eliminate GRASS entirely, accepting a re-validation pass.

## 5. Portable stack (proposed environment)

rasterio (+ rasterio.mask, rasterio.sample, rasterio.warp), geopandas, shapely, pyproj, numpy, networkx; plus either grass-session (path a) or pysheds/richdem/whitebox (path b). The existing compute modules already rely only on numpy + GDAL (osgeo) + rasterio.

## 6. Validation gate

The rebuilt standalone must reproduce qgis_CA's DHSVM inputs on the CA test case before it is trusted. Under path (a) the rasters (DEM, slope, soildepth, flow_acc) and the derived bins should match byte-for-byte where the algorithms are identical; the stream files should match exactly. Under path (b), expect agreement within tolerance and document the differences. Reuse the existing diagnostics: `compare_soildepth.py` for soil depth, plus simple raster and text diffs for the others.

## 7. Suggested rebuild order

1. Scaffold the standalone package and environment; decide path a vs b.
2. Port the trivial stages first: clip (rasterio.mask), slope (numpy + slope_fill), and the already-portable compute modules (soildepth, bins, states). Validate each against qgis_CA.
3. Hydrology core (r.watershed / r.stream.extract / r.to.vect) per the chosen path; validate flow_acc and the stream network against qgis_CA.
4. Vector I/O layer: geopandas rewrite of the field/topology/sampling stages and channelclass; validate stream.class.dat / stream.network.dat / stream.map.dat against qgis_CA.
5. dhsvm_area_config on rasterio + pyproj.
6. Veg-class scripts (make_2class / make_4class / rdnbr_to_veg_bs) with hardcoded paths parameterized.
7. plot_mainstem last (diagnostic), or leave QGIS-only.
