# Cross-engine comparison, evidence (2026-09-23)

Outputs of `scripts/diagnostics/cross_engine/compare_engines.py` for CA
and AR, the WW-DHSVM stream files it wrote, the Toolkit diagnostic run
on the WW-DHSVM network, and the DHSVM evaluation tables of the
`wA_` reruns. The record is `docs/audit/cross_engine_comparison_2026_09_23.md`.

- `CA/`, `AR/`: `comparison_<case>.md` and `.json` (network summaries,
  channel-cell overlap, WW-DHSVM checks, input sha256), `drop_sweep_ww_<case>.csv`
  (the ported constant-drop sweep, 10 to 600 cells), `orientation_<case>_ww.txt`
  (`check_network_orientation.py` on the WW-DHSVM rasters), `ww_streams/`
  (`stream.class.dat`, `stream.network.dat`, `stream.map.dat`,
  `Channel.State.01.01.2016.00.00.00`; the map header carries a generation
  timestamp, the body is reproducible).
- `dhsvm_eval_CA.csv`, `dhsvm_eval_AR.csv`: the metrics of the April control
  (`oA`), Tier E (`nA`) and WW-DHSVM (`wA`) runs from `tierE_eval_CA.py` and
  `tierE_eval_AR.py`.

The WW-DHSVM rasters (`ww_rasters/`, GeoTIFF) are not versioned; the script
regenerates them.
