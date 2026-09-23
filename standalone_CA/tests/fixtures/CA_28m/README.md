# CA 28 m regression fixture

The five rasters are the GRASS-stage outputs of the standalone pipeline
on Camp Branch at 28 m (74 x 82 cells, EPSG:32617, 28.158 m), taken from
the Tier E run `/work/ys451/dhsvm_ca/tierE/fixed_CA_28m` on DCC
(2026-09-22, feature branch `feat-tierE-raster-network` at 71b8a4c, merged
as main 8173bee). They are what `segments_from_raster.py` reads:

| file | stage | sha256 |
|---|---|---|
| elev_clipped.tif | clip.py | 181539c9.. |
| slope_filled.tif | r.slope.aspect + slope_fill.py | c6c55348.. |
| flow_acc.tif | r.watershed | f0126f0d.. |
| stream_raster.tif | r.stream.extract | 6a6d112c.. |
| stream_dir.tif | r.stream.extract (direction) | 7e1c7acb.. |

`expected/` holds what the three network stages wrote from them on DCC:
`segments.csv` and `stream_cells.csv` from the same run, and
`stream.class.dat`, `stream.map.dat`, `stream.network.dat`, which are the
`DEM_CA_tierE` stream files the audit's DHSVM reruns used (sha256
26983e53.., a4f0f02e.., 6a01bd2f..; see
`docs/audit/tier_e_network_orientation_2026_09_22.md`).

`test_ca28m_regression.py` reruns the stages on the rasters and requires
the CSVs to agree value by value and the three stream files byte for
byte. A deliberate change to the network stage therefore means
regenerating `expected/` and recording why in
`standalone_CA/docs/validation_log.md`; an unexplained difference is a
regression.
