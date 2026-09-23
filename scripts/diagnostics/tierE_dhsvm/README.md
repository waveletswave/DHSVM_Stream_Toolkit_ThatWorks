# Tier E DHSVM reruns (CA, 2026-09-22)

The three scripts that produced the DHSVM-side numbers in
`docs/audit/tier_e_network_orientation_2026_09_22.md`. They ran on the
Mac from `DHSVM-PNNL-2025/TestCase/CA/`, against the DHSVM3.2 build of
2025-01-13; every path is a default that `TIERE_*` environment variables
override.

- `tierE_rerun_CA.py`: writes one configuration per run and network
  from the manuscript templates (`CA_prefire_S4h.dhs`, `CA_LAI70.dhs`,
  `CA_LAI20.dhs`), changing only the lines listed in its header, checks
  the sha256 of every input it points at, runs DHSVM and reports the
  Mass.Final.Balance totals. `--base jun` (today's `DEM_CA_0406`) writes
  `old_<tag>` and `new_<tag>`; `--base apr` (the April manuscript set,
  `DEM_CA_apr`) writes `oA_<tag>` and `nA_<tag>`; `--smoke` runs nine
  days. The Output Directory value is kept at or under 78 characters,
  the limit of this DHSVM build.
- `tierE_compare_CA.py`: manuscript versus control (reproducibility),
  control versus Tier E (network effect), column by column for
  Aggregated.Values and Mass.Balance, the Stream.Flow totals, R12 (the
  SAVE outlets against the routed total) and the channel lateral inflow,
  which is the manuscript's simulated Q.
- `tierE_eval_CA.py`: the manuscript's streamflow metrics with the
  manuscript's own code, imported by path from
  `SpongeBurn/10_CA_Calib/DHSVM_CA_Stitch_Optimal_LAI_V2.py`.
