# AR: Toolkit versus WW-DHSVM, A_c = 47571.5 m2 (60 cells)

Grid 55 x 72 at 28.158 m, 2870 basin cells (2.276 km2). Same grid, mask, DEM and support area for both engines; no stream burning.

| | Toolkit (GRASS MFD, r.stream.extract) | WW-DHSVM (pyflwdir D8) |
|---|---|---|
| channel cells | 170 | 203 |
| segments | 19 | 29 |
| map records | 170 | 203 |
| cells under two segments | 0 | 0 |
| total length (m) | 5498.3 | 6602.4 |
| mean segment length (m) | 289.4 | 227.7 |
| drainage density (km/km2) | 2.416 | 2.902 |
| outlet segments | [19] | [29] |
| outlet tail cells (row, col, z) | {19: (53, 41, 796.61)} | {29: (53, 41, 796.61)} |
| SAVE rows | 1 | 1 |
| max routing rank | 7 | 8 |
| rank histogram | {1: 10, 2: 4, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1} | {1: 15, 2: 4, 3: 3, 4: 2, 5: 2, 6: 1, 7: 1, 8: 1} |
| max in-degree | 2 | 2 |
| in-degree histogram | {0: 10, 2: 9} | {0: 15, 2: 14} |
| classes used | {13: 16, 14: 3} | {13: 12, 14: 4, 15: 4, 16: 5, 17: 1, 18: 3} |
| max-|acc| cell | (53, 41) | (53, 41) |
| its |acc| (cells) | 2602.6 | 2870.0 |
| max-|acc| cell in an outlet segment | True | True |
| lowest channel cell | (53, 41) | (53, 41) |
| its elevation (m) | 796.61 | 796.61 |
| lowest channel cell in an outlet segment | True | True |

Channel-cell overlap: 151 common, 19 Toolkit only, 52 WW-DHSVM only, Jaccard 0.68. By Toolkit |acc| band (cells): 0-60: 2 of 3 / 36; 60-120: 45 of 56 / 57; 120-300: 49 of 55 / 55; 300-1000: 27 of 28 / 27; 1000-max: 28 of 28 / 28.
The Toolkit's max-|acc| cell (53, 41) is a WW-DHSVM channel cell: True; WW-DHSVM's max-area cell (53, 41) is a Toolkit channel cell: True.

Class tables (id: width, depth, Manning n):
- Toolkit: 13: 0.50, 0.100, 0.045; 14: 1.00, 0.100, 0.045
- WW-DHSVM: 13: 0.74, 0.150, 0.100; 14: 0.98, 0.150, 0.100; 15: 1.30, 0.179, 0.100; 16: 1.72, 0.216, 0.100; 17: 2.27, 0.260, 0.100; 18: 3.01, 0.314, 0.100

WW-DHSVM checks: topology ok True, outlet invariants ok True, errors [], warnings []; depression fill raised 0 cells; drop analysis objective 320 cells (0.25371466472630033 km2), band (320, 600), giving 71 channel cells in 7 segments.

WW-DHSVM stream files (sha256): class 696ebf8beab15040, network 3ce3cac9365a9682, map 0413270a3e800c4a; Channel.State written for 29 segments at 0.05 m initial depth.
