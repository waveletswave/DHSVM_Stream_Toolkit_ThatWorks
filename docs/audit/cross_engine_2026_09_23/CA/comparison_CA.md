# CA: Toolkit versus WW-DHSVM, A_c = 47571.5 m2 (60 cells)

Grid 74 x 82 at 28.158 m, 4334 basin cells (3.436 km2). Same grid, mask, DEM and support area for both engines; no stream burning.

| | Toolkit (GRASS MFD, r.stream.extract) | WW-DHSVM (pyflwdir D8) |
|---|---|---|
| channel cells | 231 | 266 |
| segments | 22 | 30 |
| map records | 231 | 266 |
| cells under two segments | 0 | 0 |
| total length (m) | 7600.8 | 8901.2 |
| mean segment length (m) | 345.5 | 296.7 |
| drainage density (km/km2) | 2.212 | 2.59 |
| outlet segments | [12, 22] | [29] |
| outlet tail cells (row, col, z) | {12: (67, 68, 829.61), 22: (65, 68, 838.0)} | {29: (68, 68, 826.81)} |
| SAVE rows | 2 | 1 |
| max routing rank | 8 | 10 |
| rank histogram | {1: 12, 2: 2, 3: 2, 4: 2, 5: 1, 6: 1, 7: 1, 8: 1} | {1: 16, 2: 4, 3: 2, 4: 2, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1} |
| max in-degree | 2 | 3 |
| in-degree histogram | {0: 12, 2: 10} | {0: 16, 2: 13, 3: 1} |
| classes used | {13: 17, 14: 5} | {6: 1, 13: 11, 14: 5, 15: 3, 16: 4, 17: 2, 18: 4} |
| max-|acc| cell | (63, 67) | (68, 68) |
| its |acc| (cells) | 3387.5 | 4334.0 |
| max-|acc| cell in an outlet segment | True | True |
| lowest channel cell | (67, 68) | (68, 68) |
| its elevation (m) | 829.61 | 826.81 |
| lowest channel cell in an outlet segment | True | True |

Channel-cell overlap: 153 common, 78 Toolkit only, 113 WW-DHSVM only, Jaccard 0.445. By Toolkit |acc| band (cells): 0-60: 1 of 5 / 77; 60-120: 31 of 56 / 52; 120-300: 37 of 73 / 46; 300-1000: 36 of 44 / 43; 1000-max: 48 of 53 / 48.
The Toolkit's max-|acc| cell (63, 67) is a WW-DHSVM channel cell: True; WW-DHSVM's max-area cell (68, 68) is a Toolkit channel cell: False.

Class tables (id: width, depth, Manning n):
- Toolkit: 13: 0.50, 0.100, 0.045; 14: 1.00, 0.100, 0.045
- WW-DHSVM: 6: 3.57, 0.352, 0.030; 13: 0.72, 0.150, 0.100; 14: 1.00, 0.150, 0.100; 15: 1.37, 0.186, 0.100; 16: 1.89, 0.230, 0.100; 17: 2.59, 0.284, 0.100; 18: 3.57, 0.352, 0.100

WW-DHSVM checks: topology ok True, outlet invariants ok True, errors [], warnings []; depression fill raised 0 cells; drop analysis objective 120 cells (0.09514299927236262 km2), band (120, 600), giving 191 channel cells in 17 segments.

WW-DHSVM stream files (sha256): class 8e99b33ae558a55b, network 890970f541bfbc55, map 22d4dddb3edebefa; Channel.State written for 30 segments at 0.05 m initial depth.
