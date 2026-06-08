# QGIS Scripts (Reference Implementation)

## ⚠️ Important note on rowcolmap

This directory contains TWO scripts that both compute (Row, Col) for stream segments — they were originally adapted from PNNL ArcGIS code at different times and the older one had a coordinate-origin bug:

| File                                                            | Origin      | Status                                                                       |
| --------------------------------------------------------------- | ----------- | ---------------------------------------------------------------------------- |
| `rowcolmap_robust.py` (formerly `roadaspect.py`)            | top-left    | ✅ Use this                                                                  |
| `rowcolmap_legacy_bottom_left.py` (formerly `rowcolmap.py`) | bottom-left | ❌ Bug: produces vertically-flipped Row indices, breaks DHSVM stream.map.dat |

The **standalone** equivalent (`../standalone/rowcolmap.py`) matches `rowcolmap_robust.py` byte-for-byte.

If you have legacy DHSVM runs whose stream.map.dat was generated with the bottom-left version, the basin will appear vertically flipped from DHSVM's internal viewpoint. Re-run preprocessing with the robust version.
