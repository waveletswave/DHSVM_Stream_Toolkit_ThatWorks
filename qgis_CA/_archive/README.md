# Archived scripts

These scripts are NOT called by `prep_dhsvm_inputs.py` and are retained only for historical reference.

## Files

### `roadaspect.py`
A misnamed file: despite its filename, the body is actually a "robust" rewrite of `rowcolmap` — i.e. it computes (Row, Col) for stream segments using the correct **top-left** origin. It was at one point a candidate for replacing `rowcolmap.py`, but the master pipeline (`prep_dhsvm_inputs.py`) ended up implementing Row/Col logic inline (see `_ensure_fields_rowcol_len` and `_write_stream_map`), so this standalone module is no longer imported anywhere.

### `rowcolmap.py`
The original ArcGIS-port version. Has a known bug: it computes the row index from the **bottom-left** origin of the DEM:

```python
row = int((centroid.y() - dem_raster.extent().yMinimum()) / cell_size_y)
```

DHSVM expects top-left origin (row 0 at the top, increasing downward). If this script is ever re-imported and used, it will produce a vertically-flipped Row index and break `stream.map.dat`. **Do not use.**

## Why archived rather than deleted

Keeping these files (a) preserves the design history that explains why the master pipeline embeds Row/Col logic inline rather than importing a helper module, and (b) makes the bottom-left bug a discoverable cautionary example for anyone maintaining DHSVM preprocessing in QGIS. Both files are inert as long as nothing imports them.
