# standalone_CA network stage before Tier E

`vector_attrs.py` and `stream_network.py` as they were at commit 1aaee3d.
They built the stream network from r.to.vect lines re-linked by geometry,
taking the endpoint with the larger raw flow-accumulation value as
downstream. r.watershed writes negative accumulation on a basin-clipped
DEM, so the rule reversed every segment; the sink merge then made a
headwater the outlet (Tier E audit, docs/audit/). Kept for the record;
the replacements are `segments_from_raster.py` and `network_files.py`.
