"""Validation helpers: byte-level and per-pixel raster comparison against qgis_CA."""
import hashlib
import numpy as np
import rasterio


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_rasters(out_path, ref_path, atol=0.0):
    """Compare a produced raster against a reference. Returns True if they match.

    Reports, in order: byte-identity, then metadata (shape/dtype/nodata/transform/
    crs), then per-pixel agreement including where nodata sits. The per-pixel
    breakdown is what tells us *how* a mismatch fails (boundary ring, nodata value,
    or grid offset), which points at the specific clip parameter to change.
    """
    print(f"OUT: {out_path}")
    print(f"REF: {ref_path}")

    h_out, h_ref = sha256(out_path), sha256(ref_path)
    byte_identical = h_out == h_ref
    print(f"  byte-identical: {byte_identical}")
    print(f"    sha256 out: {h_out}")
    print(f"    sha256 ref: {h_ref}")

    with rasterio.open(out_path) as o, rasterio.open(ref_path) as r:
        ao, ar = o.read(1), r.read(1)
        print(f"  shape   out/ref: {o.shape} / {r.shape}")
        print(f"  dtype   out/ref: {ao.dtype} / {ar.dtype}")
        print(f"  nodata  out/ref: {o.nodata} / {r.nodata}")
        print(f"  crs     out/ref: {o.crs} / {r.crs}")
        to = tuple(round(v, 6) for v in o.transform)[:6]
        tr = tuple(round(v, 6) for v in r.transform)[:6]
        print(f"  transform match: {to == tr}")
        if to != tr:
            print(f"    out: {to}")
            print(f"    ref: {tr}")

        if o.shape != r.shape:
            print("  -> shape differs; cannot do per-pixel comparison")
            return False

        # nodata masks (handle both an explicit nodata tag and NaN)
        def nodata_mask(arr, nod):
            m = np.isnan(arr)
            if nod is not None and not np.isnan(nod):
                m = m | (arr == nod)
            return m

        mo, mr = nodata_mask(ao, o.nodata), nodata_mask(ar, r.nodata)
        print(f"  nodata count out/ref: {int(mo.sum())} / {int(mr.sum())}")
        nodata_layout_match = np.array_equal(mo, mr)
        print(f"  nodata layout identical: {nodata_layout_match}")

        both_data = (~mo) & (~mr)
        if both_data.any():
            diff = np.abs(ao[both_data].astype("float64") - ar[both_data].astype("float64"))
            print(f"  data pixels compared: {int(both_data.sum())}")
            print(f"  max abs diff (data): {float(diff.max())}")
            print(f"  n pixels diff > {atol}: {int((diff > atol).sum())}")

        # pixels where one is nodata and the other is not -> boundary/all_touched tell
        only_out = mo & (~mr)
        only_ref = mr & (~mo)
        if only_out.any() or only_ref.any():
            print(f"  nodata-only-in-out: {int(only_out.sum())} | "
                  f"nodata-only-in-ref: {int(only_ref.sum())}")
            print("    (nonzero here => boundary ring mismatch, check ALL_TOUCHED)")

        matched = byte_identical or (
            nodata_layout_match and both_data.any()
            and float(np.abs(ao[both_data].astype('float64')
                             - ar[both_data].astype('float64')).max()) <= atol
        )
        return matched
