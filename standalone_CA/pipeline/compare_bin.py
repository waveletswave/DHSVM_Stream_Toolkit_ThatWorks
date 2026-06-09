"""Byte-level comparison for raw .bin files (no header, so sha256 is exact)."""
import sys
sys.path.insert(0, ".")
import hashlib
import numpy as np
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_bin(out_path, ref_path, dtype):
    out_path, ref_path = Path(out_path), Path(ref_path)
    so, sr = out_path.stat().st_size, ref_path.stat().st_size
    ho, hr = sha256(out_path), sha256(ref_path)
    identical = ho == hr
    flag = "OK " if identical else "XX "
    line = (f"{flag}{out_path.name:32s} size {so}/{sr} "
            f"{'identical' if identical else 'DIFFER'}")
    print(line)
    if not identical:
        if so == sr:
            a = np.fromfile(out_path, dtype=dtype)
            b = np.fromfile(ref_path, dtype=dtype)
            d = np.abs(a.astype("float64") - b.astype("float64"))
            print(f"     same size; n differing values: {int((d != 0).sum())} "
                  f"/ {a.size}; max abs diff: {float(d.max())}")
        else:
            print(f"     SIZE MISMATCH — {so - sr} bytes; check dtype/shape/order")
    return identical


if __name__ == "__main__":
    from paths import OUT, REF
    BIN_OUT = OUT / "bins"
    REF_BIN = REF / "DHSVM_input_binaries"

    checks = [
        ("dem.bin", "float32"),
        ("mask.bin", "int8"),
        ("soil.bin", "int8"),
        ("veg.bin", "int8"),
        ("soildepth_uniform_2.0m.bin", "float32"),
        ("soildepth_uniform_2.5m.bin", "float32"),
        ("soildepth_uniform_3.0m.bin", "float32"),
        ("soildepth_uniform_3.5m.bin", "float32"),
        ("soildepth_uniform_4.0m.bin", "float32"),
    ]
    all_ok = True
    print("=== bins byte comparison vs qgis_CA ===")
    for name, dt in checks:
        ok = compare_bin(BIN_OUT / name, REF_BIN / name, dt)
        all_ok = all_ok and ok
    print(f"\nBINS ALL byte-identical: {all_ok}")
