#!/usr/bin/env python3
"""
check_slope_consumers.py -- does the slope conditioning touch any stream cell?

The conditioning filled the boundary-ring + interior-void cells (where the OLD
slope was NoData). The channel files (stream.class.dat, stream.network.dat) are
built by sampling slope ALONG the stream lines. So those files can only change
if a conditioned cell coincides with (or sits next to) a stream cell. This
checks that overlap. If it is zero, the channel files are unchanged and need
no comparison. Read-only.

    python check_slope_consumers.py            # auto-resolves Intermediate_GIS
"""
import argparse, sys
from pathlib import Path
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--old-slope', default=None)
    ap.add_argument('--dem', default=None)
    ap.add_argument('--stream', default=None)
    a = ap.parse_args()
    try:
        import rasterio
    except ModuleNotFoundError:
        sys.exit("[error] rasterio not importable. conda activate dhsvm_rs")

    WS = Path(__file__).resolve().parent.parent
    I = WS / "Intermediate_GIS"
    old_slope = a.old_slope or str(I / "stream_slope_OLD.tif")
    dem       = a.dem       or str(I / "elev_clipped.tif")
    stream    = a.stream    or str(I / "stream_raster.tif")

    def rd(fp):
        with rasterio.open(fp) as s:
            return s.read(1).astype(float), s.nodata

    try:
        elev, e_nd = rd(dem)
        sl, s_nd = rd(old_slope)
        st, t_nd = rd(stream)
    except Exception as ex:
        sys.exit(f"[error] could not read a raster ({ex}). "
                 f"Need stream_slope_OLD.tif, elev_clipped.tif, stream_raster.tif.")

    valid = np.isfinite(elev) & ((elev != e_nd) if e_nd is not None else True)
    cond  = valid & (~np.isfinite(sl) | ((sl == s_nd) if s_nd is not None else False))
    strm  = np.isfinite(st) & ((st != t_nd) if t_nd is not None else True) & (st != 0)

    # dilate streams by 1 cell (along-line sampling can land on a neighbour)
    dil = np.zeros_like(strm)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            sh = np.zeros_like(strm)
            ys = slice(max(0, dy), strm.shape[0] + min(0, dy))
            xs = slice(max(0, dx), strm.shape[1] + min(0, dx))
            yd = slice(max(0, -dy), strm.shape[0] + min(0, -dy))
            xd = slice(max(0, -dx), strm.shape[1] + min(0, -dx))
            sh[yd, xd] = strm[ys, xs]
            dil |= sh

    n_cond = int(cond.sum()); n_strm = int(strm.sum())
    ov_exact = int((cond & strm).sum())
    ov_dil   = int((cond & dil).sum())
    print(f"conditioned cells (the filled 318): {n_cond}")
    print(f"stream-raster cells               : {n_strm}")
    print(f"conditioned on a stream cell      : {ov_exact}")
    print(f"conditioned on stream cell +/-1   : {ov_dil}")
    print()
    if ov_dil == 0:
        print("-> stream.class.dat, stream.network.dat, and Channel.State are "
              "UNCHANGED by the conditioning (no stream cell, or its neighbour, "
              "was touched). No need to compare them.")
    else:
        print(f"-> up to {ov_dil} stream-adjacent cell(s) were conditioned; the "
              "channel slope sampling for those segments may differ slightly. "
              "diff stream.network.dat (and stream.class.dat) old-vs-new to see "
              "the magnitude.")


if __name__ == "__main__":
    main()
