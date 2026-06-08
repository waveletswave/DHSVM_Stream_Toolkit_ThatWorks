#!/usr/bin/env python3
"""
diagnose_slope_gaps.py -- characterize the interior slope-NoData cells.

Read-only. Answers: are the interior missing cells one patch or scattered?
do they ring interior DEM voids (NoData elev)? do they sit on stream cells?
are the slope and DEM grids aligned? Nothing is written.

    python diagnose_slope_gaps.py            # auto-resolves Intermediate_GIS
    python diagnose_slope_gaps.py --slope ../Intermediate_GIS/stream_slope.tif \
                                  --dem   ../Intermediate_GIS/elev_clipped.tif
"""
import argparse, sys
from pathlib import Path
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slope', default=None)
    ap.add_argument('--dem',   default=None)
    ap.add_argument('--stream', default=None)
    a = ap.parse_args()
    try:
        import rasterio
    except ModuleNotFoundError:
        sys.exit("[error] rasterio not importable. conda activate dhsvm_rs")

    WS = Path(__file__).resolve().parent.parent
    slope_fp = a.slope or str(WS/"Intermediate_GIS"/"stream_slope.tif")
    dem_fp   = a.dem   or str(WS/"Intermediate_GIS"/"elev_clipped.tif")
    stream_fp = a.stream or str(WS/"Intermediate_GIS"/"stream_raster.tif")

    with rasterio.open(slope_fp) as s:
        slope = s.read(1).astype(float); s_nd = s.nodata
        s_tf, s_bounds = s.transform, s.bounds
    with rasterio.open(dem_fp) as d:
        elev = d.read(1).astype(float); e_nd = d.nodata
        d_tf, d_bounds = d.transform, d.bounds

    # grid alignment
    print("grid alignment (slope vs dem):")
    print(f"  transform equal: {s_tf == d_tf}")
    print(f"  bounds equal:    {s_bounds == d_bounds}")
    if s_tf != d_tf or s_bounds != d_bounds:
        print(f"  slope transform: {s_tf}")
        print(f"  dem   transform: {d_tf}")
        print(f"  slope bounds: {s_bounds}")
        print(f"  dem   bounds: {d_bounds}")
        print("  -> grids differ; some interior 'gaps' may be misalignment, not real NoData.\n")
    else:
        print("  -> grids identical; gaps are real slope NoData.\n")

    valid_dem = np.isfinite(elev) & ((elev != e_nd) if e_nd is not None else True)
    slope_valid = np.isfinite(slope) & ((slope != s_nd) if s_nd is not None else True)
    missing = valid_dem & ~slope_valid

    outside = ~valid_dem
    touch = np.zeros_like(missing)
    for dy in (-1,0,1):
        for dx in (-1,0,1):
            if dy==0 and dx==0: continue
            sh = np.zeros_like(outside)
            ys=slice(max(0,dy),outside.shape[0]+min(0,dy)); xs=slice(max(0,dx),outside.shape[1]+min(0,dx))
            yd=slice(max(0,-dy),outside.shape[0]+min(0,-dy)); xd=slice(max(0,-dx),outside.shape[1]+min(0,-dx))
            sh[yd,xd]=outside[ys,xs]; touch|=sh
    interior = missing & ~touch
    print(f"total slope-NoData (DEM-valid): {int(missing.sum())}  "
          f"| boundary {int((missing&touch).sum())}  interior {int(interior.sum())}")

    ij = np.argwhere(interior)
    if len(ij)==0:
        print("no interior gaps."); return

    # connected components (8-conn) over interior cells, simple flood fill
    seen = np.zeros_like(interior); comps=[]
    iset = {tuple(p) for p in ij}
    for p in ij:
        p=tuple(p)
        if seen[p]: continue
        stack=[p]; seen[p]=True; size=0
        while stack:
            y,x=stack.pop(); size+=1
            for dy in(-1,0,1):
                for dx in(-1,0,1):
                    q=(y+dy,x+dx)
                    if q in iset and not seen[q]:
                        seen[q]=True; stack.append(q)
        comps.append(size)
    comps.sort(reverse=True)
    print(f"interior cells form {len(comps)} connected patch(es); sizes: {comps}")

    # do interior cells ring an elev-void? count NoData-elev neighbours
    nodata_elev = ~valid_dem
    ring_void = 0
    for (y,x) in ij:
        nb=False
        for dy in(-1,0,1):
            for dx in(-1,0,1):
                if dy==0 and dx==0: continue
                yy,xx=y+dy,x+dx
                if 0<=yy<elev.shape[0] and 0<=xx<elev.shape[1] and nodata_elev[yy,xx]:
                    nb=True
        ring_void += int(nb)
    print(f"interior cells adjacent to a NoData-elev cell: {ring_void} "
          f"(if >0, they ring DEM voids)")

    # overlap with stream raster
    try:
        with rasterio.open(stream_fp) as st:
            stream = st.read(1)
            stv = np.isfinite(stream) & (stream != (st.nodata if st.nodata is not None else -9999))
            stv &= (stream != 0)
        n_on_stream = int((interior & stv).sum())
        print(f"interior cells on a stream-raster cell: {n_on_stream} of {len(ij)}")
    except Exception as ex:
        print(f"(stream raster not checked: {ex})")

    # elev at interior cells vs basin
    print(f"elev at interior cells: mean {elev[interior].mean():.1f}  "
          f"min {elev[interior].min():.1f}  max {elev[interior].max():.1f}  "
          f"| basin elev mean {elev[valid_dem].mean():.1f}")


if __name__ == "__main__":
    main()
