#!/usr/bin/env python3
"""Combined quick-look diagnostic for standalone DHSVM pipeline outputs.

One figure, six panels, so a user sees the whole pipeline at a glance:
  1 DEM with stream network and watershed boundary
  2 soil depth
  3 flow accumulation (log scale)
  4 slope
  5 flow direction (D8 classes)
  6 basin location (lon/lat, with a basemap if cartopy is available)

The four field panels and the flow-direction panel use DHSVM grid row/column
indices, the native coordinate system for debugging DHSVM inputs (cell
(row, col) maps straight to the .bin files). The location panel uses lon/lat,
which anyone can read to see where the basin sits. The watershed boundary is
the contour of the valid-data mask, so it needs no separate boundary file.

This is a diagnostic for catching gross geometric errors, not a validation
gate. The byte-level regression test remains the source of truth.
"""

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless: write files, no display server needed
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LightSource, LogNorm, ListedColormap, BoundaryNorm
import rasterio
from rasterio import features
from rasterio.transform import Affine
import geopandas as gpd
from shapely.geometry import shape

# The box_aspect + fixed-limits combination makes matplotlib's axes logger emit
# "Ignoring fixed y limits ..." on the cartopy panel. The figure renders
# correctly; silence just that logger for clean output.
logging.getLogger("matplotlib.axes._base").setLevel(logging.ERROR)


# ---- io ----

def read_raster(path):
    """Read single-band raster. Return (masked_array, affine_transform, crs)."""
    with rasterio.open(path) as src:
        return src.read(1, masked=True), src.transform, src.crs


def stream_polylines(stream_path, transform, ref_crs):
    """Read stream vector, return list of (cols, rows) polylines in pixel space.

    The shapefile may ship with an empty .prj (crs is None); it comes from the
    same GRASS region as the rasters, so it is plotted in their pixel frame.
    """
    streams = gpd.read_file(stream_path)
    if ref_crs is not None and streams.crs is not None and streams.crs != ref_crs:
        streams = streams.to_crs(ref_crs)
    inv = ~transform  # UTM (x, y) -> fractional (col, row), corner-referenced
    out = []
    for geom in streams.geometry:
        if geom is None:
            continue
        parts = [geom] if geom.geom_type == "LineString" else list(geom.geoms)
        for line in parts:
            xs, ys = line.xy
            cols, rows = [], []
            for x, y in zip(xs, ys):
                c, r = inv * (x, y)
                # inv maps a cell center to (col+0.5, row+0.5), but imshow
                # centers cell (row, col) at axis (col, row). Subtract 0.5 from
                # each so the stream line lands in the raster and boundary
                # frame; the boundary applies the same -0.5 via
                # Affine.translation(-0.5, -0.5) in draw_boundary.
                cols.append(c - 0.5)
                rows.append(r - 0.5)
            out.append((cols, rows))
    return out


def boundary_polylines(boundary_path, transform, ref_crs):
    """Read a watershed-boundary polygon, return its exterior and interior
    rings as (cols, rows) polylines in pixel space.

    Same frame convention as stream_polylines: inv maps a cell center to
    (col+0.5, row+0.5), so subtract 0.5 to share the raster frame. Drawing the
    actual polygon (instead of the valid-data mask contour) gives a smooth line
    and, as a side benefit, makes any polygon-vs-rasterized-mask mismatch
    visible: the smooth line will not hug the colored cells if the two disagree.
    """
    gdf = gpd.read_file(boundary_path)
    if ref_crs is not None and gdf.crs is not None and gdf.crs != ref_crs:
        gdf = gdf.to_crs(ref_crs)
    inv = ~transform

    def ring_to_pix(coords):
        cols, rows = [], []
        for x, y in coords:
            c, r = inv * (x, y)
            cols.append(c - 0.5)
            rows.append(r - 0.5)
        return cols, rows

    out = []
    for geom in gdf.geometry:
        if geom is None:
            continue
        polys = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
        for poly in polys:
            out.append(ring_to_pix(poly.exterior.coords))
            for ring in poly.interiors:
                out.append(ring_to_pix(ring.coords))
    return out


def basin_footprint_wgs84(arr, transform, crs):
    """Largest valid-data polygon, reprojected to WGS84 lon/lat."""
    if np.ma.is_masked(arr):
        m = (~np.ma.getmaskarray(arr)).astype("uint8")
    else:
        m = np.ones(arr.shape, "uint8")
    polys = [shape(g) for g, v in features.shapes(m, mask=m.astype(bool),
                                                  transform=transform) if v == 1]
    if not polys:
        return None
    biggest = max(polys, key=lambda p: p.area)
    return gpd.GeoDataFrame(geometry=[biggest], crs=crs).to_crs(4326)


# ---- drawing helpers ----

def hillshade(dem):
    ls = LightSource(azdeg=315, altdeg=45)
    return ls.hillshade(np.ma.filled(dem.astype(float), np.nan), vert_exag=2)


def data_mask(arr):
    if np.ma.is_masked(arr):
        return (~np.ma.getmaskarray(arr)).astype(float)
    return np.ones(arr.shape, float)


def data_bbox(arr):
    """(rmin, rmax, cmin, cmax) grid indices of the valid-data bounding box.

    Inclusive indices of the first and last rows and cols holding any valid
    cell. Panels crop their view to this box so the basin fills the frame
    instead of floating in the padded-grid nodata margin.
    """
    if np.ma.is_masked(arr):
        valid = ~np.ma.getmaskarray(arr)
    else:
        valid = np.ones(arr.shape, bool)
    rows = np.where(np.any(valid, axis=1))[0]
    cols = np.where(np.any(valid, axis=0))[0]
    if rows.size == 0 or cols.size == 0:        # all nodata: fall back to full grid
        return 0, arr.shape[0] - 1, 0, arr.shape[1] - 1
    return int(rows[0]), int(rows[-1]), int(cols[0]), int(cols[-1])


def _frame_aspect(bbox):
    """box_aspect (height / width) of the cropped view, buffer included.

    Shared by rowcol_axes and the location panel so all six frames keep
    identical proportions after cropping.
    """
    rmin, rmax, cmin, cmax = bbox
    m = max(rmax - rmin, cmax - cmin) * 0.04
    h = (rmax - rmin + 1) + 2 * m
    w = (cmax - cmin + 1) + 2 * m
    return h / w


def draw_boundary(ax, arr, boundary_lines=None):
    """Watershed boundary.

    If boundary_lines (rings from the boundary shapefile, via boundary_polylines)
    are supplied, draw them as a single smooth line. Otherwise fall back to the
    exact cell-edge polygon of the valid-data mask.

    Mask-contour fallback uses features.shapes rather than contour. contour
    draws the 0.5 isoline between cell centers, so the true edge along the
    outermost row/col gets clipped and the boundary looks open. The polygon edge
    closes correctly and hugs the real basin outline. Same method as the
    location footprint.
    """
    if boundary_lines is not None:
        for cols, rows in boundary_lines:
            ax.plot(cols, rows, color="black", linewidth=1.6, zorder=5)
        return
    m = data_mask(arr).astype("uint8")
    tr = Affine.translation(-0.5, -0.5)  # cell corners -> imshow pixel coords
    for geom, val in features.shapes(m, mask=m.astype(bool), transform=tr):
        if val != 1:
            continue
        poly = shape(geom)
        xs, ys = poly.exterior.xy
        ax.plot(xs, ys, color="black", linewidth=1.6, zorder=5)
        for ring in poly.interiors:
            xr, yr = ring.xy
            ax.plot(xr, yr, color="black", linewidth=1.0, zorder=5)


def draw_streams(ax, polylines):
    for cols, rows in polylines:
        ax.plot(cols, rows, color="#0d3b66", linewidth=2.0, zorder=4,
                solid_capstyle="round",
                path_effects=[pe.Stroke(linewidth=3.0, foreground="white"),
                              pe.Normal()])


def rowcol_axes(ax, bbox, grid=False):
    rmin, rmax, cmin, cmax = bbox
    ax.set_xlabel("Column index")
    ax.set_ylabel("Row index")
    # Crop the view to the basin's valid-data bbox, but keep the tick *labels*
    # as true DHSVM grid indices: stream_network.py writes col = int((x - xmin)
    # / px), so cell (row, col) maps straight to the .bin files. (0, 0) sits
    # off-frame whenever the basin does not reach the padded-grid corner, the
    # common case. Do not re-base labels to 0; that would break the
    # cell-to-.bin correspondence these panels exist to check.
    # Small buffer so the basin does not sit flush against the frame.
    m = max(rmax - rmin, cmax - cmin) * 0.04
    x0, x1 = cmin - 0.5 - m, cmax + 0.5 + m
    y_top, y_bot = rmin - 0.5 - m, rmax + 0.5 + m
    ax.set_xlim(x0, x1)
    ax.set_ylim(y_bot, y_top)  # descending: imshow origin is upper
    # Drive the frame shape with box_aspect (not aspect='equal'); this keeps
    # cells square, avoids the limit-vs-aspect warning, and makes every panel
    # share the same frame proportions.
    ax.set_aspect("auto")
    ax.set_box_aspect((y_bot - y_top) / (x1 - x0))
    if grid:
        ax.set_xticks(np.arange(cmin - 0.5, cmax + 1, 1), minor=True)
        ax.set_yticks(np.arange(rmin - 0.5, rmax + 1, 1), minor=True)
        ax.grid(which="minor", color="0.6", linewidth=0.25, alpha=0.4)


# ---- panels ----

def panel_dem(ax, fig, dem, polylines, bbox, boundary_lines=None):
    nrows, ncols = dem.shape
    valid = int(np.ma.count(dem))          # basin cells; nodata margin excluded
    total = nrows * ncols
    ax.imshow(hillshade(dem), cmap="gray", origin="upper", zorder=0)
    im = ax.imshow(dem, cmap="terrain", origin="upper", alpha=0.85, zorder=1)
    fig.colorbar(im, ax=ax, shrink=0.85, label="Elevation (m)")
    draw_streams(ax, polylines)
    draw_boundary(ax, dem, boundary_lines)
    ax.set_title(f"DEM + stream network    {nrows} rows x {ncols} cols\n"
                 f"basin: {valid} of {total} cells ({100 * valid / total:.0f}% data)")
    rowcol_axes(ax, bbox, grid=True)


def panel_continuous(ax, fig, arr, label, cmap, bbox, boundary_lines=None):
    im = ax.imshow(arr, cmap=cmap, origin="upper")
    fig.colorbar(im, ax=ax, shrink=0.85, label=label)
    draw_boundary(ax, arr, boundary_lines)
    vmin, vmax = float(np.ma.min(arr)), float(np.ma.max(arr))
    ax.set_title(f"{label.split(' (')[0]}    min {vmin:.3g}  max {vmax:.3g}")
    rowcol_axes(ax, bbox)


def panel_flow_accum(ax, fig, acc, bbox, boundary_lines=None):
    mag = np.ma.masked_less_equal(np.ma.abs(acc), 0)
    im = ax.imshow(mag, cmap="cubehelix_r", origin="upper",
                   norm=LogNorm(vmin=1, vmax=float(np.ma.max(mag))))
    fig.colorbar(im, ax=ax, shrink=0.85, label="Flow accumulation (cells, log)")
    draw_boundary(ax, acc, boundary_lines)
    ax.set_title("Flow accumulation")
    rowcol_axes(ax, bbox)


def panel_flow_dir(ax, fig, arr, bbox, boundary_lines=None):
    # GRASS r.watershed drainage: 1-8 are D8 directions; negatives mark cells
    # that drain out of the region (basin edge). Collapse all negatives into a
    # single "out" class in neutral gray so the eight real directions read
    # cleanly and the edge cells are self-explanatory.
    a = np.ma.where(arr < 0, -1, arr).astype(float)
    vals = np.unique(a.compressed())          # e.g. [-1, 1, 2, ... 8]
    pos = [v for v in vals if v >= 0]
    posmap = {v: plt.cm.tab10(i % 10) for i, v in enumerate(pos)}
    # "out" gets bright cyan, deliberately outside the tab10 direction colors so
    # it cannot collide with any of the eight D8 directions.
    colors = [(0.0, 0.9, 1.0, 1.0) if v < 0 else posmap[v] for v in vals]
    cmap = ListedColormap(colors)
    bounds = np.concatenate([vals - 0.5, [vals[-1] + 0.5]])
    norm = BoundaryNorm(bounds, cmap.N)
    im = ax.imshow(a, cmap=cmap, norm=norm, origin="upper")
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, ticks=vals,
                        label="Flow direction (D8)")
    cbar.ax.set_yticklabels(["out" if v < 0 else f"{int(v)}" for v in vals])
    draw_boundary(ax, a, boundary_lines)
    ax.set_title("Flow direction (D8)")
    rowcol_axes(ax, bbox)


def _load_basemap_features():
    """Return cartopy Natural Earth features if all are reachable, else None.

    cartopy downloads feature data lazily at draw time, so reachability is
    verified up front by forcing each feature's geometries to load. If any
    download fails (e.g. offline cluster node), the basemap is skipped.
    """
    try:
        import cartopy.feature as cfeature
    except Exception:
        return None
    feats = [
        cfeature.NaturalEarthFeature("physical", "land", "50m",
                                     facecolor="#efeadb", edgecolor="none"),
        cfeature.NaturalEarthFeature("physical", "ocean", "50m",
                                     facecolor="#cfe8f3", edgecolor="none"),
        cfeature.NaturalEarthFeature("cultural", "admin_1_states_provinces_lines",
                                     "50m", facecolor="none", edgecolor="0.5"),
        cfeature.NaturalEarthFeature("physical", "coastline", "50m",
                                     facecolor="none", edgecolor="0.3"),
    ]
    try:
        for f in feats:
            list(f.geometries())  # force download now; raises if unreachable
        return feats
    except Exception:
        return None


def panel_location(fig, pos, footprint, bbox):
    """Watershed footprint on a lon/lat map. The view range adapts to the
    watershed's own extent, so a small one zooms in and a large one pulls back.
    Adds a cartopy basemap when reachable, otherwise footprint + graticule."""
    minx, miny, maxx, maxy = footprint.total_bounds
    lon_c, lat_c = (minx + maxx) / 2.0, (miny + maxy) / 2.0
    span = max(maxx - minx, maxy - miny, 1e-3)  # extent in degrees
    try:
        import cartopy.crs as ccrs
        have_cartopy = True
    except Exception:
        have_cartopy = False

    if have_cartopy:
        ax = fig.add_subplot(*pos, projection=ccrs.PlateCarree())
        feats = _load_basemap_features()
        if feats:
            for f in feats:
                ax.add_feature(f, zorder=0)
            half = span / 2.0 + max(span * 2.0, 0.4)   # room for geographic context
        else:
            half = span / 2.0 + max(span * 0.5, 0.02)  # focus on the footprint
        footprint.plot(ax=ax, facecolor="#d33", edgecolor="black",
                       alpha=0.85, zorder=3, transform=ccrs.PlateCarree())
        ax.set_extent([lon_c - half, lon_c + half, lat_c - half, lat_c + half],
                      crs=ccrs.PlateCarree())
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.7")
        gl.top_labels = gl.right_labels = False
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    else:
        ax = fig.add_subplot(*pos)
        half = span / 2.0 + max(span * 0.5, 0.02)
        footprint.plot(ax=ax, facecolor="#d33", edgecolor="black", alpha=0.85)
        ax.set_xlim(lon_c - half, lon_c + half)
        ax.set_ylim(lat_c - half, lat_c + half)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        # Equidistant lon/lat inside the frame: stretch the lat axis by
        # 1/cos(lat) so a degree of longitude and latitude have equal ground
        # length, keeping the footprint undistorted.
        ax.set_aspect(1.0 / math.cos(math.radians(lat_c)))
    # Outer frame matches the data panels (same cropped-basin box_aspect) so
    # the panel is no larger than the others in either dimension; the
    # equidistant map sits inside, with the leftover space as margin.
    ax.set_box_aspect(_frame_aspect(bbox))
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.8)
        s.set_edgecolor("black")
    ax.set_title(f"Watershed location    {lat_c:.3f} N, {abs(lon_c):.3f} W")
    return ax


# ---- driver ----

def main():
    p = argparse.ArgumentParser(
        description="Combined quick-look for DHSVM pipeline outputs")
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--out", type=Path, default=None,
                   help="output PNG (default: run-dir/quicklook/quicklook.png)")
    p.add_argument("--dem", default="elev_clipped.tif")
    p.add_argument("--stream", default="streamfile.shp")
    p.add_argument("--boundary", default=None,
                   help="watershed-boundary polygon (shp). If given, drawn as a "
                        "smooth line; otherwise the valid-data mask contour is used")
    p.add_argument("--soil-depth", default="soildepth.tif")
    p.add_argument("--flow-accum", default="flow_acc.tif")
    p.add_argument("--flow-dir", default="flow_dir.tif")
    p.add_argument("--slope", default="slope_filled.tif")
    p.add_argument("--slope-unit", default="deg",
                   help="label for the slope colorbar: deg, %%, or fraction")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--title", default="DHSVM pipeline quick-look")
    args = p.parse_args()

    rd = args.run_dir
    dem, transform, crs = read_raster(rd / args.dem)
    bbox = data_bbox(dem)  # crop every panel to the basin, not the padded grid
    polylines = stream_polylines(rd / args.stream, transform, crs)
    boundary_lines = (boundary_polylines(Path(args.boundary), transform, crs)
                      if args.boundary else None)
    footprint = basin_footprint_wgs84(dem, transform, crs)

    fig = plt.figure(figsize=(19, 11.5), constrained_layout=True)
    fig.suptitle(args.title, fontsize=17, fontweight="bold")

    # Top row: DEM+stream, soil depth, slope. Bottom row: flow accumulation,
    # flow direction, basin location (flow_acc and flow_dir sit adjacent).
    panel_dem(fig.add_subplot(2, 3, 1), fig, dem, polylines, bbox, boundary_lines)

    soil, _, _ = read_raster(rd / args.soil_depth)
    panel_continuous(fig.add_subplot(2, 3, 2), fig, soil,
                     "Soil depth (m)", "YlOrBr", bbox, boundary_lines)

    slope, _, _ = read_raster(rd / args.slope)
    panel_continuous(fig.add_subplot(2, 3, 3), fig, slope,
                     f"Slope ({args.slope_unit})", "magma", bbox, boundary_lines)

    acc, _, _ = read_raster(rd / args.flow_accum)
    panel_flow_accum(fig.add_subplot(2, 3, 4), fig, acc, bbox, boundary_lines)

    fdir, _, _ = read_raster(rd / args.flow_dir)
    panel_flow_dir(fig.add_subplot(2, 3, 5), fig, fdir, bbox, boundary_lines)

    if footprint is not None:
        panel_location(fig, (2, 3, 6), footprint, bbox)

    out = args.out or (rd / "quicklook" / "quicklook.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}  ({args.dpi} dpi)")


if __name__ == "__main__":
    main()
