# Stream Initiation Threshold

How the stream network threshold is set in the standalone pipeline, the methods
considered, and the reasoning. Recorded so the choice is not lost and can be
revisited.

## What the threshold controls

The pipeline extracts the stream network with r.stream.extract, which marks a
cell as channel once its upslope contributing area exceeds a threshold. This
threshold is the critical support area at the channel head: how much drainage
area must accumulate before a hillslope becomes a channel. It is a geomorphic
quantity, not a numerical knob. It sets drainage density, and therefore how many
channel segments DHSVM routes.

Only one threshold is involved. In this pipeline r.watershed runs without a
threshold; it produces only flow accumulation and drainage direction.
r.stream.extract holds the single threshold that defines the network.

## Current value and its problems

The threshold is 60 cells, selected by visual inspection in QGIS. Two problems.

First, a cell count is not a physical quantity. At 30 m, 60 cells is about
0.054 km square. At 10 m the same 60 cells is about 0.006 km square, nine times
smaller, so the identical setting produces a far denser network on a finer DEM.
The threshold has to be expressed as area, not cells, for the network to mean
the same thing across resolutions.

Second, visual selection is not objective or reproducible. Geomorphology has
established methods to set this threshold on a principled basis.

## Methods considered

Constant support area (O'Callaghan and Mark, 1984). All channel heads begin at
one contributing area A_c; the cell threshold is A_c divided by cell area. This
is what r.stream.extract implements. Simple and native to GRASS, but a single
A_c cannot capture the way drainage density varies across a basin.

Slope-area threshold (Montgomery and Dietrich, 1988, 1989, 1992). Channel
initiation depends on both contributing area and local slope, A times S to a
power above a constant, because steeper slopes channelize at smaller areas. This
inverse area-slope relationship holds in humid terrain, which describes the
Southern Appalachian study area. More physical, but r.stream.extract is
constant-area; a slope-area rule has to be built by hand in raster algebra.
Likely more than DHSVM needs here.

Constant drop analysis (Broscoe, 1959; Tarboton, Bras, and Rodriguez-Iturbe,
1991). An objective statistical method. The constant stream drop law holds that
the mean elevation drop along streams of each Strahler order is statistically
the same. The method sweeps a range of thresholds and tests whether the mean
drop of first-order streams differs from the higher orders by a t-test. The
smallest threshold for which the difference is not significant (absolute t below
2) gives the finest network consistent with the geomorphic law, and is taken as
the objective threshold. TauDEM implements this, and it is available in the QGIS
TauDEM tools. It uses only the DEM, depends on no regional data, and replaces
visual selection with a reproducible criterion.

Drainage density calibration. Adjust the threshold so the extracted network
matches an observed one, such as the USGS NHD high-resolution flowlines. Most
accurate against the real network, but limited to areas where such data exist.

## Decision

Two layers, kept separate.

Mechanism. The threshold is expressed as a physical support area A_c. The cell
count is computed from A_c and the DEM cell size at run time, so it scales with
resolution. A 30 m to 10 m change rescales the cell count automatically and the
network keeps the same physical drainage density. A_c is also a parameter that
transfers across watersheds, which suits a pipeline meant to run on any basin
from its boundary.

Value. A_c is set by constant drop analysis, not by eye. This is a one-time
analysis per watershed and DEM, run in QGIS or TauDEM, that yields the objective
A_c; the pipeline then consumes that area. NHD serves as a sanity check where
available, confirming the extracted channel heads sit near the observed ones.

The drop analysis is not part of the pipeline. It is the offline step that
determines the parameter. The pipeline takes the resulting area and computes the
cell threshold for whatever DEM resolution it is given.

## Resolution dependence

Rescaling A_c to a cell count assumes the channel-head support area does not
change with resolution, which is the constant-area assumption. A finer DEM
resolves smaller valleys and finer convergence, so the objective threshold need
not equal the rescaled coarse one. The rigorous check is to re-run drop analysis
on the finer DEM and compare its A_c to the rescaled value. Agreement means A_c
is stable across the scale and rescaling is sound. Disagreement means the finer
DEM reveals network structure the coarse one could not, which is itself a
finding. The quick-look figure shows the difference between the two networks
directly.

## References

Broscoe, A. J. (1959). Quantitative analysis of longitudinal stream profiles of
small watersheds. Office of Naval Research Project NR 389-042, Technical Report
18, Department of Geology, Columbia University.

Montgomery, D. R., and Dietrich, W. E. (1988). Where do channels begin? Nature,
336, 232-234.

Montgomery, D. R., and Dietrich, W. E. (1989). Source areas, drainage density,
and channel initiation. Water Resources Research, 25(8), 1907-1918.

Montgomery, D. R., and Dietrich, W. E. (1992). Channel initiation and the
problem of landscape scale. Science, 255, 826-830.

O'Callaghan, J. F., and Mark, D. M. (1984). The extraction of drainage networks
from digital elevation data. Computer Vision, Graphics, and Image Processing,
28, 323-344.

Tarboton, D. G., Bras, R. L., and Rodriguez-Iturbe, I. (1991). On the extraction
of channel networks from digital elevation data. Hydrological Processes, 5(1),
81-100.
