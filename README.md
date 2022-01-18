## R code for different methods of measuring extent of occurrence (EOO) and area of occupancy (AOO)

A set of R scripts for various methods of measuring geographic range size for IUCN Red List assessments (a different script for each method) used in the publication (in review):

**The effect of sampling effort and methodology on range size estimates of poorly-recorded species for IUCN Red List assessments.**

*Charles J. Marsh, Mindy M. Syfert, Elina Aletrari, Yoni Gavish, William E. Kunin, Neil Brummitt*

Scripts can be found in the EOO and AOO folders, and uses data collected for the pteridophyte IUCN Sampled Red List Index (SLRI).  All methods generate range size in km<sup>2</sup>. The scripts use an incomplete subset of the data for South and Central American endemics and should not be used for generating real IUCN Red List assessments or other products.

Links for each method are:

### EOO measures

- [MCP](https://github.com/charliem2003/AOOvEOO/blob/master/EOO_measures/MCP.R) (IUCN-recommended method; IUCN 2019)
- [Alpha hull](https://github.com/charliem2003/AOOvEOO/blob/master/EOO_measures/alpha_hull.R) (explores alpha values between 1-6; Burgman and Fox, 2003)
- [LoCoH](https://github.com/charliem2003/AOOvEOO/blob/master/EOO_measures/LoCoH.R) using the 3 potential methods (Getz et al., 2007; Getz and Wilmers, 2004)

### AOO measures

- [Grid overlay](https://github.com/charliem2003/AOOvEOO/blob/master/AOO_measures/grid_overlay.R) at multiple grain sizes (IUCN-recommended method at 2×2 km; IUCN 2019)
- [Variables grain](https://github.com/charliem2003/AOOvEOO/blob/master/AOO_measures/variable_grain.R) (Willis et al., 2003)
- [Circular buffers](https://github.com/charliem2003/AOOvEOO/blob/master/AOO_measures/circular_buffers.R) (Breiner and Bergamini, 2018)
- [CMC](https://github.com/charliem2003/AOOvEOO/blob/master/AOO_measures/CMC.R) (Cartographic Method of Conglomerates; Hernández and Navarro, 2007)
- [Occupancy downscaling](https://github.com/charliem2003/AOOvEOO/blob/master/AOO_measures/occupancy_downscaling.R) (Groom et al., 2018; Marsh et al., 2019) using the downscale package (Marsh et al., 2018)

**References where you can find the original description of the methods**

Breiner, F.T., Bergamini, A., **2018**. Improving the estimation of area of occupancy for IUCN Red List assessments by using a circular buffer approach. *Biodiversity and Conservation* 27, 2443–2448. https://doi.org/10.1007/s10531-018-1555-5

Burgman, M.A., Fox, J.C., **2003**. Bias in species range estimates from minimum convex polygons: implications for conservation and options for improved planning. *Animal Conservation* 6, 19–28. https://doi.org/10.1017/S1367943003003044

Getz, W.M., Fortmann-Roe, S., Cross, P.C., Lyons, A.J., Ryan, S.J., Wilmers, C.C., **2007**. LoCoH: nonparameteric kernel methods for constructing home ranges and utilization distributions. *PLoS ONE* 2, e207. https://doi.org/10.1371/journal.pone.0000207

Getz, W.M., Wilmers, C.C., **2004**. A local nearest-neighbor convex-hull construction of home ranges and utilization distributions. *Ecography* 27, 489–505. https://doi.org/10.1111/j.0906-7590.2004.03835.x

Groom, Q.J., Marsh, C.J., Gavish, Y., Kunin, W.E., **2018**. How to predict fine resolution occupancy from coarse occupancy data. *Methods in Ecology and Evolution* 9, 2273–2284. https://doi.org/10.1111/2041-210X.13078

Hernández, H.M., Navarro, M., **2007**. A new method to estimate areas of occupancy using herbarium data. *Biodiversity and Conservation* 16, 2457–2470. https://doi.org/10.1007/s10531-006-9134-6

IUCN, **2019**. Guidelines for using the IUCN Red List categories and criteria (No. Version 14). Prepared by the Standards and Petitions Subcommittee.

Marsh, C.J., Barwell, L.J., Gavish, Y., Kunin, W.E., **2018**. downscale: an R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes. *Journal of Statistical Software* 86, 1–20. https://doi.org/10.18637/jss.v086.c03

Marsh, C.J., Gavish, Y., Kunin, W.E., Brummitt, N.A., **2019**. Mind the gap: Can downscaling Area of Occupancy overcome sampling gaps when assessing IUCN Red List status? *Diversity and Distributions* 025, 1832–1845. https://doi.org/10.1111/ddi.12983

Willis, F., Moat, J., Paton, A., **2003**. Deﬁning a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. *Biodiversity and Conservation* 12, 1537–1552.
