# Data provenance and licences

Every input to the analysis is third-party. This file records their provenance and license.
The analysis runs without the missing files, since `grids.rds` in the repository root
caches the per-cell covariates (`npp`, `bio_15`, `tri`, `river_density`) for all cell
areas. The river network and environmental rasters are needed only to rebuild that
cache from scratch.

## Included in the repository

| File | Source | Licence |
|---|---|---|
| `languoid.csv` | Glottolog 5.3 (Hammarström, Forkel, Haspelmath & Bank 2026), <https://glottolog.org> | CC BY 4.0 |
| `shp/glottography.*` | Glottography Consortium 2025, `asher2007world` v1.0.0, <https://doi.org/10.5281/zenodo.15287258>, digitised from Asher, R. E. & Christopher J. Moseley (eds.) 2007, *Atlas of the World's Languages*, 2nd edn, Routledge, and described in Ranacher et al. 2025, *Scientific Data* 12:1466, <https://doi.org/10.1038/s41597-025-05828-6> | CC BY 4.0 |
| `projected_shp/south_america.*` | Natural Earth, Admin 0 countries, dissolved to a single South America polygon and reprojected to EPSG:6933. <https://www.naturalearthdata.com> | Public domain |
| `c14dates.csv` | South American subset of the global food-production transition database in Gregorio de Souza, Ruiz-Pérez, Ruiz-Giralt, Lancelotti & Madella 2025, *Scientific Reports* 15:8301. <https://doi.org/10.1038/s41598-025-92782-3> | CC BY 4.0 |

## Not included

These are excluded because their licences do not permit redistribution.

| File | Source |
|---|---|
| `projected_shp/rivers.*` | MERIT Hydro–Vector, Lin, Pan, Wood, Yamazaki & Allen 2021, *Scientific Data* 8:28. <https://doi.org/10.1038/s41597-021-00819-9> |
| `projected_rasters/bio_15.tif` | WorldClim 2.1, precipitation seasonality (BIO15), Fick & Hijmans 2017. <https://worldclim.org> |
| `projected_rasters/tri.tif` | Terrain ruggedness index, derived from WorldClim 2.1 elevation with `gdaldem TRI` (GDAL 3.8.4). |
| `projected_rasters/npp.tif` | Net primary productivity, Global Agro-Ecological Zones (GAEZ), Fischer et al. 2008, IIASA/FAO. <https://gaez.fao.org> |

## Rebuilding grids.rds

Only needed if you want to regenerate the per-cell covariates rather than use
the cached ones.

1. Download the four files above.
2. Reproject each to **EPSG:6933** (WGS 84 / NSIDC EASE-Grid 2.0 Global) and
   place them at the paths in the table.
   The river layer needs an `ID` field; the rasters must each carry a single
   band named `bio_15`, `npp` or `tri` respectively, since `extract_raster_values()`
   takes the grid column name from the layer name.
3. Delete `grids.rds` and run `Rscript main.R`. The grids are rebuilt and
   re-cached.

## Citing the data

Please cite the original sources and this repository in case of reuse.