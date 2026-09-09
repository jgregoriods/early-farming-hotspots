# Early farming and hotspots of linguistic diversity in South America

Analysis code for Gregorio de Souza, Ruiz-Pérez & van Gijn, *Early Farming and
Hotspots of Linguistic Diversity in South America*.

## Running it

Requires R (developed under 4.6.1) and:

```r
install.packages(c("sf", "terra", "spdep", "ggplot2", "patchwork",
                   "glue", "here", "ggrepel", "RColorBrewer"))
```

Then, from the repository root:

```sh
Rscript main.R
```

Writes:

| Output | Contents |
|---|---|
| `plots/Fig1.tif`, `Fig2.tif`, `Fig3.tif` | Main figures |
| `plots/SI/S1_Fig.tif`, `S2_Fig.tif`, `S3_Fig.tif` | Supporting figures |
| `out/table_1.txt`, `out/table_2.txt` | Tables 1 and 2, at the reference cell area |
| `out/S1_Text.pdf` | S1 Text: the full regression sweep, all 20 fits |
| `out/S1_Table.pdf` | S1 Table: farming transition age in hotspots vs. other cells, all cell areas |
| `out/*.csv` | Source coverage, source agreement, richness correlations |

## Data

See `data/PROVENANCE.md` for the source and licence of every input.

## Licence

Code is MIT licensed (see `LICENSE`). The data are under the separate licences
recorded in `data/PROVENANCE.md`. Please cite the original sources and this
repository in case of reuse.