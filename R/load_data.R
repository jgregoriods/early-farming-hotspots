suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(here)
})

# Base vectors for South America and rivers

south_america <- vect(here("data/projected_shp/south_america.shp"))

require_data <- function(path, source) {
    if (!file.exists(path))
        stop("Missing restricted input: ", path, "\n  Source: ", source,
             "\n  This file is not redistributable and is not in the repository.",
             "\n  It is only needed to rebuild grids.rds; see data/PROVENANCE.md.",
             call. = FALSE)
    path
}

get_rivers <- function()
    vect(require_data(here("data/projected_shp/rivers.shp"),
                      "MERIT Hydro-Vector, Lin et al. 2021 (CC BY-NC)"))

get_rasters <- function() {
    dir <- here("data/projected_rasters")
    files <- file.path(dir, c("bio_15.tif", "npp.tif", "tri.tif"))
    for (f in files)
        require_data(f, "WorldClim 2.1 (bio_15, tri) / GAEZ (npp)")
    lapply(files, rast)
}

# Remove unattested, unclassifiable, bookkeeping ...

EXCLUDED_FAMILY_IDS <- c(
    "unat1236", "uncl1493", "book1242", "sign1238",
    "pidg1258", "mixe1287", "spee1234", "arti1236"
)

# Remove non-indigenous languages

INDO_EUROPEAN_FAMILY_ID <- "indo1319"
NON_INDIGENOUS_GLOTTOCODES <- c("cari1276")

languoid <- read.csv(here("data/languoid.csv"))
languoid <- languoid[languoid$level == "language" & !is.na(languoid$latitude) & !is.na(languoid$longitude), ]
glottolog <- st_as_sf(languoid, coords = c("longitude", "latitude"), crs = 4326) |> st_transform(crs(south_america))
glottolog <- glottolog[lengths(st_intersects(glottolog, st_as_sf(south_america))) > 0, ]

glottolog <- glottolog[
    !glottolog$family_id %in% EXCLUDED_FAMILY_IDS &
        !glottolog$family_id %in% INDO_EUROPEAN_FAMILY_ID &
        !glottolog$id %in% NON_INDIGENOUS_GLOTTOCODES,
]

names(glottolog)[names(glottolog) == "id"] <- "glottocode"
names(glottolog)[names(glottolog) == "family_id"] <- "family"

glottolog$isolate <- glottolog$family == ""

# Fill isolates' family column with their own glottocode

glottolog$family[glottolog$family == ""] <- glottolog$glottocode[glottolog$family == ""]

glottography_raw <- st_read(here("data/shp/glottography.shp"), quiet = TRUE)
glottography <- aggregate(
    glottography_raw[c("name", "cldf_langu", "family_id")],
    by = list(id = glottography_raw$id),
    FUN = function(x) x[1]
) |> st_transform(crs(south_america))

glottography <- glottography[lengths(st_intersects(glottography, st_as_sf(south_america))) > 0, ]
glottography <- glottography[!(glottography$name == "No label" & is.na(glottography$cldf_langu)), ]

# Fill in missing family_id in Glottography from Glottolog

languoid_full <- read.csv(here("data/languoid.csv"))
backfill_family_id <- languoid_full$family_id[match(glottography$cldf_langu, languoid_full$id)]
needs_backfill <- is.na(glottography$family_id) & !is.na(backfill_family_id) & backfill_family_id != ""
glottography$family_id[needs_backfill] <- backfill_family_id[needs_backfill]

family_rows <- languoid_full[languoid_full$level == "family", ]
family_name <- function(codes) family_rows$name[match(codes, family_rows$id)]

glottography <- glottography[
    !glottography$family_id %in% EXCLUDED_FAMILY_IDS &
        !glottography$family_id %in% INDO_EUROPEAN_FAMILY_ID &
        !glottography$id %in% NON_INDIGENOUS_GLOTTOCODES,
]

glottography$language_id <- ifelse(is.na(glottography$cldf_langu), glottography$name, glottography$cldf_langu)
glottography$isolate <- is.na(glottography$family_id)
glottography$family_id[glottography$isolate] <- glottography$language_id[glottography$isolate]

# Function to load the C14 dates for transition to cultivation from Souza et al. 2025

load_farming_dates <- function(path = here("data/c14dates.csv")) {
    dates <- read.csv(path)
    dates <- st_as_sf(dates, coords = c("Longitude", "Latitude"), crs = 4326) |>
        st_transform(crs(south_america))
}

# Compare Glottolog and Glottography coverage and agreement

source_coverage <- function() {
    families <- unique(glottolog$family)
    data.frame(
        level = c("language (per point)", "family (per distinct family)"),
        n_glottolog = c(nrow(glottolog), length(families)),
        n_matched = c(
            sum(glottolog$glottocode %in% unique(na.omit(glottography$cldf_langu))),
            sum(families %in% unique(glottography$family_id))
        )
    ) |> transform(pct_matched = round(100 * n_matched / n_glottolog, 1))
}

nearest_polygon_km <- function(point_key, polygon_key) {
    out <- rep(NA_real_, nrow(glottolog))
    for (key in unique(glottolog[[point_key]])) {
        polys <- glottography[!is.na(glottography[[polygon_key]]) &
                                  glottography[[polygon_key]] == key, ]
        if (!nrow(polys)) next
        pts <- which(glottolog[[point_key]] == key)
        nearest <- st_nearest_feature(glottolog[pts, ], polys)
        out[pts] <- as.numeric(st_distance(glottolog[pts, ], polys[nearest, ],
                                           by_element = TRUE)) / 1000
    }
    return(out)
}

source_agreement <- function() {
    do.call(rbind, lapply(
        list(c("language", "glottocode", "cldf_langu"), c("family", "family", "family_id")),
        function(spec) {
            d <- nearest_polygon_km(spec[2], spec[3])
            matched <- !is.na(d)
            outside <- matched & d > 0
            data.frame(
                level = spec[1],
                n_matched = sum(matched),
                n_within = sum(matched & d == 0),
                pct_within = round(100 * sum(matched & d == 0) / sum(matched), 1),
                n_outside = sum(outside),
                mean_dist_km = round(mean(d[outside]), 1),
                median_dist_km = round(median(d[outside]), 1)
            )
        }
    ))
}
