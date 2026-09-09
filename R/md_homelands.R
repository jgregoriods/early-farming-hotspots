# Wichmann's "minimal distance" homeland estimate.

suppressPackageStartupMessages({
    library(here)
    library(readr)
})

R_EARTH <- 6371  # km

# Great-circle distances between every pair of points, in km

haversine_matrix <- function(lon, lat) {
    lon <- lon * pi / 180
    lat <- lat * pi / 180
    a <- sin(outer(lat, lat, "-") / 2)^2 +
        outer(cos(lat), cos(lat), "*") * sin(outer(lon, lon, "-") / 2)^2
    a[a > 1] <- 1
    2 * R_EARTH * asin(sqrt(a))
}

# The language with the smallest mean distance to all the others

md_homeland <- function(languages) {
    d <- haversine_matrix(languages$longitude, languages$latitude)
    mean_d <- colSums(d) / (nrow(d) - 1)
    i <- which.min(mean_d)
    data.frame(md_n = nrow(languages),
               md_language = languages$name[i],
               md_glottocode = languages$id[i],
               md_lon = round(languages$longitude[i], 4),
               md_lat = round(languages$latitude[i], 4),
               md_mean_dist_km = round(mean_d[i], 1))
}

languoid <- read.csv(here("data", "languoid.csv"))
languoid <- languoid[languoid$level == "language" &
                     !is.na(languoid$latitude) & !is.na(languoid$longitude), ]

homelands <- read.csv(here("data", "family_homelands.csv"))
homelands[grep("^md_", names(homelands))] <- NULL

md <- do.call(rbind, lapply(homelands$family_code, function(code)
    md_homeland(languoid[languoid$family_id == code, ])))

homelands <- cbind(homelands, md)
write_csv(homelands, here("data", "family_homelands.csv"), na = "")

for (k in seq_len(nrow(homelands)))
    cat(sprintf("%-18s n = %3d  %-36s %8.3f %8.3f\n",
                homelands$family_name[k], homelands$md_n[k],
                homelands$md_language[k], homelands$md_lon[k], homelands$md_lat[k]))
