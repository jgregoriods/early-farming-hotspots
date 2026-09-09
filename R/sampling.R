suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(here)
    library(glue)
})

GRID_AREAS <- c(25000, 50000, 100000, 175000, 250000)
REFERENCE_AREA <- 100000

# Keep only cells with > 10% land
MIN_OVERLAP <- 0.1

area_label <- function(area_km2)
    paste0(format(area_km2, scientific = FALSE, trim = TRUE), "km2")

area_pretty <- function(area_km2)
    paste0(format(area_km2, big.mark = ",", scientific = FALSE, trim = TRUE), " km²")

hex_width_for_area <- function(area_km2) sqrt(2 * area_km2 / sqrt(3))

make_grid <- function(area_km2, boundary) {
    grid <- st_make_grid(boundary, cellsize = hex_width_for_area(area_km2) * 1000,
                         square = FALSE) |> st_as_sf()

    intersection <- st_intersects(grid, st_as_sf(boundary))
    grid <- grid[lengths(intersection) > 0, ]
    grid$CELL_ID <- 1:nrow(grid)

    grid$overlap <- sapply(1:nrow(grid), function(i) {
        cell <- grid[i, ]
        intersection_area <- st_area(st_intersection(cell, st_as_sf(boundary)))
        cell_area <- st_area(cell)
        return(as.numeric(intersection_area / cell_area))
    })

    return(grid)
}

extract_raster_values <- function(grid, rasters) {
    cat("Extracting raster values...\n")
    for (i in seq_along(rasters)) {
        raster <- rasters[[i]]
        raster_name <- names(raster)

        agg_fun <- median

        grid[[raster_name]] <- extract(raster, vect(grid), fun=agg_fun, na.rm=T)[, 2]
    }
    return(grid)
}

extract_river_density <- function(grid, rivers, boundary = south_america) {
    cat("Computing drainage density...\n")
    land <- st_union(st_make_valid(st_as_sf(boundary)))
    cell_land <- st_intersection(grid[, "CELL_ID"], land)
    cell_land$land_km2 <- as.numeric(st_area(cell_land)) / 1e6

    clipped <- st_intersection(st_as_sf(rivers)["ID"], cell_land[, "CELL_ID"])
    clipped$len_km <- as.numeric(st_length(clipped)) / 1000
    len_by_cell <- aggregate(len_km ~ CELL_ID, data = st_drop_geometry(clipped), FUN = sum)

    total_len <- len_by_cell$len_km[match(grid$CELL_ID, len_by_cell$CELL_ID)]
    total_len[is.na(total_len)] <- 0
    land_km2 <- cell_land$land_km2[match(grid$CELL_ID, cell_land$CELL_ID)]
    return(total_len / land_km2)
}

farming_age_by_cell <- function(grid, dates) {
    joined <- st_drop_geometry(st_join(grid, dates))
    age_by_cell <- aggregate(AgeCalBP ~ CELL_ID, data = joined, FUN = max)
    return(age_by_cell$AgeCalBP[match(grid$CELL_ID, age_by_cell$CELL_ID)])
}

extract_shp_values <- function(grid, shp, col_name, agg_fun) {
    i <- suppressWarnings(st_intersection(grid, st_as_sf(shp)))
    i_summary <- aggregate(i[[col_name]], by = list(i$CELL_ID), FUN = agg_fun)
    res <- i_summary$x[match(grid$CELL_ID, i_summary$Group.1)]
    return(res)
}

.grid_cache <- NULL

get_grid <- function(area_km2) {
    path <- here("grids.rds")
    if (is.null(.grid_cache)) {
        .grid_cache <<- if (file.exists(path)) readRDS(path) else list()
    }
    key <- area_label(area_km2)
    if (key %in% names(.grid_cache)) {
        kv(area_pretty(area_km2), glue("{nrow(.grid_cache[[key]])} cells"), "cached")
        return(.grid_cache[[key]])
    }
    kv(area_pretty(area_km2), "building ...")
    grid <- make_grid(area_km2, south_america)
    grid <- extract_raster_values(grid, get_rasters())
    grid$river_density <- extract_river_density(grid, get_rivers())
    .grid_cache[[key]] <<- grid
    saveRDS(.grid_cache, path)
    kv(area_pretty(area_km2), glue("{nrow(grid)} cells"), "built and cached")
    return(grid)
}
