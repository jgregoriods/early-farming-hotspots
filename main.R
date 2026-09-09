suppressPackageStartupMessages({
    library(glue)
    library(here)
    library(sf)
    library(spdep)
    library(ggplot2)
    library(patchwork)
})

# Console formatting

REPORT_WIDTH <- 74

section <- function(title) {
    cat("\n", title, "\n", strrep("-", REPORT_WIDTH), "\n", sep = "")
}

kv <- function(label, value, note = "") {
    cat(format(label, width = 38),
        format(format(value), width = 12, justify = "right"),
        if (nzchar(note)) paste0("  ", note) else "", "\n", sep = "")
}

item <- function(...) cat("  ", ..., "\n", sep = "")

wrapped <- function(text) {
    cat(strwrap(text, width = REPORT_WIDTH, prefix = "  "), sep = "\n")
    cat("\n")
}

source(here("R/si.R"))
source(here("R/load_data.R"))
source(here("R/sampling.R"))
source(here("R/modeling.R"))
source(here("R/plotting.R"))

dir.create(here("out"), showWarnings = FALSE)
dir.create(here("plots"), showWarnings = FALSE)
dir.create(here("plots", "SI"), showWarnings = FALSE)

# Dates from Souza et al. 2025

farming_dates <- load_farming_dates()

section("Data")
kv("Languages (Glottolog, South America)", nrow(glottolog))
kv("Families, excluding isolates", length(unique(glottolog$family[!glottolog$isolate])))
kv("Isolates", sum(glottolog$isolate))
kv("Glottography range polygons", nrow(glottography))
kv("Farming transition dates", nrow(farming_dates))

targets <- list(
    n_family = "Genealogical (points)",
    n_family_poly = "Genealogical (polygons)",
    n_language = "Language (points)",
    n_language_poly = "Language (polygons)"
)

# Compare Glottolog and Glottography

coverage <- source_coverage()
agreement <- source_agreement()

write.csv(coverage, here("out", "source_coverage.csv"), row.names = FALSE)
write.csv(agreement, here("out", "source_agreement.csv"), row.names = FALSE)

section("Glottolog vs Glottography")
kv("Languages with a polygon", glue("{coverage$n_matched[1]}/{coverage$n_glottolog[1]}"),
   glue("{coverage$pct_matched[1]}%"))
kv("Families with a polygon", glue("{coverage$n_matched[2]}/{coverage$n_glottolog[2]}"),
   glue("{coverage$pct_matched[2]}%"))
kv("Points inside their own polygon", glue("{agreement$n_within[1]}/{agreement$n_matched[1]}"),
   glue("{agreement$pct_within[1]}%, rest mean {agreement$mean_dist_km[1]} km"))
kv("Points inside a same-family polygon", glue("{agreement$n_within[2]}/{agreement$n_matched[2]}"),
   glue("{agreement$pct_within[2]}%, rest mean {agreement$mean_dist_km[2]} km"))

# Precalculate grids for all resolutions

prepare_grid <- function(cell_size) {
    grid <- get_grid(cell_size)
    grid$n_family <- extract_shp_values(grid, glottolog, "family", function(x) length(unique(x)))
    grid$n_language <- extract_shp_values(grid, glottolog, "glottocode", function(x) length(unique(x)))
    grid$n_family_poly <- extract_shp_values(grid, glottography, "family_id", function(x) length(unique(x)))
    grid$n_language_poly <- extract_shp_values(grid, glottography, "language_id", function(x) length(unique(x)))
    grid <- grid[grid$overlap >= MIN_OVERLAP, ]

    grid$farming_age <- farming_age_by_cell(grid, farming_dates)

    for (col in names(targets)) {
        stats <- hotspot_stats(grid, col)
        grid[[glue("gi_{col}")]] <- stats$gi
        grid[[glue("hotspot_{col}")]] <- stats$hotspot
        grid[[glue("autocov_{col}")]] <- stats$autocov
    }
    return(grid)
}

section("Analysis grids")
grids <- lapply(GRID_AREAS, prepare_grid)
names(grids) <- area_label(GRID_AREAS)

# Correlation between genealogical and language richness

richness_cor <- do.call(rbind, lapply(GRID_AREAS, function(cell_size) {
    g <- st_drop_geometry(grids[[area_label(cell_size)]])
    data.frame(
        area_km2 = cell_size,
        dataset = c("glottolog (points)", "glottography (polygons)"),
        n_cells = c(sum(!is.na(g$n_family) & !is.na(g$n_language)),
                    sum(!is.na(g$n_family_poly) & !is.na(g$n_language_poly))),
        cor_genealogical_language = c(cor(g$n_family, g$n_language, use = "complete.obs"),
                                      cor(g$n_family_poly, g$n_language_poly, use = "complete.obs"))
    )
}))

write.csv(richness_cor, here("out", "richness_correlation.csv"), row.names = FALSE)

r2 <- function(x) formatC(x, format = "f", digits = 2)
r_ref <- richness_cor$cor_genealogical_language[richness_cor$area_km2 == REFERENCE_AREA &
                                                richness_cor$dataset == "glottography (polygons)"]

section("Genealogical x language richness")
kv("Across the sweep, both datasets",
   glue("r = {r2(min(richness_cor$cor_genealogical_language))}-{r2(max(richness_cor$cor_genealogical_language))}"))
kv(glue("Polygons at {area_pretty(REFERENCE_AREA)}"), glue("r = {r2(r_ref)}"))

# Figure 1
# Language families of South America

glottolog$family_name <- family_name(glottolog$family)
glottography$family_name <- family_name(glottography$family_id)

N_MAJOR_FAMILIES <- 8
MAJOR_FAMILIES <- names(sort(table(glottolog$family_name), decreasing = TRUE))[1:N_MAJOR_FAMILIES]

glottolog$family_group <- family_group(glottolog, MAJOR_FAMILIES)
glottography$family_group <- family_group(glottography, MAJOR_FAMILIES)

glottography_sa <- suppressWarnings(st_intersection(glottography, sa_map))

family_strong <- c(family_colours(N_MAJOR_FAMILIES), NEUTRAL_STRONG)
family_soft <- c(lighten(family_colours(N_MAJOR_FAMILIES)), NEUTRAL_SOFT)

p_family <- ggplot() +
    land_base() +
    geom_sf(data = glottography_sa, aes(fill = family_group), color = "grey65", linewidth = 0.1) +
    geom_sf(data = glottolog, fill = family_strong[as.integer(glottolog$family_group)],
            shape = 21, size = 2.1, color = "grey20", stroke = 0.22) +
    coastline() +
    scale_fill_manual(name = "Language family", values = family_soft,
                      guide = guide_legend(override.aes = list(fill = family_strong))) +
    map_theme() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.84, 0.19),
        legend.justification = c(0.5, 0.5),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(10, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.margin = margin(4, 6, 4, 4)
    )
section("Figure 1: language families")
kv("Families drawn by name", N_MAJOR_FAMILIES)
wrapped(paste(MAJOR_FAMILIES, collapse = ", "))
save_figure(here("plots", "Fig1.tif"), p_family, width = 6.5, height = 7.5)

# Figure 2
# Richness counts and hotspots

section("Figure 2: richness counts and hotspots")

grid_ref <- grids[[area_label(REFERENCE_AREA)]]

row_label <- function(p, text) p + labs(y = text) +
    theme(axis.title.y = element_text(size = 9, angle = 90))

panels <- c(
    lapply(seq_along(targets), function(i) {
        p <- count_panel(grid_ref, names(targets)[i], targets[[i]], base_size = 9, coords = FALSE)
        if (i == 1) row_label(p, "Richness") else p
    }),
    lapply(seq_along(targets), function(i) {
        p <- hotspot_panel(grid_ref, names(targets)[i], NULL, base_size = 9, coords = FALSE)
        if (i == 1) row_label(p, "Gi* hot spots") else p
    })
)

save_figure(here("plots", "Fig2.tif"),
            wrap_plots(panels, nrow = 2) &
                theme(legend.key.width = unit(7, "pt"), legend.text = element_text(size = 6.5),
                      legend.title = element_text(size = 8)),
            width = 7.5, height = 5.6)

# S1 Fig
# Hotspot sweep across cell sizes

sweep_panel <- function(grid, col, area, first_row, first_col) {
    hotspot_panel(grid, col, if (first_row) sub(" (", "\n(", targets[[col]], fixed = TRUE) else NULL,
                  base_size = 9, coords = FALSE) +
        labs(y = if (first_col) area else NULL) +
        theme(axis.title.y = element_text(size = 8, angle = 90),
              plot.title = element_text(size = 7.5, hjust = 0.5, lineheight = 0.95))
}

sweep <- list()
for (i in seq_along(GRID_AREAS)) {
    grid <- grids[[area_label(GRID_AREAS[i])]]
    for (j in seq_along(targets)) {
        sweep[[length(sweep) + 1]] <- sweep_panel(
            grid, names(targets)[j], area_pretty(GRID_AREAS[i]), i == 1, j == 1)
    }
}

save_figure(here("plots", "SI", "S1_Fig.tif"),
            wrap_plots(sweep, nrow = length(GRID_AREAS), byrow = TRUE),
            width = 5.6, height = 8.7)

# Farming age in hotspots vs. rest of the cells

dates_rows <- list()
raw_rows <- list()

for (cell_size in GRID_AREAS) {
    grid <- grids[[area_label(cell_size)]]

    for (col in names(targets)) {
        dated <- dated_cells(grid, col, glue("hotspot_{col}"))
        raw_rows[[glue("{cell_size}_{col}")]] <- data.frame(
            area_km2 = cell_size, target = col,
            farming_age = dated$farming_age, hotspot = dated$hotspot
        )

        row <- wilcoxon_dates(dated, col)
        if (!is.null(row)) {
            row$area_km2 <- cell_size
            dates_rows[[glue("{cell_size}_{col}")]] <- row
        }
    }
}

dates_table <- do.call(rbind, dates_rows)
dates_table <- dates_table[, c("area_km2", "target", "n_dated", "n_hotspot_dated",
                               "median_age_hotspot", "median_age_rest", "p")]
raw_df <- do.call(rbind, raw_rows)
raw_df$target_label <- factor(unlist(targets[raw_df$target]), levels = unlist(targets))
raw_df$hotspot_label <- factor(ifelse(raw_df$hotspot == 1, "Hot spot", "Rest"), levels = c("Rest", "Hot spot"))
raw_df$area_label <- factor(area_pretty(raw_df$area_km2), levels = area_pretty(GRID_AREAS))

set.seed(42)
p_dates <- ggplot(raw_df, aes(x = hotspot_label, y = farming_age, fill = hotspot_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.6) +
    scale_fill_manual(values = c("Rest" = REST_FILL, "Hot spot" = HOTSPOT_FILL), guide = "none") +
    facet_grid(target_label ~ area_label) +
    labs(x = NULL, y = "Cultivation age (cal BP)") +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0, size = 7), strip.text.x = element_text(size = 7),
          axis.text = element_text(size = 6.5), axis.title = element_text(size = 8))

ref_dates <- dates_table[dates_table$area_km2 == REFERENCE_AREA &
                         dates_table$target == "n_family_poly", ]

section(glue("Cultivation age in hotspots, {area_pretty(REFERENCE_AREA)}"))
kv("Dated cells, of which hotspots", glue("{ref_dates$n_dated} / {ref_dates$n_hotspot_dated}"))
kv("Median age in hotspot cells", glue("{round(ref_dates$median_age_hotspot)} cal BP"))
kv("Median age in the rest", glue("{round(ref_dates$median_age_rest)} cal BP"))
kv("Wilcoxon rank-sum", glue("p = {signif(ref_dates$p, 2)}"))

save_figure(here("plots", "SI", "S2_Fig.tif"), p_dates, width = 7.5, height = 6.5)

si_dates <- data.frame(
    `Cell area` = area_pretty(dates_table$area_km2),
    `Richness measure` = unlist(targets[dates_table$target]),
    `n dated` = dates_table$n_dated,
    `n hotspot` = dates_table$n_hotspot_dated,
    `Median hotspot` = round(dates_table$median_age_hotspot),
    `Median rest` = round(dates_table$median_age_rest),
    p = ifelse(dates_table$p >= 0.001,
               formatC(dates_table$p, format = "f", digits = 3),
               formatC(dates_table$p, format = "e", digits = 2)),
    check.names = FALSE
)

s2 <- text_table(si_dates, align = c("l", "l", "r", "r", "r", "r", "r"))

write_pdf(s2, here("out", "S1_Table.pdf"))
item("out/S1_Table.pdf")

# Figure 3
# Genealogical hotspots and farming transition ages

p_richness <- ggplot() +
    land_base() +
    geom_sf(data = grid_ref, aes(fill = n_family_poly), color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
    geom_sf(data = glottography, fill = NA, color = "black", linewidth = 0.25) +
    coastline() +
    richness_scale("Genealogical richness\n(polygons)") +
    map_theme()

date_levels <- c("All dates", "\u2265 5,000 cal BP", "\u2265 8,000 cal BP")

old_dates <- farming_dates[!is.na(farming_dates$AgeCalBP) & farming_dates$AgeCalBP >= 5000, ]
oldest_dates <- farming_dates[!is.na(farming_dates$AgeCalBP) & farming_dates$AgeCalBP >= 8000, ]
farming_dates$date_class <- factor("All dates", levels = date_levels)
old_dates$date_class <- factor("\u2265 5,000 cal BP", levels = date_levels)
oldest_dates$date_class <- factor("\u2265 8,000 cal BP", levels = date_levels)

DATE_FILL <- "white"
DATE_SIZES <- c(0.9, 2.4, 4.0)

p_overlay <- ggplot() +
    land_base() +
    geom_sf(data = grid_ref[grid_ref$hotspot_n_family_poly, ], fill = HOTSPOT_FILL,
            color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
    geom_sf(data = farming_dates, aes(shape = date_class), size = DATE_SIZES[1], color = "grey20") +
    coastline() +
    geom_sf(data = old_dates, aes(shape = date_class, fill = date_class),
            size = DATE_SIZES[2], color = "black", stroke = 0.6) +
    geom_sf(data = oldest_dates, aes(shape = date_class, fill = date_class),
            size = DATE_SIZES[3], color = "black", stroke = 0.6) +
    scale_shape_manual(name = NULL, breaks = date_levels,
                       values = setNames(c(16, 21, 21), date_levels)) +
    scale_fill_manual(name = NULL, breaks = date_levels,
                      values = setNames(c(NA, DATE_FILL, DATE_FILL), date_levels), guide = "none") +
    guides(shape = guide_legend(override.aes = list(
        color = c("grey20", "black", "black"),
        fill = c(NA, DATE_FILL, DATE_FILL),
        size = DATE_SIZES
    ))) +
    map_theme()

section("Figure 3: hotspots and farming dates")
save_figure(here("plots", "Fig3.tif"),
            (p_richness + p_overlay) &
                theme(legend.key.width = unit(8, "pt"), legend.key.height = unit(14, "pt"),
                      legend.text = element_text(size = 7), legend.title = element_text(size = 8)),
            width = 7.5, height = 5.0)

# S3 Fig
# Alternative homelands against genealogical hotspots

homelands <- read.csv(here("data", "family_homelands.csv"))
homelands <- homelands[order(-homelands$n_languages_sa), ]
has_alt <- !is.na(homelands$alt_lon)

SOURCE_LEVELS <- c("Published", "Alternative", "Minimal distance")

as_points <- function(lon, lat, abbrev, source) {
    st_sf(abbrev = abbrev, source = factor(source, levels = SOURCE_LEVELS),
          geometry = st_geometry(st_transform(
              st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326),
              st_crs(grid_ref))))
}

published <- as_points(homelands$homeland_lon, homelands$homeland_lat,
                       homelands$abbrev, "Published")
alternative <- as_points(homelands$alt_lon[has_alt], homelands$alt_lat[has_alt],
                         homelands$abbrev[has_alt], "Alternative")
minimal <- as_points(homelands$md_lon, homelands$md_lat,
                     homelands$abbrev, "Minimal distance")

# Families are keyed by abbreviation, expanded in the caption

homeland_pts <- rbind(published, alternative, minimal)
labelled <- homeland_pts
labelled$label <- paste0(labelled$abbrev,
                         c(Published = "", Alternative = "-alt",
                           `Minimal distance` = "-md")[as.character(labelled$source)])

p_homelands <- ggplot() +
    land_base() +
    geom_sf(data = grid_ref[grid_ref$hotspot_n_family_poly, ], fill = HOTSPOT_FILL,
            color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
    coastline() +
    geom_sf(data = homeland_pts, aes(fill = source, shape = source),
            size = 2.6, color = "grey15", stroke = 0.5) +
    ggrepel::geom_text_repel(data = labelled, aes(geometry = geometry, label = label),
                             stat = "sf_coordinates", size = 3, color = "grey10", seed = 42,
                             bg.color = "white", bg.r = 0.16,
                             point.padding = 0.5, box.padding = 0.3, force = 1,
                             max.overlaps = Inf, min.segment.length = 0,
                             segment.color = "grey40", segment.size = 0.3) +
    scale_fill_manual(name = NULL, values = c(Published = "grey15", Alternative = "white",
                                              `Minimal distance` = MD_FILL)) +
    scale_shape_manual(name = NULL, values = c(Published = 21, Alternative = 21,
                                               `Minimal distance` = 24)) +
    map_theme() +
    theme(legend.position = "inside", legend.position.inside = c(0.82, 0.12),
          legend.text = element_text(size = 8),
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

section("S3 Fig: family homelands")
kv("Homelands with an alternative", sum(has_alt), paste(homelands$family_name[has_alt], collapse = ", "))
kv("Minimal-distance estimates", nrow(homelands), "md_* columns, R/md_homelands.R")
wrapped(glue("Caption key: {paste(homelands$abbrev, homelands$family_name, \
             sep = ' = ', collapse = '; ')}. Suffix -alt marks the alternative \
             homeland, -md the minimal-distance estimate."))
for (k in seq_len(nrow(homelands)))
    item(format(homelands$family_name[k], width = 18),
         format(glue("n = {homelands$md_n[k]}"), width = 8),
         format(homelands$md_language[k], width = 36),
         glue("{homelands$md_mean_dist_km[k]} km"))
save_figure(here("plots", "SI", "S3_Fig.tif"), p_homelands, width = 6.5, height = 7.5)

# Logistic and autologistic regression

section("Logistic regression")

fits <- list()

for (cell_size in GRID_AREAS) {
    grid <- grids[[area_label(cell_size)]]

    for (col in names(targets)) {
        res <- fit_pair(grid, col, cell_size)
        fits[[glue("{cell_size}_{col}")]] <- list(
            res = res, cell_size = cell_size, area = area_pretty(cell_size),
            label = targets[[col]]
        )

        kv(glue("{area_pretty(cell_size)}, {targets[[col]]}"),
           glue("{format(res$n, width = 3)} cells, {format(res$n_hotspot, width = 2)} hot"))
        for (w in res$warnings) item("! ", w)
    }
}

MAIN_TABLES <- c(table_1 = "n_family_poly", table_2 = "n_language_poly")

main_fits <- fits[glue("{REFERENCE_AREA}_{MAIN_TABLES}")]
names(main_fits) <- names(MAIN_TABLES)

for (nm in names(main_fits)) {
    f <- main_fits[[nm]]
    writeLines(fit_block(f$res, f$area, f$label), here("out", glue("{nm}.txt")))
}

si <- list()
for (cell_size in GRID_AREAS) {
    area_fits <- Filter(function(f) f$cell_size == cell_size, fits)
    si <- c(si, unlist(lapply(area_fits, function(f) fit_block(f$res, f$area, f$label))))
}

write_pdf(si, here("out", "S1_Text.pdf"))

for (nm in names(main_fits)) {
    f <- main_fits[[nm]]
    section(glue("{sub('table_', 'Table ', nm)}: {f$label} at {f$area}"))
    cat("\nPlain\n")
    print(summary(f$res$model))
    print(f$res$moran)
    cat("\nAutologistic\n")
    print(summary(f$res$model_auto))
    print(f$res$moran_auto)
}

section("Output")
kv("Tables", "table_1.txt, table_2.txt, S1_Text.pdf, S1_Table.pdf", "out/")
kv("Figures", glue("{length(list.files(here('plots'), pattern = '[.]tif$', recursive = TRUE))} files"), "plots/")
