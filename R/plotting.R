# Map style

LAND_FILL <- "white"
COAST_COLOR <- "grey20"
COAST_WIDTH <- 0.45
CELL_BORDER <- "grey70"
CELL_BORDER_WIDTH <- 0.1
NODATA_FILL <- "grey93"
COLD_FILL <- "grey84"
RICHNESS_PALETTE <- "Blues"
RICHNESS_DESAT <- 0.45
HOTSPOT_FILL <- "#B5433E"
MD_FILL <- "white"
REST_FILL <- "grey65"

sa_map <- st_as_sf(south_america)

land_base <- function() geom_sf(data = sa_map, fill = LAND_FILL, color = COAST_COLOR, linewidth = COAST_WIDTH)
coastline <- function() geom_sf(data = sa_map, fill = NA, color = COAST_COLOR, linewidth = COAST_WIDTH)

map_theme <- function(base_size = 11, coords = TRUE) {
    theme_minimal(base_size = base_size) +
        theme(
            axis.title = element_blank(),
            axis.text = if (coords) element_text(size = base_size - 4, color = "grey45")
                        else element_blank(),
            panel.grid = if (coords) element_line(color = "grey90", linewidth = 0.3)
                         else element_blank(),
            plot.title = element_text(size = base_size - 1, hjust = 0.5),
            plot.margin = margin(2, 2, 2, 2)
        )
}

desaturate <- function(cols, amount) {
    m <- col2rgb(cols)
    g <- matrix(rep(0.299 * m[1, ] + 0.587 * m[2, ] + 0.114 * m[3, ], each = 3), nrow = 3)
    rgb(t(m + (g - m) * amount), maxColorValue = 255)
}

richness_scale <- function(name) {
    scale_fill_gradientn(name = name, na.value = NODATA_FILL,
                         colours = desaturate(RColorBrewer::brewer.pal(9, RICHNESS_PALETTE), RICHNESS_DESAT))
}

# Figure 2 panels

count_panel <- function(grid, richness_col, title = NULL, base_size = 11, coords = TRUE) {
    ggplot() +
        land_base() +
        geom_sf(data = grid, aes(fill = .data[[richness_col]]), color = CELL_BORDER,
                linewidth = CELL_BORDER_WIDTH) +
        coastline() +
        richness_scale("count") +
        labs(title = title) +
        map_theme(base_size, coords)
}

hotspot_panel <- function(grid, richness_col, title = NULL, base_size = 11, coords = TRUE) {
    ggplot() +
        land_base() +
        geom_sf(data = grid, fill = NODATA_FILL, color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
        geom_sf(data = grid[!is.na(grid[[richness_col]]), ], fill = COLD_FILL,
                color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
        geom_sf(data = grid[grid[[glue("hotspot_{richness_col}")]], ], fill = HOTSPOT_FILL,
                color = CELL_BORDER, linewidth = CELL_BORDER_WIDTH) +
        coastline() +
        labs(title = title) +
        map_theme(base_size, coords)
}

# Language families (Figure 1)

family_group <- function(x, majors) {
    factor(ifelse(x$isolate, "Isolate",
                  ifelse(x$family_name %in% majors, x$family_name, "Other")),
           levels = c(majors, "Other", "Isolate"))
}

# Saving

FIGURE_DPI <- 300
MAX_IN <- c(width = 7.5, height = 8.75)
MAX_MB <- 10

save_figure <- function(file, plot, width, height) {
    ggsave(file, plot, device = "tiff", width = width, height = height,
           dpi = FIGURE_DPI, compression = "lzw", bg = "white")
    mb <- file.size(file) / 1024^2
    flags <- c(if (width > MAX_IN[["width"]] || height > MAX_IN[["height"]])
                   glue("over {MAX_IN[['width']]} x {MAX_IN[['height']]} in"),
               if (mb > MAX_MB) glue("over {MAX_MB} MB"))
    item(sub(paste0(here(), "/"), "", file, fixed = TRUE),
         glue("  {width} x {height} in, {formatC(mb, format = 'f', digits = 1)} MB"),
         if (length(flags)) glue("  [{paste(flags, collapse = ', ')}]") else "")
}

# Language family fills (Figure 1)

FAMILY_PALETTE <- "Classic Tableau"
FILL_LIGHTEN <- 0.6
NEUTRAL_STRONG <- c("grey72", "grey52")
NEUTRAL_SOFT <- c("grey93", "grey84")

lighten <- function(cols, amount = FILL_LIGHTEN) {
    m <- col2rgb(cols)
    rgb(t(m + (255 - m) * amount), maxColorValue = 255)
}

family_colours <- function(n, palette = FAMILY_PALETTE) {
    cols <- unname(palette.colors(NULL, palette))
    cols[apply(col2rgb(cols), 2, function(x) diff(range(x)) >= 12)][seq_len(n)]
}
