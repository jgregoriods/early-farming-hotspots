# Significance threshold for Gi* z-scores
GI_THRESHOLD <- 1.96

hotspot_stats <- function(grid, richness_col) {
    analysis <- !is.na(grid[[richness_col]])
    nb <- suppressWarnings(poly2nb(grid[analysis, ]))
    listw_gi <- nb2listw(include.self(nb), zero.policy = TRUE)
    listw_nb <- nb2listw(nb, zero.policy = TRUE)

    gi <- rep(NA_real_, nrow(grid))
    gi[analysis] <- as.numeric(localG(grid[[richness_col]][analysis], listw_gi, zero.policy = TRUE))
    hotspot <- !is.na(gi) & gi > GI_THRESHOLD

    autocov <- rep(NA_real_, nrow(grid))
    autocov[analysis] <- lag.listw(listw_nb, as.numeric(hotspot[analysis]), zero.policy = TRUE)

    return(data.frame(gi = gi, hotspot = hotspot, autocov = autocov))
}

dated_cells <- function(grid, richness_col, hotspot_col) {
    dated <- st_drop_geometry(grid[!is.na(grid[[richness_col]]) & !is.na(grid$farming_age), ])
    dated$hotspot <- as.integer(dated[[hotspot_col]])
    return(dated)
}

wilcoxon_dates <- function(dated, richness_col) {
    if (sum(dated$hotspot) < 3 || sum(dated$hotspot) > nrow(dated) - 3) {
        return(NULL)
    }
    wt <- wilcox.test(farming_age ~ hotspot, data = dated)
    return(data.frame(
        target = richness_col,
        n_dated = nrow(dated),
        n_hotspot_dated = sum(dated$hotspot),
        median_age_hotspot = median(dated$farming_age[dated$hotspot == 1]),
        median_age_rest = median(dated$farming_age[dated$hotspot == 0]),
        p = wt$p.value
    ))
}

# Logistic and autologistic regression

ENV_VARS <- c("npp", "bio_15", "tri", "river_density")
MIN_FIT_CELLS <- 5

test_spatial_autocorrelation <- function(model, grid) {
    resid <- residuals(model, type = "pearson")
    nb <- suppressWarnings(poly2nb(grid))
    listw <- nb2listw(nb, zero.policy = TRUE)
    moran_test <- moran.test(resid, listw, zero.policy = TRUE)
    return(moran_test)
}

coef_rows <- function(model, target, model_label, cell_size, n, n_hotspot, moran) {
    tab <- summary(model)$coefficients
    data.frame(
        target = target, area_km2 = cell_size, model = model_label,
        n = n, n_hotspot = n_hotspot, term = rownames(tab),
        estimate = tab[, "Estimate"], std_error = tab[, "Std. Error"],
        z_value = tab[, "z value"], p_value = tab[, "Pr(>|z|)"],
        moran_i = moran$estimate[["Moran I statistic"]], moran_p = moran$p.value,
        row.names = NULL
    )
}

fit_quietly <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}

fit_pair <- function(grid, richness_col, cell_size) {
    grid_t <- grid[!is.na(grid[[richness_col]]), ]
    df <- data.frame(
        hotspot = as.integer(grid_t[[glue("hotspot_{richness_col}")]]),
        farming_age = grid_t$farming_age,
        st_drop_geometry(grid_t[, ENV_VARS, drop = FALSE]),
        autocov = grid_t[[glue("autocov_{richness_col}")]]
    )

    keep <- complete.cases(df)
    grid_fit <- grid_t[keep, ]
    df <- df[keep, ]

    fmla <- paste("hotspot ~ farming_age +", paste(ENV_VARS, collapse = " + "))
    plain <- fit_quietly(glm(as.formula(fmla), family = binomial(), data = df))
    auto <- fit_quietly(glm(as.formula(paste(fmla, "+ autocov")), family = binomial(), data = df))
    model <- plain$value
    model_auto <- auto$value

    moran <- test_spatial_autocorrelation(model, grid_fit)
    moran_auto <- test_spatial_autocorrelation(model_auto, grid_fit)

    return(list(
        model = model, model_auto = model_auto,
        moran = moran, moran_auto = moran_auto,
        n = nrow(df), n_hotspot = sum(df$hotspot),
        warnings = unique(c(plain$warnings, auto$warnings))
    ))
}
