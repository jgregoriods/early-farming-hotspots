fit_block <- function(res, area, label) {
    c(strrep("-", REPORT_WIDTH),
      glue("{area}  |  {label}"),
      strrep("-", REPORT_WIDTH),
      glue("n = {res$n} dated cells, {res$n_hotspot} hotspots"),
      if (length(res$warnings)) paste("WARNING:", res$warnings),
      "", "Plain logistic regression", "",
      capture.output(summary(res$model)),
      "", "Moran's I on the residuals", "",
      capture.output(res$moran),
      "", "Autologistic regression", "",
      capture.output(summary(res$model_auto)),
      "", "Moran's I on the residuals", "",
      capture.output(res$moran_auto),
      "", "")
}

text_table <- function(df, align = NULL) {
    if (is.null(align)) align <- ifelse(vapply(df, is.numeric, logical(1)), "r", "l")
    body <- lapply(df, function(x) format(x, trim = TRUE, justify = "none", big.mark = ","))
    w <- mapply(function(v, h) max(nchar(c(v, h), type = "width")), body, names(df))
    line <- function(vals) paste(vapply(seq_along(vals), function(i)
        format(vals[i], width = w[i],
               justify = if (align[i] == "r") "right" else "left"), character(1)),
        collapse = "  ")
    c(line(names(df)),
      strrep("-", sum(w) + 2 * (length(w) - 1)),
      vapply(seq_len(nrow(df)), function(r)
          line(vapply(body, function(x) x[r], character(1))), character(1)))
}

write_pdf <- function(lines, file, max_cex = 0.8) {
    pdf(file, width = 8.5, height = 11, family = "mono", title = basename(file))
    on.exit(grDevices::dev.off())
    par(mar = c(2, 2, 2, 2))

    plot.new()
    cex <- min(max_cex, 1 / max(strwidth(lines, units = "user", cex = 1)))
    step <- 1.75 * strheight("Mg", units = "user", cex = cex)
    per_page <- max(1, floor(1 / step))

    pages <- split(lines, ceiling(seq_along(lines) / per_page))
    for (i in seq_along(pages)) {
        if (i > 1) plot.new()
        text(0, 1 - step * (seq_along(pages[[i]]) - 1), pages[[i]],
             adj = c(0, 1), cex = cex)
    }
    file
}
