# Base-R reimplementation with Rust-backed tree utilities.

isa <- function(x, class) inherits(x, class)

#' Build control settings for convergence checks
#'
#' @param tracer Logical flag to compute ESS tracer output.
#' @param burnin Burn-in fraction (0-1) or percent (> 1).
#' @param precision ESS precision threshold.
#' @param namesToExclude Regex of column names to ignore.
#' @param emitLogs Logical flag to emit Rust backend logs.
#' @param threads Integer thread count for Rust parallelism.
#' @param fastSplits Logical flag to use a faster split-frequency backend.
#' @return A list of control settings.
#' @export
makeControl <- function(tracer = NULL, burnin = NULL, precision = NULL, namesToExclude = NULL, emitLogs = NULL, threads = NULL, fastSplits = NULL) {
  list(
    tracer = tracer,
    burnin = burnin,
    precision = precision,
    namesToExclude = namesToExclude,
    emitLogs = emitLogs,
    threads = threads,
    fastSplits = fastSplits
  )
}

#' Run MCMC convergence checks
#'
#' @param path Optional directory containing input files.
#' @param list_files Optional vector of file paths.
#' @param format Input format string, e.g. \code{"revbayes"}.
#' @param control A control list created by \code{makeControl()}.
#' @return A list with convergence flags, burn-in, messages, and parameter tables.
#' @importFrom graphics abline hist layout legend lines par points rect
#' @export
checkConvergence <- function(path = NULL, list_files = NULL, format = "revbayes", control = makeControl()) {
  if (is.null(control$tracer)) tracer <- TRUE else tracer <- control$tracer
  if (is.null(control$precision)) precision <- 0.01 else precision <- control$precision
  if (is.null(control$burnin)) burnin <- 0.0 else burnin <- control$burnin
  if (is.null(control$namesToExclude)) {
    namesToExclude <- "br_lens|bl|Iteration|Likelihood|Posterior|Prior|Gen|LnL|LnPr|state|joint|prior|likelihood|time|loglik|iter|topo|Replicate_ID|Sample|posterior|it"
  } else {
    namesToExclude <- control$namesToExclude
  }
  if (is.null(control$emitLogs)) {
    emitLogs <- TRUE
  } else {
    emitLogs <- control$emitLogs
  }

  if (!is.null(path)) {
    list_files <- list.files(path, full.names = TRUE)
  }
  if (is.null(list_files)) {
    stop("Provide path or list_files")
  }

  if (is.null(control$threads)) threads <- NULL else threads <- control$threads
  if (is.null(control$fastSplits)) fast_splits <- NULL else fast_splits <- control$fastSplits

  check_convergence_r(
    list_files = list_files,
    format = format,
    tracer = tracer,
    burnin = burnin,
    precision = precision,
    names_to_exclude = namesToExclude,
    emit_logs = emitLogs,
    threads = threads,
    fast_splits = fast_splits
  )
}

#' @export
#' @method print convenience.diag
print.convenience.diag <- function(x, ...) {
  print(x$message)
}

printTableSplits <- function(output, splits_per_run = FALSE, filename = NULL) {
  df_splits <- table_splits(output, splits_per_run)
  if (is.null(filename)) {
    return(df_splits)
  }
  utils::write.csv(df_splits, file = filename)
}

printTableContinuous <- function(output, filename = NULL) {
  df_continuous <- table_continuous(output)
  if (is.null(filename)) {
    return(df_continuous)
  }
  utils::write.csv(df_continuous, file = filename)
}

#' @export
#' @method print list.fails
print.list.fails <- function(x, ...) {
  cat(x)
}

printConvergenceTable <- function(output) {
  list(
    Continuous = printTableContinuous(output),
    Splits = printTableSplits(output)
  )
}

printConvergenceDiag <- function(output) {
  print(output$message)
}

alpha_color <- function(cols, alpha = 1) {
  cols <- grDevices::col2rgb(cols, alpha = FALSE)
  grDevices::rgb(cols[1, ] / 255, cols[2, ] / 255, cols[3, ] / 255, alpha = alpha)
}

calcRelativeDiff <- function(dataframe1, dataframe2, stats) {
  .Call(wrap__calc_relative_diff, dataframe1, dataframe2, stats)
}

check <- function(df1, df2) {
  .Call(wrap__check_diff, df1, df2)
}

get.format <- function(format) {
  .Call(wrap__get_format, format)
}

isTree <- function(x) {
  .Call(wrap__is_tree, x)
}

read.revbayestrees <- function(file) {
  .Call(wrap__read_revbayes_trees, file)
}

essTracerC <- function(x) {
  .Call(wrap__ess_tracer_c, x)
}

essTracer <- function(input, stepSize = 1) {
  ess_tracer(input)
}

#' @export
#' @method print convenience.table
print.convenience.table <- function(x, ...) {
  print(x)
}

plotKS <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {
  col_threshold <- "gray69"
  if (is.null(fill_color)) fill_color <- "cadetblue4"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  prep <- plot_ks_data(x, precision)
  KS_values <- prep$ks_values
  minimumKS <- prep$minimum_ks
  y_topLim <- prep$y_top
  xlim <- prep$xlim
  ylim <- prep$ylim
  breaks <- prep$breaks
  par(mar = c(4.1, 3.9, 2.1, 1.0))
  if (is.null(xlab)) xlab <- "Kolmogorov-Smirnov score"
  if (is.null(ylab)) ylab <- "Counts"
  if (is.null(title)) title <- "Kolmogorov-Smirnov test"
  plot(NA,
    xlab = xlab,
    ylab = ylab,
    main = title,
    cex.main = 0.9,
    xlim = xlim,
    ylim = ylim,
    las = 1,
    bty = "l"
  )
  rect(xleft = xlim[1], ybottom = 0, xright = minimumKS, ytop = ylim[2], border = NA, col = "gray89")
  lines(x = c(minimumKS, minimumKS), y = c(0, ylim[2]), col = col_threshold, lwd = 2, lty = 2)
  hist(KS_values,
    xlim = xlim,
    ylim = ylim,
    breaks = breaks,
    col = fill_color,
    border = FALSE,
    las = 1,
    add = TRUE,
    ...
  )
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

plotKSPooled <- function(x, precision = 0.01, filename = NULL, ...) {
  col_threshold <- "gray69"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 6, height = 6)
  }
  prep <- plot_ks_pooled_data(x, precision)
  ks_values <- prep$ks_values
  minimumKS <- prep$minimum_ks
  break_values <- prep$breaks
  colors_hist <- grDevices::hcl.colors(length(ks_values), palette = "Dark 3")
  par(mar = c(4.1, 3.9, 3.1, 5.5), xpd = FALSE)
  plot(NA,
    xlab = "Kolmogorov-Smirnov score",
    ylab = "Counts",
    main = "Kolmogorov-Smirnov test",
    cex.main = 0.9,
    xlim = prep$xlim,
    ylim = prep$ylim,
    las = 1,
    bty = "l"
  )
  lines(x = c(minimumKS, minimumKS), y = c(0, prep$ylim[2]), col = col_threshold, lwd = 2, lty = 2)
  rect(xleft = prep$xlim[1], ybottom = 0, xright = minimumKS, ytop = prep$ylim[2], border = NA, col = "gray89")
  hist(ks_values[[1]],
    xlim = prep$xlim,
    ylim = prep$ylim,
    breaks = break_values,
    col = alpha_color(colors_hist[1], 0.7),
    border = TRUE,
    yaxs = "i",
    las = 1,
    add = TRUE,
    ...
  )
  for (i in 2:length(ks_values)) {
    hist(ks_values[[i]],
      add = TRUE,
      col = alpha_color(colors_hist[i], 0.7),
      breaks = break_values,
      border = TRUE,
      yaxs = "i",
      las = 1,
      ...
    )
  }
  abline(v = minimumKS, col = "antiquewhite4", lwd = 2, lty = 2)
  legend("topright",
    legend = prep$labels,
    fill = alpha_color(colors_hist, 0.7),
    box.lty = 0,
    inset = c(-0.25, 0),
    cex = 0.9,
    xpd = TRUE
  )
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

plotEssSplits <- function(x, per_run = FALSE, breaks = NULL, precision = 0.01, fill_color = NULL, filename = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {
  col_threshold <- "gray69"
  if (is.null(fill_color)) fill_color <- "seagreen4"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  prep <- plot_ess_splits_data(x, per_run, precision, breaks)
  minimumESS <- prep$minimum_ess
  ESS_values <- vector()
  if (per_run == TRUE) {
    runs <- prep$runs
    n_runs <- length(runs)
    par(mfrow = c(n_runs / 2, n_runs / 2), oma = c(1.5, 1.5, 1.5, 0) + 0.1, mar = c(2, 2, 2.3, 2))
    layout(matrix(c(1:n_runs), nrow = n_runs / 2, ncol = n_runs / 2, byrow = TRUE))
    for (i in 1:n_runs) {
      run <- runs[[i]]
      ESS_values <- run$values
      y_topLim <- run$y_top
      x_topLim <- run$x_top
      plot(NA,
        xlab = NA,
        ylab = NA,
        main = run$name,
        xlim = c(0, x_topLim),
        ylim = c(0, y_topLim + 1),
        las = 1,
        bty = "l"
      )
      rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim + 1, border = NA, col = "gray89")
      lines(x = c(minimumESS, minimumESS), y = c(0, y_topLim + 1), col = col_threshold, lwd = 2, lty = 2)
      hist(ESS_values,
        xlab = NA,
        ylab = NA,
        xlim = c(0, x_topLim),
        ylim = c(0, y_topLim + 1),
        breaks = prep$breaks,
        border = FALSE,
        col = fill_color,
        add = TRUE,
        ...
      )
    }
    if (is.null(xlab)) xlab <- "Effective Sample Size"
    if (is.null(ylab)) ylab <- "Counts"
    if (is.null(title)) title <- "Effective Sample Size for splits per run"
    title(main = title, cex.main = 0.9, xlab = xlab, ylab = ylab, outer = TRUE, line = 0.5)
  } else {
    ESS_values <- prep$values
    x_topLim <- prep$x_top
    y_topLim <- prep$y_top
    par(mar = c(4.1, 3.9, 2.1, 1.0))
    if (is.null(xlab)) xlab <- "Effective Sample Size"
    if (is.null(ylab)) ylab <- "Counts"
    if (is.null(title)) title <- "Effective Sample Size for splits"
    plot(NA,
      xlab = xlab,
      ylab = ylab,
      main = title,
      cex.main = 0.9,
      xlim = c(0, x_topLim),
      ylim = c(0, y_topLim + 1),
      las = 1,
      bty = "l"
    )
    rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim + 1, border = NA, col = "gray89")
    lines(x = c(minimumESS, minimumESS), y = c(0, y_topLim + 1), col = col_threshold, lwd = 2, lty = 2)
    hist(ESS_values,
      xlim = c(0, x_topLim),
      ylim = c(0, y_topLim + 1),
      breaks = prep$breaks,
      border = FALSE,
      col = fill_color,
      add = TRUE,
      ...
    )
  }
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

plotEssContinuous <- function(x, per_run = FALSE, breaks = NULL, precision = 0.01, fill_color = NULL, filename = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {
  col_threshold <- "gray69"
  if (is.null(fill_color)) fill_color <- "seagreen4"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  prep <- plot_ess_continuous_data(x, per_run, precision, breaks)
  minimumESS <- prep$minimum_ess
  ESS_values <- vector()
  if (per_run == TRUE) {
    runs <- prep$runs
    n_runs <- length(runs)
    par(mfrow = c(n_runs / 2, n_runs / 2), oma = c(1.5, 1.5, 1.5, 0) + 0.1, mar = c(2, 2, 2.3, 2))
    layout(matrix(c(1:n_runs), nrow = n_runs / 2, ncol = n_runs / 2, byrow = TRUE))
    for (i in 1:n_runs) {
      run <- runs[[i]]
      ESS_values <- run$values
      y_topLim <- run$y_top
      x_topLim <- run$x_top
      plot(NA,
        xlab = NA,
        ylab = NA,
        main = run$name,
        xlim = c(0, x_topLim),
        ylim = c(0, y_topLim + 1),
        las = 1,
        bty = "l"
      )
      rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim + 1, border = NA, col = "gray89")
      lines(x = c(minimumESS, minimumESS), y = c(0, y_topLim + 1), col = col_threshold, lwd = 2, lty = 2)
      hist(ESS_values,
        xlab = NA,
        ylab = NA,
        xlim = c(0, x_topLim),
        ylim = c(0, y_topLim + 1),
        breaks = prep$breaks,
        border = FALSE,
        col = fill_color,
        add = TRUE,
        ...
      )
    }
    if (is.null(xlab)) xlab <- "Effective Sample Size"
    if (is.null(ylab)) ylab <- "Counts"
    if (is.null(title)) title <- "Effective Sample Size for continuous parameters per run"
    title(main = title, cex.main = 0.9, xlab = xlab, ylab = ylab, outer = TRUE, line = 0.5)
  } else {
    ESS_values <- prep$values
    x_topLim <- prep$x_top
    y_topLim <- prep$y_top
    par(mar = c(4.1, 3.9, 2.1, 1.0))
    if (is.null(xlab)) xlab <- "Effective Sample Size"
    if (is.null(ylab)) ylab <- "Counts"
    if (is.null(title)) title <- "Effective Sample Size for continuous parameters"
    plot(NA,
      xlab = xlab,
      ylab = ylab,
      main = title,
      cex.main = 0.9,
      xlim = c(0, x_topLim),
      ylim = c(0, y_topLim + 1),
      las = 1,
      bty = "l"
    )
    rect(xleft = minimumESS, ybottom = 0, xright = x_topLim, ytop = y_topLim + 1, border = NA, col = "gray89")
    lines(x = c(minimumESS, minimumESS), y = c(0, y_topLim + 1), col = col_threshold, lwd = 2, lty = 2)
    hist(ESS_values,
      xlim = c(0, x_topLim),
      ylim = c(0, y_topLim + 1),
      breaks = prep$breaks,
      border = FALSE,
      col = fill_color,
      add = TRUE,
      ...
    )
  }
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

plotDiffSplits <- function(output, minimumESS = 625, fill_color = NULL, filename = NULL, per_run = FALSE, xlab = NULL, ylab = NULL, title = NULL, ...) {
  col_threshold <- "gray69"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  prep <- plot_diff_splits_data(output, minimumESS)
  exp_diff_runs <- prep$expected
  diff_points <- prep$points
  plot(NA,
    xlab = if (is.null(xlab)) "Split frequency" else xlab,
    ylab = if (is.null(ylab)) "Difference" else ylab,
    main = if (is.null(title)) "Split difference between runs" else title,
    cex.main = 0.9,
    xlim = c(0, 1),
    ylim = c(0, 0.5),
    las = 1,
    bty = "l"
  )
  for (i in seq_along(diff_points)) {
    if (length(diff_points[[i]]$x) > 0) {
      points(exp_diff_runs[1, ], exp_diff_runs[2, ], type = "l", col = col_threshold, lwd = 2, lty = 2)
      points(diff_points[[i]]$x, diff_points[[i]]$y, pch = 19, cex = 0.6)
    }
  }
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

readTrace <- function(paths, format = "simple", delim, burnin = 0, check.names = FALSE, skip, ...) {
  .Call(wrap__read_trace, paths, format, delim, burnin, check.names, skip)
}

loadTrees <- function(file, type = NA, format = "mb", gens.per.tree = NA, trim = 1, logfile = NA, skip = NA) {
  if (is.na(gens.per.tree)) gens.per.tree <- NaN
  if (is.na(type)) type <- NULL
  if (is.na(logfile)) logfile <- NULL
  if (is.na(skip)) skip <- -1L
  .Call(wrap__load_trees_with_type, file, type, format, gens.per.tree, trim, logfile, skip)
}

loadMulti <- function(path = NULL, tree_files = NULL, log_files = NULL, format = "revbayes", labels = NA, ...) {
  if (length(labels) == 1 && is.na(labels)) labels <- NULL
  .Call(wrap__load_multi, path, tree_files, log_files, format, labels, -1L)
}

loadFiles <- function(path = NULL, list_files = NULL, format, tree_name = "psi") {
  .Call(wrap__load_files, path, list_files, format, tree_name)
}

getInfo <- function(all_runs, run, namesToExclude, trees = FALSE, splitWindows = FALSE) {
  .Call(wrap__get_info, all_runs, run, namesToExclude, trees, splitWindows)
}

removeBurnin <- function(output, burnin) .Call(wrap__remove_burnin, output, burnin)

clade.freq.named <- function(x, start, end, rooted = FALSE, ...) .Call(wrap__clade_freq_named, x, start, end)

clade.freq.tree <- function(x, rooted = FALSE, ...) .Call(wrap__clade_freq_tree, x)

clade.freq.trees <- function(x, start, end, rooted = FALSE, ...) .Call(wrap__clade_freq_trees, x, start, end)

check.clades.freq <- function(runs, freq) .Call(wrap__check_clades_freq, runs, freq)

essContParam <- function(runs, windows = FALSE, namesToExclude, tracer) {
  .Call(wrap__ess_cont_param, runs, windows, namesToExclude, tracer)
}

essSplitFreq <- function(runs, tracer) .Call(wrap__ess_split_freq, runs, tracer)

splitFreq <- function(runs, windows = FALSE) .Call(wrap__split_freq, runs, windows)

meanContParam <- function(runs, namesToExclude) .Call(wrap__mean_cont_param, runs, namesToExclude)

ksTest <- function(runs, windows = FALSE, namesToExclude) {
  .Call(wrap__ks_test, runs, windows, namesToExclude)
}

ksThreshold <- function(alpha, ess) .Call(wrap__ks_threshold, alpha, ess)

minESS <- function(per) .Call(wrap__min_ess, per)

expectedDiffSplits <- function(ess) .Call(wrap__expected_diff_splits, ess)

se <- function(x) .Call(wrap__se, x)

quants <- function(x) .Call(wrap__quants, x)
