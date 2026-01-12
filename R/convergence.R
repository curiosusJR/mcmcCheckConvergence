# Base-R reimplementation with Rust-backed tree utilities.

isa <- function(x, class) inherits(x, class)

#' Build control settings for convergence checks
#'
#' @param tracer Logical flag to compute ESS tracer output.
#' @param burnin Burn-in fraction (0-1) or percent (> 1).
#' @param precision ESS precision threshold.
#' @param namesToExclude Regex of column names to ignore.
#' @param emitLogs Logical flag to emit Rust backend logs.
#' @return A list of control settings.
#' @export
makeControl <- function(tracer = NULL, burnin = NULL, precision = NULL, namesToExclude = NULL, emitLogs = NULL) {
  control <- vector(mode = "list", length = 5)
  names(control) <- c("tracer", "burnin", "precision", "namesToExclude", "emitLogs")
  control$tracer <- tracer
  control$burnin <- burnin
  control$precision <- precision
  control$namesToExclude <- namesToExclude
  control$emitLogs <- emitLogs
  control
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

  result <- check_convergence_r(
    list_files = list_files,
    format = format,
    tracer = tracer,
    burnin = burnin,
    precision = precision,
    names_to_exclude = namesToExclude,
    emit_logs = emitLogs
  )

  output <- list()
  output$burnin <- result$burnin
  output$message <- result$message
  output$message_complete <- result$message_complete
  output$converged <- result$converged

  if (length(result$tree_freqs) > 0 && is.null(names(result$tree_freqs))) {
    names(result$tree_freqs) <- paste0("Run_", seq_along(result$tree_freqs))
  }
  if (length(result$tree_ess) > 0 && is.null(names(result$tree_ess))) {
    names(result$tree_ess) <- paste0("Run_", seq_along(result$tree_ess))
  }
  if (length(result$cont_ess) > 0 && is.null(names(result$cont_ess))) {
    names(result$cont_ess) <- paste0("Run_", seq_along(result$cont_ess))
  }

  compar_names <- result$compar_names
  if (is.null(compar_names) || length(compar_names) == 0) {
    compar_names <- character(0)
    if (length(result$tree_freqs) > 1) {
      for (r1 in 1:(length(result$tree_freqs) - 1)) {
        for (r2 in (r1 + 1):length(result$tree_freqs)) {
          compar_names <- c(compar_names, paste("Run_", r1, "_Run_", r2, sep = ""))
        }
      }
    }
  }
  if (length(result$tree_compare) == length(compar_names) && is.null(names(result$tree_compare))) {
    names(result$tree_compare) <- compar_names
  }
  if (length(result$tree_compare_freqs) == length(compar_names) && is.null(names(result$tree_compare_freqs))) {
    names(result$tree_compare_freqs) <- compar_names
  }
  if (length(result$cont_compare) == length(compar_names) && is.null(names(result$cont_compare))) {
    names(result$cont_compare) <- compar_names
  }

  tree_ess_df <- result$tree_ess_df

  freq_per_run <- result$freq_per_run

  tree_frequencies <- result$tree_freqs
  if (length(result$tree_compare_freqs) > 0) {
    tree_frequencies <- result$tree_compare_freqs
  }
  output$tree_parameters <- list(
    frequencies = tree_frequencies,
    ess = tree_ess_df,
    compare_runs = result$tree_compare,
    freq_per_run = freq_per_run
  )
  class(output$tree_parameters) <- "convenience.table"

  cont_ess_df <- result$cont_ess_df

  cont_compare_df <- result$cont_compare_df

  output$continuous_parameters <- list(
    means = result$cont_means,
    ess = cont_ess_df,
    compare_runs = cont_compare_df
  )
  class(output$continuous_parameters) <- "convenience.table"

  if (length(result$failed) > 0) {
    output$failed <- result$failed
    class(output$failed) <- "list.fails"
  }

  if (length(result$failed_names) > 0) {
    output$failed_names <- result$failed_names
  }

  class(output$message) <- "list.fails"
  class(output$message_complete) <- "list.fails"
  class(output) <- "convenience.diag"
  output
}

#' @export
#' @method print convenience.diag
print.convenience.diag <- function(x, ...) {
  print(x$message)
}

printTableSplits <- function(output, splits_per_run = FALSE, filename = NULL) {
  if (splits_per_run == TRUE) {
    df_splits <- output$tree_parameters$freq_per_run
  } else {
    if (length(output$tree_parameters$frequencies) == 1) {
      names_spits <- names(output$tree_parameters$frequencies[[1]])
      freq_splits <- as.vector(unlist(output$tree_parameters$frequencies[[1]]))
      df_splits <- data.frame(row.names = names_spits, frequencies = freq_splits, stringsAsFactors = FALSE)
      df_splits <- as.data.frame(df_splits[row.names(output$tree_parameters$ess), , drop = FALSE])
      row.names(df_splits) <- row.names(output$tree_parameters$ess)
    } else {
      df_3 <- align_named_vectors(output$tree_parameters$frequencies)
      rownames(df_3) <- names(output$tree_parameters$frequencies)
      df_3 <- t(df_3)
      df_3 <- (rowSums(as.data.frame(df_3)) * 2) / (length(output$tree_parameters$ess) - 1)
      df_3 <- (df_3) / 4
      df_splits <- as.data.frame(df_3, stringsAsFactors = FALSE)
      df_splits <- as.data.frame(df_splits[row.names(output$tree_parameters$ess), , drop = FALSE])
      row.names(df_splits) <- row.names(output$tree_parameters$ess)
    }
    colnames(df_splits) <- "frequencies"
    df_splits$ESS <- rowSums(output$tree_parameters$ess)
    df_splits <- df_splits[order(df_splits$frequencies), ]
  }
  if (is.null(filename)) {
    return(df_splits)
  }
  utils::write.csv(df_splits, file = filename)
}

printTableContinuous <- function(output, filename = NULL) {
  df_continuous <- as.data.frame(output$continuous_parameters$means)
  colnames(df_continuous) <- "means"
  df_continuous$ESS <- rowSums(output$continuous_parameters$ess)
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

plotKS <- function(x, precision = 0.01, fill_color = NULL, filename = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {
  col_threshold <- "gray69"
  if (is.null(fill_color)) fill_color <- "cadetblue4"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01, minimumESS)
  KS_values <- vector()
  for (i in seq_len(ncol(x$continuous_parameters$compare_runs))) {
    KS_values <- c(KS_values, x$continuous_parameters$compare_runs[[i]])
  }
  y_topLim <- max(hist(KS_values, plot = FALSE)$counts)
  par(mar = c(4.1, 3.9, 2.1, 1.0))
  if (is.null(xlab)) xlab <- "Kolmogorov-Smirnov score"
  if (is.null(ylab)) ylab <- "Counts"
  if (is.null(title)) title <- "Kolmogorov-Smirnov test"
  plot(NA,
    xlab = xlab,
    ylab = ylab,
    main = title,
    cex.main = 0.9,
    xlim = c((min(minimumKS, KS_values) - 0.01), (max(minimumKS, KS_values) + 0.01)),
    ylim = c(0, y_topLim + 1),
    las = 1,
    bty = "l"
  )
  rect(xleft = min(minimumKS, KS_values) - 0.01, ybottom = 0, xright = minimumKS, ytop = y_topLim + 1, border = NA, col = "gray89")
  lines(x = c(minimumKS, minimumKS), y = c(0, y_topLim + 1), col = col_threshold, lwd = 2, lty = 2)
  hist(KS_values,
    xlim = c((min(minimumKS, KS_values) - 0.01), (max(minimumKS, KS_values) + 0.01)),
    ylim = c(0, y_topLim + 1),
    breaks = seq(0, max(minimumKS, KS_values) + 0.01, 0.0023),
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
  minimumESS <- minESS(precision)
  minimumKS <- ksThreshold(0.01, minimumESS)
  ks_values <- list()
  for (i in seq_len(nrow(x$continuous_parameters$compare_runs))) {
    ks_values[[i]] <- as.numeric(x$continuous_parameters$compare_runs[i, ])
  }
  break_values <- seq(0, max(minimumKS, unlist(ks_values)) + 0.01, 0.0023)
  colors_hist <- grDevices::hcl.colors(length(ks_values), palette = "Dark 3")
  par(mar = c(4.1, 3.9, 3.1, 5.5), xpd = FALSE)
  plot(NA,
    xlab = "Kolmogorov-Smirnov score",
    ylab = "Counts",
    main = "Kolmogorov-Smirnov test",
    cex.main = 0.9,
    xlim = c((min(minimumKS, unlist(ks_values)) - 0.01), (max(minimumKS, unlist(ks_values)) + 0.01)),
    ylim = c(0, (length(ks_values[[1]]) - 2)),
    las = 1,
    bty = "l"
  )
  lines(x = c(minimumKS, minimumKS), y = c(0, length(ks_values[[1]]) - 2), col = col_threshold, lwd = 2, lty = 2)
  rect(xleft = min(minimumKS, unlist(ks_values)) - 0.01, ybottom = 0, xright = minimumKS, ytop = length(ks_values[[1]]) - 2, border = NA, col = "gray89")
  hist(ks_values[[1]],
    xlim = c((min(minimumKS, unlist(ks_values)) - 0.01), (max(minimumKS, unlist(ks_values)) + 0.01)),
    ylim = c(0, (length(ks_values[[1]]) - 2)),
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
    legend = rownames(x$continuous_parameters$compare_runs),
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
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  if (per_run == TRUE) {
    n_runs <- ncol(x$tree_parameters$ess)
    par(mfrow = c(n_runs / 2, n_runs / 2), oma = c(1.5, 1.5, 1.5, 0) + 0.1, mar = c(2, 2, 2.3, 2))
    layout(matrix(c(1:n_runs), nrow = n_runs / 2, ncol = n_runs / 2, byrow = TRUE))
    for (i in 1:n_runs) {
      ESS_values <- x$tree_parameters$ess[, i]
      ESS_values <- ESS_values[!is.na(ESS_values)]
      if (is.null(breaks)) {
        breaks <- seq(0, (max(minimumESS, ESS_values)) + 50, 25)
      }
      y_topLim <- max(hist(ESS_values, plot = FALSE)$counts)
      x_topLim <- max(minimumESS, ESS_values) + (max(minimumESS, ESS_values)) / 10
      plot(NA,
        xlab = NA,
        ylab = NA,
        main = colnames(x$tree_parameters$ess)[i],
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
        breaks = breaks,
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
    for (i in 1:length(x$tree_parameters$ess)) {
      ESS_values <- c(ESS_values, x$tree_parameters$ess[[i]])
    }
    ESS_values <- ESS_values[!is.na(ESS_values)]
    if (is.null(breaks)) {
      breaks <- seq(0, (max(minimumESS, ESS_values)) + 50, 25)
    }
    x_topLim <- max(minimumESS, ESS_values) + (max(minimumESS, ESS_values)) / 10
    y_topLim <- max(hist(ESS_values, breaks = breaks, plot = FALSE)$counts)
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
      breaks = breaks,
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
  minimumESS <- minESS(precision)
  ESS_values <- vector()
  if (per_run == TRUE) {
    n_runs <- ncol(x$continuous_parameters$ess)
    par(mfrow = c(n_runs / 2, n_runs / 2), oma = c(1.5, 1.5, 1.5, 0) + 0.1, mar = c(2, 2, 2.3, 2))
    layout(matrix(c(1:n_runs), nrow = n_runs / 2, ncol = n_runs / 2, byrow = TRUE))
    for (i in 1:n_runs) {
      ESS_values <- x$continuous_parameters$ess[, i]
      ESS_values <- ESS_values[!is.na(ESS_values)]
      if (is.null(breaks)) {
        breaks <- seq(0, (max(minimumESS, ESS_values)) + 50, 25)
      }
      y_topLim <- max(hist(ESS_values, plot = FALSE)$counts)
      x_topLim <- max(minimumESS, ESS_values) + (max(minimumESS, ESS_values)) / 10
      plot(NA,
        xlab = NA,
        ylab = NA,
        main = colnames(x$continuous_parameters$ess)[i],
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
        breaks = breaks,
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
    for (i in 1:length(x$continuous_parameters$ess)) {
      ESS_values <- c(ESS_values, x$continuous_parameters$ess[[i]])
    }
    ESS_values <- ESS_values[!is.na(ESS_values)]
    if (is.null(breaks)) {
      breaks <- seq(0, (max(minimumESS, ESS_values)) + 50, 25)
    }
    x_topLim <- max(minimumESS, ESS_values) + (max(minimumESS, ESS_values)) / 10
    y_topLim <- max(hist(ESS_values, breaks = breaks, plot = FALSE)$counts)
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
      breaks = breaks,
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

plotDiffSplits <- function(x, precision = 0.01, filename = NULL, ...) {
  col_threshold <- "gray69"
  if (!is.null(filename)) {
    grDevices::pdf(file = filename, width = 4.5, height = 4.5)
  }
  minimumESS <- minESS(precision)
  if (minimumESS == 625) {
    fdir <- system.file("thresholds/expectedDiff_625.rds", package = "mcmcCheckConvergence")
    exp_diff_runs <- readRDS(fdir)
  } else {
    exp_diff_runs <- expectedDiffSplits(minimumESS)
  }
  diff_splits <- x$tree_parameters$compare_runs
  plot(NA,
    xlab = "Split frequency",
    ylab = "Difference",
    main = "Split difference between runs",
    cex.main = 0.9,
    xlim = c(0, 1),
    ylim = c(0, 0.5),
    las = 1,
    bty = "l"
  )
  for (i in seq_along(diff_splits)) {
    if (length(diff_splits[[i]]) > 0) {
      points(exp_diff_runs[1, ], exp_diff_runs[2, ], type = "l", col = col_threshold, lwd = 2, lty = 2)
      points(as.numeric(names(diff_splits[[i]])), diff_splits[[i]], pch = 19, cex = 0.6)
    }
  }
  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}

# Rust-backed wrappers (override R implementations).
align_named_vectors <- function(vec_list) .Call(wrap__align_named_vectors, vec_list)

readTrace <- function(paths, format = "simple", delim, burnin = 0, check.names = FALSE, skip, ...) {
  .Call(wrap__read_trace, paths, format, delim, burnin, check.names, skip)
}

loadTrees <- function(file, type = NA, format = "revbayes", gens.per.tree = NA, trim = 1, logfile = NA, skip = NA) {
  if (is.na(gens.per.tree)) gens.per.tree <- NaN
  if (is.na(logfile)) logfile <- NULL
  if (is.na(skip)) skip <- -1L
  .Call(wrap__load_trees, file, format, gens.per.tree, trim, logfile, skip)
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
