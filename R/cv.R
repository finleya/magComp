magcomp_rebuild_data <- function(data, rows)
{
  fold_data <- data$data[rows,,drop = FALSE]
  for(nm in data$category_names){
    fold_data[[nm]] <- factor(fold_data[[nm]],
                              levels = data$category_levels[[nm]],
                              ordered = isTRUE(data$category_ordered[nm]))
  }

  args <- list(
    data = fold_data,
    unit = data$unit_name,
    design = data$design
  )
  if(length(data$category_names) == 1L)
    args$category <- data$category_names
  if(length(data$category_names) > 1L)
    args$categories <- data$category_names
  if(length(data$coord_names) > 0L)
    args$coords <- data$coord_names
  if(length(data$tree_factor_name) > 0L)
    args$tree_factor <- data$tree_factor_name
  if(length(data$tree_count_name) > 0L)
    args$tree_count <- data$tree_count_name
  if(length(data$attribute_names) > 0L)
    args$attributes <- data$attribute_names
  if(length(data$category_ordered) > 0L)
    args$ordered <- data$category_ordered

  do.call(magCompData, args)
}

magcomp_cv_folds <- function(units, n_fold, fold_id)
{
  n <- length(units)
  if(is.null(fold_id)){
    if(n_fold > n)
      stop("error: n_fold cannot exceed the number of sampling units")
    fold_id <- sample(rep(seq_len(n_fold), length.out = n))
    names(fold_id) <- units
    return(fold_id)
  }

  if(is.null(names(fold_id))){
    stop("error: fold_id must be named by sampling unit when supplied")
  }
  if(!all(units %in% names(fold_id))){
    stop("error: fold_id names must include all sampling units")
  }
  fold_id <- fold_id[units]
  if(anyNA(fold_id))
    stop("error: fold_id contains missing values")
  fold_id
}

magcomp_cv_observed <- function(data, type, composition = TRUE)
{
  counts <- data$counts
  if(!composition)
    counts <- matrix(rowSums(counts), ncol = 1L, dimnames = list(data$unit_levels, "total"))
  exposure <- data$sampler_area
  total <- rowSums(counts)
  values <- lapply(type, function(quantity){
    switch(quantity,
      count = counts,
      density = counts / exposure,
      magnitude = matrix(total / exposure, ncol = 1L,
                         dimnames = list(data$unit_levels, "total")),
      proportion = {
        p <- counts / total
        p[is.nan(p)] <- NA_real_
        p
      })
  })
  stats::setNames(values, type)
}

# Keep only user-supplied entries, taking their validated array names from
# preparation. Covariance settings retain partial overrides and fixed markers;
# omitted values must still receive each training fold's own defaults.
magcomp_cv_supplied_starts <- function(x, prepared)
{
  out <- list()
  for(nm in names(x)){
    if(is.null(x[[nm]])) next
    out[[nm]] <- if(nm == "covariance") x[[nm]] else
      if(is.list(x[[nm]])) magcomp_cv_supplied_starts(x[[nm]], prepared[[nm]]) else prepared[[nm]]
  }
  out
}

magcomp_cv_starting <- function(starting, model, train_data)
{
  subset_term <- function(value, term, field){
    if(is.null(value[[field]])) return(value)
    levels <- if(term$term_type %in% c("iid", "factor_iid"))
      unique(as.character(train_data$unit_data[[term$args[1L]]])) else
      train_data$unit_levels
    value[[field]] <- if(field == "scores")
      value[[field]][levels,,drop = FALSE] else value[[field]][levels]
    value
  }
  for(term in model$magnitude$process_terms){
    nm <- term$name
    if(!is.null(starting$magnitude$terms[[nm]]))
      starting$magnitude$terms[[nm]] <- subset_term(
        starting$magnitude$terms[[nm]], term, "effects")
  }
  for(b in names(model$composition$allocation_components)){
    for(term in model$composition$allocation_components[[b]]$process_terms){
      nm <- term$name
      if(!is.null(starting$composition$blocks[[b]]$terms[[nm]]))
        starting$composition$blocks[[b]]$terms[[nm]] <- subset_term(
          starting$composition$blocks[[b]]$terms[[nm]], term, "scores")
    }
  }
  starting
}

# Total count/density intervals must be formed by summing categories within
# each draw before summarizing. There is no useful total-proportion output.
magcomp_cv_totals <- function(prediction, summary, probs)
{
  if(prediction$metadata$count_conditioning != "none") return(list())
  types <- intersect(prediction$metadata$type, c("count", "density"))
  stats::setNames(lapply(types, function(type)
    magcomp_summarize_prediction(magcomp_prediction_total(prediction[[type]]),
      summary, probs, magnitude = TRUE)), types)
}

magCompCV <- function(data,
                      family = fixedAreaStems(),
                      magnitude = ~ 1,
                      composition = ~ 1,
                      smooth = magCompSmooth(),
                      priors = magCompPriors(),
                      starting = list(),
                      control = magCompControl(),
                      tuning = list(),
                      n_fold = 10,
                      fold_id = NULL,
                      type = NULL,
                      target = c("expected", "observation"),
                      summary = c("mean", "samples", "quantile"),
                      probs = c(0.025, 0.5, 0.975),
                      keep_fits = FALSE)
{
  if(!inherits(data, "magCompData")){
    stop("error: magCompCV() currently requires data from magCompData()")
  }
  if(is.null(data$counts))
    stop("error: magCompCV() requires categorical count data")

  if(is.null(magnitude) && is.null(composition))
    stop("error: at least one of magnitude or composition must be specified")
  specification <- list(family = family, magnitude = magnitude, composition = composition)
  type <- magcomp_prediction_types(type, specification)
  target <- match.arg(target)
  summary <- match.arg(summary)
  probs <- magcomp_prediction_probs(probs)
  magcomp_check_prediction_family(specification)
  n_fold <- as.integer(n_fold)[1L]
  if(is.na(n_fold) || n_fold < 2L)
    stop("error: n_fold must be at least 2")
  control$run_sampler <- TRUE
  data$sampler_area <- magcomp_sampler_area(data, family)
  if("density" %in% type && is.null(data$sampler_area))
    stop("error: density predictions require plot area or BAF information")

  units <- data$unit_levels
  fold_id <- magcomp_cv_folds(units, n_fold, fold_id)
  fold_levels <- unique(fold_id)
  if(length(fold_levels) < 2L)
    stop("error: fold_id must identify at least two folds")
  start_model <- NULL
  if(length(starting)){
    prepare_control <- control
    prepare_control$run_sampler <- FALSE
    start_model <- magComp(data, family, magnitude, composition, smooth,
                           priors, starting, prepare_control, tuning = tuning)
    # Validate and name only supplied starts. Unspecified values must still
    # be initialized using each training fold, without held-out responses.
    starting <- magcomp_cv_supplied_starts(starting, start_model$starting)
  }
  predictions <- totals <- stats::setNames(vector("list", length(fold_levels)),
                                           as.character(fold_levels))
  fits <- if(isTRUE(keep_fits)) vector("list", length(fold_levels)) else NULL
  if(isTRUE(keep_fits))
    names(fits) <- as.character(fold_levels)

  for(f in seq_along(fold_levels)){
    fold <- fold_levels[f]
    holdout_units <- units[fold_id == fold]
    train_rows <- !(as.character(data$data[[data$unit_name]]) %in%
                      holdout_units)
    train_data <- magcomp_rebuild_data(data, train_rows)
    fit <- magComp(
      train_data,
      family = family,
      magnitude = magnitude,
      composition = composition,
      smooth = smooth,
      priors = priors,
      starting = if(is.null(start_model)) starting else
        magcomp_cv_starting(starting, start_model, train_data),
      tuning = tuning,
      control = control
    )

    newdata <- data$unit_data[holdout_units,,drop = FALSE]
    holdout_area <- if(!is.null(data$sampler_area))
      data$sampler_area[match(holdout_units, data$unit_levels)] else NULL
    holdout_total <- if(is.null(magnitude) && (target == "observation" ||
      any(c("count", "density") %in% type)))
      rowSums(data$counts)[match(holdout_units, data$unit_levels)] else NULL
    # One prediction call per fold keeps all quantities and totals aligned.
    # The shared prediction summary helper applies the requested summary.
    draws <- predict(fit, newdata, type = type, target = target,
                     summary = "samples", probs = probs, area = holdout_area,
                     total_count = holdout_total)
    totals[[f]] <- magcomp_cv_totals(draws, summary, probs)
    predictions[[f]] <- magcomp_summarize_prediction_result(draws, summary, probs)
    if(isTRUE(keep_fits))
      fits[[f]] <- fit
  }

  out <- list(
    call = match.call(),
    fold_id = fold_id,
    observed = magcomp_cv_observed(data, type, !is.null(composition)),
    predictions = predictions,
    totals = totals,
    fits = fits,
    type = type,
    target = target,
    count_conditioning = if(is.null(magnitude)) "observed_total" else "none",
    summary = summary,
    probs = probs
  )
  class(out) <- "magCompCV"
  out
}

print.magCompCV <- function(x, ...)
{
  cat("magComp cross-validation\n")
  cat("  folds:      ", length(x$predictions), "\n", sep = "")
  cat("  target:     ", x$target, "\n", sep = "")
  if(x$count_conditioning == "observed_total")
    cat("  composition assessment: conditional on observed totals\n")
  cat("  summary:    ", x$summary, "\n", sep = "")
  cat("  quantities: ", paste(x$type, collapse = ", "), "\n", sep = "")
  invisible(x)
}

plot.magCompCV <- function(x, type = x$type[1L], level = c("cell", "total"),
                           file = NULL, main = NULL, ...)
{
  if(length(type) != 1L || !type %in% x$type)
    stop("error: type must select one quantity returned by magCompCV()")
  level <- match.arg(level)
  if(level == "total" && type == "proportion")
    stop("error: total proportions are not meaningful; use level = 'cell'")
  if(level == "total" && x$count_conditioning == "observed_total")
    stop("error: totals are conditioned on in composition-only CV; use level = 'cell'")
  scalar <- type == "magnitude" || level == "total"
  tab <- lapply(seq_along(x$predictions), function(f){
    prediction <- x$predictions[[f]]
    units <- prediction$metadata$unit_ids
    value <- if(level == "total" && type != "magnitude")
      x$totals[[f]][[type]] else prediction[[type]]
    observed <- x$observed[[type]][units,,drop = FALSE]
    if(level == "total") observed <- rowSums(observed)
    na_rm <- x$target == "observation" && type == "proportion"
    mean_value <- quantile_value <- NULL
    if(x$summary == "samples"){
      mean_value <- magcomp_summarize_prediction(value, "mean", x$probs,
                                                 scalar, na_rm)
      quantile_value <- magcomp_summarize_prediction(value, "quantile", x$probs,
                                                     scalar, na_rm)
    }else if(x$summary == "mean"){
      mean_value <- value
    }else{
      quantile_value <- value
    }
    n <- length(observed)
    q <- if(!is.null(quantile_value)) matrix(quantile_value, nrow = n) else NULL
    center <- if(!is.null(mean_value)) as.numeric(mean_value) else
      q[,which.min(abs(x$probs - 0.5))]
    data.frame(observed = as.numeric(observed), center = center,
      lower = if(!is.null(q) && ncol(q) > 1L) q[,which.min(x$probs)] else NA_real_,
      upper = if(!is.null(q) && ncol(q) > 1L) q[,which.max(x$probs)] else NA_real_)
  })
  tab <- do.call(rbind, tab)
  if(!any(is.finite(tab$observed) & is.finite(tab$center)))
    stop("error: no finite observed/predicted pairs to plot")
  center_label <- if(x$summary == "quantile")
    paste0(100 * x$probs[which.min(abs(x$probs - 0.5))], "% quantile") else "mean"
  title <- main %||% paste("magComp CV", x$target, level, type)
  if(!is.null(file)){
    grDevices::png(file, width = 950, height = 800)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  lim <- range(unlist(tab), finite = TRUE)
  graphics::plot(tab$observed, tab$center, pch = 19,
                 col = grDevices::adjustcolor("gray25", 0.55),
                 xlab = "Observed", ylab = paste("Prediction", center_label),
                 main = title, xlim = lim, ylim = lim, ...)
  graphics::segments(tab$observed, tab$lower, tab$observed, tab$upper,
                     col = grDevices::adjustcolor("steelblue", 0.25))
  graphics::abline(0, 1, lwd = 2)
  invisible(x)
}
