standSummary <- function(object, ...)
{
  UseMethod("standSummary")
}

standSummary.default <- function(object, ...)
{
  stop("error: standSummary() requires output from fitted() or predict() with summary = 'samples'")
}

# Column selectors follow magCompData(): bare names or quoted column names.
magcomp_stand_column <- function(data, expr, field)
{
  name <- magcomp_name_vector(expr, field)
  if(length(name) != 1L || is.na(name) || !nzchar(name))
    stop("error: ", field, " must identify one column in data")
  if(!name %in% names(data))
    stop("error: ", field, " column not found: ", name)
  value <- data[[name]]
  if(!is.atomic(value) || !is.null(dim(value)) || length(value) != nrow(data))
    stop("error: ", field, " column must be a vector")
  value
}

standSummary.magCompPrediction <- function(object, data, unit, group = NULL,
                                         represented_area,
                                         type = c("density", "total", "proportion"),
                                         summary = c("mean", "samples", "quantile"),
                                         probs = c(0.025, 0.5, 0.975), ...)
{
  magcomp_check_prediction_args(list(...), "standSummary")
  type <- match.arg(type)
  summary <- match.arg(summary)
  probs <- magcomp_prediction_probs(probs)
  meta <- object$metadata
  if(!identical(meta$summary, "samples"))
    stop("error: standSummary() requires retained draws; use fitted() or predict() with summary = 'samples'")
  if(!meta$source %in% c("fitted", "predict"))
    stop("error: standSummary() requires location-level output from fitted() or predict()")

  # Magnitude is an overall density. Give it a single category named total.
  samples <- object$density
  if(is.null(samples)){
    if(is.null(object$magnitude))
      stop("error: standSummary() requires density or magnitude samples; proportions alone lack abundance information")
    if(type == "proportion")
      stop("error: stand proportions require category density samples")
    d <- dim(object$magnitude)
    samples <- array(object$magnitude, dim = c(d[1L], 1L, d[-1L]),
      dimnames = append(dimnames(object$magnitude), list("total"), after = 1L))
  }
  ids <- meta$unit_ids
  d <- dim(samples)
  if(!is.numeric(samples) || !length(ids) || !length(d) %in% c(3L, 4L) ||
     !identical(dimnames(samples)[[1L]], ids) || anyNA(ids) ||
     anyDuplicated(ids) || any(!nzchar(ids)) ||
     any(!is.finite(samples)) || any(samples < 0))
    stop("error: standSummary() requires finite nonnegative density samples with prediction unit IDs")

  if(missing(data) || !is.data.frame(data) || nrow(data) == 0L)
    stop("error: data must be a nonempty data frame of locations and represented areas")
  if(missing(unit))
    stop("error: unit is required; identify the prediction ID column in data")
  if(missing(represented_area))
    stop("error: represented_area is required; identify the land-area column in data")
  location <- magcomp_stand_column(data, substitute(unit), "unit")
  if(anyNA(location))
    stop("error: unit must not contain missing IDs")
  location <- as.character(location)
  if(any(!nzchar(location)))
    stop("error: unit must not contain empty IDs")
  location_rows <- match(location, ids)
  if(anyNA(location_rows))
    stop("error: unit IDs in data must match prediction unit IDs")

  group_expr <- substitute(group)
  group <- if(is.null(group_expr)) rep("all", nrow(data)) else
    magcomp_stand_column(data, group_expr, "group")
  if(anyNA(group))
    stop("error: group must not contain missing or empty labels")
  group <- as.character(group)
  if(any(!nzchar(group)))
    stop("error: group must not contain missing or empty labels")
  area <- magcomp_stand_column(data, substitute(represented_area), "represented_area")
  if(!is.numeric(area))
    stop("error: represented_area column must be numeric")
  magcomp_positive_finite(area, "represented_area")

  # Map cells to groups without copying their draws for every intersection.
  # Repeated cell/group entries add areas; draw and chain columns stay aligned.
  stand_ids <- unique(group)
  weights <- Matrix::sparseMatrix(
    i = match(group, stand_ids), j = location_rows, x = area,
    dims = c(length(stand_ids), length(ids)), dimnames = list(stand_ids, ids))
  totals <- as.matrix(weights %*% matrix(samples, nrow = d[1L]))
  stand_area <- as.numeric(Matrix::rowSums(weights))
  if(any(!is.finite(totals)) || any(!is.finite(stand_area)))
    stop("error: aggregation overflowed; check densities and represented areas")
  stand_ids <- rownames(totals)
  totals <- array(totals, dim = c(length(stand_ids), d[-1L]),
    dimnames = c(list(stand_ids), dimnames(samples)[-1L]))
  total <- magcomp_prediction_total(totals)
  if(any(!is.finite(total)))
    stop("error: aggregation overflowed; check densities and represented areas")
  positive <- matrix(total > 0, nrow = length(stand_ids))
  value <- switch(type,
    total = totals,
    density = sweep(totals, 1L, stand_area, "/"),
    proportion = sweep(totals, c(1L, seq.int(3L, length(dim(totals)))), total, "/"))
  if(type == "proportion") value[is.nan(value)] <- NA_real_

  point <- meta$family == "point_basal_area"
  units <- c(density = if(point) "basal area per unit area" else "stems per unit area",
             total = if(point) "basal area" else "stems", proportion = "fraction")
  result <- stats::setNames(list(value), type)
  result$metadata <- list(
    source = "standSummary", prediction_source = meta$source,
    target = meta$target, type = type, summary = "samples",
    family = meta$family, units = units[type], unit_ids = stand_ids,
    represented_area = stats::setNames(stand_area, stand_ids),
    components = meta$components, count_conditioning = meta$count_conditioning,
    n_samples = meta$n_samples, n_chains = meta$n_chains,
    sample_ids = meta$sample_ids, chain_ids = meta$chain_ids,
    probs = NULL, proportion_conditioning = "positive_total",
    zero_total_probability = stats::setNames(rowMeans(!positive), stand_ids),
    n_positive_draws = stats::setNames(rowSums(positive), stand_ids),
    n_draws = meta$n_draws,
    aggregation = list(unit_ids = ids,
      data = data.frame(unit = location, group = group, represented_area = area,
                        stringsAsFactors = FALSE),
      exposure = meta$exposure, total_count = meta$total_count))
  class(result) <- c("magCompPrediction", "list")
  magcomp_summarize_prediction_result(result, summary, probs)
}
