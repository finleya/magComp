# The two public location paths share quantity construction and summaries.
# Observation draws always use the complete count vector, even for a subset
# of requested outputs, so output selection does not change the simulation.
magcomp_prediction_types <- function(type, object = NULL)
{
  allowed <- c("magnitude", "proportion", "density", "count")
  defaults <- allowed
  if(!is.null(object)){
    if(is.null(object$composition)) allowed <- defaults <- c("magnitude", "count")
    if(is.null(object$magnitude)){
      allowed <- c("proportion", "density", "count")
      defaults <- "proportion"
    }
  }
  if(is.null(type)) return(defaults)
  if(!is.character(type) || length(type) == 0L || anyNA(type) ||
     anyDuplicated(type) || any(!type %in% allowed))
    stop("error: type must be a nonempty, unique subset of ",
         paste(allowed, collapse = ", "), " for this model")
  type
}

magcomp_prediction_probs <- function(probs)
{
  if(!is.numeric(probs) || length(probs) == 0L || anyNA(probs) ||
     any(!is.finite(probs)) || any(probs <= 0 | probs >= 1))
    stop("error: probs must be nonempty numeric probabilities between zero and one")
  probs
}

magcomp_check_prediction_fit <- function(object)
{
  if(!isTRUE(object$fitted) || is.null(object$backend))
    stop("error: object does not contain sampler output")
  magcomp_check_prediction_family(object)
}

# Select parameter draws before constructing new latent effects or counts.
# All backend sample arrays end in draw, or draw x chain for multiple chains.
magcomp_prediction_draws <- function(object, n_draws)
{
  if(is.null(n_draws)) return(object)
  if(!is.numeric(n_draws) || length(n_draws) != 1L ||
     !is.finite(n_draws) || n_draws < 1 || n_draws != floor(n_draws))
    stop("error: n_draws must be a positive integer or NULL")

  reference <- object$backend$magnitude_samples %||% object$backend$proportion_samples
  multi <- object$backend$n_chains > 1L
  n <- dim(reference)[length(dim(reference)) - as.integer(multi)]
  if(n_draws > n)
    stop("error: n_draws exceeds the retained draws per chain (", n, ")")
  # Midpoints of equally sized intervals through the retained chain.
  ids <- as.integer(floor((seq_len(n_draws) - 0.5) * n / n_draws) + 1L)
  select <- function(x){
    if(is.null(x)) return(NULL)
    if(is.list(x)) return(lapply(x, select))
    indices <- lapply(dim(x), seq_len)
    indices[[length(indices) - as.integer(multi)]] <- ids
    do.call(`[`, c(list(x), indices, list(drop = FALSE)))
  }
  fields <- grep("_samples$", names(object$backend), value = TRUE)
  object$backend[fields] <- lapply(object$backend[fields], select)
  object$backend$sample_ids <- ids
  object
}

# Keep the same latent draws used to construct predictions, with only the
# model indexing needed by samples(). No training arrays or fitted object
# are retained in this context. The environment is local to one predict call.
magcomp_prediction_latent <- function(object, newdata, unit_ids)
{
  n <- length(unit_ids)
  terms <- function(x, rows){
    lapply(x, function(t) list(name = t$name, term_type = t$term_type,
                              coef_name = t$coef_name,
                              group = seq_len(rows)))
  }
  context <- list(
    data = list(unit_data = newdata, unit_levels = unit_ids,
                allocation_blocks = list()),
    magnitude = if(is.null(object$magnitude)) NULL else
      list(process_terms = terms(object$magnitude$process_terms, n)),
    composition = if(is.null(object$composition)) NULL else
      list(allocation_components = list()),
    backend = object$backend[c("n_chains", "chain_ids", "sample_ids")],
    control = object$control[c("n_burn", "n_thin")]
  )
  for(b in names(object$composition$allocation_components)){
    block <- object$data$allocation_blocks[[b]]
    np <- if(identical(block$type, "conditional")) length(block$parent_levels) else 1L
    context$data$allocation_blocks[[b]] <- list(
      counts = matrix(integer(0L), 0L, ncol(block$counts)),
      unit_index = rep(seq_len(n), np),
      parent_index = rep(seq_len(np), each = n),
      parent_specific_coef = block$parent_specific_coef,
      parent_levels = block$parent_levels
    )
    context$composition$allocation_components[[b]] <- list(
      process_terms = terms(object$composition$allocation_components[[b]]$process_terms,
                            n * np))
  }
  latent <- new.env(parent = emptyenv())
  latent$model <- context
  latent$draws <- magcomp_raw_samples(object)
  latent
}

magcomp_check_prediction_family <- function(object)
{
  if(!magcomp_prediction_family_supported(object))
    stop("error: prediction supports fixedAreaStems() and pointBasalArea() fits")
}

magcomp_prediction_total_count <- function(total_count, newdata, unit_ids)
{
  if(is.null(total_count))
    stop("error: composition-only count predictions and observation draws require total_count")
  if(is.character(total_count) && length(total_count) == 1L){
    if(!total_count %in% names(newdata))
      stop("error: total_count column not found: ", total_count)
    total_count <- newdata[[total_count]]
  }
  if(!is.numeric(total_count) || length(total_count) != length(unit_ids))
    stop("error: total_count must be a numeric vector with one value per prediction unit")
  if(!is.null(names(total_count))){
    if(anyDuplicated(names(total_count)) || !setequal(names(total_count), unit_ids))
      stop("error: total_count names must match prediction unit IDs")
    total_count <- total_count[unit_ids]
  }
  if(anyNA(total_count) || any(!is.finite(total_count)) || any(total_count < 0) ||
     any(total_count != floor(total_count)) || any(total_count > .Machine$integer.max))
    stop("error: total_count must contain nonnegative integer counts within the supported integer limit")
  stats::setNames(as.numeric(total_count), unit_ids)
}

# Retain a category dimension of length one for magnitude-only count output.
magcomp_total_count_array <- function(x)
{
  array(x, dim = append(dim(x), 1L, after = 1L),
        dimnames = append(dimnames(x), list("total"), after = 1L))
}

magcomp_multinomial_draws <- function(proportion, total_count)
{
  d <- dim(proportion)
  n <- d[1L]
  J <- d[2L]
  probabilities <- array(proportion, c(n, J, prod(d[-(1:2)])))
  out <- array(0L, dim = dim(probabilities))
  for(s in seq_len(dim(probabilities)[3L])) for(i in seq_len(n))
    out[i,,s] <- as.vector(stats::rmultinom(1L, total_count[i], probabilities[i,,s]))
  array(out, dim = d, dimnames = dimnames(proportion))
}

# A spatial process has one value at a coordinate in each posterior draw.
# Reuse effects at exact coordinate ties before conditioning on distinct new
# locations; repeated coordinates do not require another spatial draw.
magcomp_spatial_prediction_samples <- function(effect_samples, coords_obs,
                                                coords_new, covariance, variance,
                                                m, method, is_factor, cov_model, st_scale)
{
  n_obs <- nrow(coords_obs)
  n_new <- nrow(coords_new)
  coords <- rbind(coords_obs, coords_new)
  ord <- do.call(order, unname(as.data.frame(coords)))
  distinct <- c(TRUE, rowSums(coords[ord[-1L],,drop = FALSE] !=
    coords[ord[-length(ord)],,drop = FALSE]) > 0L)
  # Stable ordering puts observed rows before new rows at an identical point.
  first <- integer(nrow(coords))
  first[ord] <- ord[distinct][cumsum(distinct)]
  first <- first[n_obs + seq_len(n_new)]
  observed <- first <= n_obs
  unique_new <- which(first == n_obs + seq_len(n_new))
  draws <- if(length(unique_new))
    .Call("magcomp_spatial_predict_cpp", effect_samples, coords_obs,
          coords_new[unique_new,,drop = FALSE], covariance, variance,
          as.integer(m), method, is_factor, cov_model, as.double(st_scale), PACKAGE = "magComp") else NULL
  if(length(unique_new) == n_new)
    return(draws)
  d <- dim(effect_samples)
  out <- matrix(NA_real_, nrow = n_new, ncol = prod(d[-1L]))
  if(any(observed))
    out[observed,] <- matrix(effect_samples, nrow = n_obs)[first[observed],,drop = FALSE]
  if(any(!observed))
    out[!observed,] <- matrix(draws, nrow = length(unique_new))[
      match(first[!observed] - n_obs, unique_new),,drop = FALSE]
  array(out, dim = c(n_new, d[-1L]))
}

magcomp_summarize_prediction <- function(samples, summary, probs,
                                         magnitude = FALSE, na_rm = FALSE,
                                         unit_levels = dimnames(samples)[[1L]])
{
  if(summary == "samples")
    return(samples)
  d <- dim(samples)
  n <- d[1L]
  J <- if(magnitude) 1L else d[2L]
  draws <- matrix(samples, nrow = n * J)
  categories <- if(magnitude) "total" else dimnames(samples)[[2L]]
  if(summary == "mean"){
    value <- rowMeans(draws, na.rm = na_rm)
    value[is.nan(value)] <- NA_real_
    return(matrix(value, nrow = n, ncol = J,
                  dimnames = list(unit_levels, categories)))
  }
  # Give quantile and unit/category dimensions explicitly: apply() simplifies
  # a single probability and would otherwise break singleton summaries.
  value <- apply(draws, 1L, function(x){
    if(na_rm) x <- x[!is.na(x)]
    if(length(x) == 0L) return(rep(NA_real_, length(probs)))
    stats::quantile(x, probs = probs, names = FALSE)
  })
  q <- matrix(value, nrow = length(probs), ncol = n * J)
  labels <- paste0(100 * probs, "%")
  if(magnitude){
    out <- t(q)
    dimnames(out) <- list(unit_levels, labels)
    return(out)
  }
  out <- aperm(array(q, dim = c(length(probs), n, J)), c(2L, 3L, 1L))
  dimnames(out) <- list(unit_levels, categories, labels)
  out
}

magcomp_prediction_total <- function(count)
{
  d <- dim(count)
  total <- apply(count, c(1L, seq.int(3L, length(d))), sum)
  array(total, dim = d[-2L], dimnames = dimnames(count)[-2L])
}

magcomp_summarize_prediction_result <- function(result, summary, probs)
{
  meta <- result$metadata
  for(quantity in meta$type)
    result[[quantity]] <- magcomp_summarize_prediction(
      result[[quantity]], summary, probs, magnitude = quantity == "magnitude",
      na_rm = meta$proportion_conditioning == "positive_total" && quantity == "proportion")
  result$metadata$summary <- summary
  result$metadata$probs <- if(summary == "quantile") probs else NULL
  result
}

magcomp_prediction_result <- function(object, magnitude, proportion, exposure,
                                       target, type, summary, probs, source,
                                       total_count = NULL)
{
  reference <- magnitude %||% proportion
  sample_axis <- if(is.null(magnitude)) 3L else 2L
  n <- dim(reference)[1L]
  n_samples <- dim(reference)[sample_axis]
  n_chains <- if(length(dim(reference)) > sample_axis) dim(reference)[sample_axis + 1L] else 1L
  unit_ids <- dimnames(reference)[[1L]]
  samples <- list()
  zero_probability <- positive_draws <- NULL
  if(is.null(object$magnitude)){
    if(target == "expected"){
      if("proportion" %in% type) samples$proportion <- proportion
      if(any(c("count", "density") %in% type)){
        count <- sweep(proportion, 1L, total_count, "*")
        if("count" %in% type) samples$count <- count
        if("density" %in% type) samples$density <- sweep(count, 1L, exposure, "/")
      }
    }else{
      count <- magcomp_multinomial_draws(proportion, total_count)
      if("count" %in% type) samples$count <- count
      if("density" %in% type) samples$density <- sweep(count, 1L, exposure, "/")
      if("proportion" %in% type){
        p <- sweep(count, 1L, total_count, "/")
        p[is.nan(p)] <- NA_real_
        samples$proportion <- p
      }
      zero_probability <- stats::setNames(as.numeric(total_count == 0), unit_ids)
      positive_draws <- stats::setNames(as.numeric(total_count > 0) * n_samples * n_chains, unit_ids)
    }
  }else if(is.null(object$composition)){
    if(target == "expected"){
      if("magnitude" %in% type) samples$magnitude <- magnitude
      if("count" %in% type)
        samples$count <- magcomp_total_count_array(sweep(magnitude, 1L, exposure, "*"))
    }else{
      count_mean <- magcomp_total_count_array(sweep(magnitude, 1L, exposure, "*"))
      size <- object$backend$nb_size_samples %||% object$backend$nb_size
      count <- magcomp_posterior_predict_nb(count_mean, size)
      total <- magcomp_prediction_total(count)
      if("count" %in% type) samples$count <- count
      if("magnitude" %in% type) samples$magnitude <- sweep(total, 1L, exposure, "/")
      positive <- matrix(total > 0, nrow = n)
      zero_probability <- stats::setNames(rowMeans(!positive), unit_ids)
      positive_draws <- stats::setNames(rowSums(positive), unit_ids)
    }
  }else if(target == "observation"){
    count_mean <- magcomp_predict_category_samples("count", magnitude, proportion, exposure)
    size <- object$backend$nb_size_samples %||% object$backend$nb_size
    count <- magcomp_posterior_predict_nb(count_mean, size)
    total <- magcomp_prediction_total(count)
    positive <- matrix(total > 0, nrow = n)
    zero_probability <- stats::setNames(rowMeans(!positive), unit_ids)
    positive_draws <- stats::setNames(rowSums(positive), unit_ids)
    if("count" %in% type) samples$count <- count
    if("magnitude" %in% type)
      samples$magnitude <- sweep(total, 1L, exposure, "/")
    if("density" %in% type)
      samples$density <- sweep(count, 1L, exposure, "/")
    if("proportion" %in% type){
      p <- sweep(count, c(1L, seq.int(3L, length(dim(count)))), total, "/")
      p[is.nan(p)] <- NA_real_
      samples$proportion <- p
    }
  }else{
    if("magnitude" %in% type) samples$magnitude <- magnitude
    if("proportion" %in% type) samples$proportion <- proportion
    if(any(c("density", "count") %in% type)){
      density <- magcomp_predict_category_samples("density", magnitude, proportion, NULL)
      if("density" %in% type) samples$density <- density
      if("count" %in% type) samples$count <- sweep(density, 1L, exposure, "*")
    }
  }
  result <- samples[type]
  point <- object$family$name == "point_basal_area"
  density_unit <- if(point) "basal area per unit area" else "stems per unit area"
  units <- c(magnitude = density_unit, proportion = "fraction",
             density = density_unit, count = if(point) "in-tree tally" else "stems")
  result$metadata <- list(
    source = source, target = target, type = type, summary = "samples",
    family = object$family$name, units = units[type], unit_ids = unit_ids,
    exposure = if(is.null(exposure)) NULL else stats::setNames(exposure, unit_ids),
    components = c(if(!is.null(object$magnitude)) "magnitude", if(!is.null(object$composition)) "composition"),
    count_conditioning = if(is.null(object$magnitude) && !is.null(total_count))
      if(source == "fitted") "observed_total" else "supplied_total" else "none",
    total_count = total_count,
    n_samples = n_samples, n_chains = n_chains,
    sample_ids = object$backend$sample_ids %||% seq_len(n_samples),
    chain_ids = object$backend$chain_ids %||% paste0("chain", seq_len(n_chains)),
    probs = NULL,
    proportion_conditioning = if(target == "observation") "positive_total" else "none",
    zero_total_probability = zero_probability,
    n_positive_draws = positive_draws,
    n_draws = n_samples * n_chains
  )
  class(result) <- c("magCompPrediction", "list")
  magcomp_summarize_prediction_result(result, summary, probs)
}

summary.magCompPrediction <- function(object, probs = c(0.025, 0.5, 0.975), ...)
{
  if(object$metadata$summary != "samples")
    stop("error: summary() requires retained prediction draws; use fitted() or predict() with summary = 'samples'")
  probs <- magcomp_prediction_probs(probs)
  meta <- object$metadata
  rows <- lapply(meta$type, function(quantity){
    samples <- object[[quantity]]
    magnitude <- quantity == "magnitude"
    categories <- if(magnitude) "total" else dimnames(samples)[[2L]]
    draws <- matrix(samples, nrow = length(meta$unit_ids) * length(categories))
    na_rm <- meta$proportion_conditioning == "positive_total" && quantity == "proportion"
    means <- magcomp_summarize_prediction(samples, "mean", probs,
                                          magnitude, na_rm)
    quantiles <- magcomp_summarize_prediction(samples, "quantile", probs,
                                              magnitude, na_rm)
    quantiles <- matrix(quantiles, nrow = nrow(draws))
    colnames(quantiles) <- paste0("q", gsub("\\.", "_", as.character(100 * probs)))
    labels <- expand.grid(unit = meta$unit_ids, category = categories,
                          stringsAsFactors = FALSE)
    data.frame(unit = labels$unit, quantity = quantity,
               category = labels$category, mean = as.numeric(means),
               sd = apply(draws, 1L, stats::sd, na.rm = na_rm), quantiles,
               n_draws = if(na_rm) rowSums(!is.na(draws)) else ncol(draws),
               row.names = NULL, check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  class(out) <- c("magCompPredictionSummary", "data.frame")
  attr(out, "metadata") <- meta
  attr(out, "probs") <- probs
  out
}

print.magCompPredictionSummary <- function(x, digits = 4, ...)
{
  meta <- attr(x, "metadata")
  # Data-frame column subsetting can drop the prediction metadata.
  if(is.null(meta)) return(print.magCompSummary(x, digits = digits, ...))
  cat("magComp ", meta$source, " summary (target = '", meta$target, "')\n", sep = "")
  if(meta$count_conditioning %in% c("observed_total", "supplied_total"))
    cat("Counts condition on ",
        if(meta$count_conditioning == "observed_total") "observed" else "supplied",
        " totals.\n", sep = "")
  if(meta$proportion_conditioning == "positive_total" && "proportion" %in% meta$type)
    cat("Proportions condition on a positive total.\n")
  print.magCompSummary(x, digits = digits, ...)
}

print.magCompPrediction <- function(x, ...)
{
  meta <- x$metadata
  cat("magComp ", meta$source, " result\n", sep = "")
  cat("  target:     ", meta$target, "\n", sep = "")
  cat("  summary:    ", meta$summary, "\n", sep = "")
  cat("  quantities: ", paste(meta$type, collapse = ", "), "\n", sep = "")
  cat("  units:      ", length(meta$unit_ids), "\n", sep = "")
  cat("  draws:      ", meta$n_samples, " per chain; ", meta$n_chains,
      " chain(s)\n", sep = "")
  if(meta$count_conditioning %in% c("observed_total", "supplied_total"))
    cat("  count predictions: conditional on ",
        if(meta$count_conditioning == "observed_total") "observed" else "supplied",
        " totals\n", sep = "")
  if(meta$proportion_conditioning == "positive_total" &&
     "proportion" %in% meta$type)
    cat("  proportion summaries condition on a positive total;",
        " zero-total information is in $metadata\n", sep = "")
  invisible(x)
}
