# Public draws use the same component/block/term paths as starting values.
magcomp_raw_samples <- function(object)
{
  if(!inherits(object, "magComp") || !isTRUE(object$fitted) ||
     is.null(object$backend))
    stop("error: object must be a fitted magComp model")

  backend <- object$backend
  out <- list()
  if(!is.null(object$magnitude))
    out$magnitude <- list(coefficients = backend$magnitude_coefficient_samples, terms = list())
  if(!is.null(object$composition)) out$composition <- list(blocks = list())
  if(!is.null(backend$nb_size_samples)) out$magnitude$size <- backend$nb_size_samples
  for(term in object$magnitude$process_terms){
    t <- term$name
    part <- list(effects = backend$magnitude_process_effect_samples[[t]],
                 variance = backend$magnitude_process_variance_samples[[t]])
    if(!is.null(term$cov_model))
      part$covariance <- magcomp_split_covariance_samples(backend$magnitude_process_covariance_samples[[t]], term)
    out$magnitude$terms[[t]] <- part
  }
  for(b in names(object$composition$allocation_components)){
    block <- list(coefficients = backend$composition_block_coefficient_samples[[b]],
                  terms = list())
    for(term in object$composition$allocation_components[[b]]$process_terms){
      part <- list(loadings = magcomp_block_loading_samples(object, b, term),
                   scores = magcomp_block_factor_samples(object, b, term))
      if(!is.null(term$cov_model))
        part$covariance <- magcomp_split_covariance_samples(magcomp_block_covariance_samples(object, b, term), term)
      block$terms[[term$name]] <- part
    }
    out$composition$blocks[[b]] <- block
  }
  out
}

# ASCII symbols keep printed tables and plot labels portable across devices.
magcomp_parameter_symbol <- function(component, parameter, term_type = "")
{
  covariance_names <- c("phi", "nu", "lambda", "alpha", "phi_1", "lambda_1", "phi_2", "lambda_2",
                        "a", "c", "beta", "gamma", "delta")
  if(parameter %in% covariance_names)
    return(paste0(parameter, if(component == "magnitude") "_v" else "_w"))
  if(component == "magnitude"){
    return(switch(parameter,
      coefficients = "beta", size = "r_M",
      effects = if(term_type == "iid") "delta" else "v",
      variance = if(term_type == "iid") "sigma_delta^2" else "sigma_v^2",
      ""))
  }
  switch(parameter,
    coefficients = "alpha",
    loadings = if(term_type == "factor_iid") "gamma" else "lambda",
    scores = if(term_type == "factor_iid") "u" else "w",
    factor_scale = "1",
    intercept_smoothing_variance = "tau_0^2",
    slope_smoothing_variance = "tau_l^2", "")
}

# Selection and derived quantities share the retained-draw and chain indexing.
samples <- function(object, ...)
{
  UseMethod("samples")
}

samples.magComp <- function(object, type = c("parameters", "effects", "combined"),
                            start = 1, thin = 1, chains = NULL,
                            combine_chains = FALSE, group = NULL, ...)
{
  magcomp_check_prediction_args(list(...), "samples")
  type <- match.arg(type)
  out <- magcomp_raw_samples(object)
  magcomp_select_samples(object, out, type, start, thin, chains,
                         combine_chains, group)
}

samples.magCompPrediction <- function(object, type = c("parameters", "effects", "combined"),
                                      start = 1, thin = 1, chains = NULL,
                                      combine_chains = FALSE, group = NULL, ...)
{
  magcomp_check_prediction_args(list(...), "samples")
  type <- match.arg(type)
  if(is.null(object$latent))
    stop("error: samples() requires predict(..., keep_effects = TRUE)")
  magcomp_select_samples(object$latent$model, object$latent$draws, type,
                         start, thin, chains, combine_chains, group)
}

# Fit and prediction extraction share selection, grouping, and combination.
magcomp_select_samples <- function(object, out, type, start, thin, chains,
                                   combine_chains, group)
{
  if(!is.logical(combine_chains) || length(combine_chains) != 1L || is.na(combine_chains))
    stop("error: combine_chains must be TRUE or FALSE")
  for(value in list(start = start, thin = thin)){
    if(!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
       value < 1 || value != floor(value))
      stop("error: start and thin must be positive integers")
  }
  nc <- object$backend$n_chains %||% 1L
  reference <- out$magnitude$coefficients %||% out$composition$blocks[[1L]]$coefficients
  nd <- dim(reference)[length(dim(reference)) - as.integer(nc > 1L)]
  if(start > nd) stop("error: start exceeds the retained draws per chain")
  ids <- seq.int(start, nd, by = thin)
  chain_names <- object$backend$chain_ids %||% magcomp_chain_names(nc)
  if(is.null(chains)) chains <- seq_len(nc)
  if(is.character(chains)) chains <- match(chains, chain_names)
  if(!is.numeric(chains) || !length(chains) || anyNA(chains) ||
     any(!is.finite(chains)) || any(chains != floor(chains)) ||
     any(chains < 1 | chains > nc) || anyDuplicated(chains))
    stop("error: chains must select distinct existing chains by number or name")
  if(!is.null(group) && type != "combined")
    stop("error: group is only used with type = 'combined'")

  select <- function(x){
    if(is.null(x)) return(NULL)
    if(is.list(x)) return(lapply(x, select))
    # Default extraction shares the fitted arrays instead of copying every
    # latent draw just to calculate a scalar parameter summary.
    if(length(ids) == nd && identical(as.integer(chains), seq_len(nc))) return(x)
    d <- dim(x)
    axis <- length(d) - as.integer(nc > 1L)
    index <- lapply(d, seq_len)
    index[[axis]] <- ids
    if(nc > 1L) index[[axis + 1L]] <- chains
    x <- do.call(`[`, c(list(x), index, list(drop = FALSE)))
    dn <- dimnames(x)
    if(nc > 1L) dn[[axis + 1L]] <- chain_names[chains]
    dimnames(x) <- dn
    if(nc > 1L && length(chains) == 1L)
      x <- array(x, dim = dim(x)[-(axis + 1L)], dimnames = dn[-(axis + 1L)])
    x
  }
  out <- lapply(out, select)
  original_ids <- (object$backend$sample_ids %||% seq_len(nd))[ids]
  meta <- list(type = type, n_samples = length(ids), n_chains = length(chains),
               draw_ids = original_ids,
               iterations = object$control$n_burn + original_ids * object$control$n_thin,
               chain_ids = chain_names[chains], combine_chains = combine_chains,
               group = group)
  describe_terms <- function(terms) stats::setNames(
    vapply(terms, `[[`, character(1L), "term_type"),
    vapply(terms, `[[`, character(1L), "name"))
  term_types <- list(magnitude = describe_terms(object$magnitude$process_terms),
                     composition = lapply(object$composition$allocation_components,
                       function(b) describe_terms(b$process_terms)))
  if(type == "combined"){
    out <- magcomp_combined_samples(object, out, meta)
    if(combine_chains) out <- structure(lapply(out, as.numeric), info = attr(out, "info"))
  }else{
    if(type == "effects"){
      if(!is.null(out$magnitude)){
        out$magnitude <- list(terms = lapply(out$magnitude$terms, function(t) t["effects"]))
      }
      if(!is.null(out$composition)){
        out$composition$blocks <- lapply(out$composition$blocks, function(b)
          list(terms = lapply(b$terms, function(t) t["scores"])))
      }
    }
    if(combine_chains && length(chains) > 1L){
      pool <- function(x){
        if(is.null(x)) return(NULL)
        if(is.list(x)) return(lapply(x, pool))
        d <- dim(x)
        leading <- seq_len(length(d) - 2L)
        dn <- c(dimnames(x)[leading], list(NULL))
        array(x, dim = c(d[leading], prod(tail(d, 2L))), dimnames = dn)
      }
      out <- lapply(out, pool)
    }
  }
  attr(out, "sample_info") <- meta
  attr(out, "term_types") <- term_types
  class(out) <- "magCompSamples"
  out
}

# Average coefficients over fitted units within each requested group. These are
# not area-weighted response summaries; those belong to standSummary().
magcomp_combined_samples <- function(object, raw, meta)
{
  data <- object$data$unit_data
  units <- object$data$unit_levels
  group <- meta$group
  if(is.null(group)){
    values <- units
  }else if(identical(group, ".all")){
    values <- rep("all", length(units))
  }else{
    if(!is.character(group) || length(group) != 1L || !group %in% names(data))
      stop("error: group must name a fitted-unit column or be '.all'")
    values <- data[[group]]
    if(anyNA(values) || is.list(values)) stop("error: group must contain nonmissing scalar values")
  }
  levels <- unique(values)
  groups <- lapply(levels, function(g) which(values == g))
  ns <- meta$n_samples
  nc <- meta$n_chains
  draw_matrix <- function(x) matrix(x, nrow = ns, ncol = nc,
                                     dimnames = list(as.character(meta$draw_ids), meta$chain_ids))
  coefficient_name <- function(t) t$coef_name %||% "(Intercept)"
  # Flatten each array once. Rebuilding a location-by-draw matrix for every
  # location would repeatedly copy the complete prediction field.
  magnitude_effects <- lapply(raw$magnitude$terms, function(t)
    matrix(t$effects, nrow = dim(t$effects)[1L]))
  composition_effects <- lapply(raw$composition$blocks, function(b)
    lapply(b$terms, function(t) list(
      scores = matrix(t$scores, nrow = dim(t$scores)[1L]),
      loadings = matrix(t$loadings, nrow = dim(t$loadings)[1L]),
      K = dim(t$scores)[2L])))
  weighted_effect <- function(x, index){
    # One row per group/location, then draw and chain; no covariate multiplier
    # is applied because this is a varying coefficient, not a linear predictor.
    rows <- x[index,,drop = FALSE]
    draw_matrix(colMeans(rows))
  }
  draws <- list()
  info <- list()
  add <- function(value, component, block, coefficient, stick, parent, g){
    index <- paste(c(if(nzchar(stick)) stick, coefficient), collapse = ":")
    info[[length(info) + 1L]] <<- data.frame(
      component = component, block = block, term = "", parameter = "combined_coefficients",
      index = index, symbol = paste0(if(component == "magnitude") "beta" else "alpha",
        "[", index, "] + random effects"),
      unit = if(is.null(group)) as.character(levels[g]) else "",
      group = if(is.null(group)) NA_character_ else levels[g],
      coefficient = coefficient, stick = stick, parent = parent,
      stringsAsFactors = FALSE)
    draws[[length(draws) + 1L]] <<- value
  }
  terms <- object$magnitude$process_terms
  coefficients <- unique(vapply(terms, coefficient_name, character(1L)))
  for(coefficient in coefficients){
    selected <- terms[vapply(terms, coefficient_name, character(1L)) == coefficient]
    beta <- raw$magnitude$coefficients
    base <- if(coefficient %in% rownames(beta))
      draw_matrix(magcomp_coefficient_draws(beta, match(coefficient, rownames(beta)))) else draw_matrix(0)
    for(g in seq_along(groups)){
      value <- base
      for(term in selected){
        index <- if(term$term_type == "iid") term$group[groups[[g]]] else groups[[g]]
        value <- value + weighted_effect(magnitude_effects[[term$name]], index)
      }
      add(value, "magnitude", "", coefficient, "", "", g)
    }
  }
  for(b in names(object$composition$allocation_components)){
    component <- object$composition$allocation_components[[b]]
    block <- object$data$allocation_blocks[[b]]
    terms <- component$process_terms
    coefficients <- unique(vapply(terms, coefficient_name, character(1L)))
    beta <- raw$composition$blocks[[b]]$coefficients
    for(coefficient in coefficients){
      selected <- terms[vapply(terms, coefficient_name, character(1L)) == coefficient]
      for(j in seq_len(dim(beta)[1L])){
        parent_index <- if(isTRUE(block$parent_specific_coef))
          (j - 1L) %/% (ncol(block$counts) - 1L) + 1L else NULL
        block_rows <- if(is.null(parent_index)) seq_along(block$unit_index) else
          which(block$parent_index == parent_index)
        base <- if(coefficient %in% colnames(beta))
          draw_matrix(magcomp_composition_draws(beta, j, match(coefficient, colnames(beta)))) else draw_matrix(0)
        for(g in seq_along(groups)){
          rows <- block_rows[block$unit_index[block_rows] %in% groups[[g]]]
          value <- base
          for(term in selected){
            part <- composition_effects[[b]][[term$name]]
            index <- if(term$term_type == "factor_iid") term$group[rows] else block$unit_index[rows]
            K <- part$K
            score_rows <- part$scores[index,,drop = FALSE]
            average_scores <- matrix(colMeans(score_rows), nrow = K)
            loading <- part$loadings[j,]
            value <- value + draw_matrix(colSums(matrix(loading, nrow = K) * average_scores))
          }
          parent <- if(is.null(parent_index)) "" else block$parent_levels[parent_index]
          add(value, "composition", b, coefficient, rownames(beta)[j], parent, g)
        }
      }
    }
  }
  if(!length(draws)) stop("error: no varying coefficients; type = 'combined' requires random-effect terms")
  info <- do.call(rbind, info)
  names(draws) <- make.unique(paste(info$component, info$block, info$index,
                                  if(is.null(group)) info$unit else info$group, sep = ":"))
  attr(draws, "info") <- info
  draws
}

print.magCompSamples <- function(x, ...)
{
  meta <- attr(x, "sample_info")
  cat("magComp posterior samples\n")
  cat("  type: ", meta$type, "\n", sep = "")
  cat("  draws: ", meta$n_samples, " per chain; ", meta$n_chains, " chain(s)\n", sep = "")
  if(meta$combine_chains) cat("  chains pooled; original chain and draw identities retained\n")
  if(!is.null(meta$group)) cat("  averaged within: ", meta$group, "\n", sep = "")
  invisible(x)
}
