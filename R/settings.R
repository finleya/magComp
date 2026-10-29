# Public model settings are component/block/term lists. The sampler's flat
# arrays are an internal representation, constructed only after validation.
magcomp_check_names <- function(x, allowed, path)
{
  if(!is.list(x) || (length(x) && (is.null(names(x)) ||
     anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))))
    stop("error: ", path, " must be a list with unique, nonempty names")
  bad <- setdiff(names(x), allowed)
  if(length(bad))
    stop("error: unknown setting in ", path, ": ", paste(bad, collapse = ", "))
  invisible(x)
}

magcomp_positive <- function(x, n = 1L, path)
{
  if(!is.numeric(x) || length(x) != n || anyNA(x) ||
     any(!is.finite(x)) || any(x <= 0))
    stop("error: ", path, " must contain ", n, " positive finite value(s)")
  as.numeric(x)
}

magcomp_resolve_priors <- function(priors, magnitude, composition, starting = list(), tuning = list())
{
  # Revalidate even when a caller edits a settings template directly.
  priors <- do.call(magCompPriors, unclass(priors))
  defaults <- magCompPriors()
  magcomp_check_names(tuning, c("magnitude", "composition"), "tuning")
  magcomp_check_names(tuning$magnitude %||% list(), "terms", "tuning$magnitude")
  magcomp_check_names(tuning$composition %||% list(), "blocks", "tuning$composition")
  for(component in c("magnitude", "composition")){
    active <- if(component == "magnitude") !is.null(magnitude) else !is.null(composition)
    if(!active && !identical(priors[[component]], defaults[[component]]))
      stop("error: priors$", component, " is unused for the disabled component")
  }
  resolved_tuning <- list()
  mag <- priors$magnitude
  comp <- priors$composition
  if(identical(mag$size, "auto"))
    mag["size"] <- list(if(length(magnitude$process_terms)) NULL else
      c(shape = 1, rate = 0.01))
  if(!is.null(mag$size) && length(magnitude$process_terms))
    stop("error: estimating nb_size is currently supported only for ",
         "magnitude models without random-effect process terms")
  if(any(!vapply(magnitude$process_terms, function(t)
                 t$term_type %in% c("iid", "gp", "nngp"), logical(1L))))
    stop("error: magnitude formulas require iid(), gp(), or nngp() process terms")
  mag_names <- vapply(magnitude$process_terms, "[[", character(1L), "name")
  magcomp_check_names(mag$terms, mag_names, "priors$magnitude$terms")
  magcomp_check_names(tuning$magnitude$terms %||% list(), mag_names, "tuning$magnitude$terms")
  for(i in seq_along(magnitude$process_terms)){
    term <- magnitude$process_terms[[i]]
    p <- mag$terms[[term$name]] %||% list()
    spatial <- term$term_type %in% c("gp", "nngp")
    tune <- tuning$magnitude$terms[[term$name]] %||% list()
    magcomp_check_names(tune, if(spatial) "covariance" else character(), "tuning$magnitude$terms")
    if(!spatial && length(starting$magnitude$terms[[term$name]]$covariance))
      stop("error: covariance starting values supplied for nonspatial term ", term$name)
    if(!spatial && length(p$covariance))
      stop("error: covariance prior supplied for nonspatial term ", term$name)
    term$variance_prior <- p$variance %||%
      if(spatial) mag$spatial_variance else mag$iid_variance
    if(spatial) term <- magcomp_resolve_covariance(term, p$covariance %||% list(),
      starting$magnitude$terms[[term$name]]$covariance %||% list(),
      tuning$magnitude$terms[[term$name]]$covariance %||% list(),
      paste0("magnitude$terms$", term$name, "$covariance"))
    if(spatial && any(!term$covariance_fixed))
      resolved_tuning$magnitude$terms[[term$name]] <- list(covariance = magcomp_covariance_tuning_list(term))
    magnitude$process_terms[[i]] <- term
  }
  blocks <- composition$allocation_components
  magcomp_check_names(comp$blocks, names(blocks), "priors$composition$blocks")
  magcomp_check_names(tuning$composition$blocks %||% list(), names(blocks), "tuning$composition$blocks")
  for(b in names(blocks)){
    bp <- comp$blocks[[b]] %||% comp[setdiff(names(comp), "blocks")]
    terms <- blocks[[b]]$process_terms
    if(any(!vapply(terms, function(t)
                   t$term_type %in% c("factor_iid", "factor_gp", "factor_nngp"), logical(1L))))
      stop("error: composition formulas require factor_iid(), factor_gp(), or factor_nngp() terms")
    term_names <- vapply(terms, "[[", character(1L), "name")
    bt <- tuning$composition$blocks[[b]] %||% list()
    magcomp_check_names(bt, "terms", "tuning$composition$blocks")
    magcomp_check_names(bt$terms %||% list(), term_names, "tuning$composition$blocks$terms")
    magcomp_check_names(bp$terms %||% list(), term_names,
                        paste0("priors$composition$blocks$", b, "$terms"))
    for(i in seq_along(terms)){
      term <- terms[[i]]
      p <- bp$terms[[term$name]] %||% list()
      term$loading_variance <- p$loading_variance %||% bp$loading_variance
      spatial <- term$term_type %in% c("factor_gp", "factor_nngp")
      magcomp_check_names(bt$terms[[term$name]] %||% list(), if(spatial) "covariance" else character(),
                          "tuning$composition$blocks$terms")
      if(!spatial && length(starting$composition$blocks[[b]]$terms[[term$name]]$covariance))
        stop("error: covariance starting values supplied for nonspatial term ", term$name)
      if(!spatial && length(p$covariance))
        stop("error: covariance prior supplied for nonspatial term ", term$name)
      if(spatial) term <- magcomp_resolve_covariance(term, p$covariance %||% list(),
        starting$composition$blocks[[b]]$terms[[term$name]]$covariance %||% list(),
        tuning$composition$blocks[[b]]$terms[[term$name]]$covariance %||% list(),
        paste0("composition$blocks$", b, "$terms$", term$name, "$covariance"))
      if(spatial && any(!term$covariance_fixed))
        resolved_tuning$composition$blocks[[b]]$terms[[term$name]] <- list(covariance = magcomp_covariance_tuning_list(term))
      terms[[i]] <- term
    }
    blocks[[b]]$process_terms <- terms
    blocks[[b]]$coefficient_priors <- list(intercept_variance = bp$intercept_variance,
      slope_variance = bp$slope_variance, intercept_smoothing_variance = bp$intercept_smoothing_variance,
      slope_smoothing_variance = bp$slope_smoothing_variance)
  }
  if(!is.null(composition)){
    composition$allocation_components <- blocks
    if(length(blocks) == 1L) composition$process_terms <- blocks[[1L]]$process_terms
  }
  # Only magnitude priors remain global in the native kernel.
  kernel <- if(is.null(magnitude)) list() else
    list(coefficient_variance = mag$coefficient_variance, nb_size = mag$size)
  priors$magnitude <- if(is.null(magnitude)) NULL else mag
  if(is.null(composition)) priors$composition <- NULL
  list(priors = priors, magnitude = magnitude, composition = composition, tuning = resolved_tuning, kernel = kernel)
}

magcomp_start_value <- function(value, expected, path, positive = FALSE)
{
  if(is.null(value)) return(NULL)
  if(!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
     length(value) != length(expected))
    stop("error: ", path, " must contain ", length(expected), " finite numeric value(s)")
  if(!identical(dim(value), dim(expected)))
    stop("error: ", path, " must have dimensions ",
         if(is.null(dim(expected))) "of a vector" else paste(dim(expected), collapse = " x "))
  if(positive && any(value <= 0))
    stop("error: ", path, " must be positive and finite")
  # Named inputs are matched, never silently interpreted in the wrong order.
  if(is.null(dim(expected))){
    if(!is.null(names(value))){
      if(anyDuplicated(names(value)) || !setequal(names(value), names(expected)))
        stop("error: ", path, " names must match the model template")
      value <- value[names(expected)]
    }
    return(stats::setNames(as.numeric(value), names(expected)))
  }
  for(axis in seq_along(dim(expected))){
    supplied <- dimnames(value)[[axis]]
    wanted <- dimnames(expected)[[axis]]
    if(!is.null(supplied) && (anyDuplicated(supplied) || !setequal(supplied, wanted)))
      stop("error: ", path, " dimnames must match the model template: ", paste(wanted, collapse = ", "))
  }
  if(!is.null(rownames(value))) value <- value[rownames(expected),,drop = FALSE]
  if(!is.null(colnames(value))) value <- value[,colnames(expected),drop = FALSE]
  storage.mode(value) <- "double"
  dimnames(value) <- dimnames(expected)
  value
}

magcomp_start_shapes <- function(data, magnitude, composition)
{
  out <- list()
  if(!is.null(magnitude))
    out$magnitude <- list(coefficients = stats::setNames(rep(0, ncol(magnitude$X)),
                                                        colnames(magnitude$X)), terms = list())
  if(!is.null(composition)) out$composition <- list(blocks = list())
  for(term in magnitude$process_terms){
    levels <- if(term$term_type == "iid") term$group_levels else data$unit_levels
    entry <- list(effects = stats::setNames(rep(0, term$n_group), levels), variance = 1)

    out$magnitude$terms[[term$name]] <- entry
  }
  for(b in names(composition$allocation_components)){
    block <- data$allocation_blocks[[b]]
    cmp <- composition$allocation_components[[b]]
    rows <- block$stick_names
    if(is.null(rows)) rows <- colnames(block$counts)[seq_len(block$coef_stick_count)]
    entry <- list(coefficients = matrix(0, block$coef_stick_count, ncol(cmp$X),
                                       dimnames = list(rows, colnames(cmp$X))), terms = list())
    for(term in cmp$process_terms){
      loading <- matrix(0, block$coef_stick_count, term$K,
                        dimnames = list(rows, term$factor_names))
      loading[cbind(term$anchor_index, seq_len(term$K))] <- 1
      scores <- matrix(0, term$n_group, term$K,
                       dimnames = list(magcomp_factor_group_names(term), term$factor_names))
      part <- list(loadings = loading, scores = scores)

      entry$terms[[term$name]] <- part
    }
    out$composition$blocks[[b]] <- entry
  }
  out
}

magcomp_start_to_kernel <- function(starting)
{
  out <- list(magnitude_coefficients = starting$magnitude$coefficients, nb_size = starting$magnitude$size,
              magnitude_process = list(), magnitude_process_variance = list(),
              composition_block_coefficient = list(),
              composition_block_loading = list(), composition_block_score = list())
  for(t in names(starting$magnitude$terms)){
    term <- starting$magnitude$terms[[t]]
    if(!is.null(term$effects)) out$magnitude_process[[t]] <- term$effects
    if(!is.null(term$variance)) out$magnitude_process_variance[[t]] <- term$variance
  }
  for(b in names(starting$composition$blocks)){
    block <- starting$composition$blocks[[b]]
    if(!is.null(block$coefficients)) out$composition_block_coefficient[[b]] <- block$coefficients
    for(t in names(block$terms)){
      term <- block$terms[[t]]
      if(!is.null(term$loadings)) out$composition_block_loading[[b]][[t]] <- term$loadings
      if(!is.null(term$scores)) out$composition_block_score[[b]][[t]] <- term$scores
    }
  }
  out
}

magcomp_start_from_kernel <- function(x, magnitude, composition)
{
  out <- list()
  if(!is.null(magnitude)) out$magnitude <- list(coefficients = x$magnitude_coefficients, terms = list())
  if(!is.null(composition)) out$composition <- list(blocks = list())
  if(!is.null(x$nb_size)) out$magnitude$size <- x$nb_size
  for(term in magnitude$process_terms){
    t <- term$name
    entry <- list(effects = x$magnitude_process[[t]],
                  variance = unname(x$magnitude_process_variance[[t]]))
    if(!is.null(term$cov_model)) entry$covariance <- term$covariance_starting
    out$magnitude$terms[[t]] <- entry
  }
  for(b in names(composition$allocation_components)){
    entry <- list(coefficients = x$composition_block_coefficient[[b]], terms = list())
    for(term in composition$allocation_components[[b]]$process_terms){
      t <- term$name
      part <- list(loadings = x$composition_block_loading[[b]][[t]],
                   scores = x$composition_block_score[[b]][[t]])
      if(!is.null(term$cov_model)) part$covariance <- term$covariance_starting
      entry$terms[[t]] <- part
    }
    out$composition$blocks[[b]] <- entry
  }
  out
}

magcomp_validate_starts <- function(x, shapes, path = "starting")
{
  magcomp_check_names(x, names(shapes), path)
  out <- list()
  for(nm in names(x)){
    value <- x[[nm]]
    if(is.null(value)) next
    here <- paste0(path, "$", nm)
    out[[nm]] <- if(is.list(shapes[[nm]]))
      magcomp_validate_starts(value, shapes[[nm]], here) else
      magcomp_start_value(value, shapes[[nm]], here, nm %in% c("variance", "size"))
  }
  out
}

magcomp_resolve_starting <- function(data, magnitude, composition, starting, priors, control)
{
  shapes <- magcomp_start_shapes(data, magnitude, composition)
  if(!is.null(priors$magnitude$size)) shapes$magnitude$size <- control$nb_size
  regular <- starting
  for(t in names(regular$magnitude$terms)) regular$magnitude$terms[[t]]$covariance <- NULL
  for(b in names(regular$composition$blocks))
    for(t in names(regular$composition$blocks[[b]]$terms))
      regular$composition$blocks[[b]]$terms[[t]]$covariance <- NULL
  public <- magcomp_validate_starts(regular, shapes)
  for(b in names(composition$allocation_components)){
    for(term in composition$allocation_components[[b]]$process_terms){
      loading <- public$composition$blocks[[b]]$terms[[term$name]]$loadings
      if(is.null(loading)) next
      for(k in seq_len(term$K)){
        if(loading[term$anchor_index[k], k] <= 0 ||
           (k > 1L && any(loading[term$anchor_index[seq_len(k - 1L)], k] != 0)))
          stop("error: starting loadings for ", b, "/", term$name,
               " must respect anchor constraints (positive diagonal, upper entries zero)")
      }
    }
  }
  kernel <- magcomp_start_to_kernel(public)
  kernel <- magcomp_prepare_starting(data, magnitude, composition, kernel, control$composition_link)
  # Complete all entries, including non-factor and conditional blocks, so the
  # displayed template is exactly the state passed to the sampler.
  if(is.null(kernel$magnitude_coefficients)) kernel$magnitude_coefficients <- shapes$magnitude$coefficients
  if(!is.null(priors$magnitude$size) && is.null(kernel$nb_size)) kernel$nb_size <- control$nb_size
  for(b in names(composition$allocation_components)){
    block <- data$allocation_blocks[[b]]
    cmp <- composition$allocation_components[[b]]
    if(is.null(kernel$composition_block_coefficient[[b]])){
      coef <- shapes$composition$blocks[[b]]$coefficients
      J <- ncol(block$counts) - 1L
      for(c in seq_len(nrow(coef))){
        j <- (c - 1L) %% J + 1L
        keep <- if(isTRUE(block$parent_specific_coef))
          block$parent_index == (c - 1L) %/% J + 1L else rep(TRUE, nrow(block$counts))
        success <- sum(block$counts[keep,j])
        trials <- sum(block$counts[keep,j:ncol(block$counts),drop = FALSE])
        p <- max(1e-12, min(1 - 1e-12, (success + 0.5) / (trials + 1)))
        coef[c,1L] <- if(control$composition_link == "probit") stats::qnorm(p) else stats::qlogis(p)
      }
      kernel$composition_block_coefficient[[b]] <- coef
    }
    for(term in cmp$process_terms){
      t <- term$name
      shape <- shapes$composition$blocks[[b]]$terms[[t]]
      if(is.null(kernel$composition_block_loading[[b]][[t]]))
        kernel$composition_block_loading[[b]][[t]] <- shape$loadings
      if(is.null(kernel$composition_block_score[[b]][[t]]))
        kernel$composition_block_score[[b]][[t]] <- shape$scores

    }
  }
  for(nm in names(kernel))
    if(is.list(kernel[[nm]]) && !length(kernel[[nm]])) kernel[nm] <- list(NULL)
  list(public = magcomp_start_from_kernel(kernel, magnitude, composition), kernel = kernel)
}

# Inspect a prepared or fitted model without running or advancing its sampler.
magCompSettings <- function(object)
{
  if(!inherits(object, "magComp"))
    stop("error: object must be a prepared or fitted magComp model")
  rows <- list()
  terms <- list()
  priors <- object$priors
  add <- function(component, block, term, parameter, status, prior, value, constraint = ""){
    initial <- if(length(value) <= 4L) paste(signif(value, 5), collapse = ", ") else
      if(is.null(dim(value))) paste0("vector [", length(value), "]") else
        paste0("matrix [", paste(dim(value), collapse = " x "), "]")
    term_type <- ""
    if(nzchar(term)){
      terms <- if(component == "magnitude") object$magnitude$process_terms else
        object$composition$allocation_components[[block]]$process_terms
      term_type <- Filter(function(x) x$name == term, terms)[[1L]]$term_type
    }
    symbol <- magcomp_parameter_symbol(component, parameter, term_type)
    rows[[length(rows) + 1L]] <<- data.frame(component, block, term, parameter,
                                            symbol, status, prior, initial, constraint)
  }
  normal <- function(v) paste0("Normal(0, variance=", v, ")")
  ig <- function(v) paste0("IG(shape=", v[1L], ", scale=", v[2L], ")")
  uniform <- function(v) paste0("Uniform(", paste(signif(v, 5), collapse = ", "), ")")
  mp <- priors$magnitude
  starts <- object$starting
  if(!is.null(object$magnitude)) {
    add("magnitude", "", "", "coefficients", "estimated",
        normal(mp$coefficient_variance), starts$magnitude$coefficients)
    add("magnitude", "", "", "size", if(is.null(mp$size)) "fixed" else "estimated",
        if(is.null(mp$size)) "" else paste0("Gamma(shape=", mp$size[1L], ", rate=", mp$size[2L], ")"),
        starts$magnitude$size %||% object$control$nb_size,
        "size estimation requires no magnitude process terms")
  }
  for(term in object$magnitude$process_terms){
    t <- term$name
    p <- list(variance = term$variance_prior)
    value <- starts$magnitude$terms[[t]]
    terms[[length(terms) + 1L]] <- data.frame(component = "magnitude", block = "",
                                             term = t, formula = term$label)
    add("magnitude", "", t, "effects", "estimated", "Gaussian process/group prior", value$effects)
    add("magnitude", "", t, "variance", "estimated", ig(p$variance), value$variance)
    if(!is.null(term$cov_model)){
      p$covariance <- term$covariance_priors
      for(j in seq_along(term$covariance_parameters)) for(k in seq_len(nrow(term$covariance_values))){
        estimated <- !term$covariance_fixed[k,j]
        add("magnitude", "", t, term$covariance_parameters[j],
            if(estimated) "estimated" else "fixed",
            if(estimated) uniform(c(term$covariance_lower[k,j], term$covariance_upper[k,j])) else "",
            term$covariance_values[k,j],
            paste(term$cov_model, rownames(term$covariance_values)[k]))
      }
    }
    priors$magnitude$terms[[t]] <- p
  }
  cp <- priors$composition
  for(b in names(object$composition$allocation_components)){
    block <- object$composition$allocation_components[[b]]
    bp <- cp$blocks[[b]] %||% cp[setdiff(names(cp), "blocks")]
    value <- starts$composition$blocks[[b]]
    add("composition", b, "", "coefficients", "estimated",
        paste0("intercept ", normal(bp$intercept_variance), "; slopes ", normal(bp$slope_variance)),
        value$coefficients)
    add("composition", b, "", "intercept_smoothing_variance", "fixed",
        "", bp$intercept_smoothing_variance, object$data$allocation_blocks[[b]]$smooth$intercept)
    add("composition", b, "", "slope_smoothing_variance", "fixed",
        "", bp$slope_smoothing_variance, object$data$allocation_blocks[[b]]$smooth$slope)
    bp$terms <- list()
    for(term in block$process_terms){
      t <- term$name
      part <- value$terms[[t]]
      p <- list(loading_variance = term$loading_variance)
      terms[[length(terms) + 1L]] <- data.frame(component = "composition", block = b,
                                               term = t, formula = term$label)
      add("composition", b, t, "loadings", "estimated with constraints",
          paste0(normal(term$loading_variance),
                 "; anchor diagonal: same normal restricted to (0, Inf)"), part$loadings,
          paste0("anchor diagonal > 0; upper anchor entries = 0; anchors: ",
                 paste(term$anchors, collapse = ", ")))
      add("composition", b, t, "scores", "estimated", "unit-scale Gaussian factors", part$scores)
      add("composition", b, t, "factor_scale", "fixed", "", 1)
      if(!is.null(term$cov_model)){
        p$covariance <- term$covariance_priors
        for(j in seq_along(term$covariance_parameters)) for(k in seq_len(nrow(term$covariance_values))){
          estimated <- !term$covariance_fixed[k,j]
          add("composition", b, t, term$covariance_parameters[j],
              if(estimated) "estimated" else "fixed",
              if(estimated) uniform(c(term$covariance_lower[k,j], term$covariance_upper[k,j])) else "",
              term$covariance_values[k,j],
              paste(term$cov_model, rownames(term$covariance_values)[k]))
        }
      }
      bp$terms[[t]] <- p
    }
    priors$composition$blocks[[b]] <- bp
  }
  structure(list(parameters = do.call(rbind, rows),
                 terms = if(length(terms)) do.call(rbind, terms) else
                   data.frame(component = character(), block = character(),
                              term = character(), formula = character()),
                 priors = priors, starting = starts, tuning = object$tuning), class = "magCompSettings")
}

print.magCompSettings <- function(x, ...)
{
  cat("magComp parameter settings\n")
  if(nrow(x$terms)){
    cat("\nFormula term names\n")
    print(x$terms, row.names = FALSE)
  }
  cat("\nResolved priors, initial values, and constraints\n")
  print(x$parameters, row.names = FALSE)
  cat("\nEdit $priors or $starting and pass them to magComp().\n")
  invisible(x)
}
