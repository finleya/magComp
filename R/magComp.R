magcomp_observed_magnitude_start <- function(data, magnitude)
{
  if(is.null(data$counts))
    return(NULL)

  y <- data$counts
  area <- data$sampler_area %||% data$area
  if(is.null(area))
    return(NULL)

  eta <- log((rowSums(y) + 0.5) / area)
  if(any(!is.finite(eta)))
    return(NULL)

  X <- magnitude$X
  fit <- tryCatch(stats::lm.fit(X, eta), error = function(e) NULL)
  magnitude_coefficients <- rep(0, ncol(X))
  if(!is.null(fit)){
    keep <- is.finite(fit$coefficients)
    magnitude_coefficients[keep] <- fit$coefficients[keep]
  }
  names(magnitude_coefficients) <- colnames(X)

  list(eta = eta, magnitude_coefficients = magnitude_coefficients)
}

magcomp_sampler_area <- function(data, family)
{
  if(identical(family$name, "fixed_area_stems"))
    return(data$area)

  if(identical(family$name, "point_basal_area")){
    baf <- data$baf
    if(is.null(baf))
      return(NULL)
    baf <- as.numeric(baf)
    if(anyNA(baf) || any(!is.finite(baf)) || any(baf <= 0))
      return(NULL)
    out <- 1 / baf
    names(out) <- data$unit_levels
    return(out)
  }

  NULL
}

magcomp_start_coefficient_vector <- function(magnitude_coefficients, default, p)
{
  out <- default
  if(is.null(magnitude_coefficients))
    return(out)

  magnitude_coefficients <- as.numeric(magnitude_coefficients)
  if(length(magnitude_coefficients) == 1L && is.finite(magnitude_coefficients)){
    out[1L] <- magnitude_coefficients
  }else if(length(magnitude_coefficients) == p && all(is.finite(magnitude_coefficients))){
    out <- magnitude_coefficients
  }
  names(out) <- names(default)
  out
}

magcomp_magnitude_residual_start <- function(data, magnitude, starting)
{
  observed <- magcomp_observed_magnitude_start(data, magnitude)
  if(is.null(observed))
    return(starting)

  if(is.null(starting$magnitude_coefficients))
    starting$magnitude_coefficients <- observed$magnitude_coefficients

  if(length(magnitude$process_terms) < 1L)
    return(starting)

  valid_terms <- Filter(function(term)
    term$term_type %in% c("iid", "nngp", "gp"),
    magnitude$process_terms)
  if(length(valid_terms) < 1L)
    return(starting)

  X <- magnitude$X
  magnitude_coefficients <- magcomp_start_coefficient_vector(starting$magnitude_coefficients, observed$magnitude_coefficients,
                                      ncol(X))
  residual <- observed$eta - drop(X %*% magnitude_coefficients)
  residual[!is.finite(residual)] <- 0

  missing_terms <- Filter(function(term)
    is.null(starting$magnitude_process[[term$name]]),
    valid_terms)
  residual_share <- if(length(missing_terms) > 0L)
    residual / length(missing_terms) else residual

  if(is.null(starting$magnitude_process))
    starting$magnitude_process <- list()

  for(term in valid_terms){
    covariate <- term$covariate %||% rep(1, data$n_unit)

    if(identical(term$term_type, "iid")){
      if(is.null(starting$magnitude_process[[term$name]])){
        value <- numeric(term$n_group)
        for(g in seq_len(term$n_group)){
          keep <- term$group == g
          denom <- sum(covariate[keep]^2)
          if(is.finite(denom) && denom > sqrt(.Machine$double.eps)){
            value[g] <- sum(covariate[keep] * residual_share[keep]) /
              denom
          }
        }
        names(value) <- term$group_levels
        starting$magnitude_process[[term$name]] <- value
      }
    }

    if(is.null(starting$magnitude_process[[term$name]])){
      effect <- residual_share
      if(!is.null(term$covariate)){
        keep <- abs(covariate) > sqrt(.Machine$double.eps)
        effect[keep] <- residual_share[keep] / covariate[keep]
        effect[!keep] <- 0
      }
      names(effect) <- data$unit_levels
      starting$magnitude_process[[term$name]] <- effect
    }

    if(is.null(starting$magnitude_process_variance[[term$name]])){
      variance <- stats::var(as.numeric(
        starting$magnitude_process[[term$name]]
      ))
      if(!is.finite(variance) || variance <= 0)
        variance <- 1
      if(is.null(starting$magnitude_process_variance))
        starting$magnitude_process_variance <- list()
      starting$magnitude_process_variance[[term$name]] <- max(variance, 1e-4)
    }
  }

  starting
}

magcomp_prepare_starting <- function(data, magnitude, composition, starting,
                                     composition_link = "probit")
{
  if(!is.null(magnitude))
    starting <- magcomp_magnitude_residual_start(data, magnitude, starting)
  starting <- magcomp_factor_svd_start(data, composition, starting,
                                       composition_link = composition_link)

  magcomp_constrain_factor_starts(composition, starting)
}

magcomp_block_smooth_spec <- function(block, smooth)
{
  spec <- list(intercept = smooth$intercept,
               slope = smooth$slope,
               loading = smooth$loading,
               ordered = smooth$ordered)
  block_specs <- smooth$blocks %||% list()
  if(length(block_specs) > 0L){
    block_spec <- NULL
    if(!is.null(block_specs[[block$name]])){
      block_spec <- block_specs[[block$name]]
    }else if(!is.null(block_specs[[block$child]])){
      block_spec <- block_specs[[block$child]]
    }
    if(!is.null(block_spec)){
      provided <- attr(block_spec, "provided") %||%
        intersect(names(block_spec), c("intercept", "slope", "loading"))
      for(nm in intersect(provided, c("intercept", "slope", "loading"))){
        spec[[nm]] <- block_spec[[nm]] %||% "none"
      }
    }
  }

  spec$ordered <- isTRUE(block$ordered)
  if(!spec$ordered &&
     (spec$intercept != "none" || spec$slope != "none" ||
      spec$loading != "none")){
    stop("error: smoothing requested for nominal allocation block '",
         block$name, "'")
  }
  spec
}

magcomp_attach_block_smoothing <- function(data, smooth)
{
  if(length(data$allocation_blocks) == 0L)
    return(data)

  block_names <- names(data$allocation_blocks)
  response_names <- vapply(data$allocation_blocks, `[[`, character(1L),
                           "child")
  valid_names <- unique(c(block_names, response_names))
  unknown <- setdiff(names(smooth$blocks %||% list()), valid_names)
  if(length(unknown) > 0L){
    stop("error: unknown smoothing block name(s): ",
         paste(unknown, collapse = ", "))
  }

  for(block_name in block_names){
    block <- data$allocation_blocks[[block_name]]
    block$smooth <- magcomp_block_smooth_spec(block, smooth)
    data$allocation_blocks[[block_name]] <- block
  }
  data
}

magComp <- function(data,
                    family = fixedAreaStems(),
                    magnitude = ~ 1,
                    composition = ~ 1,
                    smooth = magCompSmooth(),
                    priors = magCompPriors(),
                    starting = list(),
                    control = magCompControl(),
                    tuning = list(),
                    ...)
{
  if(!inherits(data, "magCompData"))
    data <- magCompData(data, ...)

  if(!inherits(family, "magCompFamily"))
    stop("error: family must be created by fixedAreaStems() ",
         "or pointBasalArea()")
  if(identical(family$name, "continuous_gaussian"))
    stop("error: continuousGaussian() is a reserved family and is not implemented")
  if(!inherits(smooth, "magCompSmooth"))
    stop("error: smooth must be created by magCompSmooth()")
  if(!inherits(priors, "magCompPriors"))
    stop("error: priors must be created by magCompPriors()")
  if(!inherits(control, "magCompControl"))
    stop("error: control must be created by magCompControl()")

  if(!identical(data$design$type, "unknown") &&
     !identical(family$design, "weighted_tree_list") &&
     !identical(data$design$type, family$design)){
    stop("error: data design '", data$design$type,
         "' is not compatible with family '", family$name, "'")
  }

  if(is.null(magnitude) && is.null(composition))
    stop("error: at least one of magnitude or composition must be specified")
  if(!is.null(composition) && !length(data$allocation_blocks))
    stop("error: composition requires category counts; use composition = NULL for totals only")
  if(is.null(composition) && !identical(smooth, magCompSmooth()))
    stop("error: smoothing is unused when composition = NULL")

  magnitude_component <- if(is.null(magnitude)) NULL else magcomp_build_component(
    formula = magnitude,
    data = data$unit_data,
    component = "magnitude"
  )
  composition_component <- if(is.null(composition)) NULL else magcomp_build_composition_component(
    composition = composition,
    data = data
  )
  for(component in Filter(Negate(is.null), c(list(magnitude_component),
                     composition_component$allocation_components))){
    columns <- c(all.vars(component$fixed_formula),
      unlist(lapply(component$process_terms, function(term)
        c(term$args, term$coef_name)), use.names = FALSE))
    magcomp_check_unit_columns(data$data, data$data[[data$unit_name]], columns)
  }
  if(!is.null(composition_component)){
    composition_component <- magcomp_finalize_component_process_terms(composition_component, data)
    data <- magcomp_attach_block_smoothing(data, smooth)
  }
  data$sampler_area <- magcomp_sampler_area(data, family)
  if(isTRUE(control$run_sampler))
    magcomp_check_sampler_support(data, family, magnitude_component,
                                 composition_component, smooth, priors)
  resolved <- magcomp_resolve_priors(priors, magnitude_component, composition_component, starting, tuning)
  priors <- resolved$priors
  tuning <- resolved$tuning
  magnitude_component <- resolved$magnitude
  composition_component <- resolved$composition
  starts <- magcomp_resolve_starting(data, magnitude_component, composition_component,
                                     starting, priors, control)
  starting <- starts$public

  config <- list(
    y = data$counts,
    allocation_blocks = if(is.null(composition_component)) list() else data$allocation_blocks,
    area = data$sampler_area,
    baf = data$baf,
    X_magnitude = magnitude_component$X,
    X_composition = composition_component$X,
    coords = data$coords,
    family = family,
    smooth = smooth,
    priors = resolved$kernel,
    starting = starts$kernel,
    control = control,
    magnitude_process = magnitude_component$process_terms,
    composition_process = composition_component$process_terms,
    composition_blocks = composition_component$allocation_components
  )

  backend <- NULL
  if(isTRUE(control$run_sampler)){
    backend <- magcomp_run_backend_chains(
      config = config,
      data = data,
      magnitude = magnitude_component,
      composition = composition_component
    )
  }

  out <- list(
    call = match.call(),
    data = data,
    family = family,
    magnitude = magnitude_component,
    composition = composition_component,
    smooth = smooth,
    priors = priors,
    starting = starting,
    tuning = tuning,
    control = control,
    config = config,
    backend = backend,
    fitted = isTRUE(control$run_sampler)
  )

  class(out) <- "magComp"
  out
}
