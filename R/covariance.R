# Covariance metadata comes from the same registry used by both native kernels.
magCompCovModels <- function()
{
  .Call("magcomp_cov_models", PACKAGE = "magComp")
}

uniform <- function(lower, upper)
{
  if(!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1L ||
     length(upper) != 1L || !is.finite(lower) || !is.finite(upper) || lower >= upper)
    stop("error: uniform() requires finite scalar bounds with lower < upper")
  structure(c(lower = lower, upper = upper), class = "magCompUniform")
}

fixed <- function(value)
{
  if(!is.numeric(value) || !length(value) || anyNA(value) || any(!is.finite(value)))
    stop("error: fixed() requires finite numeric values")
  structure(value, class = "magCompFixed")
}

magcomp_prepare_covariance <- function(term, coords, label)
{
  model <- term$params$cov_model %||% "exp"
  registry <- magCompCovModels()
  if(!is.character(model) || length(model) != 1L || !model %in% names(registry))
    stop("error: unknown cov_model; choose ", paste(names(registry), collapse = ", "))
  term$cov_model <- model
  term$space_time <- registry[[model]]$space_time
  term$covariance_parameters <- registry[[model]]$parameters
  term$st_scale <- term$params$st_scale %||% 1
  if(!is.numeric(term$st_scale) || length(term$st_scale) != 1L ||
     !is.finite(term$st_scale) || term$st_scale <= 0)
    stop("error: st_scale must be a positive finite number")
  if(!term$space_time && term$st_scale != 1)
    stop("error: st_scale applies only to space-time covariance models")
  term
}

magcomp_search_coords <- function(term, coords = term$coords)
{
  if(isTRUE(term$space_time)) coords[,ncol(coords)] <- coords[,ncol(coords)] * term$st_scale
  coords
}

magcomp_covariance_defaults <- function(term)
{
  coords <- unique(term$coords)
  space <- if(term$space_time) coords[,-ncol(coords),drop = FALSE] else coords
  d <- as.matrix(stats::dist(unique(space)))
  D <- max(d)
  spatial <- NULL
  multi <- NULL
  if(is.finite(D) && D > 0){
    diag(d) <- Inf
    nearest <- apply(d, 1L, min)
    nearest <- nearest[is.finite(nearest) & nearest > 0]
    rmin <- max(0.5 * stats::median(nearest), sqrt(.Machine$double.eps) * D)
    spatial <- uniform(3 / (2 * D), 3 / min(rmin, 0.5 * D))
    multi <- list(phi_1 = uniform(3 / D, 30 / D),
                  phi_2 = uniform(30 / D, 3 / min(rmin, 0.05 * D)))
  }
  temporal <- NULL
  T <- tmin <- NULL
  if(term$space_time){
    times <- sort(unique(coords[,ncol(coords)]))
    T <- diff(range(times))
    if(T > 0){
      tmin <- 0.5 * min(diff(times))
      temporal <- uniform(3 / (2 * T), 3 / tmin)
    }
  }
  priors <- switch(term$cov_model,
    exp = list(phi = spatial),
    matern = list(phi = spatial, nu = uniform(0.05, 2)),
    sep_exp = list(phi = spatial, lambda = temporal),
    multi_res_sep_exp = list(alpha = uniform(0, 1), phi_1 = multi$phi_1,
      lambda_1 = temporal, phi_2 = multi$phi_2, lambda_2 = temporal),
    gneiting = list(a = if(!is.null(temporal)) {
      power <- ncol(space) / 2
      uniform((0.05^(-1 / power) - 1) / (2 * T)^2,
              (0.05^(-1 / power) - 1) / tmin^2)
    }, c = spatial, alpha = uniform(0, 1), beta = uniform(0, 1),
      gamma = uniform(0, 1), delta = uniform(0, 2)))
  starts <- lapply(priors, function(x) if(!is.null(x)) mean(x) else NULL)
  if(term$cov_model %in% c("exp", "matern", "sep_exp") && D > 0) starts$phi <- 6 / D
  if(term$cov_model == "gneiting"){
    starts$alpha <- fixed(1)
    starts$gamma <- fixed(0.5)
    starts$delta <- fixed(0)
  }
  list(priors = priors, starting = starts)
}

# Common scalar settings apply to all factors; named lists/vectors override
# individual factors. Missing entries inherit their automatic defaults.
magcomp_by_factor <- function(value, factors, prior = FALSE)
{
  out <- rep(list(NULL), length(factors))
  names(out) <- factors
  if(is.null(value)) return(out)
  if(prior && inherits(value, "magCompUniform")) return(stats::setNames(rep(list(value), length(factors)), factors))
  if(!prior && is.numeric(value)){
    is_fixed <- inherits(value, "magCompFixed")
    if(length(value) == 1L && is.null(names(value))) value <- rep(value, length(factors))
    entries <- lapply(seq_along(value), function(i)
      if(is_fixed) fixed(as.numeric(value[i])) else as.numeric(value[i]))
    if(!is.null(names(value))) names(entries) <- names(value)
    value <- entries
  }
  if(!is.list(value)) stop("error: covariance settings require a value or a named list by factor")
  if(is.null(names(value)) && length(value) == length(factors)) names(value) <- factors
  magcomp_check_names(value, factors, "covariance factors")
  out[names(value)] <- value
  out
}

magcomp_covariance_list <- function(values, fixed_values, factors)
{
  out <- lapply(seq_len(ncol(values)), function(p){
    x <- stats::setNames(values[,p], factors)
    flags <- fixed_values[,p]
    if(all(flags)) return(fixed(if(length(x) == 1L) unname(x) else x))
    if(!any(flags)) return(if(length(x) == 1L) unname(x) else x)
    stats::setNames(lapply(seq_along(x), function(k)
      if(flags[k]) fixed(unname(x[k])) else unname(x[k])), factors)
  })
  stats::setNames(out, colnames(values))
}

magcomp_resolve_covariance <- function(term, priors = list(), starting = list(), tuning = list(), path)
{
  factors <- term$factor_names %||% "process"
  K <- length(factors)
  model <- magCompCovModels()[[term$cov_model]]
  pars <- model$parameters
  for(x in list(priors, starting, tuning)) magcomp_check_names(x, pars, path)
  defaults <- magcomp_covariance_defaults(term)
  values <- lower <- upper <- step <- matrix(NA_real_, K, length(pars), dimnames = list(factors, pars))
  held <- matrix(FALSE, K, length(pars), dimnames = dimnames(values))
  prior_list <- list()
  for(p in seq_along(pars)){
    name <- pars[p]
    pr <- magcomp_by_factor(priors[[name]], factors, prior = TRUE)
    init <- magcomp_by_factor(starting[[name]], factors)
    tune <- magcomp_by_factor(tuning[[name]], factors)
    prior_list[[name]] <- stats::setNames(vector("list", K), factors)
    for(k in seq_len(K)){
      where <- paste0(path, "$", name, "[", factors[k], "]")
      explicit_start <- !is.null(init[[k]])
      value <- init[[k]] %||% if(is.null(pr[[k]])) defaults$starting[[name]] else NULL
      if(!is.null(value) && (!is.numeric(value) || length(value) != 1L || !is.finite(value)))
        stop("error: covariance starting values must be finite scalars for ", where)
      if(!is.null(tune[[k]]) && (!is.numeric(tune[[k]]) || length(tune[[k]]) != 1L ||
                               !is.finite(tune[[k]]) || tune[[k]] <= 0 || tune[[k]] > sqrt(.Machine$double.xmax)))
        stop("error: covariance tuning must be a positive finite standard deviation for ", where)
      held[k,p] <- inherits(value, "magCompFixed")
      bounds <- pr[[k]] %||% defaults$priors[[name]]
      domain <- model$domains[p]
      if(held[k,p]){
        if(!is.null(pr[[k]])) stop("error: prior supplied for fixed covariance parameter ", where)
        if(!is.null(tune[[k]])) stop("error: tuning supplied for fixed covariance parameter ", where)
        lower[k,p] <- upper[k,p] <- as.numeric(value)
        step[k,p] <- 0
      }else{
        if(is.null(bounds)) stop("error: cannot choose automatic support for ", where,
          "; supply a prior or fix the parameter when coordinates have no spatial or temporal span")
        if(!inherits(bounds, "magCompUniform")) stop("error: ", where, " prior must use uniform(lower, upper)")
        if((domain == "positive" && bounds[1] < 0) ||
           (domain %in% c("unit", "positive_unit") && (bounds[1] < 0 || bounds[2] > 1)) ||
           (domain == "nonnegative" && bounds[1] < 0))
          stop("error: prior support outside the covariance parameter domain for ", where)
        lower[k,p] <- bounds[1]
        upper[k,p] <- bounds[2]
        # An explicit ordinary start releases a default-fixed shape parameter.
        if(is.null(value)) value <- mean(bounds)
        if(!explicit_start && (value <= bounds[1] || value >= bounds[2])) value <- mean(bounds)
        if(value <= bounds[1] || value >= bounds[2])
          stop("error: starting value must be strictly inside prior support for ", where)
        step[k,p] <- tune[[k]] %||% 0.1
        if(!is.finite(step[k,p]) || step[k,p] <= 0) stop("error: covariance tuning must be positive for ", where)
        prior_list[[name]][[k]] <- bounds
      }
      value <- as.numeric(value)
      if(length(value) != 1L || !is.finite(value) ||
         (domain %in% c("positive", "positive_unit") && value <= 0) ||
         (domain %in% c("unit", "positive_unit") && (value < 0 || value > 1)) ||
         (domain == "nonnegative" && value < 0))
        stop("error: invalid covariance starting/fixed value for ", where)
      values[k,p] <- value
    }
    prior_list[[name]] <- Filter(Negate(is.null), prior_list[[name]])
    if(!length(prior_list[[name]])) prior_list[name] <- list(NULL)
    if(K == 1L && length(prior_list[[name]])) prior_list[[name]] <- prior_list[[name]][[1L]]
  }
  term$covariance_values <- values
  term$covariance_lower <- lower
  term$covariance_upper <- upper
  term$covariance_fixed <- held
  term$covariance_tuning <- step
  term$covariance_priors <- Filter(Negate(is.null), prior_list)
  term$covariance_starting <- magcomp_covariance_list(values, held, factors)
  term
}

magcomp_correlation <- function(coords, cov_model = "exp", parameters)
{
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if(any(!is.finite(coords))) stop("error: coordinates must be finite")
  model <- magCompCovModels()[[cov_model]]
  if(is.null(model)) stop("error: unknown covariance model")
  if(!is.null(names(parameters))) parameters <- parameters[model$parameters]
  .Call("magcomp_cor_matrix", coords, cov_model, as.double(parameters), PACKAGE = "magComp")
}

# Public arrays retain only estimated entries; the backend retains all values
# so prediction can use fixed and estimated parameters in the same kernel.
magcomp_split_covariance_samples <- function(x, term)
{
  out <- list()
  K <- nrow(term$covariance_values)
  for(p in seq_along(term$covariance_parameters)){
    keep <- which(!term$covariance_fixed[,p])
    if(!length(keep)) next
    index <- c(list(x, keep + K * (p - 1L)), rep(list(TRUE), length(dim(x)) - 1L), list(drop = FALSE))
    value <- do.call(`[`, index)
    dimnames(value)[[1L]] <- rownames(term$covariance_values)[keep]
    out[[term$covariance_parameters[p]]] <- value
  }
  out
}

magcomp_covariance_tuning_list <- function(term)
{
  out <- list()
  for(p in seq_along(term$covariance_parameters)){
    keep <- which(!term$covariance_fixed[,p])
    if(!length(keep)) next
    value <- term$covariance_tuning[keep,p]
    if(nrow(term$covariance_tuning) == 1L) value <- unname(value)
    out[[term$covariance_parameters[p]]] <- value
  }
  out
}
