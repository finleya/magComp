magcomp_observed_stick_eta <- function(y, link = c("probit", "logit"))
{
  link <- match.arg(link)
  n <- nrow(y)
  J <- ncol(y)
  Jstick <- J - 1L
  eta <- matrix(NA_real_, nrow = n, ncol = Jstick,
                dimnames = list(rownames(y), colnames(y)[seq_len(Jstick)]))
  remaining <- rowSums(y)

  for(j in seq_len(Jstick)){
    keep <- remaining > 0
    prob <- (y[keep,j] + 0.5) / (remaining[keep] + 1)
    eta[keep,j] <- if(identical(link, "probit")) stats::qnorm(prob) else
      stats::qlogis(prob)
    remaining <- remaining - y[,j]
  }

  eta
}

magcomp_constrain_factor_start <- function(loading, factor, anchor_index,
                                           preserve_surface = FALSE)
{
  K <- ncol(loading)
  if(isTRUE(preserve_surface) && !is.null(factor) &&
     qr(loading[anchor_index,,drop = FALSE])$rank == K){
    loading_names <- dimnames(loading)
    factor_names <- dimnames(factor)
    # LQ decomposition of the anchor block: a common orthogonal rotation
    # preserves factor %*% t(loading), without fixing any loading scale.
    rotation <- qr.Q(qr(t(loading[anchor_index,,drop = FALSE])))
    loading <- loading %*% rotation
    factor <- factor %*% rotation
    signs <- sign(diag(loading[anchor_index,,drop = FALSE]))
    loading <- sweep(loading, 2L, signs, "*")
    factor <- sweep(factor, 2L, signs, "*")
    # Allocate the empirical amplitude to the loadings, matching the unit
    # factor prior scale. This only sets starts; sampled scores are not scaled.
    scale <- sqrt(colMeans(factor^2))
    scale[!is.finite(scale) | scale <= 0] <- 1
    loading <- sweep(loading, 2L, scale, "*")
    factor <- sweep(factor, 2L, scale, "/")
    dimnames(loading) <- loading_names
    dimnames(factor) <- factor_names
  }

  for(k in seq_len(K)){
    anchor <- anchor_index[k]
    if(k > 1L)
      loading[anchor_index[seq_len(k - 1L)], k] <- 0
    # Supplied valid loadings are unchanged. Automatic degenerate starts
    # need a positive initial value, not a lower bound on posterior draws.
    if(!is.finite(loading[anchor,k]) || loading[anchor,k] == 0)
      loading[anchor,k] <- 1
    else if(loading[anchor,k] < 0)
      loading[anchor,k] <- -loading[anchor,k]
  }

  list(loading = loading, factor = factor)
}

magcomp_constrain_factor_starts <- function(composition, starting)
{
  if(is.null(starting$composition_block_loading))
    return(starting)

  factor_terms <- c("factor_iid", "factor_nngp", "factor_gp")
  for(block_name in names(composition$allocation_components)){
    terms <- composition$allocation_components[[block_name]]$process_terms
    if(length(terms) == 0L)
      next

    for(term in terms){
      if(!term$term_type %in% factor_terms)
        next
      loading <- starting$composition_block_loading[[block_name]][[term$name]]
      if(is.null(loading))
        next

      factor <- NULL
      if(!is.null(starting$composition_block_score[[block_name]]))
        factor <- starting$composition_block_score[[block_name]][[term$name]]

      constrained <- magcomp_constrain_factor_start(
        loading,
        factor,
        term$anchor_index
      )
      starting$composition_block_loading[[block_name]][[term$name]] <-
        constrained$loading
      if(!is.null(factor)){
        starting$composition_block_score[[block_name]][[term$name]] <-
          constrained$factor
      }
    }
  }

  starting
}

magcomp_factor_fixed_start <- function(eta, X, coef_start = NULL)
{
  residual <- matrix(0, nrow = nrow(eta), ncol = ncol(eta),
                     dimnames = dimnames(eta))
  if(is.null(X))
    X <- matrix(1, nrow = nrow(eta), ncol = 1,
                dimnames = list(rownames(eta), "(Intercept)"))
  if(nrow(X) != nrow(eta))
    return(list(coef = NULL, residual = residual))

  coef <- NULL
  if(!is.null(coef_start)){
    if(length(coef_start) != ncol(eta) * ncol(X)){
      stop("error: starting composition coefficients must have length ",
           "coef_stick_count * ncol(X_composition)")
    }
    coef <- matrix(as.numeric(coef_start),
                   nrow = ncol(eta),
                   ncol = ncol(X),
                   dimnames = list(colnames(eta), colnames(X)))
    if(anyNA(coef) || any(!is.finite(coef))){
      stop("error: starting composition coefficients must contain finite ",
           "values")
    }
  }else{
    coef <- matrix(0, nrow = ncol(eta), ncol = ncol(X),
                   dimnames = list(colnames(eta), colnames(X)))
  }

  if(is.null(coef_start) && ncol(X) > 0L){
    intercept <- match("(Intercept)", colnames(X), nomatch = 0L)
    complete_x <- stats::complete.cases(X)
    for(j in seq_len(ncol(eta))){
      keep <- is.finite(eta[,j]) & complete_x
      if(sum(keep) > ncol(X)){
        fit <- stats::lm.fit(X[keep,,drop = FALSE], eta[keep,j])
        keep_coef <- is.finite(fit$coefficients)
        coef[j,keep_coef] <- fit$coefficients[keep_coef]
      }else if(any(keep) && intercept > 0L){
        value <- mean(eta[keep,j])
        if(is.finite(value))
          coef[j,intercept] <- value
      }
    }
  }

  if(ncol(X) > 0L){
    residual <- eta - X %*% t(coef)
  }else{
    residual <- eta
  }

  residual[!is.finite(residual)] <- 0
  dimnames(residual) <- dimnames(eta)

  list(coef = coef, residual = residual)
}

magcomp_factor_eta_residual <- function(eta, X)
{
  magcomp_factor_fixed_start(eta, X)$residual
}

magcomp_factor_block_fixed_start <- function(eta, X, coef_start, block)
{
  if(!isTRUE(block$parent_specific_coef))
    return(magcomp_factor_fixed_start(eta, X, coef_start))

  Jstick <- ncol(eta)
  coef <- matrix(0, nrow = block$coef_stick_count, ncol = ncol(X),
                  dimnames = list(block$stick_names, colnames(X)))
  residual <- eta
  for(a in seq_along(block$parent_levels)){
    rows <- which(block$parent_index == a)
    sticks <- (a - 1L) * Jstick + seq_len(Jstick)
    supplied <- if(is.null(coef_start)) NULL else
      coef_start[sticks,,drop = FALSE]
    fixed <- magcomp_factor_fixed_start(eta[rows,,drop = FALSE],
                                        X[rows,,drop = FALSE], supplied)
    coef[sticks,] <- fixed$coef
    residual[rows,] <- fixed$residual
  }
  # Missing continuation ratios provide no information to group projections.
  residual[!is.finite(eta)] <- NA_real_
  list(coef = coef, residual = residual)
}

magcomp_factor_start_terms <- function(terms)
{
  factor_terms <- Filter(
    function(term) term$term_type %in% c("factor_iid", "factor_nngp",
                                         "factor_gp"),
    terms
  )
  if(length(factor_terms) == 0L)
    return(factor_terms)

  c(
    Filter(function(term) identical(term$term_type, "factor_iid"),
           factor_terms),
    Filter(function(term) !identical(term$term_type, "factor_iid"),
           factor_terms)
  )
}

magcomp_factor_group_index <- function(term, n)
{
  if(identical(term$term_type, "factor_iid"))
    return(as.integer(term$group))

  as.integer(term$spatial_group %||% seq_len(n))
}

magcomp_factor_group_names <- function(term)
{
  if(identical(term$term_type, "factor_iid"))
    return(term$group_levels)

  term$location_levels %||% as.character(seq_len(term$n_group))
}

magcomp_factor_start_matrix <- function(value, n_row, n_col, row_names,
                                        col_names, label)
{
  if(is.null(value))
    return(NULL)

  if(length(value) != n_row * n_col)
    stop("error: starting ", label, " must have length ", n_row * n_col)

  out <- matrix(as.numeric(value),
                nrow = n_row,
                ncol = n_col,
                dimnames = list(row_names, col_names))
  if(anyNA(out) || any(!is.finite(out)))
    stop("error: starting ", label, " must contain finite values")

  out
}

magcomp_factor_project_residual <- function(residual, term, block)
{
  n <- nrow(residual)
  group <- magcomp_factor_group_index(term, n)
  n_group <- term$n_group
  covariate <- term$covariate %||% rep(1, n)
  Jstick <- ncol(residual)
  conditional <- isTRUE(block$parent_specific_coef)
  surface <- matrix(0, nrow = n_group, ncol = block$coef_stick_count,
                    dimnames = list(magcomp_factor_group_names(term),
                                    block$stick_names))

  for(g in seq_len(n_group)){
    keep_g <- group == g & is.finite(covariate)
    if(!any(keep_g))
      next
    for(c in seq_len(ncol(surface))){
      j <- (c - 1L) %% Jstick + 1L
      keep <- keep_g & is.finite(residual[,j])
      if(conditional)
        keep <- keep & block$parent_index == (c - 1L) %/% Jstick + 1L
      denom <- sum(covariate[keep]^2)
      if(is.finite(denom) && denom > .Machine$double.eps){
        surface[g,c] <- sum(covariate[keep] * residual[keep,j]) / denom
      }
    }
  }

  surface[!is.finite(surface)] <- 0
  surface
}

magcomp_factor_surface_svd_start <- function(surface, term)
{
  K <- term$K
  loading <- matrix(0, nrow = ncol(surface), ncol = K,
                    dimnames = list(colnames(surface), term$factor_names))
  factor <- matrix(0, nrow = nrow(surface), ncol = K,
                   dimnames = list(rownames(surface), term$factor_names))

  if(any(abs(surface) > 0)){
    sv <- svd(surface, nu = K, nv = K)
    r <- min(K, length(sv$d), ncol(sv$u), ncol(sv$v))
    if(r > 0L){
      scale <- sqrt(sv$d[seq_len(r)])
      factor[,seq_len(r)] <- sweep(sv$u[,seq_len(r),drop = FALSE],
                                   2L, scale, "*")
      loading[,seq_len(r)] <- sweep(sv$v[,seq_len(r),drop = FALSE],
                                    2L, scale, "*")
    }
  }

  magcomp_constrain_factor_start(loading, factor, term$anchor_index,
                                 preserve_surface = TRUE)
}

magcomp_factor_from_loading <- function(surface, loading, factor_names)
{
  fit <- tryCatch(stats::lm.fit(loading, t(surface)),
                  error = function(e) NULL)
  factor <- matrix(0, nrow = nrow(surface), ncol = ncol(loading),
                   dimnames = list(rownames(surface), factor_names))
  if(!is.null(fit)){
    value <- matrix(fit$coefficients,
                    nrow = ncol(loading),
                    ncol = nrow(surface))
    value <- t(value)
    value[!is.finite(value)] <- 0
    factor[,] <- value
  }

  factor
}

magcomp_loading_from_factor <- function(surface, factor)
{
  fit <- tryCatch(stats::lm.fit(factor, surface),
                  error = function(e) NULL)
  loading <- matrix(0, nrow = ncol(surface), ncol = ncol(factor),
                    dimnames = list(colnames(surface), colnames(factor)))
  if(!is.null(fit)){
    value <- matrix(fit$coefficients,
                    nrow = ncol(factor),
                    ncol = ncol(surface))
    value <- t(value)
    value[!is.finite(value)] <- 0
    loading[,] <- value
  }

  loading
}

magcomp_factor_term_contribution <- function(term, loading, factor, residual,
                                             block)
{
  group <- magcomp_factor_group_index(term, nrow(residual))
  covariate <- term$covariate %||% rep(1, nrow(residual))
  surface <- factor %*% t(loading)
  if(isTRUE(block$parent_specific_coef)){
    Jstick <- ncol(residual)
    out <- matrix(0, nrow = nrow(residual), ncol = Jstick)
    for(j in seq_len(Jstick)){
      sticks <- (block$parent_index - 1L) * Jstick + j
      out[,j] <- surface[cbind(group, sticks)]
    }
  }else{
    out <- surface[group,,drop = FALSE]
  }
  out <- sweep(out, 1L, covariate, "*")
  out[!is.finite(out)] <- 0
  dimnames(out) <- dimnames(residual)
  out
}

magcomp_set_factor_start <- function(starting, block_name, term, loading,
                                     factor)
{
  if(is.null(starting$composition_block_loading))
    starting$composition_block_loading <- list()
  if(is.null(starting$composition_block_loading[[block_name]]))
    starting$composition_block_loading[[block_name]] <- list()
  starting$composition_block_loading[[block_name]][[term$name]] <- loading

  if(is.null(starting$composition_block_score))
    starting$composition_block_score <- list()
  if(is.null(starting$composition_block_score[[block_name]]))
    starting$composition_block_score[[block_name]] <- list()
  starting$composition_block_score[[block_name]][[term$name]] <- factor

  starting
}

magcomp_set_composition_coef_start <- function(starting, block_name, coef)
{
  if(is.null(coef))
    return(starting)
  if(is.null(starting$composition_block_coefficient))
    starting$composition_block_coefficient <- list()
  if(is.null(starting$composition_block_coefficient[[block_name]]))
    starting$composition_block_coefficient[[block_name]] <- coef

  starting
}

magcomp_starting_factor_term <- function(starting, block_name, term,
                                         residual, block)
{
  loading_start <- NULL
  if(!is.null(starting$composition_block_loading[[block_name]])){
    loading_start <-
      starting$composition_block_loading[[block_name]][[term$name]]
  }
  factor_start <- NULL
  if(!is.null(starting$composition_block_score[[block_name]])){
    factor_start <-
      starting$composition_block_score[[block_name]][[term$name]]
  }

  loading_start <- magcomp_factor_start_matrix(
    loading_start,
    n_row = block$coef_stick_count,
    n_col = term$K,
    row_names = block$stick_names,
    col_names = term$factor_names,
    label = paste("loading values for", term$name)
  )
  factor_start <- magcomp_factor_start_matrix(
    factor_start,
    n_row = term$n_group,
    n_col = term$K,
    row_names = magcomp_factor_group_names(term),
    col_names = term$factor_names,
    label = paste("factor values for", term$name)
  )

  if(!is.null(loading_start) && !is.null(factor_start)){
    return(list(loading = loading_start, factor = factor_start))
  }

  surface <- magcomp_factor_project_residual(residual, term, block)
  if(!is.null(loading_start)){
    loading <- magcomp_constrain_factor_start(
      loading_start,
      NULL,
      term$anchor_index
    )$loading
    factor <- magcomp_factor_from_loading(surface, loading,
                                          term$factor_names)
    return(list(loading = loading, factor = factor))
  }

  if(!is.null(factor_start)){
    loading <- magcomp_loading_from_factor(surface, factor_start)
    loading <- magcomp_constrain_factor_start(loading, NULL, term$anchor_index)$loading
    return(list(loading = loading, factor = factor_start))
  }

  magcomp_factor_surface_svd_start(surface, term)
}

magcomp_factor_svd_start <- function(data, composition, starting,
                                     block_name = NULL,
                                     composition_link = "probit")
{
  if(is.null(block_name)){
    for(name in names(composition$allocation_components))
      starting <- magcomp_factor_svd_start(
        data, composition, starting, name, composition_link
      )
    return(starting)
  }

  component <- composition$allocation_components[[block_name]]
  block <- data$allocation_blocks[[block_name]]
  if(is.null(component) || is.null(block))
    return(starting)
  terms <- magcomp_factor_start_terms(component$process_terms)
  if(length(terms) == 0L)
    return(starting)
  if(is.null(data$counts))
    return(starting)

  y <- block$counts
  eta <- magcomp_observed_stick_eta(y, composition_link)
  coef_start <- NULL
  if(!is.null(starting$composition_block_coefficient[[block_name]]))
    coef_start <- starting$composition_block_coefficient[[block_name]]
  fixed_start <- magcomp_factor_block_fixed_start(eta, component$X,
                                                  coef_start, block)
  starting <- magcomp_set_composition_coef_start(starting, block_name,
                                                fixed_start$coef)
  residual <- fixed_start$residual
  available <- is.finite(residual)

  for(term in terms){
    term_start <- magcomp_starting_factor_term(starting, block_name, term,
                                               residual, block)
    starting <- magcomp_set_factor_start(starting, block_name, term,
                                        term_start$loading,
                                        term_start$factor)
    residual <- residual - magcomp_factor_term_contribution(
      term,
      term_start$loading,
      term_start$factor,
      residual,
      block
    )
    residual[!is.finite(residual)] <- 0
    residual[!available] <- NA_real_
  }

  starting
}
