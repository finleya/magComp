print.magComp <- function(x, ...)
{
  cat("magComp model\n")
  cat("  family:       ", x$family$name, "\n", sep = "")
  cat("  units:        ", x$data$n_unit, "\n", sep = "")
  if(!is.null(x$data$counts))
    cat("  categories:   ", ncol(x$data$counts), "\n", sep = "")
  cat("  magnitude:    ",
      paste(deparse(x$magnitude$formula), collapse = " "), "\n", sep = "")
  if(inherits(x$composition$formula, "magCompComposition")){
    comp_text <- paste(names(x$composition$formula), collapse = ", ")
  }else{
    comp_text <- paste(deparse(x$composition$formula), collapse = " ")
  }
  cat("  composition:  ", comp_text, "\n", sep = "")
  cat("  magnitude terms: ",
      length(x$magnitude$process_terms), "\n", sep = "")
  cat("  composition terms:",
      length(x$composition$process_terms), "\n", sep = "")
  if(!is.null(x$backend) && !is.null(x$backend$n_chains))
    cat("  chains:       ", x$backend$n_chains, "\n", sep = "")
  cat("  fitted:       ", x$fitted, "\n", sep = "")
  invisible(x)
}

magcomp_has_intercept <- function(X)
{
  is.matrix(X) &&
    ncol(X) >= 1L &&
    identical(colnames(X)[1L], "(Intercept)") &&
    all(X[,1L] == 1) &&
    all(is.finite(X))
}

magcomp_check_sampler_support <- function(data, family, magnitude,
                                          composition, smooth, priors)
{
  if(!family$name %in% c("fixed_area_stems", "point_basal_area")){
    stop("error: the sampler currently supports fixedAreaStems() and ",
         "pointBasalArea() only")
  }
  if(is.null(data$counts))
    stop("error: count-based families require a category count matrix")
  if(!is.null(composition)) {
    if(length(data$category_names) < 1L || length(data$category_names) > 2L){
      stop("error: the sampler currently supports one or two ",
           "categorical response variables")
    }
    if(length(data$allocation_blocks) < 1L ||
       length(data$allocation_blocks) > 2L){
      stop("error: the sampler currently requires one or two ",
           "allocation blocks")
    }
    if(!identical(data$allocation_blocks[[1L]]$type, "marginal")){
      stop("error: the first allocation block must be marginal")
    }
    if(length(data$allocation_blocks) == 2L &&
       !identical(data$allocation_blocks[[2L]]$type, "conditional")){
      stop("error: the second allocation block must be conditional")
    }
    if(is.null(composition$allocation_components) ||
       length(composition$allocation_components) != length(data$allocation_blocks)){
      stop("error: the sampler requires one composition formula ",
           "per allocation block")
    }
    if(ncol(data$counts) < 2L)
      stop("error: count-based families require at least two categories")
  }
  if(!is.null(magnitude)) {
    if(is.null(data$sampler_area)){
      if(identical(family$name, "fixed_area_stems")){
        stop("error: fixedAreaStems() requires fixedArea(area)")
      }else{
        stop("error: pointBasalArea() requires pointSampling(baf)")
      }
    }
    if(anyNA(data$sampler_area) || any(!is.finite(data$sampler_area)) ||
       any(data$sampler_area <= 0)){
      stop("error: sampler exposure values must be positive and finite")
    }
    if(!magcomp_has_intercept(magnitude$X))
      stop("error: the sampler currently requires a magnitude ",
           "formula with an intercept")
  }
  for(block_name in names(composition$allocation_components)){
    if(!magcomp_has_intercept(composition$allocation_components[[block_name]]$X)){
      stop("error: the sampler currently requires composition ",
           "formulas with an intercept")
    }
  }
  for(term in magnitude$process_terms){
    if(!term$term_type %in% c("iid", "gp", "nngp")){
      stop("error: the sampler currently supports iid(), gp(), ",
           "and nngp() ",
           "magnitude process terms only")
    }
    if(identical(term$term_type, "iid")){
      if(is.null(term$group) || is.null(term$n_group))
        stop("error: iid() grouping information was not prepared")
      if(!is.null(term$coef_name) && is.null(term$covariate))
        stop("error: iid() covariate information was not prepared")
    }else if(identical(term$term_type, "gp")){
      if(is.null(term$coords))
        stop("error: magnitude GP coordinate information was not prepared")
      if(!is.null(term$coef_name) && is.null(term$covariate))
        stop("error: magnitude GP covariate information was not prepared")
    }else{
      if(is.null(term$coords) || is.null(term$order) ||
         is.null(term$neighbor_count) || is.null(term$neighbor_index)){
        stop("error: magnitude NNGP spatial information was not prepared")
      }
      if(!is.null(term$coef_name) && is.null(term$covariate)){
        stop("error: magnitude NNGP spatial covariate information was not ",
             "prepared")
      }
    }
  }
  for(block_name in names(composition$allocation_components)){
    terms <- composition$allocation_components[[block_name]]$process_terms
    for(term in terms){
      if(!term$term_type %in% c("factor_iid", "factor_gp",
                                "factor_nngp")){
        stop("error: the sampler currently supports factor_iid(), ",
             "factor_gp(), and factor_nngp() composition terms only")
      }
      if(term$K > data$allocation_blocks[[block_name]]$coef_stick_count)
        stop("error: composition factor K cannot exceed the number of ",
             "coefficient-stick rows")
      if(is.null(term$anchor_index) || length(term$anchor_index) != term$K)
        stop("error: composition factor anchors were not finalized")
      if(identical(term$term_type, "factor_iid")){
        if(is.null(term$group) || is.null(term$n_group))
          stop("error: factor_iid() grouping information was not prepared")
        if(!is.null(term$coef_name) && is.null(term$covariate))
          stop("error: factor_iid() covariate information was not prepared")
      }else{
        if(identical(term$term_type, "factor_gp")){
          if(is.null(term$coords))
            stop("error: factor_gp() coordinate information was not prepared")
          if(!is.null(term$coef_name) && is.null(term$covariate))
            stop("error: factor_gp() covariate information was not prepared")
        }else{
          if(is.null(term$coords) || is.null(term$order) ||
             is.null(term$neighbor_count) || is.null(term$neighbor_index)){
            stop("error: factor_nngp() spatial information was not prepared")
          }
          if(!is.null(term$coef_name) && is.null(term$covariate)){
            stop("error: factor_nngp() covariate information was not ",
                 "prepared")
          }
        }
      }
    }
  }
  if(!identical(smooth$loading, "none")){
    stop("error: the sampler currently requires ",
         "magCompSmooth(loading = 'none')")
  }
  for(block_name in names(data$allocation_blocks)){
    block_smooth <- data$allocation_blocks[[block_name]]$smooth
    if(!is.null(block_smooth) && !identical(block_smooth$loading, "none")){
      stop("error: the sampler currently requires ",
           "loading = 'none' for allocation block '", block_name, "'")
    }
  }

  invisible(TRUE)
}

magcomp_chain_names <- function(n_chains)
{
  paste0("chain", seq_len(n_chains))
}

magcomp_bind_chain_samples <- function(chains, name)
{
  first <- chains[[1L]][[name]]
  d <- dim(first)
  if(is.null(d))
    stop("error: backend object '", name, "' is not an array")

  n_chains <- length(chains)
  out <- array(NA_real_, dim = c(d, n_chains))

  for(ch in seq_len(n_chains)){
    x <- chains[[ch]][[name]]
    if(!identical(dim(x), d))
      stop("error: backend object '", name, "' has inconsistent dimensions")

    if(length(d) == 2L){
      out[,,ch] <- x
    }else if(length(d) == 3L){
      out[,,,ch] <- x
    }else{
      stop("error: backend object '", name, "' has unsupported dimensions")
    }
  }

  out
}

magcomp_average_chain_matrix <- function(chains, name)
{
  first <- chains[[1L]][[name]]
  out <- first

  if(length(chains) == 1L)
    return(out)

  for(ch in seq_len(length(chains))[-1L]){
    x <- chains[[ch]][[name]]
    if(!identical(dim(x), dim(first)))
      stop("error: backend object '", name, "' has inconsistent dimensions")
    out <- out + x
  }

  out / length(chains)
}

magcomp_bind_chain_sample_lists <- function(chains, name)
{
  first <- chains[[1L]][[name]]
  if(is.null(first))
    return(NULL)
  if(!is.list(first) || is.null(names(first)))
    stop("error: backend object '", name, "' must be a named list")

  out <- vector("list", length(first))
  names(out) <- names(first)
  n_chains <- length(chains)
  for(b in seq_along(first)){
    d <- dim(first[[b]])
    if(is.null(d))
      stop("error: backend object '", name, "' contains a non-array entry")

    out[[b]] <- array(NA_real_, dim = c(d, n_chains))
    for(ch in seq_len(n_chains)){
      x <- chains[[ch]][[name]][[b]]
      if(!identical(dim(x), d))
        stop("error: backend object '", name,
             "' has inconsistent dimensions")

      if(length(d) == 2L){
        out[[b]][,,ch] <- x
      }else if(length(d) == 3L){
        out[[b]][,,,ch] <- x
      }else{
        stop("error: backend object '", name,
             "' has unsupported dimensions")
      }
    }
  }

  out
}

magcomp_bind_chain_named_numeric <- function(chains, name, chain_ids)
{
  first <- chains[[1L]][[name]]
  if(is.null(first))
    return(NULL)
  item_names <- names(first)
  out <- matrix(NA_real_, nrow = length(first), ncol = length(chains),
                dimnames = list(item_names, chain_ids))
  for(ch in seq_along(chains)){
    x <- chains[[ch]][[name]]
    if(length(x) != length(first))
      stop("error: backend object '", name, "' has inconsistent length")
    if(!is.null(item_names) && !identical(names(x), item_names))
      stop("error: backend object '", name, "' has inconsistent names")
    out[,ch] <- as.numeric(x)
  }
  out
}

magcomp_bind_chain_named_history <- function(chains, name, chain_ids)
{
  first <- chains[[1L]][[name]]
  if(is.null(first))
    return(NULL)
  if(!is.list(first) || is.null(names(first)))
    stop("error: backend object '", name, "' must be a named list")

  out <- vector("list", length(first))
  names(out) <- names(first)
  for(block_name in names(first)){
    out[[block_name]] <- lapply(chains, function(x){
      value <- x[[name]][[block_name]]
      if(is.null(value))
        stop("error: backend object '", name, "' has inconsistent names")
      value
    })
    names(out[[block_name]]) <- chain_ids
  }
  out
}

magcomp_check_backend_chain_errors <- function(chains)
{
  for(ch in seq_along(chains)){
    if(is.null(chains[[ch]]))
      stop("error: magComp chain ", ch,
           " returned no result; its worker may have terminated", call. = FALSE)
    if(inherits(chains[[ch]], "try-error")){
      condition <- attr(chains[[ch]], "condition")
      message <- if(!is.null(condition)) conditionMessage(condition) else
        as.character(chains[[ch]])[1L]
      stop("error: magComp chain ", ch, " failed: ", message, call. = FALSE)
    }
  }

  invisible(TRUE)
}

magcomp_combine_backend_chains <- function(chains, data, magnitude, composition)
{
  magcomp_check_backend_chain_errors(chains)

  n_chains <- length(chains)
  chain_ids <- magcomp_chain_names(n_chains)

  if(n_chains == 1L){
    out <- chains[[1L]]
    out$n_chains <- 1L
    out$chain_ids <- chain_ids
    return(magcomp_label_backend(out, data, magnitude, composition))
  }

  sample_names <- c("magnitude_coefficient_samples",
                    "pi_samples",
                    "proportion_samples",
                    "magnitude_samples")
  mean_names <- c("proportion_fitted_mean",
                  "density_fitted_mean",
                  "count_fitted_mean",
                  "total_density_fitted_mean",
                  "total_count_fitted_mean")

  out <- list()
  for(name in sample_names)
    if(!is.null(chains[[1L]][[name]]))
      out[[name]] <- magcomp_bind_chain_samples(chains, name)
  if(!is.null(chains[[1L]]$composition_block_coefficient_samples)){
    out$composition_block_coefficient_samples <-
      magcomp_bind_chain_sample_lists(chains,
                                      "composition_block_coefficient_samples")
  }
  if(!is.null(chains[[1L]]$composition_block_loading_samples)){
    out$composition_block_loading_samples <-
      magcomp_bind_chain_sample_lists(chains,
                                      "composition_block_loading_samples")
  }
  if(!is.null(chains[[1L]]$composition_block_score_samples)){
    out$composition_block_score_samples <-
      magcomp_bind_chain_sample_lists(chains,
                                      "composition_block_score_samples")
  }
  if(!is.null(chains[[1L]]$composition_block_covariance_samples)){
    out$composition_block_covariance_samples <-
      magcomp_bind_chain_sample_lists(chains,
                                      "composition_block_covariance_samples")
  }
  if(!is.null(chains[[1L]]$magnitude_process_effect_samples)){
    out$magnitude_process_effect_samples <-
      magcomp_bind_chain_sample_lists(chains, "magnitude_process_effect_samples")
  }
  if(!is.null(chains[[1L]]$magnitude_process_variance_samples)){
    out$magnitude_process_variance_samples <-
      magcomp_bind_chain_sample_lists(
        chains,
        "magnitude_process_variance_samples"
      )
  }
  if(!is.null(chains[[1L]]$magnitude_process_covariance_samples)){
    out$magnitude_process_covariance_samples <-
      magcomp_bind_chain_sample_lists(chains,
                                      "magnitude_process_covariance_samples")
  }
  if(!is.null(chains[[1L]]$nb_size_samples)){
    out$nb_size_samples <- magcomp_bind_chain_samples(chains,
                                                      "nb_size_samples")
    out$nb_size_acceptance <- matrix(
      vapply(chains, function(x) as.numeric(x$nb_size_acceptance),
             numeric(1L)),
      nrow = 1L,
      dimnames = list("size", chain_ids)
    )
    out$nb_size_tuning <- matrix(
      vapply(chains, function(x) as.numeric(x$nb_size_tuning),
             numeric(1L)),
      nrow = 1L,
      dimnames = list("size", chain_ids)
    )
    out$nb_size_batch_acceptance <-
      stats::setNames(lapply(chains, `[[`, "nb_size_batch_acceptance"),
                      chain_ids)
    out$nb_size_tuning_history <-
      stats::setNames(lapply(chains, `[[`, "nb_size_tuning_history"),
                      chain_ids)
  }
  if(!is.null(chains[[1L]]$composition_block_covariance_acceptance)){
    out$composition_block_covariance_acceptance <-
      magcomp_bind_chain_named_numeric(
        chains,
        "composition_block_covariance_acceptance",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$composition_block_covariance_tuning)){
    out$composition_block_covariance_tuning <-
      magcomp_bind_chain_named_numeric(
        chains,
        "composition_block_covariance_tuning",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$composition_block_covariance_batch_acceptance)){
    out$composition_block_covariance_batch_acceptance <-
      magcomp_bind_chain_named_history(
        chains,
        "composition_block_covariance_batch_acceptance",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$composition_block_covariance_tuning_history)){
    out$composition_block_covariance_tuning_history <-
      magcomp_bind_chain_named_history(
        chains,
        "composition_block_covariance_tuning_history",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$magnitude_process_covariance_acceptance)){
    out$magnitude_process_covariance_acceptance <-
      magcomp_bind_chain_named_numeric(
        chains,
        "magnitude_process_covariance_acceptance",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$magnitude_process_covariance_tuning)){
    out$magnitude_process_covariance_tuning <-
      magcomp_bind_chain_named_numeric(
        chains,
        "magnitude_process_covariance_tuning",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$magnitude_process_covariance_batch_acceptance)){
    out$magnitude_process_covariance_batch_acceptance <-
      magcomp_bind_chain_named_history(
        chains,
        "magnitude_process_covariance_batch_acceptance",
        chain_ids
      )
  }
  if(!is.null(chains[[1L]]$magnitude_process_covariance_tuning_history)){
    out$magnitude_process_covariance_tuning_history <-
      magcomp_bind_chain_named_history(
        chains,
        "magnitude_process_covariance_tuning_history",
        chain_ids
      )
  }
  for(name in mean_names)
    if(!is.null(chains[[1L]][[name]]))
      out[[name]] <- magcomp_average_chain_matrix(chains, name)

  if(!is.null(chains[[1L]]$profile_stage_seconds)){
    out$profile_stage_seconds <- Reduce(
      "+",
      lapply(chains, `[[`, "profile_stage_seconds")
    )
  }
  if(!is.null(chains[[1L]]$profile_pg)){
    out$profile_pg <- Reduce("+", lapply(chains, `[[`, "profile_pg"))
    pg_total <- out$profile_pg[["total_draws"]]
    out$profile_pg[["seconds_per_draw"]] <- if(pg_total > 0)
      out$profile_pg[["total_seconds"]] / pg_total else NA_real_
    out$profile_pg[["fallback_fraction"]] <- if(pg_total > 0)
      out$profile_pg[["fallback_draws"]] / pg_total else NA_real_
    probit_total <- out$profile_pg[["probit_composition_draws"]]
    out$profile_pg[["probit_composition_seconds_per_draw"]] <-
      if(probit_total > 0)
        out$profile_pg[["probit_composition_seconds"]] / probit_total else
          NA_real_
  }

  out$nb_size <- chains[[1L]]$nb_size
  out$nngp_graph_count <- chains[[1L]]$nngp_graph_count
  out$composition_link <- chains[[1L]]$composition_link
  out$n_chains <- n_chains
  out$chain_ids <- chain_ids

  magcomp_label_backend(out, data, magnitude, composition)
}

magcomp_run_backend_chains <- function(config, data, magnitude, composition)
{
  control <- config$control
  n_chains <- control$n_chains
  n_cores <- min(control$n_cores %||% 1L, n_chains)
  chain_ids <- magcomp_chain_names(n_chains)

  run_chain <- function(ch)
  {
    if(!is.null(control$chain_seeds))
      set.seed(control$chain_seeds[ch])
    if(isTRUE(control$verbose) && n_chains > 1L)
      message("magComp chain ", ch, " of ", n_chains)

    config_ch <- config
    config_ch$control$n_chains <- 1L
    config_ch$control$chain_seeds <- NULL

    backend_ch <- .Call("magcomp_sampler", config_ch, PACKAGE = "magComp")
    backend_ch$n_chains <- 1L
    backend_ch$chain_ids <- chain_ids[ch]
    magcomp_label_backend(backend_ch, data, magnitude, composition)
  }

  if(n_chains == 1L || n_cores == 1L){
    chains <- vector("list", n_chains)
    for(ch in seq_len(n_chains))
      chains[[ch]] <- run_chain(ch)
  }else if(identical(.Platform$OS.type, "unix")){
    chains <- parallel::mclapply(seq_len(n_chains), run_chain,
                                 mc.cores = n_cores,
                                 mc.preschedule = FALSE,
                                 mc.set.seed = is.null(control$chain_seeds))
  }else{
    warning("parallel chains require forked processes on this platform; ",
            "running chains serially")
    chains <- vector("list", n_chains)
    for(ch in seq_len(n_chains))
      chains[[ch]] <- run_chain(ch)
  }

  magcomp_combine_backend_chains(chains, data, magnitude, composition)
}

magcomp_composition_factor_output_info <- function(composition, key)
{
  block_names <- names(composition$allocation_components)
  if(is.null(block_names))
    return(NULL)

  for(block_name in block_names){
    component <- composition$allocation_components[[block_name]]
    terms <- component$process_terms
    if(identical(key, block_name) && length(terms) == 1L){
      return(list(block_name = block_name, term = terms[[1L]]))
    }
    prefix <- paste0(block_name, ":")
    if(startsWith(key, prefix)){
      term_name <- substring(key, nchar(prefix) + 1L)
      for(term in terms){
        if(identical(term$name, term_name) ||
           identical(term$term_type, term_name)){
          return(list(block_name = block_name, term = term))
        }
      }
    }
  }

  NULL
}

magcomp_label_backend <- function(backend, data, magnitude, composition)
{
  units <- data$unit_levels
  categories <- if(is.null(composition)) "total" else colnames(data$counts)
  pi_stick_categories <- data$stick_names %||% categories[-length(categories)]
  magnitude_names <- colnames(magnitude$X)
  chain_ids <- backend$chain_ids

  if(!is.null(backend$magnitude_coefficient_samples)) {
    if(length(dim(backend$magnitude_coefficient_samples)) == 2L){
      dimnames(backend$magnitude_coefficient_samples) <- list(magnitude_names, NULL)
    }else{
      dimnames(backend$magnitude_coefficient_samples) <- list(magnitude_names, NULL,
                                              chain_ids)
    }
  }
  if(!is.null(backend$pi_samples)) {
    if(length(dim(backend$pi_samples)) == 3L){
      dimnames(backend$pi_samples) <- list(units, pi_stick_categories, NULL)
    }else{
      dimnames(backend$pi_samples) <- list(units, pi_stick_categories, NULL,
                                           chain_ids)
    }
  }
  if(!is.null(backend$proportion_samples)) {
    if(length(dim(backend$proportion_samples)) == 3L){
      dimnames(backend$proportion_samples) <- list(units, categories, NULL)
    }else{
      dimnames(backend$proportion_samples) <- list(units, categories, NULL,
                                                   chain_ids)
    }
  }
  if(!is.null(backend$magnitude_samples)) {
    if(length(dim(backend$magnitude_samples)) == 2L){
      dimnames(backend$magnitude_samples) <- list(units, NULL)
    }else{
      dimnames(backend$magnitude_samples) <- list(units, NULL, chain_ids)
    }
  }
  if(!is.null(backend$nb_size_samples)){
    if(length(dim(backend$nb_size_samples)) == 2L){
      dimnames(backend$nb_size_samples) <- list("size", NULL)
    }else{
      dimnames(backend$nb_size_samples) <- list("size", NULL, chain_ids)
    }
  }

  if(!is.null(backend$composition_block_coefficient_samples)){
    for(block_name in names(backend$composition_block_coefficient_samples)){
      block <- data$allocation_blocks[[block_name]]
      block_component <- composition$allocation_components[[block_name]]
      block_categories <- colnames(block$counts)
      block_sticks <- block$stick_names %||%
        block_categories[-length(block_categories)]
      block_names <- colnames(block_component$X)
      value <- backend$composition_block_coefficient_samples[[block_name]]
      if(length(dim(value)) == 3L){
        dimnames(value) <- list(block_sticks, block_names, NULL)
      }else{
        dimnames(value) <- list(block_sticks, block_names, NULL,
                                chain_ids)
      }
      backend$composition_block_coefficient_samples[[block_name]] <- value
    }
  }

  if(!is.null(backend$composition_block_loading_samples)){
    for(sample_name in names(backend$composition_block_loading_samples)){
      info <- magcomp_composition_factor_output_info(composition, sample_name)
      if(is.null(info))
        next
      block <- data$allocation_blocks[[info$block_name]]
      term <- info$term
      factor_names <- paste0("factor", seq_len(dim(
        backend$composition_block_loading_samples[[sample_name]])[2L]))
      if(!is.null(term$factor_names))
        factor_names <- term$factor_names
      value <- backend$composition_block_loading_samples[[sample_name]]
      if(length(dim(value)) == 3L){
        dimnames(value) <- list(block$stick_names, factor_names, NULL)
      }else{
        dimnames(value) <- list(block$stick_names, factor_names, NULL,
                                chain_ids)
      }
      backend$composition_block_loading_samples[[sample_name]] <- value
    }
  }

  if(!is.null(backend$composition_block_score_samples)){
    for(sample_name in names(backend$composition_block_score_samples)){
      info <- magcomp_composition_factor_output_info(composition, sample_name)
      if(is.null(info))
        next
      term <- info$term
      value <- backend$composition_block_score_samples[[sample_name]]
      factor_names <- paste0("factor", seq_len(dim(value)[2L]))
      group_names <- seq_len(dim(value)[1L])
      if(!is.null(term$factor_names))
        factor_names <- term$factor_names
      if(!is.null(term$group_levels))
        group_names <- term$group_levels
      if(!is.null(term$location_levels))
        group_names <- term$location_levels
      if(length(dim(value)) == 3L){
        dimnames(value) <- list(group_names, factor_names, NULL)
      }else{
        dimnames(value) <- list(group_names, factor_names, NULL,
                                chain_ids)
      }
      backend$composition_block_score_samples[[sample_name]] <- value
    }
  }

  if(!is.null(backend$composition_block_covariance_samples)){
    covariance_old_names <- character(0L)
    covariance_new_names <- character(0L)
    for(sample_name in names(backend$composition_block_covariance_samples)){
      value <- backend$composition_block_covariance_samples[[sample_name]]
      term_info <- magcomp_composition_factor_output_info(composition,
                                                         sample_name)
      factor_names <- paste0("factor", seq_len(dim(value)[1L]))
      if(!is.null(term_info) && !is.null(term_info$term$factor_names))
        factor_names <- term_info$term$factor_names
      old_names <- if(length(factor_names) == 1L) sample_name else
        paste0(sample_name, ":factor", seq_along(factor_names))
      new_names <- if(length(factor_names) == 1L) sample_name else
        paste(sample_name, factor_names, sep = ":")
      covariance_old_names <- c(covariance_old_names, old_names)
      covariance_new_names <- c(covariance_new_names, new_names)
      labels <- as.vector(outer(factor_names, term_info$term$covariance_parameters, paste, sep = ":"))
      dimnames(value) <- if(length(dim(value)) == 2L) list(labels, NULL) else list(labels, NULL, chain_ids)
      backend$composition_block_covariance_samples[[sample_name]] <- value
    }
    covariance_name_map <- stats::setNames(covariance_new_names, covariance_old_names)
    for(name in c("composition_block_covariance_acceptance",
                  "composition_block_covariance_tuning")){
      if(!is.null(backend[[name]])){
        row_names <- if(is.matrix(backend[[name]])) rownames(backend[[name]]) else names(backend[[name]])
        match_index <- match(row_names, names(covariance_name_map))
        row_names[!is.na(match_index)] <- covariance_name_map[match_index[!is.na(match_index)]]
        if(is.matrix(backend[[name]])) rownames(backend[[name]]) <- row_names else
          names(backend[[name]]) <- row_names
      }
    }
    for(name in c("composition_block_covariance_batch_acceptance",
                  "composition_block_covariance_tuning_history")){
      if(!is.null(backend[[name]]) && !is.null(names(backend[[name]]))){
        item_names <- names(backend[[name]])
        match_index <- match(item_names, names(covariance_name_map))
        item_names[!is.na(match_index)] <- covariance_name_map[match_index[!is.na(match_index)]]
        names(backend[[name]]) <- item_names
      }
    }
  }

  if(!is.null(backend$magnitude_process_effect_samples)){
    term_lookup <- stats::setNames(magnitude$process_terms,
                                   vapply(magnitude$process_terms, `[[`,
                                          character(1L), "name"))
    for(term_name in names(backend$magnitude_process_effect_samples)){
      value <- backend$magnitude_process_effect_samples[[term_name]]
      term <- term_lookup[[term_name]]
      row_names <- seq_len(dim(value)[1L])
      if(!is.null(term)){
        if(identical(term$term_type, "iid") &&
           !is.null(term$group_levels)){
          row_names <- term$group_levels
        }else if(term$term_type %in% c("gp", "nngp")){
          row_names <- term$location_levels %||% units
        }
      }
      if(length(dim(value)) == 2L){
        dimnames(value) <- list(row_names, NULL)
      }else{
        dimnames(value) <- list(row_names, NULL, chain_ids)
      }
      backend$magnitude_process_effect_samples[[term_name]] <- value
    }
  }

  if(!is.null(backend$magnitude_process_variance_samples)){
    for(term_name in names(backend$magnitude_process_variance_samples)){
      value <- backend$magnitude_process_variance_samples[[term_name]]
      if(length(dim(value)) == 2L){
        dimnames(value) <- list("variance", NULL)
      }else{
        dimnames(value) <- list("variance", NULL, chain_ids)
      }
      backend$magnitude_process_variance_samples[[term_name]] <- value
    }
  }

  if(!is.null(backend$magnitude_process_covariance_samples)){
    for(term_name in names(backend$magnitude_process_covariance_samples)){
      value <- backend$magnitude_process_covariance_samples[[term_name]]
      labels <- term_lookup[[term_name]]$covariance_parameters
      dimnames(value) <- if(length(dim(value)) == 2L) list(labels, NULL) else list(labels, NULL, chain_ids)
      backend$magnitude_process_covariance_samples[[term_name]] <- value
    }
  }

  if(!is.null(backend$proportion_fitted_mean))
    dimnames(backend$proportion_fitted_mean) <- list(units, categories)
  if(!is.null(backend$density_fitted_mean))
    dimnames(backend$density_fitted_mean) <- list(units, categories)
  if(!is.null(backend$count_fitted_mean))
    dimnames(backend$count_fitted_mean) <- list(units, categories)
  if(!is.null(backend$total_density_fitted_mean))
    dimnames(backend$total_density_fitted_mean) <- list(units, "total")
  if(!is.null(backend$total_count_fitted_mean))
    dimnames(backend$total_count_fitted_mean) <- list(units, "total")

  backend
}

magcomp_check_prediction_args <- function(dots, method)
{
  if(length(dots) == 0L)
    return(invisible(NULL))
  dot_names <- names(dots)
  if(is.null(dot_names))
    dot_names <- rep("", length(dots))
  unnamed <- !nzchar(dot_names)
  dot_names[unnamed] <- paste0("..", which(unnamed))
  stop("error: unused ", method, " argument(s): ",
       paste(dot_names, collapse = ", "))
}

fitted.magComp <- function(object,
                           type = NULL,
                           target = c("expected", "observation"),
                           summary = c("mean", "samples", "quantile"),
                           probs = c(0.025, 0.5, 0.975), ...)
{
  type <- magcomp_prediction_types(type, object)
  target <- match.arg(target)
  summary <- match.arg(summary)
  probs <- magcomp_prediction_probs(probs)
  magcomp_check_prediction_args(list(...), "fitted")
  magcomp_check_prediction_fit(object)
  exposure <- object$data$sampler_area
  if((!is.null(object$magnitude) && (target == "observation" || "count" %in% type) ||
      "density" %in% type) && is.null(exposure))
    stop("error: density predictions require plot area or BAF information")
  total <- if(is.null(object$magnitude) &&
    (target == "observation" || any(c("count", "density") %in% type))) rowSums(object$data$counts) else NULL
  magcomp_prediction_result(object, object$backend$magnitude_samples,
    object$backend$proportion_samples, exposure,
    target, type, summary, probs, "fitted", total)
}

magcomp_predict_unit_levels <- function(object, newdata)
{
  unit_name <- object$data$unit_name
  if(!is.null(unit_name) && unit_name %in% names(newdata)){
    out <- as.character(newdata[[unit_name]])
    if(anyNA(out) || any(!nzchar(out)))
      stop("error: prediction unit column contains missing values")
    return(make.unique(out))
  }

  rn <- rownames(newdata)
  if(is.null(rn) || anyNA(rn) || any(!nzchar(rn)) ||
     identical(rn, as.character(seq_len(nrow(newdata)))))
    return(paste0("new", seq_len(nrow(newdata))))
  make.unique(rn)
}

magcomp_predict_area <- function(object, newdata, area)
{
  n <- nrow(newdata)

  if(is.null(area)){
    if(identical(object$family$name, "point_basal_area")){
      baf_name <- object$data$design$baf_name
      if(!is.null(baf_name) && baf_name %in% names(newdata)){
        baf <- newdata[[baf_name]]
        if(!is.numeric(baf))
          stop("error: prediction BAF must be numeric")
        if(anyNA(baf) || any(!is.finite(baf)) || any(baf <= 0))
          stop("error: prediction BAF must be positive and finite")
        area <- 1 / baf
      }else{
        stop("error: prediction exposure is required: supply area = 1 / BAF ",
             "or the fitted BAF column in newdata")
      }
    }else{
      area_name <- object$data$design$area_name
      if(!is.null(area_name) && area_name %in% names(newdata)){
        area <- newdata[[area_name]]
      }else{
        stop("error: prediction area is required: supply area or the fitted ",
             "area column in newdata")
      }
    }
  }else if(is.character(area) && length(area) == 1L){
    if(!area %in% names(newdata))
      stop("error: prediction area column not found: ", area)
    area <- newdata[[area]]
  }

  if(!is.numeric(area))
    stop("error: prediction exposure must be numeric")
  if(length(area) == 1L)
    area <- rep(area, n)
  if(length(area) != n)
    stop("error: prediction area must be scalar or have one value per row")
  if(identical(object$family$name, "point_basal_area")){
    if(anyNA(area) || any(!is.finite(area)) || any(area <= 0))
      stop("error: prediction exposure must be positive and finite")
  }else{
    magcomp_positive_finite(area, "prediction area")
  }
  area
}

magcomp_prediction_family_supported <- function(object)
{
  object$family$name %in% c("fixed_area_stems", "point_basal_area")
}

magcomp_prediction_model_matrix <- function(formula, training_data, newdata,
                                            reference_X, component)
{
  training_mf <- stats::model.frame(formula, data = training_data,
                                    na.action = stats::na.fail)
  tt <- attr(training_mf, "terms")
  xlev <- stats::.getXlevels(tt, training_mf)
  mf <- stats::model.frame(formula, data = newdata,
                           na.action = stats::na.fail,
                           xlev = xlev)
  X <- stats::model.matrix(tt, data = mf,
                           contrasts.arg = attr(reference_X, "contrasts"))

  reference_names <- colnames(reference_X)
  missing_names <- setdiff(reference_names, colnames(X))
  extra_names <- setdiff(colnames(X), reference_names)
  if(length(extra_names) > 0L){
    stop("error: prediction ", component,
         " model matrix has unknown column(s): ",
         paste(extra_names, collapse = ", "))
  }

  out <- matrix(0, nrow = nrow(X), ncol = ncol(reference_X),
                dimnames = list(rownames(X), reference_names))
  out[,colnames(X)] <- X
  if(length(missing_names) > 0L)
    out[,missing_names] <- 0
  out
}

magcomp_prediction_covariate <- function(term, newdata, label)
{
  covariate_name <- term$covariate_name %||% term$coef_name
  if(is.null(covariate_name))
    return(rep(1.0, nrow(newdata)))
  if(!covariate_name %in% names(newdata))
    stop("error: prediction ", label, " covariate column not found: ",
         covariate_name)
  value <- as.numeric(newdata[[covariate_name]])
  if(length(value) != nrow(newdata) || anyNA(value) || any(!is.finite(value)))
    stop("error: prediction ", label, " covariate must be finite")
  value
}

magcomp_prediction_coords <- function(term, newdata, label)
{
  coord_names <- term$coord_names
  if(is.null(coord_names) || length(coord_names) < 1L)
    stop("error: ", label, " coordinate names are unavailable")
  missing_names <- setdiff(coord_names, names(newdata))
  if(length(missing_names) > 0L){
    stop("error: prediction ", label, " coordinate column(s) not found: ",
         paste(missing_names, collapse = ", "))
  }
  coords <- as.matrix(newdata[,coord_names,drop = FALSE])
  storage.mode(coords) <- "double"
  if(anyNA(coords) || any(!is.finite(coords)))
    stop("error: prediction ", label, " coordinates must be finite")
  coords
}

magcomp_training_process_coords <- function(object, term)
{
  coords <- term$coords
  if(is.null(coords))
    return(NULL)

  if(!is.null(term$spatial_group) &&
     length(term$spatial_group) == nrow(coords) &&
     !is.null(term$n_group) &&
     term$n_group == object$data$n_unit){
    first <- match(seq_len(object$data$n_unit), term$spatial_group)
    if(anyNA(first))
      stop("error: spatial factor coordinates are missing a sampling unit")
    coords <- coords[first,,drop = FALSE]
    rownames(coords) <- object$data$unit_levels
  }

  coords
}

magcomp_exp_distance <- function(a, b)
{
  a <- as.matrix(a)
  b <- as.matrix(b)
  out <- matrix(0, nrow = nrow(a), ncol = nrow(b))
  for(k in seq_len(ncol(a)))
    out <- out + outer(a[,k], b[,k], "-")^2
  sqrt(out)
}

magcomp_exp_corr <- function(a, b, decay)
{
  exp(-decay * magcomp_exp_distance(a, b))
}

magcomp_gp_predict_one <- function(w_obs, coords_obs, coords_new, decay,
                                   variance)
{
  n_new <- nrow(coords_new)
  R_obs <- magcomp_exp_corr(coords_obs, coords_obs, decay)
  U_obs <- chol(R_obs)
  R_new_obs <- magcomp_exp_corr(coords_new, coords_obs, decay)
  weights <- backsolve(U_obs, forwardsolve(t(U_obs), w_obs))
  mean <- drop(R_new_obs %*% weights)
  V <- forwardsolve(t(U_obs), t(R_new_obs))
  R_new <- magcomp_exp_corr(coords_new, coords_new, decay)
  cond <- R_new - crossprod(V)
  cond <- 0.5 * (cond + t(cond))

  U_new <- chol(variance * cond)
  mean + drop(t(U_new) %*% stats::rnorm(n_new))
}

magcomp_nngp_predict_one <- function(w_obs, coords_obs, coords_new, decay,
                                     variance, m)
{
  n_obs <- nrow(coords_obs)
  n_new <- nrow(coords_new)
  m <- min(as.integer(m), n_obs + n_new - 1L)
  order_new <- do.call(order, as.data.frame(coords_new))
  coords_all <- rbind(coords_obs, coords_new)
  values <- c(as.numeric(w_obs), rep(NA_real_, n_new))

  for(pos in seq_along(order_new)){
    new_id <- order_new[pos]
    global_id <- n_obs + new_id
    previous_new <- if(pos > 1L) n_obs + order_new[seq_len(pos - 1L)] else integer(0L)
    candidates <- c(seq_len(n_obs), previous_new)
    if(length(candidates) == 0L){
      mean <- 0
      var <- 1
    }else{
      d <- sqrt(rowSums(sweep(coords_all[candidates,,drop = FALSE], 2L,
                              coords_all[global_id,], "-")^2))
      neighbor <- candidates[order(d)[seq_len(min(m, length(candidates)))]]
      R_nn <- magcomp_exp_corr(coords_all[neighbor,,drop = FALSE],
                               coords_all[neighbor,,drop = FALSE], decay)
      r_in <- drop(magcomp_exp_corr(coords_all[neighbor,,drop = FALSE],
                                    coords_all[global_id,,drop = FALSE], decay))
      U <- chol(R_nn)
      b <- backsolve(U, forwardsolve(t(U), r_in))
      mean <- sum(b * values[neighbor])
      var <- 1 - sum(b * r_in)
      if(!is.finite(var) || var <= 0)
        stop("error: NNGP prediction conditional variance is non-positive or non-finite")
    }
    values[global_id] <- mean + sqrt(variance * var) * stats::rnorm(1L)
  }

  values[n_obs + seq_len(n_new)]
}

magcomp_predict_spatial_scalar_samples <- function(term, effect_samples,
                                                   covariance_samples,
                                                   variance_samples,
                                                   newdata, unit_levels,
                                                   object,
                                                   label)
{
  coords_obs <- magcomp_training_process_coords(object, term)
  coords_new <- magcomp_prediction_coords(term, newdata, label)
  method <- if(term$term_type %in% c("gp", "factor_gp")) "gp" else "nngp"
  out <- magcomp_spatial_prediction_samples(effect_samples, coords_obs,
    coords_new, covariance_samples, variance_samples %||% 1.0,
    term$m %||% 0L, method, FALSE, term$cov_model, term$st_scale)
  if(length(dim(out)) == 2L)
    dimnames(out) <- list(unit_levels, NULL)
  else
    dimnames(out) <- list(unit_levels, NULL, dimnames(effect_samples)[[3L]])
  out
}

magcomp_predict_iid_scalar_samples <- function(effect_samples, term, newdata,
                                               unit_levels, variance, label)
{
  if(is.null(term$group_name) || !term$group_name %in% names(newdata))
    stop("error: prediction ", label, " group column not found")
  group <- as.character(newdata[[term$group_name]])
  if(anyNA(group) || any(!nzchar(group)))
    stop("error: prediction ", label, " group contains missing values")

  d <- dim(effect_samples)
  n_samples <- d[2L]
  n_chains <- if(length(d) == 3L) d[3L] else 1L
  n_new <- length(group)
  observed_index <- match(group, term$group_levels)
  new_levels <- setdiff(unique(group[is.na(observed_index)]), term$group_levels)

  new_index <- match(group, new_levels)
  out <- .Call("magcomp_predict_iid_scalar_cpp",
               effect_samples,
               as.integer(observed_index),
               as.integer(new_index),
               as.integer(length(new_levels)),
               variance,
               PACKAGE = "magComp")
  if(length(dim(out)) == 2L)
    dimnames(out) <- list(unit_levels, NULL)
  else
    dimnames(out) <- list(unit_levels, NULL, dimnames(effect_samples)[[3L]])
  out
}

magcomp_prediction_block_data <- function(object, block_name, newdata,
                                          unit_levels)
{
  block <- object$data$allocation_blocks[[block_name]]
  if(identical(block$type, "marginal")){
    out <- newdata
    rownames(out) <- unit_levels
    return(out)
  }

  if(!identical(block$type, "conditional"))
    stop("error: unknown allocation block type: ", block$type)

  n <- nrow(newdata)
  n_parent <- length(block$parent_levels)
  out <- newdata[rep(seq_len(n), times = n_parent),,drop = FALSE]
  out[[block$parent]] <- factor(rep(block$parent_levels, each = n),
                                levels = block$parent_levels)
  rownames(out) <- as.vector(vapply(
    block$parent_levels,
    function(parent) paste(unit_levels, parent, sep = ":"),
    character(n)
  ))
  out
}

magcomp_multiply_scalar_effect <- function(samples, effect, covariate)
{
  d <- dim(samples)
  n <- d[1L]
  n_samples <- d[2L]
  n_chains <- if(length(d) == 3L) d[3L] else 1L

  if(length(dim(effect)) == 2L){
    for(s in seq_len(n_samples))
      samples[,s] <- samples[,s] * exp(covariate * effect[,s])
    return(samples)
  }

  for(ch in seq_len(n_chains)){
    for(s in seq_len(n_samples)){
      samples[,s,ch] <- samples[,s,ch] *
        exp(covariate * effect[,s,ch])
    }
  }
  samples
}

magcomp_predict_magnitude_samples <- function(object, newdata, unit_levels,
                                               latent = NULL)
{
  X <- magcomp_prediction_model_matrix(
    object$magnitude$fixed_formula,
    object$data$unit_data,
    newdata,
    object$magnitude$X,
    "magnitude"
  )
  magnitude_coefficients <- object$backend$magnitude_coefficient_samples

  out <- .Call("magcomp_predict_magnitude_cpp",
               X,
               magnitude_coefficients,
               NULL,
               NULL,
               PACKAGE = "magComp")
  for(term in object$magnitude$process_terms){
    covariate <- magcomp_prediction_covariate(term, newdata, "magnitude")
    effect <- NULL
    if(identical(term$term_type, "iid")){
      effect_samples <- object$backend$magnitude_process_effect_samples[[term$name]]
      variance_samples <-
        object$backend$magnitude_process_variance_samples[[term$name]]
      effect <- magcomp_predict_iid_scalar_samples(
        effect_samples,
        term,
        newdata,
        unit_levels,
        variance_samples,
        "magnitude iid"
      )
    }else if(term$term_type %in% c("gp", "nngp")){
      effect_samples <- object$backend$magnitude_process_effect_samples[[term$name]]
      covariance_samples <- object$backend$magnitude_process_covariance_samples[[term$name]]
      variance_samples <-
        object$backend$magnitude_process_variance_samples[[term$name]]
      effect <- magcomp_predict_spatial_scalar_samples(
        term,
        effect_samples,
        covariance_samples,
        variance_samples,
        newdata,
        unit_levels,
        object,
        paste0("magnitude ", term$term_type)
      )
    }else{
      stop("error: unsupported magnitude process term in prediction")
    }
    if(!is.null(latent))
      latent$draws$magnitude$terms[[term$name]]$effects <- effect
    out <- magcomp_multiply_scalar_effect(out, effect, covariate)
  }

  if(length(dim(out)) == 2L)
    dimnames(out) <- list(unit_levels, NULL)
  else
    dimnames(out) <- list(unit_levels, NULL, object$backend$chain_ids)
  out
}

magcomp_block_coef_samples <- function(object, block_name)
{
  block_samples <- object$backend$composition_block_coefficient_samples
  if(!is.null(block_samples) && !is.null(block_samples[[block_name]]))
    return(block_samples[[block_name]])
  stop("error: composition coefficient samples are unavailable for block '",
       block_name, "'")
}

magcomp_block_factor_sample_key <- function(object, block_name, term)
{
  component <- object$composition$allocation_components[[block_name]]
  block_samples <- object$backend$composition_block_score_samples
  candidates <- character(0)
  if(!is.null(term$name))
    candidates <- c(candidates, paste(block_name, term$name, sep = ":"))
  if(!is.null(term$term_type))
    candidates <- c(candidates, paste(block_name, term$term_type, sep = ":"))
  if(length(component$process_terms) == 1L)
    candidates <- c(block_name, candidates)
  candidates <- unique(candidates)
  if(!is.null(block_samples)){
    hit <- candidates[candidates %in% names(block_samples)]
    if(length(hit) > 0L)
      return(hit[1L])
  }
  if(length(component$process_terms) == 1L)
    return(block_name)
  candidates[1L]
}

magcomp_block_loading_samples <- function(object, block_name, term = NULL)
{
  block_samples <- object$backend$composition_block_loading_samples
  if(!is.null(term)){
    key <- magcomp_block_factor_sample_key(object, block_name, term)
    if(!is.null(block_samples) && !is.null(block_samples[[key]]))
      return(block_samples[[key]])
    return(NULL)
  }
  if(!is.null(block_samples) && !is.null(block_samples[[block_name]]))
    return(block_samples[[block_name]])
  NULL
}

magcomp_block_factor_samples <- function(object, block_name, term = NULL)
{
  block_samples <- object$backend$composition_block_score_samples
  if(!is.null(term)){
    key <- magcomp_block_factor_sample_key(object, block_name, term)
    if(!is.null(block_samples) && !is.null(block_samples[[key]]))
      return(block_samples[[key]])
    return(NULL)
  }
  if(!is.null(block_samples) && !is.null(block_samples[[block_name]]))
    return(block_samples[[block_name]])
  NULL
}

magcomp_block_covariance_samples <- function(object, block_name, term = NULL)
{
  block_samples <- object$backend$composition_block_covariance_samples
  if(!is.null(term)){
    key <- magcomp_block_factor_sample_key(object, block_name, term)
    if(!is.null(block_samples) && !is.null(block_samples[[key]]))
      return(block_samples[[key]])
    return(NULL)
  }
  if(!is.null(block_samples) && !is.null(block_samples[[block_name]]))
    return(block_samples[[block_name]])
  NULL
}

magcomp_predict_iid_factor_samples <- function(effect_samples, term, block_data,
                                               row_levels)
{
  if(is.null(term$group_name) || !term$group_name %in% names(block_data))
    stop("error: prediction composition factor group column not found")
  group <- as.character(block_data[[term$group_name]])
  if(anyNA(group) || any(!nzchar(group)))
    stop("error: prediction composition factor group contains missing values")

  d <- dim(effect_samples)
  K <- d[2L]
  n_samples <- d[3L]
  n_chains <- if(length(d) == 4L) d[4L] else 1L
  n_row <- length(group)
  observed_index <- match(group, term$group_levels)
  new_levels <- setdiff(unique(group[is.na(observed_index)]),
                        term$group_levels)

  new_index <- match(group, new_levels)
  out <- .Call("magcomp_predict_iid_factor_cpp",
               effect_samples,
               as.integer(observed_index),
               as.integer(new_index),
               as.integer(length(new_levels)),
               PACKAGE = "magComp")
  if(length(dim(out)) == 3L)
    dimnames(out) <- list(row_levels, dimnames(effect_samples)[[2L]], NULL)
  else
    dimnames(out) <- list(row_levels, dimnames(effect_samples)[[2L]], NULL,
                          dimnames(effect_samples)[[4L]])
  out
}

magcomp_predict_spatial_factor_samples <- function(object, term,
                                                   effect_samples,
                                                   covariance_samples,
                                                   newdata,
                                                   unit_levels,
                                                   label)
{
  coords_obs <- magcomp_training_process_coords(object, term)
  coords_new <- magcomp_prediction_coords(term, newdata, label)
  method <- if(identical(term$term_type, "factor_gp")) "gp" else "nngp"
  out <- magcomp_spatial_prediction_samples(effect_samples, coords_obs,
    coords_new, covariance_samples, 1.0,
    term$m %||% 0L, method, TRUE, term$cov_model, term$st_scale)
  if(length(dim(out)) == 3L)
    dimnames(out) <- list(unit_levels, dimnames(effect_samples)[[2L]], NULL)
  else
    dimnames(out) <- list(unit_levels, dimnames(effect_samples)[[2L]], NULL,
                          dimnames(effect_samples)[[4L]])
  out
}

magcomp_predict_block_factor_samples <- function(object, block_name, block,
                                                 term, block_data,
                                                 newdata, unit_levels)
{
  if(!term$term_type %in% c("factor_iid", "factor_gp", "factor_nngp"))
    return(NULL)

  effect_samples <- magcomp_block_factor_samples(object, block_name, term)
  if(is.null(effect_samples))
    stop("error: composition factor samples are unavailable")

  if(identical(term$term_type, "factor_iid")){
    return(magcomp_predict_iid_factor_samples(
      effect_samples,
      term,
      block_data,
      rownames(block_data)
    ))
  }

  unit_effect <- magcomp_predict_spatial_factor_samples(
    object,
    term,
    effect_samples,
    magcomp_block_covariance_samples(object, block_name, term),
    newdata,
    unit_levels,
    paste0("composition ", term$term_type)
  )
  if(identical(block$type, "conditional")){
    unit_index <- rep(seq_len(nrow(newdata)),
                      times = length(block$parent_levels))
  }else{
    unit_index <- seq_len(nrow(newdata))
  }

  if(length(dim(unit_effect)) == 3L){
    out <- unit_effect[unit_index,,,drop = FALSE]
  }else{
    out <- unit_effect[unit_index,,,,drop = FALSE]
  }
  dimnames(out)[[1L]] <- rownames(block_data)
  out
}

magcomp_predict_block_proportion_samples <- function(object, block_name,
                                                    newdata, unit_levels,
                                                    latent = NULL)
{
  block <- object$data$allocation_blocks[[block_name]]
  component <- object$composition$allocation_components[[block_name]]
  block_data <- magcomp_prediction_block_data(object, block_name, newdata,
                                              unit_levels)
  training_data <- magcomp_allocation_block_data(object$data, block)
  X <- magcomp_prediction_model_matrix(component$fixed_formula,
                                       training_data,
                                       block_data,
                                       component$X,
                                       paste0("composition block ",
                                              block_name))
  coef_samples <- magcomp_block_coef_samples(object, block_name)
  n_samples <- dim(coef_samples)[3L]
  n_chains <- if(length(dim(coef_samples)) == 4L)
    dim(coef_samples)[4L] else 1L
  J <- ncol(block$counts)
  Jstick <- J - 1L
  parent_index <- NULL
  if(isTRUE(block$parent_specific_coef))
    parent_index <- rep(seq_along(block$parent_levels), each = nrow(newdata))

  loading_samples <- NULL
  factor_samples <- NULL
  factor_covariate <- NULL
  if(length(component$process_terms) > 0L){
    loading_list <- list()
    factor_list <- list()
    covariate_list <- list()
    for(term in component$process_terms){
      if(!term$term_type %in% c("factor_iid", "factor_gp", "factor_nngp"))
        next
      loading_value <- magcomp_block_loading_samples(object, block_name, term)
      factor_value <- magcomp_predict_block_factor_samples(
        object,
        block_name,
        block,
        term,
        block_data,
        newdata,
        unit_levels
      )
      if(is.null(loading_value)){
        stop("error: posterior loading samples for composition factor term '",
             term$name, "' in block '", block_name, "' were not found")
      }
      if(is.null(factor_value)){
        stop("error: posterior factor samples for composition factor term '",
             term$name, "' in block '", block_name, "' were not found")
      }
      covariate_value <- magcomp_prediction_covariate(
        term,
        block_data,
        paste0("composition block ", block_name)
      )
      loading_list[[term$name]] <- loading_value
      factor_list[[term$name]] <- factor_value
      covariate_list[[term$name]] <- covariate_value
      if(!is.null(latent)){
        # Spatial scores are shared across child blocks at the same location.
        # IID scores retain the expanded block rows, including parent groups.
        scores <- factor_value
        if(term$term_type != "factor_iid"){
          index <- lapply(dim(scores), seq_len)
          index[[1L]] <- seq_along(unit_levels)
          scores <- do.call(`[`, c(list(scores), index, list(drop = FALSE)))
          dimnames(scores)[[1L]] <- unit_levels
        }
        latent$draws$composition$blocks[[block_name]]$terms[[term$name]]$scores <- scores
      }
    }
    if(length(loading_list) == 1L){
      loading_samples <- loading_list[[1L]]
      factor_samples <- factor_list[[1L]]
      factor_covariate <- covariate_list[[1L]]
    }else if(length(loading_list) > 1L){
      loading_samples <- loading_list
      factor_samples <- factor_list
      factor_covariate <- covariate_list
    }
  }

  out <- .Call("magcomp_predict_block_proportion_cpp",
               X,
               coef_samples,
               loading_samples,
               factor_samples,
               factor_covariate,
               parent_index %||% integer(0L),
               as.integer(J),
               object$control$composition_link %||% "logit",
               PACKAGE = "magComp")
  if(length(dim(out)) == 3L)
    dimnames(out) <- list(rownames(block_data), colnames(block$counts), NULL)
  else
    dimnames(out) <- list(rownames(block_data), colnames(block$counts), NULL,
                          object$backend$chain_ids)
  out
}

magcomp_predict_proportion_samples <- function(object, newdata, unit_levels,
                                                latent = NULL)
{
  blocks <- object$data$allocation_blocks
  block_names <- names(blocks)
  if(length(blocks) == 1L){
    out <- magcomp_predict_block_proportion_samples(object, block_names[1L],
                                                    newdata, unit_levels, latent)
    dimnames(out)[[1L]] <- unit_levels
    dimnames(out)[[2L]] <- colnames(object$data$counts)
    return(out)
  }

  if(length(blocks) != 2L)
    stop("error: prediction currently supports one or two allocation blocks")

  parent <- magcomp_predict_block_proportion_samples(object, block_names[1L],
                                                     newdata, unit_levels, latent)
  child <- magcomp_predict_block_proportion_samples(object, block_names[2L],
                                                    newdata, unit_levels, latent)
  out <- .Call("magcomp_combine_conditional_cpp",
               parent,
               child,
               PACKAGE = "magComp")
  if(length(dim(out)) == 3L)
    dimnames(out) <- list(unit_levels, colnames(object$data$counts), NULL)
  else
    dimnames(out) <- list(unit_levels, colnames(object$data$counts), NULL,
                          object$backend$chain_ids)
  out
}

magcomp_predict_category_samples <- function(type, magnitude, proportion, area)
{
  type_code <- switch(type,
                      proportion = 0L,
                      density = 1L,
                      count = 2L,
                      stop("error: unsupported prediction type"))
  area_arg <- if(identical(type, "count")) area else numeric(0L)
  out <- .Call("magcomp_predict_category_cpp",
               magnitude,
               proportion,
               as.numeric(area_arg),
               as.integer(type_code),
               PACKAGE = "magComp")
  dimnames(out) <- dimnames(proportion)
  out
}

predict.magComp <- function(object, newdata,
                            type = NULL,
                            target = c("expected", "observation"),
                            summary = c("mean", "samples", "quantile"),
                            area = NULL, probs = c(0.025, 0.5, 0.975),
                            total_count = NULL, n_draws = NULL,
                            keep_effects = FALSE, ...)
{
  type <- magcomp_prediction_types(type, object)
  target <- match.arg(target)
  summary <- match.arg(summary)
  probs <- magcomp_prediction_probs(probs)
  magcomp_check_prediction_args(list(...), "predict")
  magcomp_check_prediction_fit(object)
  if(!is.logical(keep_effects) || length(keep_effects) != 1L || is.na(keep_effects))
    stop("error: keep_effects must be TRUE or FALSE")
  if(keep_effects && summary != "samples")
    stop("error: keep_effects requires summary = 'samples'")
  if(missing(newdata) || is.null(newdata))
    stop("error: newdata is required")
  newdata <- as.data.frame(newdata)
  if(nrow(newdata) < 1L)
    stop("error: newdata must contain at least one row")
  unit_levels <- magcomp_predict_unit_levels(object, newdata)
  if(!is.null(object$magnitude) && !is.null(total_count))
    stop("error: total_count is only used for composition-only predictions")
  if(is.null(object$magnitude) && (!is.null(total_count) || target == "observation" ||
                                  any(c("count", "density") %in% type)))
    total_count <- magcomp_prediction_total_count(total_count, newdata, unit_levels)
  needs_exposure <- if(is.null(object$magnitude)) "density" %in% type else
    target == "observation" || "count" %in% type
  exposure <- if(needs_exposure || !is.null(area))
    magcomp_predict_area(object, newdata, area) else NULL
  object <- magcomp_prediction_draws(object, n_draws)
  latent <- if(keep_effects) magcomp_prediction_latent(object, newdata, unit_levels) else NULL
  magnitude <- if(!is.null(object$magnitude))
    magcomp_predict_magnitude_samples(object, newdata, unit_levels, latent) else NULL
  proportion <- if(!is.null(object$composition) &&
                    (keep_effects || target == "observation" || any(type != "magnitude")))
    magcomp_predict_proportion_samples(object, newdata, unit_levels, latent) else NULL
  result <- magcomp_prediction_result(object, magnitude, proportion, exposure,
                            target, type, summary, probs, "predict", total_count)
  if(keep_effects) result$latent <- as.list(latent)
  result
}

magcomp_posterior_predict_nb <- function(count_mean, nb_size)
{
  d <- dim(count_mean)
  n_unit <- d[1L]
  n_samples <- d[3L]
  n_chains <- if(length(d) == 4L) d[4L] else 1L
  nb_size_dim <- dim(nb_size)

  if(is.null(nb_size_dim)){
    nb_size <- as.numeric(nb_size)[1L]
    if(!is.finite(nb_size) || nb_size <= 0)
      stop("error: nb_size must be positive and finite")
  }else if(length(nb_size_dim) == 2L){
    if(nb_size_dim[1L] != 1L || nb_size_dim[2L] != n_samples)
      stop("error: nb_size sample dimensions do not match fitted counts")
  }else if(length(nb_size_dim) == 3L){
    if(nb_size_dim[1L] != 1L || nb_size_dim[2L] != n_samples ||
       nb_size_dim[3L] != n_chains)
      stop("error: nb_size sample dimensions do not match fitted counts")
  }else{
    stop("error: nb_size must be a scalar, matrix, or array")
  }

  out <- array(0L, dim = d, dimnames = dimnames(count_mean))

  for(ch in seq_len(n_chains)){
    for(s in seq_len(n_samples)){
      if(is.null(nb_size_dim)){
        size <- nb_size
      }else if(length(nb_size_dim) == 3L){
        size <- nb_size[1L,s,ch]
      }else{
        size <- nb_size[1L,s]
      }
      if(!is.finite(size) || size <= 0)
        stop("error: nb_size draws must be positive and finite")
      for(i in seq_len(n_unit)){
        cell_mean <- if(length(d) == 4L) count_mean[i,,s,ch] else count_mean[i,,s]
        if(any(!is.finite(cell_mean)) || any(cell_mean < 0))
          stop("error: predictive category means must be finite and nonnegative",
               " (unit ", i, ", sample ", s, ", chain ", ch, ")")
        mu <- sum(cell_mean)
        if(!is.finite(mu))
          stop("error: predictive total means must be finite",
               " (unit ", i, ", sample ", s, ", chain ", ch, ")")
        if(mu == 0)
          next

        total <- stats::rnbinom(1L, size = size, mu = mu)
        prob <- as.numeric(cell_mean) / mu
        prob <- prob / sum(prob)

        value <- if(length(cell_mean) == 1L) total else
          as.vector(stats::rmultinom(1L, size = total, prob = prob))
        if(length(d) == 4L)
          out[i,,s,ch] <- value
        else
          out[i,,s] <- value
      }
    }
  }

  out
}

magcomp_ppc_quantile_names <- function(probs)
{
  paste0("q", gsub("\\.", "_", as.character(100 * probs)))
}

magcomp_ppc_rows <- function(draws, observed, type, unit, category, probs,
                             group = NULL)
{
  q <- magcomp_summarize_prediction(draws, "quantile", probs,
                                    magnitude = TRUE)
  colnames(q) <- magcomp_ppc_quantile_names(probs)

  p_lower <- rowMeans(draws <= observed)
  p_upper <- rowMeans(draws >= observed)
  p_value <- pmin(1, 2 * pmin(p_lower, p_upper))

  if(is.null(group))
    group <- rep(NA_character_, length(observed))

  out <- data.frame(type = type,
                    group = group,
                    unit = unit,
                    category = category,
                    observed = as.numeric(observed),
                    mean = rowMeans(draws),
                    q,
                    p_lower = p_lower,
                    p_upper = p_upper,
                    p_value = p_value,
                    n_draws = ncol(draws),
                    check.names = FALSE)
  rownames(out) <- NULL
  out
}

magcomp_ppc_group <- function(object, group_expr, group_value)
{
  group_name <- tryCatch(magcomp_name_vector(group_expr, "group"),
                         error = function(e) NULL)

  if(length(group_name) == 1L && group_name %in% names(object$data$unit_data)){
    value <- object$data$unit_data[[group_name]]
  }else{
    value <- tryCatch(group_value, error = function(e){
      stop("error: group must identify one unit-level column or contain one ",
           "value per sampling unit", call. = FALSE)
    })
    if(is.character(value) && length(value) == 1L &&
       value %in% names(object$data$unit_data)){
      value <- object$data$unit_data[[value]]
    }
  }

  if(length(value) != object$data$n_unit)
    stop("error: group must identify one unit-level column or contain one ",
         "value per sampling unit")
  if(anyNA(value))
    stop("error: group contains missing values")

  group_factor <- factor(value, levels = unique(value))
  list(index = as.integer(group_factor),
       levels = levels(group_factor),
       n = length(levels(group_factor)))
}

ppcSummary <- function(object, ...)
{
  UseMethod("ppcSummary")
}

ppcSummary.magComp <- function(object,
                               type = NULL,
                               group = NULL,
                               probs = c(0.025, 0.5, 0.975),
                               ...)
{
  choices <- c("unit_total", "category_total", "unit_category", "group_total", "group_category")
  if(is.null(object$magnitude)) choices <- c("unit_category", "category_total", "group_category")
  if(is.null(object$composition)) choices <- c("unit_total", "group_total")
  if(is.null(type)) type <- choices[1L]
  if(length(type) != 1L || !type %in% choices)
    stop("error: this model supports posterior predictive checks: ", paste(choices, collapse = ", "))
  probs <- magcomp_prediction_probs(probs)
  magcomp_check_prediction_args(list(...), "ppcSummary")
  y_pp <- stats::fitted(object, target = "observation", type = "count",
                        summary = "samples")$count
  y_obs <- object$data$counts
  if(is.null(object$composition))
    y_obs <- matrix(rowSums(y_obs), ncol = 1L, dimnames = list(rownames(y_obs), "total"))
  n <- nrow(y_obs)
  J <- ncol(y_obs)
  n_draws <- prod(dim(y_pp)[-(1:2)])

  draw_mat <- matrix(y_pp, nrow = n * J, ncol = n_draws)
  unit_id <- rep(seq_len(n), times = J)
  category_id <- rep(seq_len(J), each = n)

  if(type == "unit_total"){
    draws <- rowsum(draw_mat, unit_id, reorder = FALSE)
    observed <- rowSums(y_obs)
    out <- magcomp_ppc_rows(draws, observed, type,
                            rownames(y_obs), NA_character_, probs)
  }else if(type == "category_total"){
    draws <- rowsum(draw_mat, category_id, reorder = FALSE)
    observed <- colSums(y_obs)
    out <- magcomp_ppc_rows(draws, observed, type,
                            NA_character_, colnames(y_obs), probs)
  }else if(type == "group_total"){
    if(missing(group) || is.null(group))
      stop("error: group is required for group posterior predictive checks")
    group_info <- magcomp_ppc_group(object, substitute(group), group)
    group_cell_id <- rep(group_info$index, times = J)
    draws <- rowsum(draw_mat, group_cell_id, reorder = FALSE)
    observed <- rowsum(matrix(rowSums(y_obs), ncol = 1L),
                       group_info$index, reorder = FALSE)[,1L]
    out <- magcomp_ppc_rows(draws, observed, type,
                            NA_character_, NA_character_, probs,
                            group = group_info$levels)
  }else if(type == "group_category"){
    if(missing(group) || is.null(group))
      stop("error: group is required for group posterior predictive checks")
    group_info <- magcomp_ppc_group(object, substitute(group), group)
    group_cell_id <- rep(group_info$index, times = J)
    group_category_id <- interaction(group_cell_id, category_id,
                                     drop = TRUE, lex.order = TRUE)
    draws <- rowsum(draw_mat, group_category_id, reorder = FALSE)
    observed <- rowsum(matrix(as.numeric(y_obs), ncol = 1L),
                       group_category_id, reorder = FALSE)[,1L]
    out <- magcomp_ppc_rows(draws, observed, type,
                            NA_character_,
                            rep(colnames(y_obs), each = group_info$n),
                            probs,
                            group = rep(group_info$levels, times = J))
  }else{
    observed <- as.numeric(y_obs)
    unit <- rep(rownames(y_obs), times = J)
    category <- rep(colnames(y_obs), each = n)
    out <- magcomp_ppc_rows(draw_mat, observed, type, unit, category, probs)
  }

  attr(out, "count_conditioning") <- if(is.null(object$magnitude)) "observed_total" else "none"
  class(out) <- c("magCompPPCSummary", "data.frame")
  out
}

print.magCompPPCSummary <- function(x, digits = 4, ...)
{
  if(identical(attr(x, "count_conditioning"), "observed_total"))
    cat("Composition checks conditional on observed totals\n")
  out <- as.data.frame(x)
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols], signif, digits = digits)
  print(out, row.names = FALSE)
  invisible(x)
}

magcomp_rhat <- function(draws)
{
  if(!is.matrix(draws))
    draws <- matrix(draws, ncol = 1L)

  n <- nrow(draws)
  m <- ncol(draws)
  if(m < 2L || n < 2L)
    return(NA_real_)

  chain_vars <- apply(draws, 2L, stats::var)
  W <- mean(chain_vars)
  if(!is.finite(W) || W <= 0)
    return(NA_real_)

  B <- n * stats::var(colMeans(draws))
  var_hat <- ((n - 1) / n) * W + B / n
  max(1.0, sqrt(var_hat / W))
}

magcomp_draw_matrix <- function(x, n_samples, n_chains)
{
  matrix(as.numeric(x), nrow = n_samples, ncol = n_chains)
}

magcomp_coefficient_draws <- function(magnitude_coefficients, k)
{
  if(length(dim(magnitude_coefficients)) == 2L)
    return(magnitude_coefficients[k,])
  magnitude_coefficients[k,,]
}

magcomp_composition_draws <- function(coef_samples, j, k)
{
  if(length(dim(coef_samples)) == 3L)
    return(coef_samples[j,k,])
  coef_samples[j,k,,]
}

magcomp_loading_draws <- function(loading_samples, j, k)
{
  if(length(dim(loading_samples)) == 3L)
    return(loading_samples[j,k,])
  loading_samples[j,k,,]
}

magcomp_scalar_draws <- function(covariance_samples, row = 1L)
{
  if(length(dim(covariance_samples)) == 2L)
    return(covariance_samples[row,])
  covariance_samples[row,,]
}

magcomp_named_acceptance <- function(x, name)
{
  if(is.null(x))
    return(NULL)
  if(is.matrix(x)){
    if(name %in% rownames(x))
      return(x[name,])
    return(NULL)
  }
  if(!is.null(names(x)) && name %in% names(x))
    return(stats::setNames(unname(x[name]), "chain1"))

  NULL
}

magcomp_sample_draws <- function(object)
{
  samples <- object
  meta <- attr(object, "sample_info")
  n_chains <- meta$n_chains
  if(meta$type == "combined"){
    out <- lapply(object, function(x) matrix(x, meta$n_samples, n_chains))
    attr(out, "info") <- attr(object, "info")
    return(out)
  }
  types <- attr(object, "term_types")
  draws <- list()
  rows <- list()
  add <- function(component, block = "", term = "", parameter, index = "",
                  symbol, x, n_samples){
    draws[[length(draws) + 1L]] <<- magcomp_draw_matrix(x, meta$n_samples, n_chains)
    rows[[length(rows) + 1L]] <<- data.frame(component, block, term, parameter,
                                            index, symbol, stringsAsFactors = FALSE)
  }
  coef <- samples$magnitude$coefficients
  for(k in seq_len(if(is.null(coef)) 0L else dim(coef)[1L]))
    add("magnitude", parameter = "coefficients", index = rownames(coef)[k],
        symbol = paste0("beta[", if(rownames(coef)[k] == "(Intercept)") "0" else rownames(coef)[k], "]"),
        x = magcomp_coefficient_draws(coef, k), n_samples = dim(coef)[2L])
  size <- samples$magnitude$size
  if(!is.null(size))
    add("magnitude", parameter = "size", symbol = "r_M",
        x = size, n_samples = dim(size)[2L])
  for(term_name in names(samples$magnitude$terms)){
    term <- list(name = term_name, term_type = unname(types$magnitude[term_name]))
    part <- samples$magnitude$terms[[term$name]]
    effects <- part$effects
    if(meta$type == "effects") for(g in seq_len(dim(effects)[1L]))
      add("magnitude", term = term$name, parameter = "effects",
          index = rownames(effects)[g],
          symbol = paste0(magcomp_parameter_symbol("magnitude", "effects", term$term_type),
                          "[", rownames(effects)[g], "]"),
          x = magcomp_scalar_draws(effects, g), n_samples = meta$n_samples)
    for(parameter in intersect("variance", names(part)))
      add("magnitude", term = term$name, parameter = parameter,
          symbol = magcomp_parameter_symbol("magnitude", parameter, term$term_type),
          x = part[[parameter]], n_samples = dim(part[[parameter]])[2L])
    for(parameter in names(part$covariance)){
      x <- part$covariance[[parameter]]
      add("magnitude", term = term$name, parameter = parameter,
          symbol = magcomp_parameter_symbol("magnitude", parameter, term$term_type),
          x = magcomp_scalar_draws(x), n_samples = dim(x)[2L])
    }
  }
  for(b in names(samples$composition$blocks)){
    block <- samples$composition$blocks[[b]]
    coef <- block$coefficients
    for(j in seq_len(if(is.null(coef)) 0L else dim(coef)[1L])){
      for(k in seq_len(dim(coef)[2L])){
        column <- colnames(coef)[k]
        stick <- rownames(coef)[j]
        add("composition", block = b, parameter = "coefficients",
            index = paste(stick, column, sep = ":"),
            symbol = paste0("alpha[", if(column == "(Intercept)") "0" else column, ",", stick, "]"),
            x = magcomp_composition_draws(coef, j, k), n_samples = dim(coef)[3L])
      }
    }
    for(term_name in names(block$terms)){
      term <- list(name = term_name, term_type = unname(types$composition[[b]][term_name]))
      part <- block$terms[[term$name]]
      scores <- part$scores
      if(meta$type == "effects") for(g in seq_len(dim(scores)[1L])) for(k in seq_len(dim(scores)[2L]))
        add("composition", block = b, term = term$name, parameter = "scores",
            index = paste(rownames(scores)[g], colnames(scores)[k], sep = ":"),
            symbol = paste0(magcomp_parameter_symbol("composition", "scores", term$term_type),
                            "[", rownames(scores)[g], ",", colnames(scores)[k], "]"),
            x = magcomp_loading_draws(scores, g, k), n_samples = meta$n_samples)
      loading <- part$loadings
      for(j in seq_len(if(is.null(loading)) 0L else dim(loading)[1L])){
        for(k in seq_len(dim(loading)[2L])){
          x <- magcomp_loading_draws(loading, j, k)
          if(all(x == 0)) next
          add("composition", block = b, term = term$name, parameter = "loadings",
              index = paste(rownames(loading)[j], colnames(loading)[k], sep = ":"),
              symbol = paste0(magcomp_parameter_symbol("composition", "loadings", term$term_type),
                              "[", rownames(loading)[j], ",", colnames(loading)[k], "]"),
              x = x, n_samples = dim(loading)[3L])
        }
      }
      for(parameter in names(part$covariance)){
        x <- part$covariance[[parameter]]
        for(k in seq_len(dim(x)[1L]))
          add("composition", block = b, term = term$name, parameter = parameter,
              index = rownames(x)[k],
              symbol = paste0(magcomp_parameter_symbol("composition", parameter), "[", rownames(x)[k], "]"),
              x = magcomp_scalar_draws(x, k), n_samples = dim(x)[2L])
      }
    }
  }
  if(!length(rows)) stop("error: no samples available for this type")
  info <- do.call(rbind, rows)
  names(draws) <- apply(info[c("component", "block", "term", "parameter", "index")],
                        1L, function(x) paste(x[nzchar(x)], collapse = ":"))
  attr(draws, "info") <- info
  draws
}

magcomp_ess_one_chain <- function(x)
{
  x <- as.numeric(x)
  n <- length(x)
  if(n < 3L)
    return(NA_real_)
  if(!is.finite(stats::var(x)) || stats::var(x) <= 0)
    return(NA_real_)

  acf_values <- stats::acf(x, plot = FALSE, lag.max = n - 1L)$acf[-1L]
  paired_sum <- 0.0
  k <- 1L
  while(k <= length(acf_values)){
    pair <- acf_values[k]
    if(k + 1L <= length(acf_values))
      pair <- pair + acf_values[k + 1L]
    if(!is.finite(pair) || pair <= 0)
      break
    paired_sum <- paired_sum + pair
    k <- k + 2L
  }

  tau <- 1.0 + 2.0 * paired_sum
  if(!is.finite(tau) || tau <= 0)
    return(NA_real_)

  min(n, n / tau)
}

magcomp_ess <- function(draws)
{
  if(!is.matrix(draws))
    draws <- matrix(draws, ncol = 1L)

  ess <- apply(draws, 2L, magcomp_ess_one_chain)
  if(all(is.na(ess)))
    return(NA_real_)

  sum(ess, na.rm = TRUE)
}

magcomp_summary_row <- function(component, parameter, draws, probs)
{
  values <- as.numeric(draws)
  qs <- stats::quantile(values, probs = probs, names = FALSE)
  q_names <- paste0("q", gsub("\\.", "_", as.character(100 * probs)))
  names(qs) <- q_names

  out <- data.frame(component = component,
                    parameter = parameter,
                    mean = mean(values),
                    sd = stats::sd(values),
                    t(qs),
                    rhat = magcomp_rhat(draws),
                    ess = magcomp_ess(draws),
                    n_samples = nrow(draws),
                    n_chains = ncol(draws),
                    check.names = FALSE)
  rownames(out) <- NULL
  out
}

summary.magComp <- function(object, probs = c(0.025, 0.5, 0.975), ...)
{
  summary(samples(object, ...), probs = probs)
}

summary.magCompSamples <- function(object, probs = c(0.025, 0.5, 0.975), ...)
{
  magcomp_check_prediction_args(list(...), "summary.magCompSamples")
  probs <- magcomp_prediction_probs(probs)
  draws <- magcomp_sample_draws(object)
  info <- attr(draws, "info")
  rows <- lapply(seq_along(draws), function(i)
    magcomp_summary_row(info$component[i], info$parameter[i], draws[[i]], probs))
  out <- do.call(rbind, rows)
  out <- cbind(info, out[setdiff(names(out), c("component", "parameter"))])
  rownames(out) <- NULL
  class(out) <- c("magCompSummary", "data.frame")
  out
}

print.magCompSummary <- function(x, digits = 4, ...)
{
  out <- as.data.frame(x)
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols], signif, digits = digits)
  print(out, row.names = FALSE)
  invisible(x)
}

magcomp_last_draw <- function(x, chain, keep_dims)
{
  d <- dim(x)
  if(keep_dims == 1L){
    if(length(d) == 2L)
      return(x[,d[2L]])
    if(length(d) == 3L)
      return(x[,d[2L],chain])
  }
  if(keep_dims == 2L){
    if(length(d) == 3L){
      out <- matrix(x[,,d[3L]], nrow = d[1L], ncol = d[2L])
      dimnames(out) <- dimnames(x)[1:2]
      return(out)
    }
    if(length(d) == 4L){
      out <- matrix(x[,,d[3L],chain], nrow = d[1L], ncol = d[2L])
      dimnames(out) <- dimnames(x)[1:2]
      return(out)
    }
  }

  stop("error: unsupported MCMC sample dimensions")
}

magcomp_mean_draw <- function(x, keep_dims)
{
  d <- dim(x)
  if(keep_dims == 1L){
    if(length(d) == 2L)
      return(rowMeans(x))
    if(length(d) == 3L)
      return(apply(x, 1L, mean))
  }
  if(keep_dims == 2L){
    if(length(d) == 3L || length(d) == 4L){
      out <- apply(x, c(1L, 2L), mean)
      out <- matrix(out, nrow = d[1L], ncol = d[2L])
      dimnames(out) <- dimnames(x)[1:2]
      return(out)
    }
  }

  stop("error: unsupported MCMC sample dimensions")
}

magcomp_factor_start_list <- function(object, samples, use, chain, keep_dims)
{
  out <- list()
  if(is.null(samples))
    return(out)

  for(sample_name in names(samples)){
    term_info <- magcomp_composition_factor_output_info(object$composition,
                                                       sample_name)
    if(is.null(term_info))
      next
    value <- if(use == "last")
      magcomp_last_draw(samples[[sample_name]], chain, keep_dims)
    else
      magcomp_mean_draw(samples[[sample_name]], keep_dims)

    block_name <- term_info$block_name
    if(is.null(out[[block_name]]) || !is.list(out[[block_name]]))
      out[[block_name]] <- list()
    out[[block_name]][[term_info$term$name]] <- value
  }

  out
}

magCompStarting <- function(object, ...)
{
  UseMethod("magCompStarting")
}

magCompStarting.magComp <- function(object, use = c("initial", "last", "mean"),
                                    chain = 1, ...)
{
  use <- match.arg(use)
  if(use == "initial") return(object$starting)
  if(!isTRUE(object$fitted) || is.null(object$backend))
    stop("error: object does not contain sampler output")

  backend <- object$backend
  n_chains <- backend$n_chains %||% 1L
  chain <- as.integer(chain)[1L]
  if(is.na(chain) || chain < 1L || chain > n_chains)
    stop("error: chain must identify one fitted chain")

  out <- list()
  if(!is.null(object$magnitude)){
    out$magnitude_coefficients <- if(use == "last")
      magcomp_last_draw(backend$magnitude_coefficient_samples, chain, 1L) else
      magcomp_mean_draw(backend$magnitude_coefficient_samples, 1L)
  }

  if(!is.null(backend$nb_size_samples)){
    if(use == "last"){
      out$nb_size <- as.numeric(magcomp_last_draw(backend$nb_size_samples,
                                                  chain, 1L))
    }else{
      out$nb_size <- as.numeric(magcomp_mean_draw(backend$nb_size_samples,
                                                  1L))
    }
  }

  if(!is.null(backend$composition_block_coefficient_samples)){
    out$composition_block_coefficient <- lapply(
      backend$composition_block_coefficient_samples,
      function(x){
        if(use == "last")
          magcomp_last_draw(x, chain, 2L)
        else
          magcomp_mean_draw(x, 2L)
      }
    )
  }

  if(!is.null(backend$composition_block_loading_samples)){
    out$composition_block_loading <- magcomp_factor_start_list(
      object,
      backend$composition_block_loading_samples,
      use,
      chain,
      2L
    )
  }

  if(!is.null(backend$composition_block_score_samples)){
    out$composition_block_score <- magcomp_factor_start_list(
      object,
      backend$composition_block_score_samples,
      use,
      chain,
      2L
    )
  }

  if(!is.null(backend$composition_block_covariance_samples)){
    out$composition_block_covariance <- magcomp_factor_start_list(
      object,
      backend$composition_block_covariance_samples,
      use,
      chain,
      1L
    )
  }

  if(!is.null(backend$magnitude_process_effect_samples)){
    out$magnitude_process <- lapply(
      backend$magnitude_process_effect_samples,
      function(x){
        if(use == "last")
          magcomp_last_draw(x, chain, 1L)
        else
          magcomp_mean_draw(x, 1L)
      }
    )
  }

  if(!is.null(backend$magnitude_process_variance_samples)){
    out$magnitude_process_variance <- lapply(
      backend$magnitude_process_variance_samples,
      function(x){
        if(use == "last")
          magcomp_last_draw(x, chain, 1L)
        else
          magcomp_mean_draw(x, 1L)
      }
    )
  }

  if(!is.null(backend$magnitude_process_covariance_samples)){
    out$magnitude_process_covariance <- lapply(
      backend$magnitude_process_covariance_samples,
      function(x){
        if(use == "last")
          magcomp_last_draw(x, chain, 1L)
        else
          magcomp_mean_draw(x, 1L)
      }
    )
  }

  magnitude <- object$magnitude
  composition <- object$composition
  for(i in seq_along(magnitude$process_terms)){
    term <- magnitude$process_terms[[i]]
    x <- out$magnitude_process_covariance[[term$name]]
    if(!is.null(x)){
      term$covariance_values[] <- x
      term$covariance_starting <- magcomp_covariance_list(term$covariance_values, term$covariance_fixed,
                                                        rownames(term$covariance_values))
      magnitude$process_terms[[i]] <- term
    }
  }
  for(b in names(composition$allocation_components))
    for(i in seq_along(composition$allocation_components[[b]]$process_terms)){
      term <- composition$allocation_components[[b]]$process_terms[[i]]
      x <- out$composition_block_covariance[[b]][[term$name]]
      if(!is.null(x)){
        term$covariance_values[] <- x
        term$covariance_starting <- magcomp_covariance_list(term$covariance_values, term$covariance_fixed,
                                                          rownames(term$covariance_values))
        composition$allocation_components[[b]]$process_terms[[i]] <- term
      }
    }
  magcomp_start_from_kernel(out, magnitude, composition)
}

plot.magComp <- function(x, parameters = NULL, file = NULL,
                         n_col = 2, n_row = 3,
                         width = 8, height = 8, res = 150, ...)
{
  plot(samples(x, ...), parameters = parameters, file = file,
       n_col = n_col, n_row = n_row, width = width, height = height, res = res)
}

plot.magCompSamples <- function(x, parameters = NULL, file = NULL,
                                n_col = 2, n_row = 3,
                                width = 8, height = 8, res = 150, ...)
{
  magcomp_check_prediction_args(list(...), "plot.magCompSamples")
  draws <- magcomp_sample_draws(x)
  info <- attr(draws, "info")
  labels <- names(draws)

  if(!is.null(parameters)){
    keep <- labels %in% parameters | info$parameter %in% parameters | info$index %in% parameters
    if(!any(keep))
      stop("error: no monitored MCMC parameters matched parameters")
    draws <- draws[keep]
    labels <- labels[keep]
    info <- info[keep,,drop = FALSE]
  }

  if(length(draws) < 1L)
    stop("error: no monitored MCMC parameters available")

  n_col <- max(1L, as.integer(n_col)[1L])
  n_row <- max(1L, as.integer(n_row)[1L])
  plots_per_page <- n_col * n_row
  opened <- FALSE

  if(!is.null(file)){
    ext <- tolower(sub("^.*\\.([[:alnum:]]+)$", "\\1", file))
    if(identical(ext, "pdf")){
      grDevices::pdf(file, width = width, height = height)
    }else if(identical(ext, "png")){
      grDevices::png(file, width = width, height = height,
                     units = "in", res = res)
      if(length(draws) > plots_per_page){
        warning("PNG traceplot output contains only the first page; use ",
                "a PDF file or filter parameters to retain all traceplots",
                call. = FALSE)
        draws <- draws[seq_len(plots_per_page)]
      }
    }else{
      stop("error: traceplot file must end in .pdf or .png")
    }
    opened <- TRUE
  }

  labels <- paste0(labels, " (", info$symbol, ")")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    if(opened)
      grDevices::dev.off()
  }, add = TRUE)

  chain_cols <- c("black", "firebrick", "steelblue", "darkgreen",
                  "purple", "orange", "gray40")
  starts <- seq(1L, length(draws), by = plots_per_page)

  for(start in starts){
    stop <- min(start + plots_per_page - 1L, length(draws))
    graphics::par(mfrow = c(n_row, n_col), mar = c(3, 3, 2, 1))
    for(i in start:stop){
      mat <- draws[[i]]
      col <- chain_cols[((seq_len(ncol(mat)) - 1L) %% length(chain_cols)) + 1L]
      graphics::matplot(attr(x, "sample_info")$iterations, mat, type = "l", lty = 1,
                        col = col, xlab = "Iteration", ylab = "",
                        main = labels[i])
    }
    if(stop - start + 1L < plots_per_page){
      for(i in seq_len(plots_per_page - (stop - start + 1L)))
        graphics::plot.new()
    }
  }

  if(opened)
    invisible(normalizePath(file, mustWork = FALSE))
  else
    invisible(labels)
}
