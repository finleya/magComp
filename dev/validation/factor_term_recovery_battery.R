# Development recovery checks; run from the package root with a current install.
# The smoke profile checks execution and constraints, not posterior convergence.
# Spatial truth uses a dense GP; NNGP recovery includes approximation error.
dev_lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(dev_lib))
  .libPaths(c(dev_lib, .libPaths()))

library(magComp)

set.seed(1)

env_int <- function(name, default)
{
  value <- as.integer(Sys.getenv(name, unset = as.character(default)))
  if(length(value) != 1L || is.na(value) || value < 1L)
    stop(name, " must be a positive integer")
  value
}

env_csv <- function(name, default)
{
  value <- Sys.getenv(name, unset = default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

draw_summary <- function(draws)
{
  draws <- as.numeric(draws)
  draws <- draws[is.finite(draws)]
  if(length(draws) == 0L)
    return(c(mean = NA_real_, q025 = NA_real_, q975 = NA_real_))
  q <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
  c(mean = mean(draws), q025 = q[1L], q975 = q[2L])
}

sample_entry <- function(x, i, j = NULL)
{
  d <- length(dim(x))
  if(is.null(j)){
    if(d == 2L)
      return(as.numeric(x[i,]))
    return(as.numeric(x[i,,]))
  }
  if(d == 2L)
    return(as.numeric(x[i,j]))
  if(d == 3L)
    return(as.numeric(x[i,j,]))
  as.numeric(x[i,j,,])
}

array_entry <- function(x, i, j)
{
  if(length(dim(x)) == 3L)
    return(as.numeric(x[i,j,]))
  as.numeric(x[i,j,,])
}

factor_correlation_draws <- function(samples, truth, k)
{
  truth <- as.numeric(truth)
  sample_dim <- dim(samples)
  if(length(sample_dim) < 3L)
    return(numeric(0L))
  if(length(sample_dim) == 3L){
    return(vapply(seq_len(dim(samples)[3L]), function(s) {
      stats::cor(truth, samples[,k,s])
    }, numeric(1L)))
  }

  out <- numeric(dim(samples)[3L] * dim(samples)[4L])
  index <- 1L
  for(chain in seq_len(dim(samples)[4L])){
    for(sample in seq_len(dim(samples)[3L])){
      out[index] <- stats::cor(truth, samples[,k,sample,chain])
      index <- index + 1L
    }
  }
  out
}

stick_prob <- function(eta)
{
  pi_value <- stats::pnorm(eta)
  p <- numeric(length(eta) + 1L)
  remaining <- 1
  for(j in seq_along(eta)){
    p[j] <- remaining * pi_value[j]
    remaining <- remaining * (1 - pi_value[j])
  }
  p[length(p)] <- remaining
  p
}

simulate_gp_scores <- function(coords, decay, factor_names)
{
  coords <- as.matrix(coords)
  distance <- as.matrix(stats::dist(coords))
  out <- matrix(0, nrow = nrow(coords), ncol = length(factor_names),
                dimnames = list(rownames(coords), factor_names))
  for(k in seq_along(factor_names)){
    correlation <- exp(-decay[k] * distance)
    out[,k] <- drop(t(chol(correlation)) %*% stats::rnorm(nrow(coords)))
  }
  out
}

unit_iid_scores <- function(group_scores, data)
{
  out <- group_scores[data$group,,drop = FALSE]
  rownames(out) <- data$plot
  out
}

term_surface <- function(scores, loading, covariate = NULL)
{
  out <- scores %*% t(loading)
  if(!is.null(covariate))
    out <- sweep(out, 1L, covariate, "*")
  out
}

loading_is_free <- function(row, column, anchor_index)
{
  if(column == 1L)
    return(TRUE)
  !row %in% anchor_index[seq_len(column - 1L)]
}

make_scenarios <- function()
{
  list(
    iid_intercept = list(
      term_specs = list(list(role = "iid_intercept", type = "iid",
                             varying = FALSE, estimate_decay = FALSE))
    ),
    iid_varying = list(
      term_specs = list(list(role = "iid_varying", type = "iid",
                             varying = TRUE, estimate_decay = FALSE))
    ),
    nngp_fixed_decay = list(
      term_specs = list(list(role = "nngp_intercept", type = "nngp",
                             varying = FALSE, estimate_decay = FALSE))
    ),
    nngp_estimated_decay = list(
      term_specs = list(list(role = "nngp_intercept", type = "nngp",
                             varying = FALSE, estimate_decay = TRUE))
    ),
    nngp_varying_estimated_decay = list(
      term_specs = list(list(role = "nngp_varying", type = "nngp",
                             varying = TRUE, estimate_decay = TRUE))
    ),
    iid_plus_iid_varying = list(
      term_specs = list(
        list(role = "iid_intercept", type = "iid",
             varying = FALSE, estimate_decay = FALSE),
        list(role = "iid_varying", type = "iid",
             varying = TRUE, estimate_decay = FALSE)
      )
    ),
    nngp_plus_nngp_varying = list(
      term_specs = list(
        list(role = "nngp_intercept", type = "nngp",
             varying = FALSE, estimate_decay = TRUE),
        list(role = "nngp_varying", type = "nngp",
             varying = TRUE, estimate_decay = TRUE)
      )
    ),
    iid_plus_nngp = list(
      term_specs = list(
        list(role = "iid_intercept", type = "iid",
             varying = FALSE, estimate_decay = FALSE),
        list(role = "nngp_intercept", type = "nngp",
             varying = FALSE, estimate_decay = TRUE)
      )
    ),
    rich_iid_nngp_varying = list(
      term_specs = list(
        list(role = "nngp_intercept", type = "nngp",
             varying = FALSE, estimate_decay = TRUE),
        list(role = "iid_varying", type = "iid",
             varying = TRUE, estimate_decay = FALSE),
        list(role = "iid_intercept", type = "iid",
             varying = FALSE, estimate_decay = FALSE),
        list(role = "nngp_varying", type = "nngp",
             varying = TRUE, estimate_decay = TRUE)
      )
    )
  )
}

add_term_names <- function(term_specs)
{
  type_count <- c(factor_iid = 0L, factor_nngp = 0L)
  for(i in seq_along(term_specs)){
    type_name <- if(identical(term_specs[[i]]$type, "iid"))
      "factor_iid" else "factor_nngp"
    type_count[type_name] <- type_count[type_name] + 1L
    term_specs[[i]]$term_name <- paste0(type_name, "_",
                                        type_count[type_name])
  }
  term_specs
}

term_formula_fragment <- function(spec)
{
  if(identical(spec$type, "iid")){
    call <- paste0(
      "factor_iid(group, K = 2, anchors = anchors, ",
      "factor_names = factor_names)"
    )
  }else{
    call <- paste0(
      "factor_nngp(x, y, K = 2, m = n_neighbor, ",
      "anchors = anchors, factor_names = factor_names)"
    )
  }
  if(isTRUE(spec$varying))
    call <- paste("z", call, sep = ":")
  call
}

make_composition_formula <- function(term_specs)
{
  rhs <- paste(c("z", vapply(term_specs, term_formula_fragment,
                             character(1L))),
               collapse = " + ")
  stats::as.formula(paste("~", rhs), env = parent.frame())
}

make_unit_data <- function(n_unit, n_group)
{
  groups <- paste0("g", seq_len(n_group))
  n_per_group <- ceiling(n_unit / n_group)
  balanced_group <- rep(groups, each = n_per_group)[seq_len(n_unit)]
  balanced_z <- rep(seq(-1, 1, length.out = n_per_group),
                    times = n_group)[seq_len(n_unit)]
  row_order <- sample(seq_len(n_unit))
  data.frame(
    plot = paste0("p", seq_len(n_unit)),
    area = rep(0.08, n_unit),
    group = balanced_group[row_order],
    z = as.numeric(scale(balanced_z[row_order])),
    x = stats::runif(n_unit),
    y = stats::runif(n_unit),
    stringsAsFactors = FALSE
  )
}

make_prediction_data <- function(n_pred, groups)
{
  data.frame(
    plot = paste0("pred", seq_len(n_pred)),
    area = rep(0.08, n_pred),
    group = sample(groups, n_pred, replace = TRUE),
    z = as.numeric(scale(stats::runif(n_pred, -1, 1))),
    x = stats::runif(n_pred),
    y = stats::runif(n_pred),
    stringsAsFactors = FALSE
  )
}

make_role_values <- function(unit_data, prediction_data)
{
  groups <- sort(unique(unit_data$group))
  all_coords <- rbind(unit_data[,c("x", "y")],
                      prediction_data[,c("x", "y")])
  rownames(all_coords) <- c(unit_data$plot, prediction_data$plot)

  group_scores <- list(
    iid_intercept = matrix(stats::rnorm(length(groups) * 2L),
                           nrow = length(groups),
                           dimnames = list(groups, factor_names)),
    iid_varying = matrix(stats::rnorm(length(groups) * 2L),
                         nrow = length(groups),
                         dimnames = list(groups, factor_names))
  )
  spatial_scores <- list(
    nngp_intercept = simulate_gp_scores(
      all_coords, role_decay$nngp_intercept, factor_names
    ),
    nngp_varying = simulate_gp_scores(
      all_coords, role_decay$nngp_varying, factor_names
    )
  )

  list(group_scores = group_scores,
       spatial_scores = spatial_scores)
}

term_scores <- function(spec, role_values, unit_data, prediction_data)
{
  if(identical(spec$type, "iid")){
    group_score <- role_values$group_scores[[spec$role]]
    return(list(train = unit_iid_scores(group_score, unit_data),
                pred = unit_iid_scores(group_score, prediction_data),
                parameter = group_score))
  }

  spatial_score <- role_values$spatial_scores[[spec$role]]
  list(train = spatial_score[unit_data$plot,,drop = FALSE],
       pred = spatial_score[prediction_data$plot,,drop = FALSE],
       parameter = spatial_score[unit_data$plot,,drop = FALSE])
}

build_truth <- function(term_specs, role_values, unit_data, prediction_data)
{
  loading <- list()
  train_factor <- list()
  prediction_factor <- list()
  parameter_factor <- list()
  decay <- list()

  for(spec in term_specs){
    scores <- term_scores(spec, role_values, unit_data, prediction_data)
    loading[[spec$term_name]] <- role_loading[[spec$role]]
    train_factor[[spec$term_name]] <- scores$train
    prediction_factor[[spec$term_name]] <- scores$pred
    parameter_factor[[spec$term_name]] <- scores$parameter
    if(identical(spec$type, "nngp"))
      decay[[spec$term_name]] <- role_decay[[spec$role]]
  }

  list(loading = loading,
       train_factor = train_factor,
       prediction_factor = prediction_factor,
       parameter_factor = parameter_factor,
       decay = decay)
}

composition_eta <- function(unit_data, term_specs, factor_scores, loading)
{
  eta <- matrix(rep(composition_intercept, each = nrow(unit_data)),
                nrow = nrow(unit_data),
                dimnames = list(unit_data$plot, names(composition_intercept)))
  eta <- eta + outer(unit_data$z, composition_slope)
  for(spec in term_specs){
    contribution <- term_surface(factor_scores[[spec$term_name]],
                                 loading[[spec$term_name]])
    if(isTRUE(spec$varying))
      contribution <- sweep(contribution, 1L, unit_data$z, "*")
    eta <- eta + contribution
  }
  eta
}

composition_prob <- function(eta)
{
  t(apply(eta, 1L, stick_prob))
}

simulate_counts <- function(unit_data, term_specs, factor_scores, loading)
{
  eta <- composition_eta(unit_data, term_specs, factor_scores, loading)
  p <- composition_prob(eta)
  total_mean <- unit_data$area *
    exp(magnitude_coefficients["(Intercept)"] + magnitude_coefficients["z"] * unit_data$z)
  total <- stats::rnbinom(nrow(unit_data), mu = total_mean, size = nb_size)
  counts <- matrix(0L, nrow = nrow(unit_data), ncol = length(classes),
                   dimnames = list(unit_data$plot, classes))
  for(i in seq_len(nrow(unit_data))){
    counts[i,] <- as.vector(stats::rmultinom(1L, total[i], p[i,]))
  }
  list(counts = counts, probability = p, total_mean = total_mean)
}

count_table <- function(unit_data, counts)
{
  do.call(rbind, lapply(seq_len(nrow(unit_data)), function(i) {
    data.frame(
      plot = unit_data$plot[i],
      class = factor(classes, levels = classes),
      tree_count = as.integer(counts[i,]),
      area = unit_data$area[i],
      group = unit_data$group[i],
      z = unit_data$z[i],
      x = unit_data$x[i],
      y = unit_data$y[i],
      stringsAsFactors = FALSE
    )
  }))
}

make_truth_start <- function(term_specs, truth, role_values, unit_data)
{
  group_order <- unique(unit_data$group)
  factor_start <- list()
  for(spec in term_specs){
    if(identical(spec$type, "iid")){
      role_score <- role_values$group_scores[[spec$role]]
      factor_start[[spec$term_name]] <- role_score[group_order,,drop = FALSE]
    }else{
      factor_start[[spec$term_name]] <- truth$train_factor[[spec$term_name]]
    }
  }

  terms <- list()
  for(spec in term_specs){
    name <- spec$term_name
    terms[[name]] <- list(loadings = truth$loading[[name]],
                          scores = factor_start[[name]])
    if(isTRUE(spec$estimate_decay))
      terms[[name]]$covariance <- list(phi = truth$decay[[name]])
  }
  list(
    magnitude = list(coefficients = magnitude_coefficients),
    composition = list(blocks = setNames(list(list(
      coefficients = cbind("(Intercept)" = composition_intercept, z = composition_slope),
      terms = terms
    )), block_name))
  )
}

add_parameter_row <- function(rows, scenario, start_type, component, term,
                              factor, parameter, truth, draws,
                              fixed = FALSE)
{
  s <- draw_summary(draws)
  rows[[length(rows) + 1L]] <- data.frame(
    scenario = scenario,
    start_type = start_type,
    component = component,
    term = term,
    factor = factor,
    parameter = parameter,
    truth = truth,
    mean = s["mean"],
    q025 = s["q025"],
    q975 = s["q975"],
    fixed = fixed,
    covered = truth >= s["q025"] && truth <= s["q975"],
    stringsAsFactors = FALSE
  )
  rows
}

summarize_parameters <- function(fit, scenario, start_type, truth)
{
  rows <- list()
  samples <- samples(fit)
  block <- samples$composition$blocks[[block_name]]

  for(j in seq_along(magnitude_coefficients)){
    rows <- add_parameter_row(
      rows, scenario, start_type, "magnitude", "fixed", NA_character_,
      names(magnitude_coefficients)[j], magnitude_coefficients[j],
      sample_entry(samples$magnitude$coefficients, j)
    )
  }

  coef_samples <- block$coefficients
  for(j in seq_along(composition_intercept)){
    rows <- add_parameter_row(
      rows, scenario, start_type, "composition_intercept", "fixed",
      NA_character_, names(composition_intercept)[j], composition_intercept[j],
      sample_entry(coef_samples, j, 1L)
    )
    rows <- add_parameter_row(
      rows, scenario, start_type, "composition_z", "fixed",
      NA_character_, names(composition_slope)[j], composition_slope[j],
      sample_entry(coef_samples, j, 2L)
    )
  }

  for(term_name in names(truth$loading)){
    loading_samples <- block$terms[[term_name]]$loadings
    for(k in seq_len(ncol(truth$loading[[term_name]]))){
      for(j in seq_len(nrow(truth$loading[[term_name]]))){
        rows <- add_parameter_row(
          rows, scenario, start_type, "loading", term_name,
          colnames(truth$loading[[term_name]])[k],
          rownames(truth$loading[[term_name]])[j],
          truth$loading[[term_name]][j,k],
          sample_entry(loading_samples, j, k),
          fixed = !loading_is_free(j, k, anchor_index)
        )
      }
    }

    factor_samples <- block$terms[[term_name]]$scores
    factor_rows <- dimnames(factor_samples)[[1L]]
    if(is.null(factor_rows))
      factor_rows <- rownames(truth$parameter_factor[[term_name]])
    truth_factor <- truth$parameter_factor[[term_name]][factor_rows,,drop = FALSE]
    for(k in seq_len(ncol(truth_factor))){
      for(i in seq_len(nrow(truth_factor))){
        rows <- add_parameter_row(
          rows, scenario, start_type, "factor_score", term_name,
          colnames(truth_factor)[k], rownames(truth_factor)[i],
          truth_factor[i,k], sample_entry(factor_samples, i, k)
        )
      }
      rows <- add_parameter_row(
        rows, scenario, start_type, "factor_score_correlation", term_name,
        colnames(truth_factor)[k], "correlation", 1,
        factor_correlation_draws(factor_samples, truth_factor[,k], k)
      )
    }
  }

  for(term_name in names(truth$decay)){
    decay_samples <- block$terms[[term_name]]$covariance$phi
    for(k in seq_along(truth$decay[[term_name]])){
      rows <- add_parameter_row(
        rows, scenario, start_type, "decay", term_name, factor_names[k],
        "decay", truth$decay[[term_name]][k],
        if(is.null(decay_samples)) truth$decay[[term_name]][k] else
          sample_entry(decay_samples, k), fixed = is.null(decay_samples)
      )
    }
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

summarize_prediction <- function(samples, truth, scenario, start_type,
                                 quantity)
{
  rows <- list()
  for(i in seq_len(dim(samples)[1L])){
    for(j in seq_len(dim(samples)[2L])){
      s <- draw_summary(array_entry(samples, i, j))
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = scenario,
        start_type = start_type,
        quantity = quantity,
        plot = dimnames(samples)[[1L]][i],
        class = dimnames(samples)[[2L]][j],
        truth = truth[i,j],
        mean = s["mean"],
        q025 = s["q025"],
        q975 = s["q975"],
        covered = truth[i,j] >= s["q025"] && truth[i,j] <= s["q975"],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

assert_mechanics <- function(fit, scenario, start_type, term_specs)
{
  term_names <- vapply(term_specs, `[[`, character(1L), "term_name")
  terms <- samples(fit)$composition$blocks[[block_name]]$terms
  starts <- magCompStarting(fit)$composition$blocks[[block_name]]$terms
  if(!identical(names(terms), term_names) || !setequal(names(starts), term_names))
    stop(scenario, "/", start_type, ": factor term identities do not match")

  for(spec in term_specs){
    term <- terms[[spec$term_name]]
    if(isTRUE(spec$estimate_decay) != !is.null(term$covariance$phi))
      stop(scenario, "/", start_type, ": decay sampling status does not match")
    loading_samples <- term$loadings
    for(k in seq_along(anchors)){
      diagonal <- as.numeric(sample_entry(loading_samples, anchor_index[k], k))
      if(!all(is.finite(diagonal) & diagonal > 0) || stats::sd(diagonal) == 0)
        stop(scenario, "/", start_type, ": diagonal loadings must be positive and sampled")
      if(k > 1L){
        for(h in seq_len(k - 1L)){
          upper <- as.numeric(sample_entry(loading_samples,
                                           anchor_index[h], k))
          if(!all(upper == 0))
            stop(scenario, "/", start_type, ": loading upper entries are not 0")
        }
      }
    }
  }
}

save_recovery_plot <- function(data, file_base, color_field = "term",
                               max_points = 6000L)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    return(invisible(NULL))

  keep <- is.finite(data$truth) & is.finite(data$mean) &
    is.finite(data$q025) & is.finite(data$q975)
  data <- data[keep,,drop = FALSE]
  if(nrow(data) == 0L)
    return(invisible(NULL))
  if(nrow(data) > max_points)
    data <- data[sort(sample(seq_len(nrow(data)), max_points)),,drop = FALSE]

  if(!"factor" %in% names(data))
    data$factor <- "value"
  data$factor[is.na(data$factor)] <- "value"
  data$panel <- paste(data$scenario, data$start_type, sep = "\n")
  data$color_value <- data[[color_field]]
  range_xy <- range(data$truth, data$mean, data$q025, data$q975,
                    finite = TRUE)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = truth, y = mean, color = color_value)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey55",
                         linewidth = 0.30) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q025, ymax = q975),
                            alpha = 0.20, linewidth = 0.20) +
    ggplot2::geom_point(ggplot2::aes(shape = factor), alpha = 0.70,
                        size = 1.1) +
    ggplot2::coord_equal(xlim = range_xy, ylim = range_xy) +
    ggplot2::facet_wrap(ggplot2::vars(panel), scales = "fixed") +
    ggplot2::labs(x = "True", y = "Posterior mean",
                  color = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(paste0(file_base, ".png"), p,
                  width = 10, height = 7, dpi = 300)
  ggplot2::ggsave(paste0(file_base, ".pdf"), p,
                  width = 10, height = 7)
  invisible(p)
}

run_scenario <- function(scenario_name, scenario, start_type)
{
  message("running ", scenario_name, " with ", start_type, " starts")

  term_specs <- add_term_names(scenario$term_specs)
  unit_data <- make_unit_data(n_unit, n_group)
  prediction_data <- make_prediction_data(n_pred, sort(unique(unit_data$group)))
  role_values <<- make_role_values(unit_data, prediction_data)
  truth <- build_truth(term_specs, role_values, unit_data, prediction_data)
  simulated <- simulate_counts(unit_data, term_specs, truth$train_factor,
                               truth$loading)
  tree_data <- count_table(unit_data, simulated$counts)

  dat <- magCompData(
    tree_data,
    unit = plot,
    category = class,
    tree_count = tree_count,
    ordered = c(class = TRUE),
    coords = c(x, y),
    design = fixedArea(area)
  )

  starting <- if(identical(start_type, "truth")){
    make_truth_start(term_specs, truth, role_values, unit_data)
  }else{
    estimated_terms <- Filter(function(spec) isTRUE(spec$estimate_decay), term_specs)
    terms <- setNames(lapply(estimated_terms, function(spec) list(covariance = list(phi = decay_start))),
                      vapply(estimated_terms, `[[`, character(1L), "term_name"))
    list(composition = list(blocks = setNames(list(list(terms = terms)), block_name)))
  }

  for(spec in term_specs){
    if(spec$type != "iid" && !isTRUE(spec$estimate_decay))
      starting$composition$blocks[[block_name]]$terms[[spec$term_name]]$covariance <-
        list(phi = fixed(role_decay[[spec$role]]))
  }

  estimated_terms <- Filter(function(spec) isTRUE(spec$estimate_decay), term_specs)
  term_priors <- setNames(
    lapply(estimated_terms, function(spec) list(covariance = list(phi = uniform(decay_support[1], decay_support[2])))),
    vapply(estimated_terms, `[[`, character(1L), "term_name")
  )

  fit <- magComp(
    dat,
    family = fixedAreaStems(),
    magnitude = ~ z,
    composition = make_composition_formula(term_specs),
    starting = starting,
    priors = magCompPriors(
      magnitude = list(coefficient_variance = 25),
      composition = list(intercept_variance = 25, slope_variance = 25,
        loading_variance = 9,
        blocks = setNames(list(list(terms = term_priors)), block_name))
    ),
    smooth = magCompSmooth(intercept = "none", slope = "none", loading = "none"),
    control = magCompControl(
      n_samples = n_burn + n_samples,
      n_burn = n_burn,
      n_thin = 1,
      n_chains = n_chains,
      n_cores = min(n_chains, 2L),
      chain_seeds = chain_seed + seq_len(n_chains),
      nb_size = nb_size,
      composition_link = "probit",
      verbose = FALSE
    )
  )

  assert_mechanics(fit, scenario_name, start_type, term_specs)

  truth_pred <- simulate_counts(prediction_data, term_specs,
                                truth$prediction_factor, truth$loading)
  truth_probability <- truth_pred$probability
  dimnames(truth_probability) <- list(prediction_data$plot, classes)
  truth_count <- sweep(truth_probability, 1L, truth_pred$total_mean, "*")

  prediction <- predict(fit, prediction_data, type = c("proportion", "count"),
                         summary = "samples", area = "area")
  pred_prop <- prediction$proportion
  prop_sum <- if(length(dim(pred_prop)) == 3L){
    as.numeric(apply(pred_prop, c(1L, 3L), sum))
  }else{
    as.numeric(apply(pred_prop, c(1L, 3L, 4L), sum))
  }
  if(!all(abs(prop_sum - 1) < 1e-10)){
    stop(scenario_name, "/", start_type,
         ": predicted proportions do not sum to one")
  }
  pred_count <- prediction$count
  if(!all(is.finite(pred_count)))
    stop(scenario_name, "/", start_type, ": predicted counts are not finite")

  mcmc_summary <- summary(fit)
  mcmc_summary <- data.frame(scenario = scenario_name,
                             start_type = start_type,
                             mcmc_summary,
                             stringsAsFactors = FALSE)

  metadata <- data.frame(
    scenario = scenario_name,
    start_type = start_type,
    n_sampling_units = n_unit,
    n_prediction_units = n_pred,
    n_groups = n_group,
    n_classes = length(classes),
    total_tree_count = sum(tree_data$tree_count),
    n_samples = n_samples,
    n_burn = n_burn,
    n_chains = n_chains,
    n_neighbor = n_neighbor,
    stringsAsFactors = FALSE
  )

  list(
    parameters = summarize_parameters(fit, scenario_name, start_type, truth),
    prediction = rbind(
      summarize_prediction(pred_prop, truth_probability, scenario_name,
                           start_type, "proportion"),
      summarize_prediction(pred_count, truth_count, scenario_name,
                           start_type, "count")
    ),
    mcmc_summary = mcmc_summary,
    metadata = metadata
  )
}

profile <- Sys.getenv("MAGCOMP_FACTOR_BATTERY_PROFILE", "smoke")
profile <- match.arg(profile, c("smoke", "standard"))
if(identical(profile, "standard")){
  n_unit_default <- 600L
  n_pred_default <- 120L
  n_group_default <- 60L
  n_samples_default <- 1500L
  n_burn_default <- 1500L
}else{
  n_unit_default <- 80L
  n_pred_default <- 20L
  n_group_default <- 10L
  n_samples_default <- 30L
  n_burn_default <- 30L
}

out_dir <- file.path("dev", "validation", "outputs",
                     "factor_term_recovery")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_unit <- env_int("MAGCOMP_FACTOR_BATTERY_N", n_unit_default)
n_pred <- env_int("MAGCOMP_FACTOR_BATTERY_N_PRED", n_pred_default)
n_group <- env_int("MAGCOMP_FACTOR_BATTERY_GROUPS", n_group_default)
n_samples <- env_int("MAGCOMP_FACTOR_BATTERY_SAMPLES", n_samples_default)
n_burn <- env_int("MAGCOMP_FACTOR_BATTERY_BURN", n_burn_default)
n_chains <- env_int("MAGCOMP_FACTOR_BATTERY_CHAINS", 2L)
n_neighbor <- env_int("MAGCOMP_FACTOR_BATTERY_M", 15L)
chain_seed <- env_int("MAGCOMP_FACTOR_BATTERY_CHAIN_SEED", 8100L)
if(n_group > n_unit)
  stop("MAGCOMP_FACTOR_BATTERY_GROUPS cannot exceed MAGCOMP_FACTOR_BATTERY_N")

classes <- c("a", "b", "c", "d", "e", "f")
block_name <- "class"
stick_names <- classes[-length(classes)]
anchors <- c("b", "d")
anchor_index <- match(anchors, stick_names)
factor_names <- c("f1", "f2")
decay_support <- c(lower = 1, upper = 10)
nb_size <- 10000
decay_start <- c(f1 = 2.5, f2 = 5.0)

magnitude_coefficients <- c("(Intercept)" = log(750), z = 0.2)
composition_intercept <- c(a = -0.5, b = -0.4, c = -0.3, d = -0.6, e = -0.8)
composition_slope <- c(a = -0.2, b = -0.1, c = 0.0, d = 0.1, e = 0.2)

role_loading <- list(
  iid_intercept = rbind(a = c(0.4, -0.2),
                        b = c(1.0, 0.0),
                        c = c(-0.3, 0.3),
                        d = c(0.2, 1.0),
                        e = c(-0.2, -0.4)),
  iid_varying = rbind(a = c(0.2, 0.3),
                      b = c(1.0, 0.0),
                      c = c(0.3, -0.2),
                      d = c(-0.2, 1.0),
                      e = c(0.1, 0.4)),
  nngp_intercept = rbind(a = c(0.3, -0.4),
                         b = c(1.0, 0.0),
                         c = c(-0.2, 0.2),
                         d = c(0.4, 1.0),
                         e = c(-0.3, 0.3)),
  nngp_varying = rbind(a = c(-0.2, 0.2),
                       b = c(1.0, 0.0),
                       c = c(0.4, -0.3),
                       d = c(0.3, 1.0),
                       e = c(0.2, -0.2))
)
for(role in names(role_loading))
  colnames(role_loading[[role]]) <- factor_names

role_decay <- list(
  nngp_intercept = c(f1 = 2.0, f2 = 4.0),
  nngp_varying = c(f1 = 3.0, f2 = 6.0)
)

scenarios <- make_scenarios()
scenario_names <- env_csv("MAGCOMP_FACTOR_BATTERY_SCENARIOS", "all")
if(!identical(scenario_names, "all")){
  missing <- setdiff(scenario_names, names(scenarios))
  if(length(missing) > 0L)
    stop("unknown scenario(s): ", paste(missing, collapse = ", "))
  scenarios <- scenarios[scenario_names]
}

start_types <- env_csv("MAGCOMP_FACTOR_BATTERY_STARTS", "svd,truth")
unknown_starts <- setdiff(start_types, c("svd", "truth"))
if(length(unknown_starts) > 0L)
  stop("unknown start type(s): ", paste(unknown_starts, collapse = ", "))

results <- list()
index <- 1L
for(scenario_name in names(scenarios)){
  for(start_type in start_types){
    results[[index]] <- run_scenario(scenario_name, scenarios[[scenario_name]],
                                     start_type)
    index <- index + 1L
  }
}

parameters <- do.call(rbind, lapply(results, `[[`, "parameters"))
prediction <- do.call(rbind, lapply(results, `[[`, "prediction"))
mcmc_summary <- do.call(rbind, lapply(results, `[[`, "mcmc_summary"))
metadata <- do.call(rbind, lapply(results, `[[`, "metadata"))

write.csv(parameters, file.path(out_dir, "parameter_recovery.csv"),
          row.names = FALSE)
write.csv(prediction, file.path(out_dir, "prediction_recovery.csv"),
          row.names = FALSE)
write.csv(mcmc_summary, file.path(out_dir, "mcmc_summary.csv"),
          row.names = FALSE)
write.csv(metadata, file.path(out_dir, "metadata.csv"), row.names = FALSE)

error_summary <- aggregate(
  abs(mean - truth) ~ scenario + start_type + component + term,
  parameters[parameters$component != "factor_score_correlation",],
  max
)
names(error_summary)[5L] <- "max_abs_error"
write.csv(error_summary, file.path(out_dir, "error_summary.csv"),
          row.names = FALSE)

correlation_summary <- parameters[
  parameters$component == "factor_score_correlation",
  c("scenario", "start_type", "term", "factor", "mean", "q025", "q975")
]
write.csv(correlation_summary,
          file.path(out_dir, "factor_score_correlation.csv"),
          row.names = FALSE)

prediction_metric_data <- transform(
  prediction,
  sq_error = (mean - truth)^2,
  width = q975 - q025
)
prediction_summary <- aggregate(
  cbind(sq_error, covered, width) ~ scenario + start_type + quantity,
  prediction_metric_data,
  mean
)
prediction_summary$rmse <- sqrt(prediction_summary$sq_error)
write.csv(prediction_summary,
          file.path(out_dir, "prediction_summary.csv"),
          row.names = FALSE)

save_recovery_plot(
  parameters[parameters$component == "loading",],
  file.path(out_dir, "loading_recovery")
)
save_recovery_plot(
  parameters[parameters$component == "decay",],
  file.path(out_dir, "decay_recovery")
)
save_recovery_plot(
  parameters[parameters$component == "factor_score",],
  file.path(out_dir, "factor_score_recovery")
)
save_recovery_plot(
  prediction[prediction$quantity == "proportion",],
  file.path(out_dir, "prediction_proportion_recovery"),
  color_field = "class"
)
save_recovery_plot(
  prediction[prediction$quantity == "count",],
  file.path(out_dir, "prediction_count_recovery"),
  color_field = "class"
)

print(error_summary)
print(correlation_summary)
cat("wrote factor-term recovery battery output to ", out_dir, "\n", sep = "")
