# Recovery experiments for positive factor-loading diagonals, outside testthat.
# Run from the package root; select a fresh installation with MAGCOMP_DEV_LIB.
# Spatial data use full unit-scale GPs, including when fitting an NNGP.
dev_lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(dev_lib)) .libPaths(c(dev_lib, .libPaths()))
library(magComp)

integer_setting <- function(name, default) {
  value <- as.integer(Sys.getenv(name, as.character(default)))
  stopifnot(length(value) == 1L, is.finite(value), value > 0L)
  value
}
burn <- integer_setting("MAGCOMP_DIAG_BURN", 10000)
draws <- integer_setting("MAGCOMP_DIAG_DRAWS", 100000)
repetitions <- integer_setting("MAGCOMP_DIAG_REPS", 10)
workers <- integer_setting("MAGCOMP_DIAG_WORKERS", 1)
output <- Sys.getenv("MAGCOMP_DIAG_OUTPUT", "dev/validation/outputs/positive_diagonal")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

probabilities <- function(eta, link) {
  split <- if(link == "probit") pnorm(eta) else plogis(eta)
  p <- matrix(0, nrow(eta), ncol(eta) + 1L)
  remaining <- rep(1, nrow(eta))
  for(j in seq_len(ncol(eta))) {
    p[,j] <- remaining * split[,j]
    remaining <- remaining * (1 - split[,j])
  }
  p[,ncol(p)] <- remaining
  p
}

make_case <- function(kind, link, child = FALSE, varying = FALSE,
                      estimate_decay = TRUE, replication = 1L) {
  # Matching seeds give GP and NNGP fits the identical generated data.
  set.seed(19000 + replication + 100L * child + 1000L * varying +
             10000L * (link == "logit"))
  spatial <- kind != "iid"
  n_group <- if(spatial) 225L else 500L
  n <- if(varying) 2L * n_group else n_group
  coords <- cbind(x = runif(n_group, 0, 15), y = runif(n_group, 0, 15))
  decay <- c(0.6, 1.2)
  scores <- matrix(rnorm(n_group * 2L), n_group, 2L)
  if(spatial) {
    distance <- as.matrix(dist(coords))
    for(k in 1:2)
      scores[,k] <- drop(t(chol(exp(-decay[k] * distance))) %*% scores[,k])
  }
  group <- rep(seq_len(n_group), each = if(varying) 2L else 1L)
  z <- if(varying) rep(c(-1, 1), n_group) else rep(0, n)
  ids <- sprintf("p%04d", seq_len(n))
  units <- data.frame(plot = ids, group = sprintf("g%04d", group),
                      x = coords[group,1], y = coords[group,2], z = z, area = 1)
  loading <- if(child) rbind(c(0.5, 0), c(-0.3, 0.4), c(0.4, -0.5),
                            c(0.35, 1.6), c(0.6, -0.4), c(-0.4, 0.6),
                            c(0.4, 0.7), c(-0.5, 0.5), c(0.3, -0.6)) else
    rbind(c(0.5, 0), c(0.35, 1.6), c(-0.4, 0.55), c(0.3, -0.7))
  anchor <- if(child) c(1L, 4L) else 1:2
  if(varying) {
    loading <- cbind(c(0.45, 0.3, -0.35, 0.25), c(1.45, -0.5, 0.8, -0.4))
    eta <- scores[group,1] %o% loading[,1] +
      (z * scores[group,2]) %o% loading[,2]
  } else eta <- scores[group,,drop = FALSE] %*% t(loading)
  intercept <- if(child) rep(c(-0.65, -0.35, 0.1), 3) else c(-0.65, -0.5, -0.35, 0.1)
  eta <- sweep(eta, 2L, intercept, "+")
  if(child) {
    p <- do.call(cbind, lapply(1:3, function(a)
      probabilities(eta[,(a - 1L) * 3L + 1:3,drop = FALSE], link) / 3))
    categories <- expand.grid(dbh = c("s", "m", "l", "xl"), species = c("a", "b", "c"))
    response <- "dbh|species"
    anchors <- c("a:s", "b:s")
  } else {
    p <- probabilities(eta, link)
    categories <- data.frame(cls = letters[1:5])
    response <- "cls"
    anchors <- c("a", "b")
  }
  size <- 40
  total <- rnbinom(n, mu = if(child) 180 else 100, size = size)
  counts <- t(vapply(seq_len(n), function(i) as.vector(rmultinom(1, total[i], p[i,])),
                      numeric(ncol(p))))
  rows <- cbind(units[rep(seq_len(n), times = ncol(p)),],
                categories[rep(seq_len(nrow(categories)), each = n),,drop = FALSE])
  rows$count <- as.vector(counts)
  dat <- if(child)
    magCompData(rows, unit = plot, categories = c(species, dbh), tree_count = count,
                design = fixedArea(area), coords = c(x, y)) else
    magCompData(rows, unit = plot, category = cls, tree_count = count,
                design = fixedArea(area), coords = c(x, y))
  factor_formula <- if(varying)
    "~ z + factor_iid(group, K = 1, name = 'intercept', anchors = 'a') + z:factor_iid(group, K = 1, name = 'slope', anchors = 'a')" else
    paste0("~ factor_", kind, "(", if(spatial) "x, y" else "group",
      ", K = 2, name = 'shape', anchors = c('", paste(anchors, collapse = "', '"), "')",
      if(kind == "nngp") ", m = 15" else "",
      ")")
  composition <- if(child)
    magCompComposition(species ~ 1, as.formula(paste("dbh | species", factor_formula))) else
    as.formula(factor_formula)
  list(data = dat, composition = composition, block = response, loading = loading,
       anchors = anchor, truth = if(varying) c(0.45, 1.45) else c(0.5, 1.6),
       decay = decay, size = size, link = link, varying = varying,
       spatial = spatial, estimate_decay = estimate_decay, seed = 20000L + replication,
       scores = scores, probability = p)
}

run_case <- function(name, spec) {
  cat(format(Sys.time()), "starting", name, "\n")
  flush.console()
  case <- do.call(make_case, spec)
  control <- magCompControl(n_burn = burn, n_samples = burn + draws, n_chains = 2,
    n_cores = 2, chain_seeds = case$seed + 1:2, nb_size = case$size,
    composition_link = case$link, run_sampler = FALSE, verbose = FALSE)
  starting <- if(case$spatial && !case$estimate_decay)
    list(composition = list(blocks = setNames(list(list(terms = list(shape = list(
      covariance = list(phi = fixed(case$decay)))))), case$block))) else NULL
  prepared <- magComp(case$data, composition = case$composition, starting = starting, control = control)
  if(case$spatial && case$estimate_decay) {
    support <- magCompSettings(prepared)$priors$composition$blocks[[case$block]]$terms$shape$covariance$phi[[1]]
    stopifnot(all(case$decay > support[1L]), all(case$decay < support[2L]))
  }
  control$run_sampler <- TRUE
  elapsed <- system.time(fit <- magComp(case$data, composition = case$composition,
                                       starting = starting, control = control))[["elapsed"]]
  terms <- samples(fit)$composition$blocks[[case$block]]$terms
  parameters <- summary(fit)
  results <- list()
  diagonal_draws <- list()
  for(k in 1:2) {
    term_name <- if(case$varying) c("intercept", "slope")[k] else "shape"
    column <- if(case$varying) 1L else k
    row <- if(case$varying) 1L else case$anchors[k]
    loading <- terms[[term_name]]$loadings
    chain_draws <- loading[row,column,,]
    stopifnot(all(is.finite(chain_draws)), all(chain_draws > 0), sd(c(chain_draws)) > 0)
    if(!case$varying && k == 2L) stopifnot(all(loading[case$anchors[1L],2,,] == 0))
    interval <- quantile(chain_draws, c(0.025, 0.975), names = FALSE)
    results[[k]] <- data.frame(scenario = name, factor = k, truth = case$truth[k],
      mean = mean(chain_draws), sd = sd(c(chain_draws)), lower = interval[1], upper = interval[2],
      covered = interval[1] <= case$truth[k] && case$truth[k] <= interval[2],
      rhat = magComp:::magcomp_rhat(chain_draws),
      ess = magComp:::magcomp_ess(chain_draws),
      probability_above_one = mean(chain_draws > 1), elapsed_seconds = elapsed)
    diagonal_draws[[k]] <- chain_draws
  }
  results <- do.call(rbind, results)
  # Retain small, auditable evidence rather than hundreds of MB of latent draws.
  saveRDS(list(diagonals = diagonal_draws, results = results, parameters = parameters,
               specification = spec, truth = case$loading, control = fit$control,
               initial = magCompStarting(fit), priors = magCompSettings(fit)$priors,
               version = as.character(packageVersion("magComp"))),
           file.path(output, paste0(name, ".rds")))
  write.csv(results, file.path(output, paste0(name, ".csv")), row.names = FALSE)
  print(results, row.names = FALSE)
  flush.console()
  invisible(results)
}

cases <- list(
  iid_probit = list(kind = "iid", link = "probit"),
  iid_logit = list(kind = "iid", link = "logit"),
  gp_fixed = list(kind = "gp", link = "probit", estimate_decay = FALSE),
  gp_estimated = list(kind = "gp", link = "probit"),
  nngp_fixed = list(kind = "nngp", link = "probit", estimate_decay = FALSE),
  nngp_estimated = list(kind = "nngp", link = "probit"),
  child_iid = list(kind = "iid", link = "probit", child = TRUE),
  child_gp = list(kind = "gp", link = "probit", child = TRUE),
  child_nngp = list(kind = "nngp", link = "probit", child = TRUE),
  iid_intercept_slope = list(kind = "iid", link = "probit", varying = TRUE))
for(r in seq_len(repetitions))
  cases[[sprintf("repeat_%02d", r)]] <- list(kind = "iid", link = "probit", replication = r + 10L)
selected <- strsplit(Sys.getenv("MAGCOMP_DIAG_CASES", "all"), ",", fixed = TRUE)[[1L]]
if(!identical(selected, "all")) {
  stopifnot(all(selected %in% names(cases)))
  cases <- cases[selected]
}
run <- function(name) run_case(name, cases[[name]])
parts <- if(workers > 1L && .Platform$OS.type != "windows")
  parallel::mclapply(names(cases), run, mc.cores = workers, mc.preschedule = FALSE) else
  lapply(names(cases), run)
if(any(vapply(parts, inherits, logical(1), "try-error")))
  stop("A recovery scenario failed; inspect the run log")
results <- do.call(rbind, parts)
write.csv(results, file.path(output, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
cat("Interval coverage:", mean(results$covered), "Maximum R-hat:", max(results$rhat),
    "Minimum ESS:", min(results$ess), "\n")
# Coverage of an individual simulated truth is not a correctness assertion.
# Report mixing failures explicitly; extend affected runs before interpreting them.
if(any(results$rhat > 1.05 | results$ess < 400))
  warning("Some diagonal chains need more sampling; inspect summary.csv")
