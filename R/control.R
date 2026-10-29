magCompControl <- function(n_samples = 2000,
                           n_burn = 1000,
                           n_thin = 1,
                           n_chains = 1,
                           n_cores = 1,
                           chain_seeds = NULL,
                           n_report = 100,
                           nb_size = 1000,
                           nb_size_tuning = 0.25,
                           metropolis = list(),
                           composition_link = c("probit", "logit"),
                           run_sampler = TRUE,
                           verbose = TRUE)
{
  composition_link <- match.arg(composition_link)
  for(name in c("n_samples", "n_burn", "n_thin")){
    value <- get(name)
    if(!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
       value != floor(value) || value < (if(name == "n_burn") 0 else 1))
      stop("error: ", name, " must be a ",
           if(name == "n_burn") "non-negative" else "positive", " integer")
    if(value >= .Machine$integer.max)
      stop("error: total MCMC iterations exceed the supported integer indexing limit")
  }
  if(n_samples <= n_burn)
    stop("error: n_samples must exceed n_burn")
  if(n_thin > n_samples - n_burn)
    stop("error: n_thin must leave at least one retained sample after burn-in")
  metropolis_defaults <- list(adapt = TRUE,
                              batch_length = 25L,
                              target_accept = 0)
  if(is.null(metropolis))
    metropolis <- list()
  if(!is.list(metropolis))
    stop("error: metropolis must be a list")
  bad_metropolis <- setdiff(names(metropolis), names(metropolis_defaults))
  if(length(bad_metropolis) > 0L){
    stop("error: unknown metropolis control(s): ",
         paste(bad_metropolis, collapse = ", "))
  }
  metropolis <- utils::modifyList(metropolis_defaults, metropolis)

  out <- list(
    n_samples = as.integer(n_samples),
    n_burn = as.integer(n_burn),
    n_thin = as.integer(n_thin),
    n_chains = as.integer(n_chains),
    n_cores = as.integer(n_cores),
    chain_seeds = chain_seeds,
    n_report = as.integer(n_report),
    nb_size = as.numeric(nb_size),
    nb_size_tuning = as.numeric(nb_size_tuning),
    metropolis = list(adapt = isTRUE(metropolis$adapt),
                      batch_length = as.integer(metropolis$batch_length),
                      target_accept = as.numeric(metropolis$target_accept)),
    composition_link = composition_link,
    run_sampler = isTRUE(run_sampler),
    verbose = isTRUE(verbose)
  )

  if(length(out$n_chains) != 1L || is.na(out$n_chains) ||
     out$n_chains < 1L)
    stop("error: n_chains must be positive")
  if(length(out$n_cores) != 1L || is.na(out$n_cores) ||
     out$n_cores < 1L)
    stop("error: n_cores must be positive")
  if(!is.null(out$chain_seeds)){
    out$chain_seeds <- as.integer(out$chain_seeds)
    if(length(out$chain_seeds) != out$n_chains ||
       anyNA(out$chain_seeds)){
      stop("error: chain_seeds must be NULL or have length n_chains")
    }
  }
  if(!is.finite(out$nb_size) || out$nb_size <= 0)
    stop("error: nb_size must be positive and finite")
  if(!is.finite(out$nb_size_tuning) || out$nb_size_tuning <= 0)
    stop("error: nb_size_tuning must be positive and finite")
  if(length(out$metropolis$batch_length) != 1L ||
     is.na(out$metropolis$batch_length) ||
     out$metropolis$batch_length < 1L){
    stop("error: metropolis$batch_length must be positive")
  }
  if(length(out$metropolis$target_accept) != 1L ||
     !is.finite(out$metropolis$target_accept) ||
     out$metropolis$target_accept < 0 ||
     out$metropolis$target_accept >= 1){
    stop("error: metropolis$target_accept must be zero (automatic) or between zero and one")
  }

  class(out) <- "magCompControl"
  out
}

print.magCompControl <- function(x, ...)
{
  cat("magComp MCMC control\n")
  cat("  n_samples:   ", x$n_samples, " (total iterations per chain)\n", sep = "")
  cat("  n_burn:      ", x$n_burn, "\n", sep = "")
  cat("  n_thin:      ", x$n_thin, "\n", sep = "")
  cat("  retained:    ", (x$n_samples - x$n_burn) %/% x$n_thin,
      " per chain\n", sep = "")
  cat("  n_chains:    ", x$n_chains, "\n", sep = "")
  cat("  n_cores:     ", x$n_cores, "\n", sep = "")
  cat("  nb_size:     ", x$nb_size, "\n", sep = "")
  cat("  nb_size_tuning: ", x$nb_size_tuning, "\n", sep = "")
  cat("  composition_link: ", x$composition_link %||% "logit", "\n", sep = "")
  cat("  metropolis:  adapt=", x$metropolis$adapt,
      ", batch_length=", x$metropolis$batch_length,
      ", target_accept=", x$metropolis$target_accept, "\n", sep = "")
  cat("  run_sampler: ", x$run_sampler, "\n", sep = "")
  invisible(x)
}
