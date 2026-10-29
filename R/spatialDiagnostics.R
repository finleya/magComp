spatialDiagnostics <- function(object, ...)
{
  UseMethod("spatialDiagnostics")
}

magcomp_spatial_draw_matrix <- function(x, row = 1L)
{
  d <- dim(x)
  if(length(d) == 2L){
    out <- matrix(x[row,], ncol = 1L)
    colnames(out) <- "chain1"
    return(out)
  }
  if(length(d) == 3L){
    out <- matrix(x[row,,], nrow = d[2L], ncol = d[3L])
    chain_names <- dimnames(x)[[3L]]
    colnames(out) <- chain_names %||% paste0("chain", seq_len(d[3L]))
    return(out)
  }

  stop("error: unsupported spatial diagnostic sample dimensions")
}

magcomp_loading_norm_matrix <- function(x, k)
{
  d <- dim(x)
  if(length(d) == 3L){
    out <- numeric(d[3L])
    for(s in seq_len(d[3L]))
      out[s] <- sqrt(sum(x[,k,s]^2))
    out <- matrix(out, ncol = 1L)
    colnames(out) <- "chain1"
    return(out)
  }
  if(length(d) == 4L){
    out <- matrix(NA_real_, nrow = d[3L], ncol = d[4L])
    for(ch in seq_len(d[4L]))
      for(s in seq_len(d[3L]))
        out[s,ch] <- sqrt(sum(x[,k,s,ch]^2))
    chain_names <- dimnames(x)[[4L]]
    colnames(out) <- chain_names %||% paste0("chain", seq_len(d[4L]))
    return(out)
  }

  stop("error: unsupported loading sample dimensions")
}

magcomp_spatial_summary_row <- function(component, term, quantity, draws,
                                        probs)
{
  values <- as.numeric(draws)
  qs <- stats::quantile(values, probs = probs, names = FALSE)
  q_names <- paste0("q", gsub("\\.", "_", as.character(100 * probs)))
  names(qs) <- q_names

  out <- data.frame(component = component,
                    term = term,
                    quantity = quantity,
                    mean = mean(values),
                    sd = stats::sd(values),
                    t(qs),
                    n_samples = nrow(draws),
                    n_chains = ncol(draws),
                    check.names = FALSE)
  rownames(out) <- NULL
  out
}

magcomp_spatial_acceptance_rows <- function(component, term, acceptance)
{
  if(is.null(acceptance) || all(!is.finite(acceptance)))
    return(NULL)

  chain_names <- names(acceptance)
  if(is.null(chain_names))
    chain_names <- paste0("chain", seq_along(acceptance))

  data.frame(component = component,
             term = term,
             chain = chain_names,
             acceptance = as.numeric(acceptance),
             row.names = NULL)
}

magcomp_spatial_correlation_rows <- function(component, term, x_name, y_name,
                                             x, y)
{
  rows <- list()
  for(ch in seq_len(ncol(x))){
    chain_name <- colnames(x)[ch] %||% paste0("chain", ch)
    if(stats::sd(x[,ch]) <= 0 || stats::sd(y[,ch]) <= 0){
      value <- NA_real_
    }else{
      value <- stats::cor(x[,ch], y[,ch])
    }
    rows[[length(rows) + 1L]] <- data.frame(
      component = component,
      term = term,
      x = x_name,
      y = y_name,
      chain = chain_name,
      correlation = value,
      row.names = NULL
    )
  }

  if(stats::sd(as.numeric(x)) <= 0 || stats::sd(as.numeric(y)) <= 0){
    value <- NA_real_
  }else{
    value <- stats::cor(as.numeric(x), as.numeric(y))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    component = component,
    term = term,
    x = x_name,
    y = y_name,
    chain = "pooled",
    correlation = value,
    row.names = NULL
  )

  do.call(rbind, rows)
}

spatialDiagnostics.magComp <- function(object, probs = c(0.025, 0.5, 0.975),
                                       ...)
{
  if(!isTRUE(object$fitted) || is.null(object$backend))
    stop("error: object does not contain sampler output")
  if(anyNA(probs) || any(probs <= 0) || any(probs >= 1))
    stop("error: probs must be between zero and one")

  backend <- object$backend
  summary_rows <- acceptance_rows <- correlation_rows <- list()
  # Add the same named-term identity to scalar and derived diagnostics.
  add <- function(kind, value, block = "", factor = "", symbol = ""){
    if(is.null(value)) return(invisible(NULL))
    value$block <- block
    value$factor <- factor
    if(kind == "summary") value$symbol <- symbol
    first <- intersect(c("component", "block", "term", "factor", "quantity", "symbol"), names(value))
    value <- value[c(first, setdiff(names(value), first))]
    if(kind == "summary") summary_rows[[length(summary_rows) + 1L]] <<- value
    if(kind == "acceptance") acceptance_rows[[length(acceptance_rows) + 1L]] <<- value
    if(kind == "correlation") correlation_rows[[length(correlation_rows) + 1L]] <<- value
  }
  for(term in object$magnitude$process_terms){
    if(!term$term_type %in% c("gp", "nngp")) next
    t <- term$name
    covariance <- backend$magnitude_process_covariance_samples[[t]]
    variance <- magcomp_spatial_draw_matrix(backend$magnitude_process_variance_samples[[t]])
    add("summary", magcomp_spatial_summary_row("magnitude", t, "variance", variance, probs), symbol = "sigma_v^2")
    for(p in seq_along(term$covariance_parameters)){
      name <- term$covariance_parameters[p]
      draws <- magcomp_spatial_draw_matrix(covariance, p)
      add("summary", magcomp_spatial_summary_row("magnitude", t, name, draws, probs),
          symbol = magcomp_parameter_symbol("magnitude", name))
      if(!term$covariance_fixed[1,p])
        add("correlation", magcomp_spatial_correlation_rows("magnitude", t, name, "variance", draws, variance))
    }
    acceptance <- magcomp_named_acceptance(backend$magnitude_process_covariance_acceptance, t)
    add("acceptance", magcomp_spatial_acceptance_rows("magnitude", t, acceptance))
  }
  for(b in names(object$composition$allocation_components)){
    for(term in object$composition$allocation_components[[b]]$process_terms){
      if(!term$term_type %in% c("factor_gp", "factor_nngp")) next
      t <- term$name
      key <- magcomp_block_factor_sample_key(object, b, term)
      covariance <- magcomp_block_covariance_samples(object, b, term)
      loading <- magcomp_block_loading_samples(object, b, term)
      factor_names <- term$factor_names
      for(k in seq_along(factor_names)){
        factor <- factor_names[k]
        norm <- magcomp_loading_norm_matrix(loading, k)
        add("summary", magcomp_spatial_summary_row("composition", t, "loading_norm", norm, probs),
            b, factor, paste0("||lambda[,", factor, "]||"))
        for(p in seq_along(term$covariance_parameters)){
          name <- term$covariance_parameters[p]
          draws <- magcomp_spatial_draw_matrix(covariance, k + term$K * (p-1L))
          add("summary", magcomp_spatial_summary_row("composition", t, name, draws, probs),
              b, factor, paste0(magcomp_parameter_symbol("composition", name), "[", factor, "]"))
          if(!term$covariance_fixed[k,p])
            add("correlation", magcomp_spatial_correlation_rows("composition", t, name, "loading_norm", draws, norm), b, factor)
        }
        acceptance_key <- if(length(factor_names) == 1L) key else paste(key, factor, sep = ":")
        acceptance <- magcomp_named_acceptance(backend$composition_block_covariance_acceptance, acceptance_key)
        add("acceptance", magcomp_spatial_acceptance_rows("composition", t, acceptance), b, factor)
      }
    }
  }
  empty <- data.frame(component = character(), block = character(), term = character(), factor = character())
  summary <- if(length(summary_rows)) do.call(rbind, summary_rows) else
    cbind(empty, data.frame(quantity = character(), symbol = character(), mean = numeric(), sd = numeric()))
  acceptance <- if(length(acceptance_rows)) do.call(rbind, acceptance_rows) else
    cbind(empty, data.frame(chain = character(), acceptance = numeric()))
  correlations <- if(length(correlation_rows)) do.call(rbind, correlation_rows) else
    cbind(empty, data.frame(x = character(), y = character(), chain = character(), correlation = numeric()))
  structure(list(summary = summary, acceptance = acceptance, correlations = correlations),
             class = "magCompSpatialDiagnostics")
}

print.magCompSpatialDiagnostics <- function(x, digits = 4, ...)
{
  cat("magComp spatial diagnostics\n")

  if(nrow(x$summary) > 0L){
    cat("\nSummary:\n")
    out <- x$summary
    numeric_cols <- vapply(out, is.numeric, logical(1L))
    out[numeric_cols] <- lapply(out[numeric_cols], signif, digits = digits)
    print(out, row.names = FALSE)
  }else{
    cat("  no spatial terms found\n")
  }

  if(nrow(x$acceptance) > 0L){
    cat("\nMetropolis acceptance:\n")
    out <- x$acceptance
    numeric_cols <- vapply(out, is.numeric, logical(1L))
    out[numeric_cols] <- lapply(out[numeric_cols], signif, digits = digits)
    print(out, row.names = FALSE)
  }

  if(nrow(x$correlations) > 0L){
    cat("\nPosterior correlations:\n")
    out <- x$correlations
    numeric_cols <- vapply(out, is.numeric, logical(1L))
    out[numeric_cols] <- lapply(out[numeric_cols], signif, digits = digits)
    print(out, row.names = FALSE)
  }

  invisible(x)
}
