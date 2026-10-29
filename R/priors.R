magcomp_prior_defaults <- function()
{
  list(magnitude = list(coefficient_variance = 100,
                        iid_variance = c(shape = 2, scale = 1),
                        spatial_variance = c(shape = 2, scale = 1),
                        size = "auto", terms = list()),
       composition = list(intercept_variance = 9, slope_variance = 100,
                          intercept_smoothing_variance = 1,
                          slope_smoothing_variance = 4,
                          loading_variance = 1, blocks = list()))
}

magcomp_prior_entries <- function(x, defaults, path)
{
  magcomp_check_names(x, names(defaults), path)
  out <- defaults
  for(nm in names(x)){
    value <- x[[nm]]
    if(nm %in% c("terms", "blocks", "covariance")){
      magcomp_check_names(value, names(value), paste0(path, "$", nm))
      out[[nm]] <- value
    }else if(nm == "size" && (is.null(value) || identical(value, "auto"))){
      out[nm] <- list(value)
    }else if(is.null(defaults[[nm]]) && is.null(value)){
      out[nm] <- list(NULL)
    }else{
      pair <- nm %in% c("iid_variance", "spatial_variance", "variance", "size")
      if(pair && !is.null(names(value))){
        expected <- if(nm == "size") c("shape", "rate") else
          c("shape", "scale")
        if(anyDuplicated(names(value)) || !setequal(names(value), expected))
          stop("error: ", path, "$", nm, " names must be ", paste(expected, collapse = ", "))
        value <- value[expected]
      }
      value <- magcomp_positive(value, if(pair) 2L else 1L, paste0(path, "$", nm))
      if(pair){
        names(value) <- if(nm == "size") c("shape", "rate") else
          c("shape", "scale")
      }
      out[[nm]] <- value
    }
  }
  out
}

magCompPriors <- function(magnitude = list(), composition = list())
{
  defaults <- magcomp_prior_defaults()
  mag <- magcomp_prior_entries(magnitude, defaults$magnitude, "priors$magnitude")
  comp <- magcomp_prior_entries(composition, defaults$composition, "priors$composition")
  for(nm in names(mag$terms))
    mag$terms[[nm]] <- magcomp_prior_entries(mag$terms[[nm]],
      list(variance = NULL, covariance = list()), paste0("priors$magnitude$terms$", nm))
  for(b in names(comp$blocks)){
    block_defaults <- comp[setdiff(names(comp), "blocks")]
    block_defaults$terms <- list()
    block <- magcomp_prior_entries(comp$blocks[[b]], block_defaults,
                                   paste0("priors$composition$blocks$", b))
    for(t in names(block$terms))
      block$terms[[t]] <- magcomp_prior_entries(block$terms[[t]],
        list(loading_variance = NULL, covariance = list()),
        paste0("priors$composition$blocks$", b, "$terms$", t))
    comp$blocks[[b]] <- block
  }
  structure(list(magnitude = mag, composition = comp), class = "magCompPriors")
}

print.magCompPriors <- function(x, ...)
{
  cat("magComp priors (Gaussian variances, IG shape/scale, Gamma size shape/rate)\n")
  print(unclass(x))
  invisible(x)
}

magcomp_smooth_method <- function(x, name)
{
  choices <- c("none", "rw1", "rw2")
  x <- as.character(x)
  if(length(x) != 1L || is.na(x) || !x %in% choices)
    stop("error: ", name, " must be one of none, rw1, or rw2")
  x
}

magcomp_smooth_spec <- function(intercept = "none", slope = "none",
                                loading = "none", ordered = TRUE)
{
  intercept <- magcomp_smooth_method(intercept, "intercept")
  slope <- magcomp_smooth_method(slope, "slope")
  loading <- magcomp_smooth_method(loading, "loading")
  if(!isTRUE(ordered) &&
     (intercept != "none" || slope != "none" || loading != "none")){
    stop("error: smoothing is only allowed for ordered categories")
  }
  list(intercept = intercept,
       slope = slope,
       loading = loading,
       ordered = isTRUE(ordered))
}

magcomp_smooth_block_spec <- function(x, name)
{
  if(is.character(x) && length(x) == 1L){
    x <- list(intercept = x, slope = x, loading = "none")
  }
  if(!is.list(x))
    stop("error: smoothing entry '", name, "' must be a list")
  unknown <- setdiff(names(x), c("intercept", "slope", "loading", "ordered"))
  if(length(unknown) > 0L){
    stop("error: unknown smoothing control(s) for '", name, "': ",
         paste(unknown, collapse = ", "))
  }
  provided <- intersect(c("intercept", "slope", "loading", "ordered"), names(x))
  out <- list(intercept = NULL,
              slope = NULL,
              loading = NULL,
              ordered = x$ordered %||% TRUE)
  for(nm in c("intercept", "slope", "loading")){
    if(!is.null(x[[nm]]))
      out[[nm]] <- magcomp_smooth_method(x[[nm]], nm)
  }
  if(!isTRUE(out$ordered)){
    values <- unlist(out[c("intercept", "slope", "loading")], use.names = FALSE)
    if(length(values) > 0L && any(values != "none"))
      stop("error: smoothing is only allowed for ordered categories")
  }
  attr(out, "provided") <- provided
  out
}

magCompSmooth <- function(intercept = c("none", "rw1", "rw2"),
                          slope = c("none", "rw1", "rw2"),
                          loading = c("none", "rw1", "rw2"),
                          ordered = TRUE,
                          ...)
{
  intercept <- match.arg(intercept)
  slope <- match.arg(slope)
  loading <- match.arg(loading)
  global <- magcomp_smooth_spec(intercept, slope, loading, ordered)

  dots <- list(...)
  if(length(dots) > 0L){
    if(is.null(names(dots)) || any(!nzchar(names(dots))))
      stop("error: named smoothing entries are required")
    if(any(duplicated(names(dots))))
      stop("error: duplicated smoothing entry names are not allowed")
    dots <- mapply(magcomp_smooth_block_spec, dots, names(dots),
                   SIMPLIFY = FALSE, USE.NAMES = TRUE)
  }

  out <- list(
    intercept = global$intercept,
    slope = global$slope,
    loading = global$loading,
    ordered = global$ordered,
    blocks = dots
  )
  class(out) <- "magCompSmooth"
  out
}

print.magCompSmooth <- function(x, ...)
{
  cat("magComp smoothing\n")
  cat("  ordered categories: ", x$ordered, "\n", sep = "")
  cat("  intercept:          ", x$intercept, "\n", sep = "")
  cat("  slope:              ", x$slope, "\n", sep = "")
  cat("  loading:            ", x$loading, "\n", sep = "")
  if(length(x$blocks) > 0L)
    cat("  block entries:      ", paste(names(x$blocks), collapse = ", "),
        "\n", sep = "")
  invisible(x)
}
