fixedArea <- function(area = NULL)
{
  expr <- substitute(area)
  out <- list(type = "fixed_area")

  if(missing(area) || identical(expr, quote(NULL))){
    out$area <- NULL
    out$area_name <- NULL
  }else{
    value <- tryCatch(eval(expr, parent.frame()), error = function(e) NULL)
    out$area <- NULL
    out$area_name <- NULL
    if(is.name(expr))
      out$area_name <- as.character(expr)
    if(is.numeric(value)){
      if(!is.null(dim(value))){
        stop("error: fixedArea() area must be scalar or have one value ",
             "per sampling unit")
      }
      value_names <- names(value)
      out$area <- value
      storage.mode(out$area) <- "double"
      names(out$area) <- value_names
    }else if(is.character(value) && length(value) == 1L){
      out$area_name <- value
    }else{
      out$area_name <- out$area_name %||% magcomp_deparse1(expr)
    }
  }

  class(out) <- "magCompDesign"
  out
}

pointSampling <- function(baf = NULL)
{
  expr <- substitute(baf)
  out <- list(type = "point_sampling")

  if(missing(baf) || identical(expr, quote(NULL))){
    out$baf <- NULL
    out$baf_name <- NULL
  }else{
    value <- tryCatch(eval(expr, parent.frame()), error = function(e) NULL)
    out$baf <- NULL
    out$baf_name <- NULL
    if(is.name(expr))
      out$baf_name <- as.character(expr)
    if(is.numeric(value)){
      if(!is.null(dim(value))){
        stop("error: pointSampling() baf must be scalar or have one value ",
             "per sampling unit")
      }
      value_names <- names(value)
      out$baf <- as.numeric(value)
      names(out$baf) <- value_names
    }else if(is.character(value) && length(value) == 1L){
      out$baf_name <- value
    }else{
      out$baf_name <- out$baf_name %||% magcomp_deparse1(expr)
    }
  }

  class(out) <- "magCompDesign"
  out
}

fixedAreaStems <- function()
{
  out <- list(
    name = "fixed_area_stems",
    design = "fixed_area",
    magnitude = "stem_density",
    composition = "stem_count_proportion",
    likelihood = "negative_binomial_multinomial"
  )
  class(out) <- "magCompFamily"
  out
}

pointBasalArea <- function()
{
  out <- list(
    name = "point_basal_area",
    design = "point_sampling",
    magnitude = "basal_area",
    composition = "basal_area_proportion",
    likelihood = "negative_binomial_multinomial"
  )
  class(out) <- "magCompFamily"
  out
}

continuousGaussian <- function(link = c("log", "identity"),
                               composition = c("logratio", "gaussian_stick"))
{
  link <- match.arg(link)
  composition <- match.arg(composition)
  out <- list(
    name = "continuous_gaussian",
    design = "weighted_tree_list",
    magnitude = "continuous_total",
    composition = composition,
    likelihood = "gaussian",
    link = link
  )
  class(out) <- "magCompFamily"
  out
}

print.magCompFamily <- function(x, ...)
{
  cat("magComp family\n")
  cat("  name:        ", x$name, "\n", sep = "")
  cat("  magnitude:   ", x$magnitude, "\n", sep = "")
  cat("  composition: ", x$composition, "\n", sep = "")
  cat("  likelihood:  ", x$likelihood, "\n", sep = "")
  invisible(x)
}
