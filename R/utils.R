`%||%` <- function(x, y)
{
  if(is.null(x)) y else x
}

magcomp_param <- function(params, name)
{
  if(is.null(params) || is.null(names(params)) || !name %in% names(params))
    return(NULL)

  params[[name]]
}

magcomp_deparse1 <- function(x)
{
  paste(deparse(x), collapse = "")
}

magcomp_name_vector <- function(expr, arg)
{
  if(missing(expr))
    return(NULL)
  if(is.null(expr))
    return(NULL)

  if(is.character(expr))
    return(expr)

  if(is.name(expr))
    return(as.character(expr))

  if(is.call(expr) && identical(as.character(expr[[1L]]), "c")){
    values <- as.list(expr)[-1L]
    return(vapply(values, magcomp_deparse1, character(1L)))
  }

  stop("error: ", arg, " must be a variable name, character vector, or c(...)")
}

magcomp_get_columns <- function(data, names, arg)
{
  if(length(names) == 0L)
    return(NULL)

  missing_names <- setdiff(names, colnames(data))
  if(length(missing_names) > 0L)
    stop("error: ", arg, " column(s) not found: ",
         paste(missing_names, collapse = ", "))

  data[,names,drop = FALSE]
}

magcomp_first_rows_by_unit <- function(data, unit)
{
  keep <- !duplicated(unit)
  out <- as.data.frame(data[keep,,drop = FALSE])
  rownames(out) <- as.character(unit[keep])
  out
}

magcomp_check_unit_columns <- function(data, unit, columns)
{
  columns <- intersect(unique(columns), names(data))
  first <- match(unit, unit)
  for(nm in columns){
    value <- data[[nm]]
    reference <- if(is.null(dim(value))) value[first] else value[first,,drop = FALSE]
    same <- value == reference
    same[is.na(value) & is.na(reference)] <- TRUE
    if(anyNA(same) || any(!same))
      stop("error: unit-level column '", nm,
           "' must be constant within each sampling unit")
  }
  invisible(NULL)
}

magcomp_unit_constant_value <- function(value, unit, unit_levels, field)
{
  if(length(value) != length(unit)){
    stop("error: design field ", field,
         " must be scalar or have one value per sampling unit")
  }
  magcomp_positive_finite(value, field)

  unit_id <- as.integer(unit)
  out <- numeric(length(unit_levels))
  for(i in seq_along(unit_levels)){
    values <- unique(as.numeric(value[unit_id == i]))
    if(length(values) != 1L){
      stop("error: design field ", field,
           " must be constant within each sampling unit")
    }
    out[i] <- values
  }
  names(out) <- unit_levels

  out
}

magcomp_design_value <- function(design, field, unit_data, n_unit,
                                 data = NULL, unit = NULL)
{
  value <- design[[field]]
  name <- design[[paste0(field, "_name")]]
  unit_levels <- rownames(unit_data) %||% as.character(seq_len(n_unit))

  if(!is.null(name) && !is.null(data) && !is.null(unit) && name %in% names(data)){
    value <- as.numeric(data[[name]])
    return(magcomp_unit_constant_value(value, unit, unit_levels, field))
  }

  if(!is.null(name) && is.null(value)){
    if(!name %in% names(unit_data))
      stop("error: design field ", field, " column not found: ", name)
    value <- as.numeric(unit_data[[name]])
    if(length(value) != n_unit || anyNA(value))
      stop("error: design field ", field, " contains missing values")
    magcomp_positive_finite(value, field)
    return(value)
  }

  if(!is.null(value)){
    if(!is.null(dim(value))){
      stop("error: design field ", field,
           " must be scalar or have one value per sampling unit")
    }

    value_names <- names(value)
    value <- as.numeric(value)
    names(value) <- value_names
    magcomp_positive_finite(value, field)

    if(length(value) == 1L)
      return(rep(unname(value), n_unit))
    if(is.null(value_names)){
      stop("error: non-scalar design field ", field,
           " values must be named by sampling unit")
    }
    if(length(value) > 1L && !is.null(value_names)){
      if(all(unit_levels %in% value_names))
        return(unname(value[unit_levels]))
      stop("error: named design field ", field,
           " values must use sampling-unit names")
    }
  }

  NULL
}

magcomp_positive_finite <- function(x, field)
{
  if(anyNA(x) || any(!is.finite(x)) || any(x <= 0))
    stop("error: design field ", field,
         " must contain positive finite values")
}

magcomp_require_formula <- function(x, arg)
{
  if(missing(x) || is.null(x))
    x <- ~ 1
  if(!inherits(x, "formula"))
    stop("error: ", arg, " must be a formula")
  if(length(x) != 2L)
    stop("error: ", arg, " must be a one-sided formula")
  x
}
