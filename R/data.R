magcomp_make_count_matrix <- function(unit_factor, category_factor,
                                      unit_levels, tree_count = NULL)
{
  if(is.null(tree_count)){
    out <- stats::xtabs(~ unit_factor + category_factor)
  }else{
    out <- stats::xtabs(tree_count ~ unit_factor + category_factor)
  }
  out <- as.matrix(out)
  if(any(!is.finite(out)) || any(rowSums(out) > .Machine$integer.max))
    stop("error: counts per sampling unit exceed the supported integer limit")
  storage.mode(out) <- "integer"
  rownames(out) <- unit_levels
  out
}

magcomp_one_allocation_block <- function(name, counts, unit_levels,
                                         ordered = FALSE)
{
  out <- list(
    name = name,
    type = "marginal",
    response = name,
    parent = NULL,
    child = name,
    ordered = isTRUE(ordered),
    parent_levels = NULL,
    levels = colnames(counts),
    coef_stick_count = ncol(counts) - 1L,
    stick_names = colnames(counts)[-ncol(counts)],
    parent_specific_coef = FALSE,
    counts = counts,
    trials = rowSums(counts),
    unit_index = seq_along(unit_levels),
    parent_index = NULL,
    row_levels = unit_levels
  )
  class(out) <- "magCompAllocationBlock"
  out
}

magcomp_build_allocation_blocks <- function(unit_factor, unit_levels,
                                            category_factors,
                                            category_names,
                                            category_ordered,
                                            tree_count = NULL)
{
  if(length(category_factors) == 0L)
    return(list())

  counts_a <- magcomp_make_count_matrix(unit_factor, category_factors[[1L]],
                                        unit_levels, tree_count)
  if(length(category_factors) == 1L){
    out <- list(magcomp_one_allocation_block(category_names[1L],
                                             counts_a, unit_levels,
                                             category_ordered[1L]))
    names(out) <- category_names[1L]
    return(out)
  }

  parent <- category_factors[[1L]]
  child <- category_factors[[2L]]
  parent_levels <- levels(parent)
  child_levels <- levels(child)
  n_unit <- length(unit_levels)
  n_parent <- length(parent_levels)
  n_child <- length(child_levels)

  tab_data <- data.frame(unit_factor = unit_factor,
                         parent = parent,
                         child = child)
  if(is.null(tree_count)){
    tab <- stats::xtabs(~ unit_factor + parent + child, data = tab_data)
  }else{
    tab_data$tree_count <- tree_count
    tab <- stats::xtabs(tree_count ~ unit_factor + parent + child,
                        data = tab_data)
  }
  child_counts <- matrix(0L, nrow = n_unit * n_parent, ncol = n_child)
  child_row_levels <- character(n_unit * n_parent)
  child_unit_index <- integer(n_unit * n_parent)
  child_parent_index <- integer(n_unit * n_parent)

  for(a in seq_len(n_parent)){
    row_id <- (a - 1L) * n_unit + seq_len(n_unit)
    child_counts[row_id,] <- tab[,a,]
    child_row_levels[row_id] <- paste(unit_levels, parent_levels[a],
                                      sep = ":")
    child_unit_index[row_id] <- seq_len(n_unit)
    child_parent_index[row_id] <- a
  }

  storage.mode(child_counts) <- "integer"
  colnames(child_counts) <- child_levels
  rownames(child_counts) <- child_row_levels

  block_a <- magcomp_one_allocation_block(category_names[1L], counts_a,
                                          unit_levels, category_ordered[1L])
  block_b <- list(
    name = paste0(category_names[2L], "|", category_names[1L]),
    type = "conditional",
    response = paste0(category_names[2L], "|", category_names[1L]),
    parent = category_names[1L],
    child = category_names[2L],
    ordered = isTRUE(category_ordered[2L]),
    parent_levels = parent_levels,
    levels = child_levels,
    coef_stick_count = n_parent * (n_child - 1L),
    stick_names = as.vector(t(outer(parent_levels,
                                    child_levels[-n_child],
                                    paste,
                                    sep = ":"))),
    parent_specific_coef = TRUE,
    counts = child_counts,
    trials = rowSums(child_counts),
    unit_index = child_unit_index,
    parent_index = child_parent_index,
    row_levels = child_row_levels
  )
  class(block_b) <- "magCompAllocationBlock"

  out <- list(block_a, block_b)
  names(out) <- c(category_names[1L], block_b$name)
  out
}

magCompData <- function(data, unit, category = NULL, categories = NULL,
                        coords = NULL, tree_factor = NULL,
                        tree_count = NULL, design = NULL,
                        attributes = NULL, ordered = NULL)
{
  if(!is.data.frame(data))
    stop("error: data must be a data frame")

  unit_name <- magcomp_name_vector(substitute(unit), "unit")
  category_name <- magcomp_name_vector(substitute(category), "category")
  categories_name <- magcomp_name_vector(substitute(categories),
                                         "categories")
  coord_names <- magcomp_name_vector(substitute(coords), "coords")
  tree_factor_name <- magcomp_name_vector(substitute(tree_factor),
                                          "tree_factor")
  tree_count_name <- magcomp_name_vector(substitute(tree_count),
                                         "tree_count")
  attribute_names <- magcomp_name_vector(substitute(attributes), "attributes")

  if(length(unit_name) != 1L)
    stop("error: unit must identify one column")
  if(!unit_name %in% names(data))
    stop("error: unit column not found: ", unit_name)
  if(length(category_name) > 0L && length(categories_name) > 0L){
    stop("error: use category or categories, not both")
  }
  if(length(category_name) > 1L)
    stop("error: category must identify one column; use categories for ",
         "multiple category variables")
  category_names <- category_name
  if(length(categories_name) > 0L)
    category_names <- categories_name
  if(length(category_names) > 2L)
    stop("error: categories must identify one or two columns")

  unit <- data[[unit_name]]
  if(anyNA(unit))
    stop("error: unit contains missing values")
  unit_factor <- factor(unit, levels = unique(unit))
  unit_levels <- levels(unit_factor)
  n_unit <- length(unit_levels)
  unit_data <- magcomp_first_rows_by_unit(data, unit_factor)

  coords_matrix <- NULL
  if(length(coord_names) > 0L){
    magcomp_check_unit_columns(data, unit, coord_names)
    coords_df <- magcomp_get_columns(unit_data, coord_names, "coords")
    coords_matrix <- as.matrix(coords_df)
    storage.mode(coords_matrix) <- "double"
    rownames(coords_matrix) <- unit_levels
  }

  count_weight <- NULL
  if(length(tree_count_name) > 0L){
    if(length(tree_count_name) != 1L)
      stop("error: tree_count must identify one column")
    if(!tree_count_name %in% names(data))
      stop("error: tree_count column not found: ", tree_count_name)
    count_weight <- as.numeric(data[[tree_count_name]])
    if(anyNA(count_weight) || any(!is.finite(count_weight)) ||
       any(count_weight < 0)){
      stop("error: tree_count must contain non-negative finite values")
    }
    if(any(abs(count_weight - round(count_weight)) >
           sqrt(.Machine$double.eps))){
      stop("error: tree_count must contain integer counts")
    }
    if(any(count_weight > .Machine$integer.max))
      stop("error: tree_count exceeds the supported integer limit")
    count_weight <- as.integer(round(count_weight))
  }

  category_factor <- NULL
  category_factors <- list()
  category_levels <- list()
  category_ordered <- logical(0L)
  cell_category <- NULL
  count_matrix <- NULL
  allocation_blocks <- list()
  stick_names <- character(0L)
  if(length(category_names) > 0L){
    missing_categories <- setdiff(category_names, names(data))
    if(length(missing_categories) > 0L){
      stop("error: category column(s) not found: ",
           paste(missing_categories, collapse = ", "))
    }
    for(k in seq_along(category_names)){
      raw_category <- data[[category_names[k]]]
      f <- if(is.factor(raw_category))
        factor(raw_category, levels = levels(raw_category),
               ordered = is.ordered(raw_category)) else factor(raw_category)
      if(anyNA(f))
        stop("error: category column contains missing values: ",
             category_names[k])
      if(length(levels(f)) < 2L)
        stop("error: category column must contain at least two levels: ",
             category_names[k])
      category_factors[[k]] <- f
      category_levels[[k]] <- levels(f)
      category_ordered[k] <- is.ordered(f)
    }
    names(category_factors) <- category_names
    names(category_levels) <- category_names
    names(category_ordered) <- category_names

    if(!is.null(ordered)){
      ordered <- stats::setNames(as.logical(ordered), names(ordered))
      if(anyNA(ordered))
        stop("error: ordered must contain TRUE/FALSE values")
      if(length(ordered) == 1L && length(category_names) == 1L &&
         is.null(names(ordered))){
        names(ordered) <- category_names
      }
      if(is.null(names(ordered)) || any(!nzchar(names(ordered)))){
        if(length(ordered) != length(category_names))
          stop("error: unnamed ordered must have one value per category")
        names(ordered) <- category_names
      }
      unknown_ordered <- setdiff(names(ordered), category_names)
      if(length(unknown_ordered) > 0L){
        stop("error: ordered names are not category columns: ",
             paste(unknown_ordered, collapse = ", "))
      }
      category_ordered[names(ordered)] <- ordered
    }

    for(nm in category_names[category_ordered]){
      if(!is.factor(data[[nm]]) && !is.numeric(data[[nm]]))
        stop("error: ordered category '", nm,
             "' must be a factor with explicitly ordered levels or numeric; ",
             "use factor(..., levels = ...) to specify the category order")
    }

    category_factor <- category_factors[[1L]]
    if(length(category_factors) == 1L){
      cell_category <- category_factor
    }else{
      cell_category <- interaction(category_factors[[1L]],
                                   category_factors[[2L]],
                                   drop = FALSE,
                                   sep = ":",
                                   lex.order = TRUE)
    }
    count_matrix <- magcomp_make_count_matrix(unit_factor, cell_category,
                                             unit_levels, count_weight)
    allocation_blocks <- magcomp_build_allocation_blocks(
      unit_factor = unit_factor,
      unit_levels = unit_levels,
      category_factors = category_factors,
      category_names = category_names,
      category_ordered = category_ordered,
      tree_count = count_weight
    )
    for(block in allocation_blocks){
      stick_names <- c(stick_names, block$stick_names)
    }
  }

  # Without categories the response is a single total per sampling unit.
  if(is.null(count_matrix))
    count_matrix <- magcomp_make_count_matrix(unit_factor,
      factor(rep("total", nrow(data)), levels = "total"), unit_levels, count_weight)

  tree_factor <- NULL
  if(length(tree_factor_name) > 0L){
    if(length(tree_factor_name) != 1L)
      stop("error: tree_factor must identify one column")
    if(!tree_factor_name %in% names(data))
      stop("error: tree_factor column not found: ", tree_factor_name)
    tree_factor <- as.numeric(data[[tree_factor_name]])
    if(anyNA(tree_factor) || any(!is.finite(tree_factor)) ||
       any(tree_factor < 0))
      stop("error: tree_factor must contain non-negative finite values")
  }

  expanded <- list()
  if(length(attribute_names) > 0L){
    if(is.null(tree_factor))
      stop("error: attributes require tree_factor")
    missing_attributes <- setdiff(attribute_names, names(data))
    if(length(missing_attributes) > 0L)
      stop("error: attribute column(s) not found: ",
           paste(missing_attributes, collapse = ", "))
    if(is.null(cell_category))
      stop("error: attributes currently require category or categories")

    for(a in attribute_names){
      contribution <- tree_factor * as.numeric(data[[a]])
      if(anyNA(contribution) || any(!is.finite(contribution)))
        stop("error: expanded attribute ", a, " contains invalid values")
      tab <- stats::xtabs(contribution ~ unit_factor + cell_category)
      expanded[[a]] <- as.matrix(tab)
      rownames(expanded[[a]]) <- unit_levels
    }
  }

  if(is.null(design))
    design <- list(type = "unknown")
  if(!inherits(design, "magCompDesign") && !is.list(design))
    stop("error: design must be created by fixedArea() or pointSampling()")

  area <- NULL
  baf <- NULL
  if(identical(design$type, "fixed_area")){
    area <- magcomp_design_value(design, "area", unit_data, n_unit,
                                 data = data, unit = unit_factor)
    if(!is.null(area))
      magcomp_positive_finite(area, "area")
  }
  if(identical(design$type, "point_sampling"))
    baf <- magcomp_design_value(design, "baf", unit_data, n_unit,
                                data = data, unit = unit_factor)

  out <- list(
    data = data,
    unit_data = unit_data,
    unit_name = unit_name,
    unit = unit_factor,
    unit_levels = unit_levels,
    n_unit = n_unit,
    category_name = category_name,
    category_names = category_names,
    category = category_factor,
    category_factors = category_factors,
    category_levels = category_levels,
    category_ordered = category_ordered,
    n_category = length(category_names),
    cell_category = cell_category,
    coord_names = coord_names,
    coords = coords_matrix,
    tree_factor_name = tree_factor_name,
    tree_factor = tree_factor,
    tree_count_name = tree_count_name,
    tree_count = count_weight,
    design = design,
    area = area,
    baf = baf,
    counts = count_matrix,
    allocation_blocks = allocation_blocks,
    stick_names = stick_names,
    expanded = expanded,
    attribute_names = attribute_names
  )

  class(out) <- "magCompData"
  out
}

print.magCompData <- function(x, ...)
{
  cat("magComp data\n")
  cat("  sampling units: ", x$n_unit, "\n", sep = "")
  cat("  unit column:    ", x$unit_name, "\n", sep = "")
  cat("  design:         ", x$design$type %||% "unknown", "\n", sep = "")
  if(length(x$category_names) > 0L)
    cat("  category vars:  ", paste(x$category_names, collapse = ", "),
        "\n", sep = "")
  if(!is.null(x$counts))
    cat("  cells:          ", ncol(x$counts), "\n", sep = "")
  if(!is.null(x$coords))
    cat("  coordinates:    ", paste(x$coord_names, collapse = ", "), "\n", sep = "")
  if(length(x$tree_count_name) > 0L)
    cat("  tree counts:    ", x$tree_count_name, "\n", sep = "")
  if(length(x$attribute_names) > 0L)
    cat("  attributes:     ", paste(x$attribute_names, collapse = ", "), "\n", sep = "")
  invisible(x)
}
