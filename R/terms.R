iid <- function(...)
{
  stop("error: iid() is only used inside magComp formulas", call. = FALSE)
}

nngp <- function(...)
{
  stop("error: nngp() is only used inside magComp formulas", call. = FALSE)
}

gp <- function(...)
{
  stop("error: gp() is only used inside magComp formulas", call. = FALSE)
}

factor_iid <- function(...)
{
  stop("error: factor_iid() is only used inside magComp composition formulas",
       call. = FALSE)
}

factor_nngp <- function(...)
{
  stop("error: factor_nngp() is only used inside magComp composition formulas",
       call. = FALSE)
}

factor_gp <- function(...)
{
  stop("error: factor_gp() is only used inside magComp composition formulas",
       call. = FALSE)
}

magcomp_composition_block_name <- function(response, parent = NULL)
{
  if(is.null(parent))
    return(response)

  paste0(response, "|", parent)
}

magCompComposition <- function(...)
{
  formulas <- list(...)
  if(length(formulas) == 1L && is.list(formulas[[1L]]) &&
     !inherits(formulas[[1L]], "formula")){
    formulas <- formulas[[1L]]
  }
  if(length(formulas) < 1L)
    stop("error: magCompComposition() requires at least one formula")

  out <- vector("list", length(formulas))
  for(i in seq_along(formulas)){
    f <- formulas[[i]]
    if(!inherits(f, "formula") || length(f) != 3L){
      stop("error: magCompComposition() entries must be two-sided ",
           "formulas")
    }

    lhs <- f[[2L]]
    parent <- NULL
    if(is.call(lhs) && identical(as.character(lhs[[1L]]), "|")){
      if(length(lhs) != 3L)
        stop("error: malformed conditional composition formula")
      response <- magcomp_deparse1(lhs[[2L]])
      parent <- magcomp_deparse1(lhs[[3L]])
    }else{
      response <- magcomp_deparse1(lhs)
    }
    if(!nzchar(response))
      stop("error: composition response name is empty")
    if(!is.null(parent) && !nzchar(parent))
      stop("error: composition parent name is empty")

    rhs <- magcomp_deparse1(f[[3L]])
    one_sided <- stats::as.formula(paste("~", rhs),
                                   env = environment(f))
    name <- magcomp_composition_block_name(response, parent)
    out[[i]] <- list(
      formula = f,
      rhs_formula = one_sided,
      response = response,
      parent = parent,
      name = name
    )
  }

  names(out) <- vapply(out, `[[`, character(1L), "name")
  if(any(duplicated(names(out))))
    stop("error: composition formulas must identify unique blocks")

  class(out) <- "magCompComposition"
  out
}

print.magCompComposition <- function(x, ...)
{
  cat("magComp composition specification\n")
  for(i in seq_along(x)){
    cat("  ", x[[i]]$name, ": ",
        paste(deparse(x[[i]]$rhs_formula), collapse = " "),
        "\n", sep = "")
  }
  invisible(x)
}

magcomp_is_process_term <- function(lbl)
{
  grepl("^(iid|nngp|gp|factor_iid|factor_nngp|factor_gp)\\(", lbl) ||
    grepl(":(iid|nngp|gp|factor_iid|factor_nngp|factor_gp)\\(", lbl)
}

magcomp_split_top_level <- function(x, split = ",")
{
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L
  out <- character(0)
  start <- 1L

  for(i in seq_along(chars)){
    ch <- chars[i]
    if(ch == "(") depth <- depth + 1L
    if(ch == ")") depth <- depth - 1L
    if(ch == split && depth == 0L){
      out <- c(out, trimws(substr(x, start, i - 1L)))
      start <- i + 1L
    }
  }

  out <- c(out, trimws(substr(x, start, nchar(x))))
  out[nzchar(out)]
}

magcomp_split_top_level_colon <- function(x)
{
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L

  for(i in seq_along(chars)){
    ch <- chars[i]
    if(ch == "(") depth <- depth + 1L
    if(ch == ")") depth <- depth - 1L
    if(ch == ":" && depth == 0L){
      return(c(trimws(substr(x, 1L, i - 1L)),
               trimws(substr(x, i + 1L, nchar(x)))))
    }
  }

  NULL
}

magcomp_parse_process_call <- function(lbl, env = parent.frame())
{
  colon_parts <- magcomp_split_top_level_colon(lbl)

  if(is.null(colon_parts)){
    proc_part <- trimws(lbl)
    coef_name <- NULL
  }else{
    left_is_process <- grepl("^[A-Za-z0-9_.]+\\(", colon_parts[1L])
    right_is_process <- grepl("^[A-Za-z0-9_.]+\\(", colon_parts[2L])

    if(left_is_process && !right_is_process){
      proc_part <- colon_parts[1L]
      coef_name <- trimws(colon_parts[2L])
    }else if(!left_is_process && right_is_process){
      proc_part <- colon_parts[2L]
      coef_name <- trimws(colon_parts[1L])
    }else{
      stop("error: unable to parse process interaction term ", lbl)
    }
  }

  m <- regexec("^([A-Za-z0-9_.]+)\\((.*)\\)$", proc_part)
  hit <- regmatches(proc_part, m)[[1L]]
  if(length(hit) != 3L)
    stop("error: malformed process term ", lbl)

  fun <- hit[2L]
  inside <- trimws(hit[3L])
  pieces <- magcomp_split_top_level(inside, split = ",")

  args <- character(0)
  params <- list()
  for(p in pieces){
    if(grepl("=", p, fixed = TRUE)){
      kv <- strsplit(p, "=", fixed = TRUE)[[1L]]
      if(length(kv) != 2L)
        stop("error: malformed parameter in ", lbl)
      key <- trimws(kv[1L])
      val <- trimws(kv[2L])
      parsed <- tryCatch(eval(parse(text = val), envir = env),
                         error = function(e) val)
      if(key %in% names(params))
        stop("error: duplicated term argument ", key)
      params[[key]] <- parsed
    }else{
      args <- c(args, trimws(p))
    }
  }

  allowed <- c("name", if(fun %in% c("nngp", "factor_nngp")) "m",
               if(fun %in% c("gp", "nngp", "factor_gp", "factor_nngp"))
                 c("cov_model", "st_scale"),
               if(startsWith(fun, "factor_")) c("K", "anchors", "factor_names"))
  magcomp_check_names(params, allowed, paste0(fun, "()"))
  if(startsWith(fun, "factor_")){
    K <- params$K %||% 1L
    if(!is.numeric(K) || length(K) != 1L || is.na(K) ||
       !is.finite(K) || K < 1 || K != as.integer(K))
      stop("error: ", fun, "() requires a positive integer K")
    params$K <- as.integer(K)
  }

  list(
    label = lbl,
    term_type = fun,
    args = args,
    params = params,
    coef_name = coef_name,
    coefficient = coef_name %||% "(Intercept)"
  )
}

magcomp_default_process_names <- function(specs)
{
  if(length(specs) == 0L)
    return(character(0))
  base <- vapply(specs, `[[`, character(1L), "term_type")
  defaults <- stats::ave(base, base, FUN = function(x) paste0(x, "_", seq_along(x)))
  out <- vapply(seq_along(specs), function(i){
    name <- specs[[i]]$params$name %||% defaults[i]
    if(!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))
      stop("error: term name must be a nonempty string")
    name
  }, character(1L))
  if(anyDuplicated(out)) stop("error: term names must be unique within a component/block")
  out
}

magcomp_check_spatial_locations <- function(coords, unit_levels, term)
{
  # Conditional allocation rows repeat a unit for each parent category.
  # Check locations once per unit, retaining all distinct sampling units.
  first <- !duplicated(unit_levels)
  coords <- coords[first,,drop = FALSE]
  unit_levels <- unit_levels[first]
  tied <- duplicated(coords) | duplicated(coords, fromLast = TRUE)
  if(any(tied))
    stop("error: spatial term '", term$name, "' (", term$term_type,
         "()) has identical coordinates for different sampling units: ",
         paste(encodeString(unit_levels[tied], quote = "'"), collapse = ", "),
         ". Spatial fitting requires distinct coordinates for different ",
         "sampling units.", call. = FALSE)
  invisible(NULL)
}

magcomp_prepare_process_term <- function(term, data, unit_levels = rownames(data))
{
  if(term$term_type == "iid"){
    if(length(term$args) != 1L)
      stop("error: iid() requires one grouping variable")

    group_name <- term$args[1L]
    term_names <- group_name
    if(!is.null(term$coef_name))
      term_names <- unique(c(term$coef_name, term_names))
    missing_names <- setdiff(term_names, names(data))
    if(length(missing_names) > 0L){
      stop("error: iid() column(s) not found: ",
           paste(missing_names, collapse = ", "))
    }

    group <- factor(data[[group_name]], levels = unique(data[[group_name]]))
    if(anyNA(group))
      stop("error: iid() grouping column contains missing values")

    if(!is.null(term$coef_name)){
      covariate <- as.numeric(data[[term$coef_name]])
      if(anyNA(covariate) || any(!is.finite(covariate)))
        stop("error: iid() covariate must be finite")
      term$covariate_name <- term$coef_name
      term$covariate <- covariate
    }

    term$group_name <- group_name
    term$group <- as.integer(group)
    term$group_levels <- levels(group)
    term$n_group <- length(term$group_levels)
  }

  if(term$term_type %in% c("nngp", "gp")){
    label <- term$term_type
    if(length(term$args) < 2L)
      stop("error: ", label, "() requires at least two coordinate variables")

    coord_names <- term$args
    term_names <- coord_names
    if(!is.null(term$coef_name))
      term_names <- unique(c(term$coef_name, term_names))
    missing_names <- setdiff(term_names, names(data))
    if(length(missing_names) > 0L){
      stop("error: ", label, "() column(s) not found: ",
           paste(missing_names, collapse = ", "))
    }

    if(!is.null(term$coef_name)){
      covariate <- as.numeric(data[[term$coef_name]])
      if(anyNA(covariate) || any(!is.finite(covariate)))
        stop("error: ", label, "() covariate must be finite")
      term$covariate_name <- term$coef_name
      term$covariate <- covariate
    }

    coords <- as.matrix(data[,coord_names,drop = FALSE])
    storage.mode(coords) <- "double"
    if(anyNA(coords) || any(!is.finite(coords)))
      stop("error: ", label, "() coordinates must be finite")
    if(nrow(coords) < 2L)
      stop("error: ", label, "() requires at least two sampling units")
    magcomp_check_spatial_locations(coords, unit_levels, term)

    term$coord_names <- coord_names
    term$coords <- coords
    term$location_levels <- rownames(data)
    if(identical(term$term_type, "nngp")){
      term$m <- as.integer(magcomp_param(term$params, "m") %||% 15L)
      if(length(term$m) != 1L || is.na(term$m) || term$m < 1L)
        stop("error: nngp() m must be a positive integer")
    }

    term <- magcomp_prepare_covariance(term, coords, label)
    term$n_group <- nrow(coords)
    if(identical(term$term_type, "nngp")){
      nngp <- magcomp_prepare_nngp_neighbors(magcomp_search_coords(term, coords), term$m)
      term$order <- nngp$order
      term$neighbor_count <- nngp$neighbor_count
      term$neighbor_index <- nngp$neighbor_index
    }
  }

  if(term$term_type == "factor_iid"){
    if(length(term$args) != 1L)
      stop("error: factor_iid() requires one grouping variable")

    group_name <- term$args[1L]
    term_names <- group_name
    if(!is.null(term$coef_name))
      term_names <- unique(c(term$coef_name, term_names))
    missing_names <- setdiff(term_names, names(data))
    if(length(missing_names) > 0L){
      stop("error: factor_iid() column(s) not found: ",
           paste(missing_names, collapse = ", "))
    }

    group <- factor(data[[group_name]], levels = unique(data[[group_name]]))
    if(anyNA(group))
      stop("error: factor_iid() grouping column contains missing values")

    if(!is.null(term$coef_name)){
      covariate <- as.numeric(data[[term$coef_name]])
      if(anyNA(covariate) || any(!is.finite(covariate)))
        stop("error: factor_iid() covariate must be finite")
      term$covariate_name <- term$coef_name
      term$covariate <- covariate
    }

    term$group_name <- group_name
    term$group <- as.integer(group)
    term$group_levels <- levels(group)
    term$n_group <- length(term$group_levels)
    term$K <- magcomp_param(term$params, "K")
  }

  if(term$term_type %in% c("factor_nngp", "factor_gp")){
    label <- term$term_type
    if(length(term$args) < 2L)
      stop("error: ", label,
           "() requires at least two coordinate variables")

    coord_names <- term$args
    term_names <- coord_names
    if(!is.null(term$coef_name))
      term_names <- unique(c(term$coef_name, term_names))
    missing_names <- setdiff(term_names, names(data))
    if(length(missing_names) > 0L){
      stop("error: ", label, "() column(s) not found: ",
           paste(missing_names, collapse = ", "))
    }

    if(!is.null(term$coef_name)){
      covariate <- as.numeric(data[[term$coef_name]])
      if(anyNA(covariate) || any(!is.finite(covariate)))
        stop("error: ", label, "() covariate must be finite")
      term$covariate_name <- term$coef_name
      term$covariate <- covariate
    }

    coords <- as.matrix(data[,coord_names,drop = FALSE])
    storage.mode(coords) <- "double"
    if(anyNA(coords) || any(!is.finite(coords)))
      stop("error: ", label, "() coordinates must be finite")
    if(nrow(coords) < 2L)
      stop("error: ", label, "() requires at least two sampling units")
    magcomp_check_spatial_locations(coords, unit_levels, term)

    term$coord_names <- coord_names
    term$coords <- coords
    term$location_levels <- rownames(data)
    term$K <- magcomp_param(term$params, "K")
    if(identical(term$term_type, "factor_nngp")){
      term$m <- as.integer(magcomp_param(term$params, "m") %||% 15L)
      if(length(term$m) != 1L || is.na(term$m) || term$m < 1L)
        stop("error: factor_nngp() m must be a positive integer")
    }

    term <- magcomp_prepare_covariance(term, coords, label)
  }

  term
}

magcomp_prepare_nngp_neighbors <- function(coords, m)
{
  coords <- as.matrix(coords)
  n <- nrow(coords)
  m <- min(as.integer(m), max(0L, n - 1L))
  coord_data <- as.data.frame(coords)
  order_index <- do.call(order, coord_data)
  neighbor_count <- integer(n)
  neighbor_index <- matrix(0L, nrow = n, ncol = m)

  if(m > 0L){
    for(pos in seq_len(n)){
      if(pos == 1L)
        next

      current <- order_index[pos]
      previous <- order_index[seq_len(pos - 1L)]
      d <- sqrt(rowSums(sweep(coords[previous,,drop = FALSE], 2L,
                              coords[current,], "-")^2))
      keep <- previous[order(d)[seq_len(min(m, length(previous)))]]
      neighbor_count[pos] <- length(keep)
      neighbor_index[pos,seq_along(keep)] <- keep
    }
  }

  list(order = as.integer(order_index),
       neighbor_count = as.integer(neighbor_count),
       neighbor_index = neighbor_index)
}

magcomp_default_factor_anchor_index <- function(y, K, parent_index = NULL,
                                                coef_stick_count = NULL)
{
  Jstick <- ncol(y) - 1L
  if(is.null(parent_index)){
    score <- rep(-Inf, Jstick)
  }else{
    parent_index <- as.integer(parent_index)
    coef_stick_count <- coef_stick_count %||%
      (max(parent_index) * Jstick)
    score <- rep(-Inf, coef_stick_count)
  }
  remaining <- rowSums(y)

  for(j in seq_len(Jstick)){
    if(is.null(parent_index)){
      keep <- remaining > 0
      if(sum(keep) > 1L){
        eta <- stats::qlogis((y[keep,j] + 0.5) / (remaining[keep] + 1))
        w <- remaining[keep]
        mu <- stats::weighted.mean(eta, w)
        value <- sum(w * (eta - mu)^2) / sum(w)
        if(is.finite(value))
          score[j] <- value
      }
    }else{
      n_parent <- coef_stick_count / Jstick
      for(a in seq_len(n_parent)){
        keep <- remaining > 0 & parent_index == a
        if(sum(keep) > 1L){
          eta <- stats::qlogis((y[keep,j] + 0.5) / (remaining[keep] + 1))
          w <- remaining[keep]
          mu <- stats::weighted.mean(eta, w)
          value <- sum(w * (eta - mu)^2) / sum(w)
          if(is.finite(value))
            score[(a - 1L) * Jstick + j] <- value
        }
      }
    }
    remaining <- remaining - y[,j]
  }

  selected <- order(score, decreasing = TRUE)
  selected <- selected[is.finite(score[selected])]
  fallback <- setdiff(seq_along(score), selected)
  c(selected, fallback)[seq_len(K)]
}

magcomp_finalize_factor_anchors <- function(term, data, block = NULL,
                                            counts = data$counts)
{
  if(!is.null(block))
    counts <- block$counts
  categories <- colnames(counts)
  stick_categories <- if(!is.null(block) && !is.null(block$stick_names))
    block$stick_names else categories[-length(categories)]
  K <- term$K
  anchors <- magcomp_param(term$params, "anchors") %||%
    magcomp_param(term$params, "anchor") %||% NULL
  label <- term$term_type

  if(K > length(stick_categories))
    stop("error: ", label,
         "() K cannot exceed the number of coefficient-stick rows")

  if(is.null(anchors)){
    anchor_index <- magcomp_default_factor_anchor_index(
      counts,
      K,
      parent_index = if(!is.null(block)) block$parent_index else NULL,
      coef_stick_count = if(!is.null(block)) block$coef_stick_count else NULL
    )
  }else if(is.numeric(anchors)){
    if(length(anchors) != K || anyNA(anchors))
      stop("error: ", label, "() anchors must have length K")
    anchor_index <- as.integer(anchors)
  }else if(is.character(anchors)){
    if(length(anchors) != K || anyNA(anchors))
      stop("error: ", label, "() anchors must have length K")
    anchor_index <- match(anchors, stick_categories)
    if(anyNA(anchor_index)){
      stop("error: ", label, "() anchors must match continuation-ratio ",
           "category names")
    }
  }else{
    stop("error: ", label, "() anchors must be numeric or character")
  }

  if(any(anchor_index < 1L) || any(anchor_index > length(stick_categories)))
    stop("error: ", label,
         "() anchor indices must identify coefficient-stick rows")
  if(any(duplicated(anchor_index)))
    stop("error: ", label, "() anchors must be unique")

  factor_names <- magcomp_param(term$params, "factor_names") %||% NULL
  if(is.null(factor_names)){
    factor_names <- stick_categories[anchor_index]
  }else{
    if(!is.character(factor_names) || length(factor_names) != K ||
       anyNA(factor_names) || any(!nzchar(factor_names))){
      stop("error: ", label,
           "() factor_names must be a non-missing character vector of length K")
    }
    if(any(duplicated(factor_names)))
      stop("error: ", label, "() factor_names must be unique")
  }

  term$anchor_index <- anchor_index
  term$anchors <- stick_categories[anchor_index]
  term$factor_names <- factor_names
  if(term$term_type %in% c("factor_nngp", "factor_gp")){
    coords_for_process <- term$coords
    if(!is.null(block) && !is.null(block$unit_index)){
      first_row <- match(seq_len(data$n_unit), block$unit_index)
      if(anyNA(first_row)){
        stop("error: ", term$term_type,
             "() spatial factor is missing a sampling unit")
      }
      coords_for_process <- term$coords[first_row,,drop = FALSE]
      term$spatial_group <- as.integer(block$unit_index)
      term$location_levels <- data$unit_levels
      term$n_group <- data$n_unit
    }else{
      term$spatial_group <- seq_len(nrow(term$coords))
      term$location_levels <- rownames(term$coords)
      term$n_group <- nrow(term$coords)
    }
  }
  if(term$term_type == "factor_nngp"){
    nngp <- magcomp_prepare_nngp_neighbors(magcomp_search_coords(term, coords_for_process), term$m)
    term$order <- nngp$order
    term$neighbor_count <- nngp$neighbor_count
    term$neighbor_index <- nngp$neighbor_index
  }
  term
}

magcomp_finalize_component_process_terms <- function(component, data)
{
  if(!is.null(component$allocation_components)){
    for(b in seq_along(component$allocation_components)){
      block_name <- names(component$allocation_components)[b]
      block_counts <- data$allocation_blocks[[block_name]]$counts
      block_component <- component$allocation_components[[b]]
      for(i in seq_along(block_component$process_terms)){
        term <- block_component$process_terms[[i]]
        if(term$term_type %in% c("factor_iid", "factor_nngp", "factor_gp")){
          term <- magcomp_finalize_factor_anchors(
            term,
            data,
            block = data$allocation_blocks[[block_name]],
            counts = block_counts
          )
        }
        block_component$process_terms[[i]] <- term
      }
      component$allocation_components[[b]] <- block_component
    }
    if(length(component$allocation_components) == 1L){
      first <- component$allocation_components[[1L]]
      component$X <- first$X
      component$fixed_formula <- first$fixed_formula
      component$process_labels <- first$process_labels
      component$process_terms <- first$process_terms
    }
    return(component)
  }

  if(length(component$process_terms) == 0L)
    return(component)
  if(is.null(data$counts))
    return(component)

  for(i in seq_along(component$process_terms)){
    term <- component$process_terms[[i]]
    if(term$term_type %in% c("factor_iid", "factor_nngp", "factor_gp"))
      term <- magcomp_finalize_factor_anchors(term, data)
    component$process_terms[[i]] <- term
  }

  component
}

magcomp_build_component <- function(formula, data, component,
                                     unit_levels = rownames(data))
{
  formula <- magcomp_require_formula(formula, component)

  tt <- stats::terms(formula, data = data, keep.order = TRUE)
  term_labels <- attr(tt, "term.labels")
  process_labels <- term_labels[vapply(term_labels, magcomp_is_process_term,
                                       logical(1L))]

  process_terms <- lapply(process_labels, magcomp_parse_process_call,
                          env = environment(formula))
  process_names <- magcomp_default_process_names(process_terms)
  for(i in seq_along(process_terms)){
    process_terms[[i]]$name <- process_names[i]
    process_terms[[i]]$component <- component
    process_terms[[i]]$term_id <- paste0(component, "_process_", i)
    process_terms[[i]] <- magcomp_prepare_process_term(process_terms[[i]],
                                                       data, unit_levels)
  }

  fixed_labels <- setdiff(term_labels, process_labels)
  has_intercept <- identical(attr(tt, "intercept"), 1L)

  if(length(fixed_labels) == 0L){
    fixed_rhs <- if(has_intercept) "1" else "0"
  }else{
    fixed_rhs <- paste(fixed_labels, collapse = " + ")
    if(!has_intercept)
      fixed_rhs <- paste("0 +", fixed_rhs)
  }
  fixed_formula <- stats::as.formula(paste("~", fixed_rhs),
                                     env = environment(formula))

  mf <- stats::model.frame(fixed_formula, data = data,
                           na.action = stats::na.fail)
  X <- stats::model.matrix(attr(mf, "terms"), data = mf)

  list(
    formula = formula,
    fixed_formula = fixed_formula,
    X = X,
    process_labels = process_labels,
    process_terms = process_terms
  )
}

magcomp_allocation_block_data <- function(data, block)
{
  if(identical(block$type, "marginal"))
    return(data$unit_data)

  if(!identical(block$type, "conditional"))
    stop("error: unknown allocation block type: ", block$type)

  out <- data$unit_data[block$unit_index,,drop = FALSE]
  out[[block$parent]] <- factor(block$parent_levels[block$parent_index],
                                levels = block$parent_levels)
  rownames(out) <- block$row_levels
  out
}

magcomp_build_composition_component <- function(composition, data)
{
  if(inherits(composition, "magCompComposition")){
    specs <- composition
    expected <- names(data$allocation_blocks)
    supplied <- names(specs)
    missing_blocks <- setdiff(expected, supplied)
    extra_blocks <- setdiff(supplied, expected)
    if(length(missing_blocks) > 0L){
      stop("error: missing composition formula(s) for allocation block(s): ",
           paste(missing_blocks, collapse = ", "))
    }
    if(length(extra_blocks) > 0L){
      stop("error: composition formula(s) do not match allocation block(s): ",
           paste(extra_blocks, collapse = ", "))
    }

    allocation_components <- vector("list", length(specs))
    names(allocation_components) <- supplied
    for(i in seq_along(specs)){
      block <- data$allocation_blocks[[supplied[i]]]
      block_data <- magcomp_allocation_block_data(data, block)
      block_component <- magcomp_build_component(
        formula = specs[[i]]$rhs_formula,
        data = block_data,
        component = "composition",
        unit_levels = data$unit_levels[block$unit_index]
      )
      block_component$formula <- specs[[i]]$formula
      block_component$allocation_block <- supplied[i]
      block_component$response <- specs[[i]]$response
      block_component$parent <- specs[[i]]$parent
      block_component$block_type <- block$type
      allocation_components[[i]] <- block_component
    }

    out <- list(
      formula = composition,
      fixed_formula = NULL,
      X = allocation_components[[1L]]$X,
      process_labels = allocation_components[[1L]]$process_labels,
      process_terms = allocation_components[[1L]]$process_terms,
      allocation_components = allocation_components,
      allocation_block_names = names(allocation_components),
      block_specific = TRUE
    )
    return(out)
  }

  component <- magcomp_build_component(
    formula = composition,
    data = data$unit_data,
    component = "composition"
  )
  if(length(data$allocation_blocks) > 0L){
    block_name <- names(data$allocation_blocks)[1L]
    block_component <- component
    block_component$allocation_block <- block_name
    block_component$response <- data$allocation_blocks[[1L]]$response
    block_component$parent <- data$allocation_blocks[[1L]]$parent
    block_component$block_type <- data$allocation_blocks[[1L]]$type
    component$allocation_components <- stats::setNames(list(block_component),
                                                       block_name)
    component$allocation_block_names <- block_name
    component$block_specific <- FALSE
  }

  component
}
