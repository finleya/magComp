test_that("component parser separates fixed effects and process terms", {
  unit_data <- data.frame(
    p25 = rnorm(5),
    p95 = rnorm(5),
    stratum = letters[c(1, 1, 2, 2, 3)],
    x = runif(5),
    y = runif(5)
  )

  comp <- magComp:::magcomp_build_component(
    ~ p25 + p95 + iid(stratum) + factor_nngp(x, y, K = 2, m = 10),
    unit_data,
    "composition"
  )

  expect_true("(Intercept)" %in% colnames(comp$X))
  expect_true("p25" %in% colnames(comp$X))
  expect_equal(length(comp$process_terms), 2L)
  expect_equal(comp$process_terms[[2L]]$term_type, "factor_nngp")
  expect_equal(comp$process_terms[[2L]]$params$K, 2L)

  mag <- magComp:::magcomp_build_component(
    ~ p95 + nngp(x, y, m = 3),
    unit_data,
    "magnitude"
  )
  expect_equal(mag$process_terms[[1L]]$term_type, "nngp")
  expect_equal(mag$process_terms[[1L]]$m, 3L)
  expect_equal(dim(mag$process_terms[[1L]]$neighbor_index), c(5L, 3L))
  expect_true(identical(mag$process_terms[[1L]]$cov_model, "exp"))

  mag_gp <- magComp:::magcomp_build_component(
    ~ p95 + gp(x, y),
    unit_data,
    "magnitude"
  )
  expect_equal(mag_gp$process_terms[[1L]]$term_type, "gp")
  expect_equal(dim(mag_gp$process_terms[[1L]]$coords), c(5L, 2L))
  expect_null(mag_gp$process_terms[[1L]]$neighbor_index)
  expect_true(identical(mag_gp$process_terms[[1L]]$cov_model, "exp"))

  mag_nngp_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:nngp(x, y, m = 2),
    unit_data,
    "magnitude"
  )
  term <- mag_nngp_slope$process_terms[[1L]]
  expect_equal(term$term_type, "nngp")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$coord_names, c("x", "y"))
  expect_equal(term$m, 2L)
  expect_equal(length(term$covariate), 5L)
  expect_equal(dim(term$neighbor_index), c(5L, 2L))

  mag_iid_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:iid(stratum),
    unit_data,
    "magnitude"
  )
  term <- mag_iid_slope$process_terms[[1L]]
  expect_equal(term$term_type, "iid")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$group_name, "stratum")
  expect_equal(length(term$covariate), 5L)

  mag_gp_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:gp(x, y),
    unit_data,
    "magnitude"
  )
  term <- mag_gp_slope$process_terms[[1L]]
  expect_equal(term$term_type, "gp")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$coord_names, c("x", "y"))
  expect_equal(length(term$covariate), 5L)
  expect_equal(dim(term$coords), c(5L, 2L))
  expect_null(term$neighbor_index)

  comp_iid_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:factor_iid(stratum, K = 2),
    unit_data,
    "composition"
  )
  term <- comp_iid_slope$process_terms[[1L]]
  expect_equal(term$term_type, "factor_iid")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$group_name, "stratum")
  expect_equal(length(term$covariate), 5L)

  comp_nngp_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:factor_nngp(x, y, K = 2, m = 2),
    unit_data,
    "composition"
  )
  term <- comp_nngp_slope$process_terms[[1L]]
  expect_equal(term$term_type, "factor_nngp")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$coord_names, c("x", "y"))
  expect_equal(term$m, 2L)
  expect_equal(length(term$covariate), 5L)

  comp_gp_slope <- magComp:::magcomp_build_component(
    ~ p95 + p25:factor_gp(x, y, K = 2),
    unit_data,
    "composition"
  )
  term <- comp_gp_slope$process_terms[[1L]]
  expect_equal(term$term_type, "factor_gp")
  expect_equal(term$coefficient, "p25")
  expect_equal(term$covariate_name, "p25")
  expect_equal(term$coord_names, c("x", "y"))
  expect_equal(length(term$covariate), 5L)
  expect_equal(dim(term$coords), c(5L, 2L))
  expect_null(term$neighbor_index)
})

test_that("spatial terms use nearest-neighbor default range supports", {
  unit_data <- data.frame(
    p25 = seq(-1, 1, length.out = 5),
    p95 = seq(10, 18, length.out = 5),
    x = 0:4,
    y = 0
  )

  mag_nngp <- magComp:::magcomp_build_component(
    ~ p95 + nngp(x, y, m = 2),
    unit_data,
    "magnitude"
  )
  expect_equal(as.numeric(magComp:::magcomp_covariance_defaults(mag_nngp$process_terms[[1L]])$priors$phi), c(3/8, 6))

  mag_gp <- magComp:::magcomp_build_component(
    ~ p95 + gp(x, y),
    unit_data,
    "magnitude"
  )
  expect_equal(as.numeric(magComp:::magcomp_covariance_defaults(mag_gp$process_terms[[1L]])$priors$phi), c(3/8, 6))

  comp_nngp <- magComp:::magcomp_build_component(
    ~ factor_nngp(x, y, K = 1, m = 2),
    unit_data,
    "composition"
  )
  expect_equal(as.numeric(magComp:::magcomp_covariance_defaults(comp_nngp$process_terms[[1L]])$priors$phi), c(3/8, 6))

  comp_gp <- magComp:::magcomp_build_component(
    ~ factor_gp(x, y, K = 1),
    unit_data,
    "composition"
  )
  expect_equal(as.numeric(magComp:::magcomp_covariance_defaults(comp_gp$process_terms[[1L]])$priors$phi), c(3/8, 6))
})

test_that("repeated parent coordinates do not change default decay priors", {
  coords <- cbind(x = c(0, 1, 3, 6, 10), y = 0)
  expected <- magComp:::magcomp_covariance_defaults(list(coords = coords, space_time = FALSE, cov_model = "exp"))
  repeated <- coords[rep(1:5, c(2, 4, 1, 2, 3)),,drop = FALSE]
  expect_equal(magComp:::magcomp_covariance_defaults(list(coords = repeated, space_time = FALSE, cov_model = "exp")), expected)
  rows <- expand.grid(plot = paste0("p", 1:5), species = c("a", "b", "c"),
                      dbh = c("small", "large"))
  rows$x <- coords[as.integer(rows$plot),1]
  rows$y <- 0
  rows$count <- 2L
  rows$area <- 1
  dat <- magCompData(rows, unit = plot, categories = c(species, dbh),
                     tree_count = count, design = fixedArea(area))
  for(kind in c("gp", "nngp")) {
    term <- paste0("factor_", kind, "(x, y, K = 1, name = 'shape')")
    model <- magComp(dat, composition = magCompComposition(
      as.formula(paste("species ~", term)),
      as.formula(paste("dbh | species ~", term))),
    control = magCompControl(run_sampler = FALSE))
    blocks <- magCompSettings(model)$priors$composition$blocks
    expect_equal(as.numeric(blocks$species$terms$shape$covariance$phi), c(0.15, 3))
    expect_equal(blocks[["dbh|species"]]$terms$shape$covariance$phi,
                 blocks$species$terms$shape$covariance$phi)
  }
})

test_that("spatial fitting rejects colocated units but allows repeated unit rows", {
  rows <- expand.grid(plot = paste0("p", 1:6), species = c("a", "b", "c"),
                      dbh = c("small", "large"))
  i <- as.integer(rows$plot)
  rows$x <- c(0, 0, 1, 2, 2, 3)[i]
  rows$y <- 0
  rows$unique_x <- i
  rows$count <- 2L
  rows$area <- 1
  make_data <- function(rows) magCompData(rows, unit = plot,
    categories = c(species, dbh), tree_count = count,
    coords = c(unique_x, y), design = fixedArea(area))
  dat <- make_data(rows)
  near <- rows
  near$x <- c(0, 1e-12, 1, 2, 2 + 1e-12, 3)[i]
  near <- make_data(near)
  same <- rows
  same$x <- 0
  same <- make_data(same)

  for(kind in c("gp", "nngp")) for(component in c("magnitude", "leading", "child")) {
    for(fixed in FALSE) {
      term <- paste0(if(component != "magnitude") "factor_" else "", kind,
        "(x, y, name = 'space'", if(component != "magnitude") ", K = 1" else "",
        ")")
      args <- if(component == "magnitude") list(magnitude = as.formula(paste("~", term))) else
        list(composition = magCompComposition(
          as.formula(paste("species ~", if(component == "leading") term else "1")),
          as.formula(paste("dbh | species ~", if(component == "child") term else "1"))))
      args$control <- magCompControl(run_sampler = FALSE)
      expect_error(do.call(magComp, c(list(data = dat), args)),
        "spatial term 'space'.*identical coordinates.*'p1', 'p2', 'p4', 'p5'")
      if(!fixed) {
        expect_error(do.call(magComp, c(list(data = same), args)),
          "spatial term 'space'.*identical coordinates.*'p1', 'p2', 'p3', 'p4', 'p5', 'p6'")
        # Exact comparisons: close coordinates must not be rounded into ties.
        # Each accepted unit also has repeated tree/parent-context rows.
        expect_no_error(do.call(magComp, c(list(data = near), args)))
      }
    }
  }
  # Check each term's coordinates, not the optional data-level coordinates.
  expect_error(magComp(dat, magnitude = ~ gp(unique_x, y, name = "valid") +
    gp(x, y, name = "space"),
    control = magCompControl(run_sampler = FALSE)), "spatial term 'space'.*identical coordinates")
  expect_no_error(magComp(dat, magnitude = ~ iid(plot),
    composition = magCompComposition(species ~ factor_iid(plot, K = 1),
                                     dbh | species ~ 1),
    control = magCompControl(run_sampler = FALSE)))
})

test_that("factor_iid terms carry group indices", {
  unit_data <- data.frame(
    plot = paste0("p", 1:4),
    stratum = c("a", "a", "b", "b"),
    z = 1:4
  )

  comp <- magComp:::magcomp_build_component(
    ~ z + factor_iid(stratum, K = 2),
    unit_data,
    "composition"
  )

  term <- comp$process_terms[[1L]]
  expect_equal(term$term_type, "factor_iid")
  expect_equal(term$group_name, "stratum")
  expect_equal(term$group, c(1L, 1L, 2L, 2L))
  expect_equal(term$group_levels, c("a", "b"))
  expect_equal(term$K, 2L)
})

test_that("factor terms finalize anchor order and names", {
  trees <- data.frame(
    plot = rep(paste0("p", 1:4), each = 3),
    cls = c("a", "a", "b",
            "b", "b", "c",
            "a", "c", "c",
            "b", "c", "c"),
    area = 0.1,
    stratum = rep(c("s1", "s1", "s2", "s2"), each = 3),
    x = rep(1:4, each = 3),
    y = rep(1:4, each = 3)
  )
  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    design = fixedArea(area)
  )

  comp <- magComp:::magcomp_build_component(
    ~ factor_iid(stratum, K = 2, anchors = c("b", "a")),
    dat$unit_data,
    "composition"
  )
  comp <- magComp:::magcomp_finalize_component_process_terms(comp, dat)
  term <- comp$process_terms[[1L]]
  expect_equal(term$anchors, c("b", "a"))
  expect_equal(term$factor_names, c("b", "a"))

  comp_named <- magComp:::magcomp_build_component(
    ~ factor_iid(stratum, K = 2, anchors = c("b", "a"),
                 factor_names = c("saw", "pole")),
    dat$unit_data,
    "composition"
  )
  comp_named <- magComp:::magcomp_finalize_component_process_terms(comp_named,
                                                                  dat)
  expect_equal(comp_named$process_terms[[1L]]$factor_names,
               c("saw", "pole"))

  comp_default <- magComp:::magcomp_build_component(
    ~ factor_iid(stratum, K = 2),
    dat$unit_data,
    "composition"
  )
  comp_default <- magComp:::magcomp_finalize_component_process_terms(
    comp_default,
    dat
  )
  expect_equal(comp_default$process_terms[[1L]]$factor_names,
               comp_default$process_terms[[1L]]$anchors)
  expect_equal(length(comp_default$process_terms[[1L]]$anchors), 2L)

  comp_gp <- magComp:::magcomp_build_component(
    ~ factor_gp(x, y, K = 1, anchors = "b", factor_names = "mid"),
    dat$unit_data,
    "composition"
  )
  comp_gp <- magComp:::magcomp_finalize_component_process_terms(comp_gp, dat)
  term <- comp_gp$process_terms[[1L]]
  expect_equal(term$term_type, "factor_gp")
  expect_equal(term$anchors, "b")
  expect_equal(term$factor_names, "mid")
  expect_equal(term$n_group, dat$n_unit)
  expect_null(term$neighbor_index)
})
