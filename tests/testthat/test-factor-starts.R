conditional_start_rows <- function()
{
  rows <- expand.grid(plot = paste0("p", 1:24), species = c("a", "b", "c"),
                       dbh = c("s", "m", "l"))
  i <- as.integer(rows$plot)
  a <- as.integer(rows$species)
  j <- as.integer(rows$dbh)
  rows$count <- 5L + (i^2 * a + i * j + 7L * a * j) %% 23L
  rows$group <- paste0("g", (i - 1L) %% 6L + 1L)
  rows$z <- seq(-1, 1, length.out = 24)[i]
  rows$x <- i / 24
  rows$y <- (i * 7L %% 23L) / 23
  rows$area <- 1
  rows
}

conditional_start_data <- function(rows = conditional_start_rows())
{
  magCompData(rows, unit = plot, categories = c(species, dbh),
              tree_count = count, coords = c(x, y), design = fixedArea(area))
}

test_that("child factor starts recover sequential residual SVD surfaces", {
  dat <- conditional_start_data()
  composition <- magCompComposition(
    species ~ z,
    dbh | species ~ z +
      factor_gp(x, y, K = 2, name = "space",
                anchors = c("b:m", "c:s")) +
      z:factor_iid(group, K = 2, name = "slope", anchors = c("b:m", "c:s")) +
      factor_iid(group, K = 2, name = "group", anchors = c("b:m", "c:s")))
  for(link in c("probit", "logit")){
    fit <- magComp(dat, composition = composition,
      control = magCompControl(composition_link = link, n_samples = 5,
        n_burn = 2, verbose = FALSE))
    block <- dat$allocation_blocks[["dbh|species"]]
    component <- fit$composition$allocation_components[["dbh|species"]]
    starts <- magCompStarting(fit)$composition$blocks[["dbh|species"]]
    # Independent wide representation: one row per unit, one column per
    # parent/child stick, with separate OLS regressions within each parent.
    residual <- matrix(0, dat$n_unit, block$coef_stick_count)
    for(a in seq_along(block$parent_levels)){
      rows <- which(block$parent_index == a)
      y <- block$counts[rows,,drop = FALSE]
      X <- component$X[rows,,drop = FALSE]
      eta <- cbind((y[,1] + 0.5) / (rowSums(y) + 1),
                    (y[,2] + 0.5) / (y[,2] + y[,3] + 1))
      eta <- if(link == "probit") qnorm(eta) else qlogis(eta)
      coef <- solve(crossprod(X), crossprod(X, eta))
      columns <- (a - 1L) * 2L + 1:2
      expect_equal(unname(starts$coefficients[columns,]), unname(t(coef)))
      residual[,columns] <- eta - X %*% coef
    }

    # Iid terms precede the spatial term even though it is first in the formula.
    for(name in c("slope", "group", "space")){
      term <- Filter(function(x) x$name == name, component$process_terms)[[1L]]
      first <- match(seq_len(dat$n_unit), block$unit_index)
      group <- if(name == "space") seq_len(dat$n_unit) else term$group[first]
      multiplier <- if(name == "slope") term$covariate[first] else rep(1, dat$n_unit)
      design <- diag(term$n_group)[group,,drop = FALSE] * multiplier
      surface <- solve(crossprod(design), crossprod(design, residual))
      sv <- svd(surface)
      expected <- sweep(sv$u[,1:2,drop = FALSE], 2L, sv$d[1:2], "*") %*%
        t(sv$v[,1:2,drop = FALSE])
      part <- starts$terms[[name]]
      actual <- part$scores %*% t(part$loadings)
      expect_equal(unname(actual), unname(expected), tolerance = 1e-8)
      expect_equal(rownames(part$loadings), block$stick_names)
      expect_equal(nrow(part$scores), term$n_group)
      expect_true(all(part$loadings[cbind(term$anchor_index, 1:2)] > 0))
      expect_equal(unname(part$loadings[term$anchor_index[1], 2]), 0)
      residual <- residual - design %*% actual
    }
    expect_true(all(is.finite(unlist(samples(fit)))))
  }
})

test_that("conditional factor starts preserve supplied named values", {
  dat <- conditional_start_data()
  composition <- magCompComposition(species ~ 1,
    dbh | species ~ z + factor_nngp(x, y, K = 2, name = "shape", m = 4,
      anchors = c("a:s", "c:m"), factor_names = c("f1", "f2")))
  model <- magComp(dat, composition = composition,
    control = magCompControl(run_sampler = FALSE))
  template <- magCompStarting(model)$composition$blocks[["dbh|species"]]
  template$terms$shape$loadings["c:m", "f1"] <- 0.4
  template$terms$shape$loadings["a:s", "f1"] <- 0.45
  template$terms$shape$loadings["c:m", "f2"] <- 1.7
  coefficients <- template$coefficients + 0.1
  for(fields in list(c("loadings", "scores"), "loadings", "scores")){
    supplied <- template$terms$shape[fields]
    supplied <- lapply(supplied, function(x)
      x[rev(seq_len(nrow(x))), rev(seq_len(ncol(x))), drop = FALSE])
    starting <- list(composition = list(blocks = list("dbh|species" = list(
      coefficients = coefficients[6:1,2:1], terms = list(shape = supplied)))))
    fit <- magComp(dat, composition = composition, starting = starting,
    control = magCompControl(run_sampler = FALSE))
    actual <- magCompStarting(fit)$composition$blocks[["dbh|species"]]
    expect_equal(actual$coefficients, coefficients)
    expect_equal(actual$terms$shape[fields], template$terms$shape[fields])
    expect_true(all(is.finite(unlist(actual))))
  }
})

test_that("uninformative child rows do not dilute factor starts", {
  rows <- conditional_start_rows()
  rows$count[rows$species == "c"] <- 0L
  composition <- magCompComposition(species ~ 1,
    dbh | species ~ z + factor_iid(group, K = 2, name = "shape",
                                   anchors = c("a:s", "b:s")))
  original <- magComp(conditional_start_data(rows), composition = composition,
    control = magCompControl(run_sampler = FALSE))
  extra <- rows
  extra$plot <- paste0(as.character(extra$plot), "_empty")
  extra$count <- 0L
  expanded <- rbind(transform(rows, plot = as.character(plot)), extra)
  again <- magComp(conditional_start_data(expanded), composition = composition,
    control = magCompControl(run_sampler = FALSE))
  initial <- magCompStarting(original)$composition$blocks[["dbh|species"]]
  expect_equal(magCompStarting(again)$composition$blocks[["dbh|species"]], initial,
                tolerance = 1e-10)
  expect_true(all(initial$coefficients[grep("^c:", rownames(initial$coefficients)),] == 0))

  # No residual information, more factors than groups, and zero slope multipliers.
  rows$count <- 0L
  rows$group <- "one"
  rows$z <- 0
  empty <- magComp(conditional_start_data(rows), composition = magCompComposition(
    species ~ 1, dbh | species ~ z + z:factor_iid(group, K = 3)),
    control = magCompControl(run_sampler = FALSE))
  starts <- magCompStarting(empty)$composition$blocks[["dbh|species"]]
  expect_true(all(is.finite(unlist(starts))))
  expect_true(all(starts$terms$factor_iid_1$scores == 0))
  term <- empty$composition$allocation_components[["dbh|species"]]$process_terms[[1]]
  expect_equal(unname(starts$terms$factor_iid_1$loadings[term$anchor_index,]), diag(3))
})
