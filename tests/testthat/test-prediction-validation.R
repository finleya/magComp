test_that("predictive counts reject invalid category means and overflowing totals", {
  predict_nb <- getFromNamespace("magcomp_posterior_predict_nb", "magComp")
  for(d in list(c(2L, 2L, 3L), c(2L, 2L, 3L, 2L))){
    for(bad in c(Inf, -Inf, NA_real_, NaN, -1)){
      lambda <- array(2, d)
      lambda[length(lambda)] <- bad
      expect_error(predict_nb(lambda, 2),
                   "predictive category means must be finite and nonnegative")
    }
  }
  expect_error(predict_nb(array(c(-1, 1), c(1L, 2L, 1L)), 2),
               "predictive category means must be finite and nonnegative")
  expect_error(
    predict_nb(array(.Machine$double.xmax, c(1L, 2L, 1L)), 2),
    "predictive total means must be finite"
  )
})

test_that("predictive counts preserve zero means and sample dimensions", {
  predict_nb <- getFromNamespace("magcomp_posterior_predict_nb", "magComp")
  set.seed(917)
  for(d in list(c(2L, 3L, 8L), c(2L, 3L, 8L, 2L))){
    labels <- list(c("empty", "nonempty"), c("a", "b", "c"), NULL)
    if(length(d) == 4L)
      labels <- c(labels, list(c("chain1", "chain2")))
    lambda <- array(c(0, 0, 0, 6, 0, 4), d, dimnames = labels)
    size_draws <- if(length(d) == 3L)
      matrix(seq_len(d[3L]), nrow = 1L) else
        array(seq_len(d[3L]), c(1L, d[3L], d[4L]))
    for(size in list(2, size_draws)){
      counts <- predict_nb(lambda, size)
      expect_identical(dim(counts), d)
      expect_identical(dimnames(counts), dimnames(lambda))
      expect_true(all(counts[lambda == 0] == 0))
      expect_true(all(is.finite(counts) & counts >= 0))
      expect_gt(sum(counts), 0)
    }
  }
})

test_that("magCompCV accepts one quantile with singleton and larger holdout folds", {
  rows <- expand.grid(plot = paste0("p", 1:4), cls = c("a", "b"))
  rows$cls <- factor(rows$cls, levels = c("a", "b"))
  rows$count <- c(1L, 4L, 2L, 3L, 2L, 1L, 3L, 2L)
  rows$area <- 0.1
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                    design = fixedArea(area))
  folds <- setNames(c(1L, 2L, 2L, 2L), dat$unit_levels)
  for(n_chains in 1:2){
    set.seed(919)
    cv <- magCompCV(
      dat, fold_id = folds, type = "count", summary = "quantile", probs = 0.5,
      control = magCompControl(n_samples = 6, n_burn = 2,
                               n_chains = n_chains,
                               chain_seeds = seq_len(n_chains) + 919L,
                               verbose = FALSE)
    )
    expect_s3_class(cv, "magCompCV")
    for(f in names(cv$predictions)){
      value <- cv$predictions[[f]]$count
      expect_equal(dim(value), c(sum(folds == f), 2L, 1L))
      expect_equal(dimnames(value)[[3L]], "50%")
      expect_equal(dim(cv$totals[[f]]$count), c(sum(folds == f), 1L))
      expect_true(all(is.finite(value)))
    }
  }
})
