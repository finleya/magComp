cv_prediction_data <- function()
{
  rows <- expand.grid(plot = paste0("p", 1:6), cls = c("a", "b"))
  rows$count <- c(0, 1, 3, 2, 0, 5, 0, 2, 1, 3, 0, 1)
  rows$area <- rep(c(0.01, 0.05, 0.1, 0.2, 0.4, 0.8), 2)
  rows$x <- rep(seq(-1, 1, length.out = 6), 2)
  magCompData(rows, unit = plot, category = cls, tree_count = count,
              design = fixedArea(area))
}

test_that("CV subsets supplied latent starts and initializes missing values within folds", {
  dat <- cv_prediction_data()
  dat$data$y <- dat$data$x^2
  dat$unit_data$y <- dat$unit_data$x^2
  mag <- ~ x + iid(plot) + gp(x, y)
  comp <- ~ factor_iid(plot, anchors = "a") + factor_nngp(x, y, m = 2, anchors = "a")
  model <- magComp(dat, magnitude = mag, composition = comp,
    control = magCompControl(run_sampler = FALSE))
  starts <- magCompStarting(model)
  # Named rows must match after reordering, including groups absent in a fold.
  starts$magnitude$terms$iid_1$effects <- rev(starts$magnitude$terms$iid_1$effects)
  starts$composition$blocks$cls$terms$factor_iid_1$scores[,1] <- seq_len(6)
  ctrl <- magCompControl(n_samples = 5, n_burn = 2, chain_seeds = 712, verbose = FALSE)
  folds <- setNames(rep(1:2, 3), dat$unit_levels)
  cv <- magCompCV(dat, magnitude = mag, composition = comp, starting = starts,
                  control = ctrl, fold_id = folds, keep_fits = TRUE)
  for(f in names(cv$fits)){
    fit <- cv$fits[[f]]
    units <- fit$data$unit_levels
    expect_equal(fit$starting$magnitude$terms$iid_1$effects,
                  starts$magnitude$terms$iid_1$effects[units])
    expect_equal(fit$starting$magnitude$terms$gp_1$effects,
                  starts$magnitude$terms$gp_1$effects[units])
    for(term in c("factor_iid_1", "factor_nngp_1"))
      expect_equal(fit$starting$composition$blocks$cls$terms[[term]]$scores,
        starts$composition$blocks$cls$terms[[term]]$scores[units,,drop = FALSE])
  }
  # Unnamed arrays follow the full-data template and must gain names before
  # fold subsetting; covariance fixed markers must survive that normalization.
  unnamed <- starts
  for(t in names(unnamed$magnitude$terms))
    unnamed$magnitude$terms[[t]]$effects <- unname(starts$magnitude$terms[[t]]$effects[
      names(model$starting$magnitude$terms[[t]]$effects)])
  for(t in names(unnamed$composition$blocks$cls$terms))
    unnamed$composition$blocks$cls$terms[[t]]$scores <-
      unname(starts$composition$blocks$cls$terms[[t]]$scores)
  unnamed$magnitude$terms$gp_1$covariance$phi <- fixed(1)
  unnamed_cv <- magCompCV(dat, magnitude = mag, composition = comp, starting = unnamed,
    control = ctrl, fold_id = folds, keep_fits = TRUE)
  for(fit in unnamed_cv$fits){
    expect_equal(fit$starting$magnitude$terms$gp_1$effects,
      starts$magnitude$terms$gp_1$effects[fit$data$unit_levels])
    expect_equal(fit$starting$composition$blocks$cls$terms$factor_nngp_1$scores,
      starts$composition$blocks$cls$terms$factor_nngp_1$scores[fit$data$unit_levels,,drop = FALSE])
    expect_s3_class(fit$starting$magnitude$terms$gp_1$covariance$phi, "magCompFixed")
  }
  partial <- list(magnitude = list(terms = list(iid_1 = list(variance = 0.5))))
  cv <- magCompCV(dat, magnitude = mag, starting = partial, control = ctrl,
                  fold_id = folds, keep_fits = TRUE)
  for(fit in cv$fits){
    direct <- magComp(fit$data, magnitude = mag, starting = partial,
    control = magCompControl(run_sampler = FALSE))
    expect_equal(fit$starting, direct$starting)
  }
})

test_that("CV returns predict defaults, quantities, exposures and observed counterparts", {
  dat <- cv_prediction_data()
  folds <- setNames(c("B", "A", "A", "B", "A", "B"), dat$unit_levels)
  cv <- magCompCV(dat, fold_id = folds, keep_fits = TRUE,
    control = magCompControl(n_samples = 7, n_burn = 2, chain_seeds = 311,
                             verbose = FALSE))
  expect_equal(cv$type, c("magnitude", "proportion", "density", "count"))
  expect_equal(c(cv$target, cv$summary), c("expected", "mean"))
  expect_equal(names(cv$predictions), c("B", "A"))
  for(f in names(cv$predictions)){
    rows <- dat$unit_data[names(folds)[folds == f],,drop = FALSE]
    direct <- predict(cv$fits[[f]], rows, area = rows$area)
    expect_equal(cv$predictions[[f]], direct)
    expect_equal(cv$totals[[f]]$density, direct$magnitude)
    expect_equal(cv$totals[[f]]$count[,1], rowSums(direct$count))
  }
  expect_equal(cv$observed$count, dat$counts)
  expect_equal(cv$observed$density, dat$counts / dat$area)
  expect_equal(cv$observed$magnitude[,1], rowSums(dat$counts) / dat$area)
  expect_true(all(is.na(cv$observed$proportion[c("p1", "p5"),])))
})

test_that("CV observation draws match direct fold fitting and prediction", {
  dat <- cv_prediction_data()
  ctrl <- magCompControl(n_samples = 10, n_burn = 2, n_chains = 2,
                         chain_seeds = c(312, 313), verbose = FALSE)
  # Include an unseen spatial location and a varying exposure. This oracle
  # checks fold construction, argument forwarding, and full joint predictions.
  dat$unit_data$y <- c(0.1, 0.4, 0.2, 0.9, 0.6, 0.3)
  dat$data$y <- dat$unit_data$y[match(dat$data$plot, dat$unit_levels)]
  cv <- magCompCV(dat, magnitude = ~nngp(x, y, m = 2), fold_id = setNames(rep(1:2, each = 3), dat$unit_levels),
    target = "observation", summary = "samples", control = ctrl, starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(1)))))))
  rebuild <- getFromNamespace("magcomp_rebuild_data", "magComp")
  for(f in names(cv$predictions)){
    heldout <- names(cv$fold_id)[cv$fold_id == f]
    training <- rebuild(dat, !dat$data$plot %in% heldout)
    fit <- magComp(training, magnitude = ~nngp(x, y, m = 2), control = ctrl, starting = list(magnitude = list(terms = list(nngp_1 = list(covariance = list(phi = fixed(1)))))))
    rows <- dat$unit_data[heldout,,drop = FALSE]
    direct <- predict(fit, rows, area = rows$area, target = "observation",
                       summary = "samples")
    expect_equal(cv$predictions[[f]], direct)
    total <- apply(direct$count, c(1L, 3L, 4L), sum)
    expect_equal(cv$totals[[f]]$count, total)
    expect_equal(direct$magnitude, sweep(total, 1L, rows$area, "/"))
    expect_equal(direct$density, sweep(direct$count, 1L, rows$area, "/"))
    expect_equal(direct$metadata$zero_total_probability,
                 setNames(rowMeans(matrix(total == 0, nrow = nrow(rows))), heldout))
  }
})

test_that("CV summaries preserve zero-total conditioning and summarize joint totals", {
  dat <- cv_prediction_data()
  folds <- setNames(c(2L, 1L, rep(2L, 4)), dat$unit_levels)
  ctrl <- magCompControl(n_samples = 14, n_burn = 2, n_chains = 2,
                         chain_seeds = c(314, 315), verbose = FALSE)
  samples <- magCompCV(dat, fold_id = folds, target = "observation",
                       summary = "samples", control = ctrl)
  for(summary in c("mean", "quantile")){
    cv <- magCompCV(dat, fold_id = folds, target = "observation",
                    summary = summary, probs = 0.5, control = ctrl)
    for(f in names(cv$predictions)) for(type in cv$type){
      raw <- samples$predictions[[f]][[type]]
      n <- dim(raw)[1L]
      J <- if(type == "magnitude") 1L else dim(raw)[2L]
      draws <- matrix(raw, n * J)
      expected <- apply(draws, 1L, function(x){
        x <- x[!is.na(x)]
        if(!length(x)) return(NA_real_)
        if(summary == "mean") mean(x) else quantile(x, 0.5, names = FALSE)
      })
      expect_equal(as.numeric(cv$predictions[[f]][[type]]), expected)
      if(type %in% c("count", "density")){
        total <- matrix(samples$totals[[f]][[type]], nrow = n)
        expected_total <- if(summary == "mean") rowMeans(total) else
          apply(total, 1L, quantile, probs = 0.5, names = FALSE)
        expect_equal(as.numeric(cv$totals[[f]][[type]]), unname(expected_total))
      }
    }
    expect_equal(cv$predictions[[1]]$metadata$proportion_conditioning, "positive_total")
    expect_equal(cv$predictions[[1]]$metadata$n_positive_draws,
                 samples$predictions[[1]]$metadata$n_positive_draws)
    image_file <- tempfile(fileext = ".png")
    expect_no_error(plot(cv, type = "count", level = "total", file = image_file))
    expect_true(file.exists(image_file))
    unlink(image_file)
  }
  expect_true(any(samples$predictions[[1]]$metadata$zero_total_probability > 0))
  image_file <- tempfile(fileext = ".png")
  expect_no_error(plot(samples, type = "proportion", file = image_file))
  unlink(image_file)
  expect_error(plot(samples, type = "proportion", level = "total"), "total proportions")
})

test_that("CV validates prediction choices before fitting and retains subsets", {
  dat <- cv_prediction_data()
  for(args in list(list(type = character()), list(type = c("count", "count")),
                   list(target = "bad"), list(summary = "bad"), list(probs = NA_real_)))
    expect_error(do.call(magCompCV, c(list(data = dat), args)))
  # Fixed-area data do not supply the BAF needed for a point-sampling fit.
  expect_error(magCompCV(dat, family = pointBasalArea(), target = "observation"),
               "plot area or BAF")
  folds <- setNames(rep(1:2, 3), dat$unit_levels)
  cv <- magCompCV(dat, fold_id = folds, type = c("density", "magnitude"),
    control = magCompControl(n_samples = 5, n_burn = 2, chain_seeds = 316,
                             verbose = FALSE))
  expect_equal(names(cv$predictions[[1]]), c("density", "magnitude", "metadata"))
  expect_equal(names(cv$observed), c("density", "magnitude"))
  expect_equal(names(cv$totals[[1]]), "density")
  expect_null(cv$fits)
  expect_equal(sort(unlist(lapply(cv$predictions, function(x) x$metadata$unit_ids),
                            use.names = FALSE)), sort(dat$unit_levels))
  image_file <- tempfile(fileext = ".png")
  expect_no_error(plot(cv, type = "density", level = "cell", file = image_file))
  expect_true(file.exists(image_file))
  unlink(image_file)
  expect_output(print(cv), "expected")
  expect_error(plot(cv, type = "count"), "type must select")
})
