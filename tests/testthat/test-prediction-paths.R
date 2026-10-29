path_fit <- function(n_samples = 12L, n_chains = 2L, point = FALSE,
                     estimated_size = FALSE)
{
  rows <- expand.grid(plot = paste0("p", 1:4), cls = c("a", "b", "c"))
  rows$count <- rep(c(1L, 4L, 2L, 3L), 3L)
  rows$area <- rep(c(0.1, 0.2, 0.3, 0.4), 3L)
  rows$baf <- rep(c(10, 20, 10, 20), 3L)
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
    design = if(point) pointSampling(baf) else fixedArea(area))
  magComp(dat, family = if(point) pointBasalArea() else fixedAreaStems(),
    priors = if(estimated_size) magCompPriors(magnitude = list(size = c(2, 0.2)))
      else magCompPriors(magnitude = list(size = NULL)),
    control = magCompControl(n_samples = 5L + n_samples, n_burn = 5L,
      n_chains = n_chains, chain_seeds = 5300L + seq_len(n_chains),
      nb_size = 3, verbose = FALSE))
}

test_that("the four paths return coherent quantities and aligned draws", {
  fit <- path_fit()
  area <- fit$data$area
  for(target in c("expected", "observation")){
    set.seed(11)
    at_fit <- fitted(fit, target = target, summary = "samples")
    set.seed(11)
    at_new <- predict(fit, fit$data$unit_data, target = target,
                      summary = "samples")
    for(result in list(at_fit, at_new)){
      expect_s3_class(result, "magCompPrediction")
      expect_named(result, c("magnitude", "proportion", "density", "count", "metadata"))
      expect_equal(dim(result$count), c(4L, 3L, 12L, 2L))
      expect_equal(dim(result$magnitude), c(4L, 12L, 2L))
      expect_identical(result$metadata$chain_ids, fit$backend$chain_ids)
      expect_equal(result$metadata$sample_ids, 1:12)
      expect_equal(result$metadata$n_draws, 24L)
      total <- apply(result$count, c(1L, 3L, 4L), sum)
      expect_equal(result$magnitude, sweep(total, 1L, area, "/"))
      expect_equal(result$density, sweep(result$count, 1L, area, "/"))
      p <- sweep(result$count, c(1L, 3L, 4L), total, "/")
      p[is.nan(p)] <- NA_real_
      expect_equal(result$proportion, p)
      expect_output(print(result), paste0("target: +", target))
    }
    for(type in c("magnitude", "proportion", "density", "count")){
      expect_equal(at_fit[[type]], at_new[[type]])
      set.seed(11)
      one <- fitted(fit, target = target, type = type, summary = "samples")
      expect_named(one, c(type, "metadata"))
      expect_identical(one[[type]], at_fit[[type]])
      set.seed(11)
      one_new <- predict(fit, fit$data$unit_data, target = target, type = type,
                         summary = "samples")
      expect_identical(one_new[[type]], at_new[[type]])
    }
  }
  default <- fitted(fit)
  expect_equal(default$metadata$target, "expected")
  expect_equal(default$metadata$summary, "mean")
  expect_equal(default$count, apply(fitted(fit, summary = "samples")$count,
                                    1:2, mean))
})

test_that("observation totals use fixed or sampled NB size at both location paths", {
  for(point in c(FALSE, TRUE)) for(estimated in c(FALSE, TRUE)){
    fit <- path_fit(point = point, estimated_size = estimated)
    expected <- predict(fit, fit$data$unit_data, summary = "samples")
    mu <- sweep(expected$magnitude, 1L, fit$data$sampler_area, "*")
    set.seed(28)
    oracle <- array(0L, dim(expected$count), dimnames = dimnames(expected$count))
    for(ch in 1:2) for(s in 1:12) for(i in 1:4){
      size <- if(estimated) fit$backend$nb_size_samples[1L,s,ch] else 3
      total <- rnbinom(1, size = size, mu = mu[i,s,ch])
      oracle[i,,s,ch] <- as.vector(rmultinom(1, total, expected$proportion[i,,s,ch]))
    }
    for(source in c("fitted", "predict")){
      set.seed(28)
      result <- if(source == "fitted")
        fitted(fit, target = "observation", summary = "samples")
      else predict(fit, fit$data$unit_data, target = "observation",
                   summary = "samples")
      expect_identical(result$count, oracle)
      if(point){
        expect_equal(result$density, sweep(result$count, 1L, fit$data$baf, "*"))
        expect_equal(result$magnitude,
          sweep(apply(result$count, c(1L, 3L, 4L), sum), 1L, fit$data$baf, "*"))
      }
    }
  }
})

test_that("prediction draw limits preserve parameter alignment and the fitted object", {
  for(chains in 1:2){
    fit <- path_fit(n_chains = chains, estimated_size = TRUE)
    original <- fit$backend
    all <- predict(fit, fit$data$unit_data, summary = "samples")
    few <- predict(fit, fit$data$unit_data, n_draws = 3, summary = "samples")
    ids <- c(3L, 7L, 11L)
    expect_identical(few$metadata$sample_ids, ids)
    expect_equal(few$metadata$n_draws, 3L * chains)
    for(type in all$metadata$type){
      indices <- lapply(dim(all[[type]]), seq_len)
      indices[[if(type == "magnitude") 2L else 3L]] <- ids
      expected <- do.call(`[`, c(list(all[[type]]), indices, list(drop = FALSE)))
      expect_equal(few[[type]], expected)
    }
    # Each selected mean must use its matching NB size and allocation draws.
    oracle <- array(0L, c(4L, 3L, 3L, chains))
    mu <- array(few$count, dim(oracle))
    prob <- array(few$proportion, dim(oracle))
    size <- array(fit$backend$nb_size_samples, c(1L, 12L, chains))
    set.seed(23)
    for(ch in seq_len(chains)) for(s in 1:3) for(i in 1:4){
      total <- rnbinom(1, size = size[1L,ids[s],ch], mu = sum(mu[i,,s,ch]))
      oracle[i,,s,ch] <- rmultinom(1, total, prob[i,,s,ch])
    }
    set.seed(23)
    actual <- predict(fit, fit$data$unit_data, n_draws = 3,
                      target = "observation", summary = "samples")
    expect_equal(as.numeric(actual$count), as.numeric(oracle))
    expect_identical(fit$backend, original)
    expect_identical(predict(fit, fit$data$unit_data, n_draws = 12,
                             summary = "samples"), all)
    one <- predict(fit, fit$data$unit_data, n_draws = 1, summary = "samples")
    expect_equal(one$metadata$n_samples, 1L)
    expect_equal(one$metadata$sample_ids, 7L)
    aggregated <- standSummary(few, data.frame(id = few$metadata$unit_ids, ha = 1),
      unit = id, represented_area = ha, summary = "samples")
    expect_identical(aggregated$metadata$sample_ids, ids)
  }
  for(n in list(0, -1, 1.5, NA_real_, Inf, "3", c(1, 2)))
    expect_error(predict(fit, fit$data$unit_data, n_draws = n), "n_draws must be")
  expect_error(predict(fit, fit$data$unit_data, n_draws = 13), "exceeds")
})

test_that("zero totals give explicit conditional proportion summaries", {
  fit <- path_fit(n_samples = 200L, n_chains = 1L, point = TRUE)
  fit$backend$magnitude_samples[1L,] <- 0
  fit$backend$magnitude_samples[2L,] <- 2
  set.seed(29)
  draws <- fitted(fit, target = "observation", summary = "samples")
  set.seed(29)
  means <- fitted(fit, target = "observation")
  set.seed(29)
  quantiles <- fitted(fit, target = "observation", summary = "quantile",
                       probs = 0.5)
  expect_true(all(draws$count[1L,,] == 0))
  expect_true(all(draws$density[1L,,] == 0))
  expect_true(all(draws$magnitude[1L,] == 0))
  expect_true(all(is.na(draws$proportion[1L,,])))
  expect_true(all(is.na(means$proportion[1L,])))
  expect_true(all(is.na(quantiles$proportion[1L,,])))
  totals <- apply(draws$count, c(1L, 3L), sum)
  expect_equal(means$metadata$zero_total_probability, rowMeans(totals == 0))
  expect_equal(means$metadata$n_positive_draws, rowSums(totals > 0))
  expect_equal(means$metadata$proportion_conditioning, "positive_total")
  expect_gt(sum(totals[2L,] == 0), 0)
  expect_gt(sum(totals[2L,] > 0), 0)
  for(i in 2:4){
    p <- draws$proportion[i,,totals[i,] > 0,drop = FALSE]
    expect_equal(as.numeric(means$proportion[i,]),
                 rowMeans(matrix(p, nrow = 3L)))
    expect_equal(as.numeric(quantiles$proportion[i,,]),
                 apply(matrix(p, nrow = 3L), 1L, median))
  }
  expect_equal(dim(quantiles$count), c(4L, 3L, 1L))
  expect_equal(dim(quantiles$magnitude), c(4L, 1L))
  tab <- summary(draws, probs = 0.5)
  prop <- subset(tab, quantity == "proportion")
  expect_equal(prop$mean, as.numeric(means$proportion))
  expect_equal(prop$q50, as.numeric(quantiles$proportion))
  expect_equal(prop$n_draws, rep(as.numeric(rowSums(totals > 0)), 3L))
  expect_true(all(is.na(prop$sd[prop$unit == "p1"])))
  expect_equal(subset(tab, quantity == "count" & unit == "p1")$mean, rep(0, 3))
  expect_output(print(tab), "Proportions condition on a positive total")
})

test_that("prediction exposure is explicit for both sampling designs and targets", {
  for(point in c(FALSE, TRUE)){
    fit <- path_fit(n_samples = 2L, n_chains = 1L, point = point)
    newdata <- data.frame(plot = "new")
    expect_error(predict(fit, newdata), "prediction (area|exposure) is required")
    expected <- predict(fit, newdata, type = c("magnitude", "proportion", "density"))
    expect_null(expected$metadata$exposure)
    expect_equal(predict(fit, newdata, area = 0.4)$count, expected$density * 0.4)
    newdata$exposure <- 0.4
    expect_equal(predict(fit, newdata, area = "exposure")$count, expected$density * 0.4)
    if(point) newdata$baf <- 2.5 else newdata$area <- 0.4
    expect_equal(predict(fit, newdata)$count, expected$density * 0.4)
    expect_error(predict(fit, newdata, area = c(1, 2)), "one value per row")
    for(area in c(0, -1, NA, Inf))
      expect_error(predict(fit, newdata, area = area), "positive|finite|missing")
    newdata$bad_area <- factor("0.4")
    expect_error(predict(fit, newdata, area = "bad_area"), "must be numeric")
    for(type in c("magnitude", "proportion"))
      expect_error(predict(fit, data.frame(plot = "new"), type = type,
                           target = "observation"), "prediction (area|exposure) is required")
  }
})

test_that("both new-location targets retain shared iid and joint spatial effects", {
  rows <- expand.grid(plot = paste0("p", 1:6), cls = c("a", "b", "c"))
  rows$count <- rep(1:6, 3)
  rows$x <- rep(seq(0, 1, length.out = 6), 3)
  rows$y <- rep(c(0, 1, 0, 1, 0, 1), 3)
  rows$stand <- rep(rep(c("s1", "s2", "s3"), each = 2), 3)
  rows$area <- 0.2
  dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                     design = fixedArea(area))
  for(kind in c("gp", "nngp")){
    mag <- as.formula(paste0("~ iid(stand) + ", kind, "(x, y)"))
    comp <- as.formula(paste0("~ factor_iid(stand, K = 1) + factor_", kind,
                              "(x, y, K = 1)"))
    fit <- magComp(dat, magnitude = mag, composition = comp,
      control = magCompControl(n_samples = 17, n_burn = 5, nb_size = 3,
        chain_seeds = 5310, verbose = FALSE))
    full <- fitted(fit, summary = "samples")
    few <- predict(fit, dat$unit_data, n_draws = 3, summary = "samples")
    expect_equal(few$magnitude, full$magnitude[,c(3, 7, 11),drop = FALSE])
    expect_equal(few$proportion, full$proportion[,,c(3, 7, 11),drop = FALSE])
    for(target in c("expected", "observation")){
      set.seed(35)
      obs <- fitted(fit, target = target, summary = "samples")
      set.seed(35)
      repeated <- predict(fit, dat$unit_data, target = target,
                            summary = "samples")
      for(type in obs$metadata$type)
        expect_equal(repeated[[type]], obs[[type]])
    }
    # Identical new covariates, new group, and new coordinate must share all
    # latent effects. Replicate observations still get separate NB draws.
    newdata <- data.frame(plot = c("new1", "new2"), x = 0.4, y = 0.4,
                           stand = "new", area = 0.2)
    set.seed(36)
    expected <- predict(fit, newdata, n_draws = 3, summary = "samples")
    expect_equal(expected$magnitude[1L,], expected$magnitude[2L,])
    expect_equal(expected$proportion[1L,,], expected$proportion[2L,,])
    set.seed(36)
    replicated <- predict(fit, newdata, target = "observation",
                            n_draws = 3, summary = "samples")
    expect_false(identical(replicated$count[1L,,], replicated$count[2L,,]))
    expect_equal(replicated$magnitude,
                 apply(replicated$count, c(1L, 3L), sum) / 0.2)
    for(type in expected$metadata$type){
      set.seed(36)
      one <- predict(fit, newdata, type = type, target = "observation",
                      n_draws = 3, summary = "samples")
      expect_identical(one[[type]], replicated[[type]])
    }
  }
})

test_that("prediction requests are validated and singleton summaries retain dimensions", {
  fit <- path_fit(n_samples = 1L, n_chains = 1L)
  for(type in list(character(), c("count", "count"), "typo", NA_character_))
    expect_error(fitted(fit, type = type), "type must be")
  for(probs in list(numeric(), NA_real_, Inf, 0, 1))
    expect_error(fitted(fit, probs = probs), "probabilities")
  expect_error(fitted(fit, target = "typo"), "arg.*should be")
  expect_error(fitted(fit, area = 1), "unused fitted argument")
  expect_error(predict(fit), "newdata is required")
  expect_error(predict(fit, data.frame(plot = character())), "at least one row")
  q <- predict(fit, data.frame(plot = "new"), area = 1,
                 summary = "quantile", probs = 0.5)
  expect_equal(dim(q$magnitude), c(1L, 1L))
  expect_equal(dim(q$count), c(1L, 3L, 1L))
  expect_equal(as.numeric(q$count), as.numeric(predict(fit, data.frame(plot = "new"), area = 1)$count))
  stand <- standSummary(fitted(fit, summary = "samples"),
    data.frame(id = fit$data$unit_levels, ha = 1), unit = id, represented_area = ha,
                      summary = "quantile", probs = 0.5)$density
  expect_equal(dim(stand), c(1L, 3L, 1L))
  set.seed(31)
  first <- fitted(fit, target = "observation", summary = "samples")
  set.seed(31)
  ppc <- ppcSummary(fit, probs = 0.5)
  expect_equal(ppc$q50, as.numeric(apply(first$count, 1L, sum)))
})

test_that("mixed spatial rows preserve coordinate identity and every sample dimension", {
  spatial <- getFromNamespace("magcomp_spatial_prediction_samples", "magComp")
  observed <- matrix(c(0, 1, 0, 1), ncol = 2L)
  novel <- matrix(c(0.3, 0.3 + 1e-12, 0.5, 0.5), ncol = 2L)
  # Very close coordinates remain distinct; exactly identical coordinates
  # reuse the draw, even when observed and new rows are interleaved.
  mixed <- rbind(observed[2L,], novel[1L,], observed[1L,], novel[2L,], novel[1L,])
  for(method in c("gp", "nngp")) for(factor in c(FALSE, TRUE)){
    for(chains in 1:2){
      d <- if(factor) c(2L, 2L, 3L) else c(2L, 3L)
      if(chains > 1L) d <- c(d, chains)
      effects <- array(as.numeric(seq_len(prod(d))) / 10, d)
      set.seed(37)
      direct <- .Call("magcomp_spatial_predict_cpp", effects, observed, novel,
                      1.0, 1.0, 2L, method, factor, "exp", 1.0, PACKAGE = "magComp")
      set.seed(37)
      actual <- spatial(effects, observed, mixed, 1.0, 1.0, 2L, method, factor, "exp", 1)
      reference <- rbind(matrix(effects, nrow = 2L), matrix(direct, nrow = 2L))
      reference <- reference[c(2L, 3L, 1L, 4L, 3L),,drop = FALSE]
      expect_identical(actual, array(reference, c(5L, d[-1L])))
    }
  }
})


test_that("summary pools stored prediction draws with labels and custom quantiles", {
  fit <- path_fit()
  draws <- predict(fit, fit$data$unit_data, summary = "samples")
  tab <- summary(draws, probs = c(0.1, 0.9))
  expect_s3_class(tab, "magCompPredictionSummary")
  expect_equal(attr(tab, "metadata"), draws$metadata)
  expect_equal(attr(tab, "probs"), c(0.1, 0.9))
  expect_equal(nrow(tab), 4L * (1L + 3L * 3L))
  for(quantity in draws$metadata$type){
    part <- tab[tab$quantity == quantity, ]
    samples <- draws[[quantity]]
    categories <- if(quantity == "magnitude") "total" else dimnames(samples)[[2L]]
    reference <- matrix(samples, nrow = 4L * length(categories))
    expect_equal(part$unit, rep(draws$metadata$unit_ids, length(categories)))
    expect_equal(part$category, rep(categories, each = 4L))
    expect_equal(part$mean, as.numeric(apply(reference, 1L, mean)))
    expect_equal(part$sd, as.numeric(apply(reference, 1L, sd)))
    expect_equal(unname(as.matrix(part[, c("q10", "q90")])),
                 unname(t(apply(reference, 1L, quantile, probs = c(0.1, 0.9)))))
  }
  expect_equal(unique(tab$n_draws), 24)
  expect_output(print(tab), "target = 'expected'")
  expect_output(print(tab[, c("unit", "quantity", "mean")]), "quantity")
  expect_error(summary(predict(fit, fit$data$unit_data)), "summary = 'samples'")
  expect_error(summary(fitted(fit, summary = "quantile")), "retained prediction draws")
  expect_error(summary(draws, probs = numeric()), "probabilities")

  # One unit, one category, one draw: quantiles stay defined but SD is not.
  one <- fitted(fit, type = "count", summary = "samples")
  one$count <- one$count[1, 1, 1, 1, drop = FALSE]
  one$metadata$unit_ids <- one$metadata$unit_ids[1]
  one$metadata$n_samples <- one$metadata$n_chains <- one$metadata$n_draws <- 1L
  one_tab <- summary(one, probs = 0.5)
  expect_equal(nrow(one_tab), 1L)
  expect_equal(one_tab$q50, as.numeric(one$count))
  expect_true(is.na(one_tab$sd))
  expect_equal(one_tab$n_draws, 1)
})
