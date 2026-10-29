component_data <- function(categories = TRUE, area = TRUE, nested = FALSE)
{
  rows <- expand.grid(plot = paste0("p", 1:8), species = c("a", "b", "c"),
                      dbh = if(nested) 1:4 else 1)
  if(nested) rows$dbh <- ordered(rows$dbh, levels = 1:4)
  i <- as.integer(rows$plot)
  rows$count <- (i + as.integer(rows$species) + as.integer(rows$dbh)) %% 5L
  rows$count[i == 1L] <- 0L
  rows$area <- 0.05 + i / 100
  rows$x <- (i - 4) / 4
  rows$y <- c(1, 3, 2, 5, 4, 8, 6, 7)[i]
  rows$group <- rep(c("g1", "g2", "g3", "g4"), 2)[i]
  do.call(magCompData, list(data = rows, unit = "plot", tree_count = "count",
    categories = if(!categories) NULL else if(nested) c("species", "dbh") else "species",
    design = if(area) fixedArea(area) else NULL))
}

component_control <- function(...)
  magCompControl(n_samples = 14, n_burn = 6, n_chains = 2,
    chain_seeds = c(710, 711), verbose = FALSE, ...)

test_that("magnitude-only fits depend on totals and retain the NB observation model", {
  categorized <- component_data()
  totals <- component_data(categories = FALSE)
  ctrl <- component_control(nb_size = 3)
  fixed_size <- magCompPriors(magnitude = list(size = NULL))
  fit <- magComp(totals, magnitude = ~ x, composition = NULL,
    priors = fixed_size, control = ctrl)
  same <- magComp(categorized, magnitude = ~ x, composition = NULL,
    priors = fixed_size, control = ctrl)
  expect_equal(samples(fit), samples(same))
  expect_named(samples(fit), "magnitude")
  retained <- predict(fit, fit$data$unit_data, summary = "samples", keep_effects = TRUE)
  expect_equal(samples(retained), samples(fit))
  expect_named(magCompStarting(fit, use = "last"), "magnitude")
  expect_named(magCompSettings(fit)$priors, "magnitude")
  expect_equal(unique(summary(fit)$component), "magnitude")
  expect_error(samples(fit, type = "combined"), "requires random-effect terms")
  expect_null(fit$backend$proportion_samples)
  expect_equal(fit$backend$profile_pg[["composition_draws"]], 0)
  expect_equal(fit$backend$profile_pg[["probit_composition_draws"]], 0)

  for(estimated in c(FALSE, TRUE)){
    if(estimated) fit <- magComp(totals, magnitude = ~ x, composition = NULL,
      priors = magCompPriors(magnitude = list(size = c(2, 0.5))), control = ctrl)
    expected <- fitted(fit, summary = "samples")
    limited <- predict(fit, fit$data$unit_data, n_draws = 3, summary = "samples")
    expect_equal(limited$count, expected$count[,,c(2, 5, 7),,drop = FALSE])
    expect_named(expected, c("magnitude", "count", "metadata"))
    expect_equal(dim(expected$count), c(8L, 1L, 8L, 2L))
    expect_equal(expected$count[,1,,], sweep(expected$magnitude, 1L, totals$area, "*"))
    expect_equal(predict(fit, fit$data$unit_data)$count, fitted(fit)$count)
    set.seed(719)
    oracle <- expected$count
    for(ch in 1:2) for(s in 1:8) for(i in 1:8){
      size <- if(estimated) fit$backend$nb_size_samples[1,s,ch] else 3
      oracle[i,1,s,ch] <- rnbinom(1, size = size, mu = expected$count[i,1,s,ch])
    }
    set.seed(719)
    observed <- fitted(fit, target = "observation", summary = "samples")
    expect_equal(observed$count, oracle)
    stand <- standSummary(observed, data.frame(id = observed$metadata$unit_ids, ha = 2),
      unit = id, represented_area = ha, type = "total", summary = "samples")
    expect_equal(as.numeric(stand$total), as.numeric(apply(observed$magnitude * 2, c(2, 3), sum)))
    expect_identical(stand$metadata$target, "observation")
    expect_equal(observed$magnitude, sweep(observed$count[,1,,], 1L, totals$area, "/"))
  }
  expect_error(fitted(fit, type = "proportion"), "for this model")
  expect_error(ppcSummary(fit, type = "category_total"), "this model supports")
  expect_error(magComp(totals, magnitude = NULL, composition = NULL), "at least one")
  expect_error(magComp(totals), "composition requires category counts")
  expect_error(magComp(totals, composition = NULL,
    starting = list(composition = list())), "unknown setting")
  expect_error(magComp(totals, composition = NULL,
    priors = magCompPriors(composition = list(slope_variance = 2))), "unused")
})

test_that("composition-only predictions condition counts on observed or supplied totals", {
  dat <- component_data(area = FALSE)
  fit <- magComp(dat, magnitude = NULL, composition = ~ x, control = component_control())
  expect_named(samples(fit), "composition")
  retained <- predict(fit, fit$data$unit_data, summary = "samples", keep_effects = TRUE)
  expect_equal(samples(retained), samples(fit))
  expect_named(magCompStarting(fit, use = "mean"), "composition")
  expect_named(magCompSettings(fit)$priors, "composition")
  expect_null(fit$backend$magnitude_samples)
  expect_null(fit$backend$nb_size)
  expect_equal(fit$backend$profile_pg[["magnitude_draws"]], 0)
  expect_equal(unique(summary(fit)$component), "composition")
  expect_error(samples(fit, type = "combined"), "requires random-effect terms")
  expect_named(fitted(fit), c("proportion", "metadata"))
  expect_equal(predict(fit, dat$unit_data)$proportion, fitted(fit)$proportion)
  probability <- fitted(fit, summary = "samples")$proportion
  limited <- predict(fit, dat$unit_data, n_draws = 3, summary = "samples")
  expect_equal(limited$proportion, probability[,,c(2, 5, 7),,drop = FALSE])
  total <- rowSums(dat$counts)
  means <- fitted(fit, type = "count", summary = "samples")
  expect_equal(means$count, sweep(probability, 1L, total, "*"))
  expect_equal(means$metadata$count_conditioning, "observed_total")
  set.seed(721)
  oracle <- array(0L, dim(probability), dimnames(probability))
  for(ch in 1:2) for(s in 1:8) for(i in 1:8)
    oracle[i,,s,ch] <- as.vector(rmultinom(1, total[i], probability[i,,s,ch]))
  set.seed(721)
  draws <- fitted(fit, type = c("count", "proportion"),
    target = "observation", summary = "samples")
  expect_identical(draws$count, oracle)
  expect_true(all(is.na(draws$proportion[1,,,])))
  expect_equal(draws$metadata$zero_total_probability, setNames(as.numeric(total == 0), names(total)))
  expect_equal(draws$metadata$n_positive_draws, setNames(16 * as.numeric(total > 0), names(total)))
  expect_error(fitted(fit, type = "density"), "area or BAF")
  expect_error(predict(fit, dat$unit_data, target = "observation"), "require total_count")
  expect_error(predict(fit, dat$unit_data, type = "count", total_count = rep(1.5, 8)), "integer")
  expect_error(predict(fit, dat$unit_data, type = "count", total_count = rep(Inf, 8)), "integer")
  expect_error(predict(fit, dat$unit_data, type = "count", total_count = c(1, 2)), "one value")
  expect_error(fitted(fit, type = "magnitude"), "for this model")
  expect_error(ppcSummary(fit, type = "unit_total"), "this model supports")

  newdata <- dat$unit_data[1,,drop = FALSE]
  newdata$N <- 20
  expected <- predict(fit, newdata, type = c("count", "density", "proportion"),
    total_count = "N", area = 0.2, summary = "samples")
  expect_equal(expected$count, expected$proportion * 20)
  expect_equal(expected$density, expected$count / 0.2)
  stand <- standSummary(expected, data.frame(id = expected$metadata$unit_ids, ha = 2),
      unit = id, represented_area = ha, summary = "samples")
  expect_equal(as.numeric(stand$density), as.numeric(expected$density))
  expect_identical(stand$metadata$count_conditioning, "supplied_total")
  expect_equal(dim(expected$count), c(1L, 3L, 8L, 2L))
  expect_equal(dim(predict(fit, newdata, summary = "quantile", probs = 0.5)$proportion), c(1L, 3L, 1L))
  expect_equal(predict(fit, dat$unit_data, type = "count", total_count = rev(total))$count,
    fitted(fit, type = "count")$count)
})

test_that("standalone spatial and iid terms retain nested factors and restart settings", {
  dat <- component_data(nested = TRUE)
  for(kind in c("iid", "gp", "nngp")){
    mag <- switch(kind, iid = ~ x + iid(group), gp = ~ x + gp(x,y),
                   nngp = ~ x + nngp(x,y,m=3))
    cmp <- switch(kind,
      iid = magCompComposition(species ~ x + factor_iid(group,K=1), dbh | species ~ x + factor_iid(group,K=1)),
      gp = magCompComposition(species ~ x + factor_gp(x,y,K=1), dbh | species ~ x + factor_gp(x,y,K=1)),
      nngp = magCompComposition(species ~ x + factor_nngp(x,y,K=1,m=3), dbh | species ~ x + factor_nngp(x,y,K=1,m=3)))
    for(mode in c("magnitude", "composition")){
      fit <- magComp(dat, magnitude = if(mode == "magnitude") mag else NULL,
        composition = if(mode == "composition") cmp else NULL,
        smooth = if(mode == "composition") magCompSmooth(dbh = list(intercept = "rw2")) else magCompSmooth(),
        control = component_control())
      expect_named(samples(fit), mode)
      settings <- magCompSettings(fit)
      again <- magComp(dat, magnitude = if(mode == "magnitude") mag else NULL,
        composition = if(mode == "composition") cmp else NULL,
        smooth = fit$smooth, priors = settings$priors,
        starting = magCompStarting(fit, use = "last", chain = 2),
        control = magCompControl(n_samples=2,n_burn=0,verbose=FALSE))
      expect_true(again$fitted)
      expect_equal(unique(summary(fit)$component), mode)
      full <- fitted(fit, summary = "samples")
      few <- predict(fit, dat$unit_data, n_draws = 3, summary = "samples")
      if(mode == "composition")
        expect_equal(few$proportion, full$proportion[,,c(2, 5, 7),,drop = FALSE])
      else expect_equal(few$magnitude, full$magnitude[,c(2, 5, 7),,drop = FALSE])
      newdata <- dat$unit_data[1:2,]
      newdata$x <- newdata$x + 0.05
      newdata$group <- "new"
      result <- predict(fit, newdata, n_draws = 3, summary = "samples")
      if(mode == "composition"){
        expect_equal(dim(result$proportion), c(2L, 12L, 3L, 2L))
        expect_equal(as.numeric(apply(result$proportion,c(1,3,4),sum)), rep(1,12))
      }else expect_true(all(is.finite(result$magnitude)))
      if(kind != "iid") expect_equal(unique(spatialDiagnostics(fit)$summary$component), mode)
    }
  }
})

test_that("standalone CV and predictive checks use the appropriate response", {
  dat <- component_data()
  folds <- setNames(rep(1:2,4), dat$unit_levels)
  for(mode in c("magnitude", "composition")){
    mag <- if(mode == "magnitude") ~ x else NULL
    comp <- if(mode == "composition") ~ x else NULL
    fit <- magComp(dat,magnitude=mag,composition=comp,control=component_control())
    ppc <- ppcSummary(fit)
    expect_equal(unique(ppc$type), if(mode == "magnitude") "unit_total" else "unit_category")
    expect_output(print(ppc), if(mode == "composition") "conditional on observed totals" else "observed")
    expect_true(all(is.finite(ppc$mean)))
    cv <- magCompCV(dat,magnitude=mag,composition=comp,control=component_control(),
      fold_id=folds,type="count",target="observation",summary="samples",keep_fits=TRUE)
    expect_equal(cv$observed$count, if(mode == "composition") dat$counts else
      matrix(rowSums(dat$counts), ncol=1L, dimnames=list(dat$unit_levels,"total")))
    for(f in names(cv$fits)){
      ids <- cv$predictions[[f]]$metadata$unit_ids
      if(mode == "composition"){
        total <- apply(cv$predictions[[f]]$count, c(1,3,4), sum)
        expect_equal(as.numeric(total), rep(unname(rowSums(dat$counts)[ids]),16))
        expect_length(cv$totals[[f]], 0L)
      }else expect_equal(cv$totals[[f]]$count, cv$predictions[[f]]$count[,1,,])
    }
    if(mode == "composition") expect_error(plot(cv,level="total"), "conditioned on")
  }
  no_area <- component_data(area=FALSE)
  expect_s3_class(magCompCV(no_area,magnitude=NULL,control=component_control(),
    fold_id=folds,summary="quantile",probs=0.5), "magCompCV")
})

test_that("standalone intercept posteriors agree with numerical integration", {
  rows <- expand.grid(plot = paste0("p", 1:4), cls = c("a", "b"))
  rows$count <- c(1,4,2,3,5,2,1,4)
  rows$area <- 1
  dat <- magCompData(rows,unit=plot,category=cls,tree_count=count,design=fixedArea(area))
  grid <- seq(-6, 6, length.out=6001)
  # Evaluate the marginal likelihood and prior independently of augmentation.
  for(mode in c("magnitude", "probit", "logit")){
    is_mag <- mode == "magnitude"
    link <- if(mode == "logit") "logit" else "probit"
    fit <- magComp(dat,magnitude=if(is_mag) ~1 else NULL,
      composition=if(is_mag) NULL else ~1,
      priors=if(is_mag) magCompPriors(magnitude=list(coefficient_variance=4, size=NULL)) else
        magCompPriors(composition=list(intercept_variance=4)),
      control=magCompControl(n_samples=7000,n_burn=1000,nb_size=3,
        composition_link=link,chain_seeds=733,verbose=FALSE))
    log_density <- vapply(grid,function(theta){
      ll <- if(is_mag) sum(dnbinom(rowSums(dat$counts),mu=exp(theta),size=3,log=TRUE)) else
        sum(dbinom(dat$counts[,1],rowSums(dat$counts),
          prob=if(link=="probit") pnorm(theta) else plogis(theta),log=TRUE))
      ll + dnorm(theta,sd=2,log=TRUE)
    },numeric(1L))
    weight <- exp(log_density-max(log_density));weight <- weight/sum(weight)
    draw <- if(is_mag) as.numeric(samples(fit)$magnitude$coefficients) else
      as.numeric(samples(fit)$composition$blocks$cls$coefficients)
    reference_mean <- sum(grid*weight)
    reference_sd <- sqrt(sum((grid-reference_mean)^2*weight))
    expect_lt(abs(mean(draw)-reference_mean)/reference_sd, 0.12)
    expect_lt(abs(sd(draw)/reference_sd-1), 0.12)
  }
})

test_that("point-sampling predictions, PPC, and CV support all component choices", {
  rows <- component_data(nested = TRUE)$data
  rows$baf <- c(10, 20)[(as.integer(rows$plot) - 1L) %% 2L + 1L]
  dat <- magCompData(rows, unit = plot, categories = c(species, dbh),
                     tree_count = count, design = pointSampling(baf))
  # Each CV fold predicts with a BAF different from its training points.
  folds <- setNames(rep(1:2, 4), dat$unit_levels)
  for(mode in c("magnitude", "full", "composition")){
    mag <- if(mode == "composition") NULL else ~ x
    comp <- if(mode == "magnitude") NULL else
      magCompComposition(species ~ x, dbh | species ~ x)
    priors <- if(mode == "magnitude")
      magCompPriors(magnitude = list(size = c(2, 0.2))) else magCompPriors()
    fit <- magComp(dat, family = pointBasalArea(), magnitude = mag,
      composition = comp, priors = priors, control = component_control(nb_size = 3))
    type <- if(mode == "magnitude") c("magnitude", "count") else
      c("proportion", "density", "count")
    total <- if(mode == "composition") rowSums(dat$counts) else NULL
    newdata <- dat$unit_data
    newdata$plot <- paste0("new", seq_len(nrow(newdata)))
    newdata$baf <- rev(newdata$baf)
    for(source in c("fitted", "predict")){
      result <- if(source == "fitted")
        fitted(fit, type = type, target = "observation", summary = "samples") else
        predict(fit, newdata, type = type, total_count = unname(total),
          target = "observation", summary = "samples")
      baf <- if(source == "fitted") dat$baf else newdata$baf
      stand <- standSummary(result, data.frame(id = result$metadata$unit_ids, ha = 2),
      unit = id, represented_area = ha, type = "total", summary = "samples")
      # Area expands basal-area densities, not tallies or reciprocal BAF.
      oracle_total <- apply(sweep(result$count, 1L, 2 * baf, "*"), c(2, 3, 4), sum)
      expect_equal(as.numeric(stand$total), as.numeric(oracle_total))
      expect_equal(unname(stand$metadata$units), "basal area")
      expect_identical(stand$metadata$target, "observation")
      expect_identical(stand$metadata$count_conditioning, result$metadata$count_conditioning)
      expect_equal(as.numeric(result$metadata$exposure), 1 / as.numeric(baf))
      expect_true(all(result$count >= 0 & result$count == floor(result$count)))
      expect_equal(unname(result$metadata$units["count"]), "in-tree tally")
      if(mode == "magnitude"){
        expect_equal(result$magnitude, sweep(result$count[,1,,], 1L, baf, "*"))
      }else{
        expect_equal(result$density, sweep(result$count, 1L, baf, "*"))
        if(mode == "composition"){
          expect_equal(as.numeric(apply(result$count, c(1,3,4), sum)),
                       rep(unname(total), 16))
          expect_true(all(is.na(result$proportion[1L,,,])))
          expect_true(all(result$density[1L,,,] == 0))
        }
      }
    }
    ppc <- ppcSummary(fit)
    expect_equal(ppc$observed, if(mode == "composition") as.numeric(dat$counts) else
      as.numeric(rowSums(dat$counts)))
    expect_true(all(is.finite(ppc$mean)))
    cv_type <- if(mode == "magnitude") c("magnitude", "count") else c("density", "count")
    cv <- magCompCV(dat, family = pointBasalArea(), magnitude = mag,
      composition = comp, priors = priors, control = component_control(nb_size = 3),
      fold_id = folds, type = cv_type, target = "observation", summary = "samples")
    quantity <- cv_type[1L]
    expect_equal(cv$observed[[quantity]], cv$observed$count * dat$baf)
    for(result in cv$predictions){
      baf <- dat$baf[match(result$metadata$unit_ids, dat$unit_levels)]
      expect_equal(as.numeric(result$metadata$exposure), 1 / as.numeric(baf))
      count <- if(mode == "magnitude") result$count[,1,,] else result$count
      expect_equal(result[[quantity]], sweep(count, 1L, baf, "*"))
      if(mode == "composition")
        expect_equal(as.numeric(apply(result$count, c(1,3,4), sum)),
                     rep(unname(total[result$metadata$unit_ids]), 16))
    }
  }
})
