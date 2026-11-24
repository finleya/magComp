shared_input_data <- function()
{
  set.seed(1)
  plots <- data.frame(id = 1:12, site = rep(1:6, each = 2),
                      x = seq(-1, 1, length.out = 12), area = rep(c(0.1, 0.2), 6))
  plots$height <- c(0.2, 0.7, 0.5, 0.2, 1, 0.7)[plots$site]
  trees <- merge(plots, expand.grid(A = c("a", "b", "c"), B = c("s", "m", "l")))
  trees$count <- rpois(nrow(trees), 3)
  magCompData(trees, unit = id, categories = c(A, B), tree_count = count,
              design = fixedArea(area))
}

test_that("shared inputs propagate through factors, slopes, restarts, and prediction", {
  data <- shared_input_data()
  for(method in c("gp", "nngp")){
    mag <- as.formula(paste0("~ x + x:", method, "(height, name='v')"))
    parent <- as.formula(paste0("A ~ x + x:factor_", method, "(height, K=2, name='w')"))
    child <- as.formula(paste0("B | A ~ x + x:factor_", method, "(height, K=1, name='u')"))
    comp <- magCompComposition(parent, child)
    fit <- magComp(data, magnitude = mag, composition = comp,
      priors = magCompPriors(magnitude = list(size = c(shape=2, rate=0.1))),
      control = magCompControl(n_samples=60, n_burn=20, n_chains=2,
                               chain_seeds=1:2, verbose=FALSE))
    raw <- samples(fit)
    map <- fit$magnitude$process_terms[[1]]$spatial_group
    expect_equal(dim(raw$magnitude$terms$v$effects), c(4L,40L,2L))
    expect_gt(length(unique(as.numeric(raw$magnitude$size))), 1)
    for(b in c("A", "B|A")){
      term <- fit$composition$allocation_components[[b]]$process_terms[[1]]
      expect_equal(term$n_group, 4L)
      expect_equal(nrow(raw$composition$blocks[[b]]$terms[[1]]$scores), 4L)
    }
    combined <- samples(fit, type="combined")
    info <- attr(combined, "info")
    for(i in seq_len(data$n_unit)){
      which_draw <- which(info$component == "magnitude" & info$unit == data$unit_levels[i])
      expect_equal(unname(combined[[which_draw]]), unname(raw$magnitude$coefficients["x",,] +
        raw$magnitude$terms$v$effects[map[i],,]))
    }
    # Check a child coefficient against its loading times the mapped factor.
    term <- fit$composition$allocation_components[["B|A"]]$process_terms[[1]]
    block <- data$allocation_blocks[["B|A"]]
    part <- raw$composition$blocks[["B|A"]]
    for(i in seq_len(data$n_unit)){
      row <- which(block$unit_index == i & block$parent_index == 1L)
      g <- term$spatial_group[row]
      index <- which(info$block == "B|A" & info$index == paste0(rownames(part$coefficients)[1],":x") &
                       info$unit == data$unit_levels[i])
      expected <- part$coefficients[1,"x",,] + part$terms$u$loadings[1,1,,] * part$terms$u$scores[g,1,,]
      expect_equal(unname(combined[[index]]), unname(expected))
    }
    last <- magCompStarting(fit, use="last", chain=2)
    prepared <- magComp(data, magnitude=mag, composition=comp, starting=last,
      priors=magCompSettings(fit)$priors, control=magCompControl(run_sampler=FALSE))
    expect_equal(prepared$starting, last)
    expect_true(all(is.finite(summary(fit)$mean)))
    expect_true(all(is.finite(summary(fit, type="combined")$mean)))
    fitted_draws <- fitted(fit, summary="samples")
    predicted <- predict(fit, newdata=data$unit_data, summary="samples", keep_effects=TRUE)
    expect_equal(predicted$magnitude, fitted_draws$magnitude)
    expect_equal(predicted$proportion, fitted_draws$proportion)
    expect_equal(samples(predicted, type="combined"), combined)
    expected <- raw$magnitude$terms$v$effects[map,,,drop=FALSE]
    for(i in seq_len(data$n_unit)) expected[i,,] <- exp(raw$magnitude$coefficients["(Intercept)",,] +
      data$unit_data$x[i] * (raw$magnitude$coefficients["x",,] + expected[i,,]))
    expect_equal(as.numeric(fitted_draws$magnitude), as.numeric(expected))
    # Identical new inputs share latent values; exposures and response draws
    # remain observation-specific. Include a fitted input and a new input.
    new <- data.frame(id=21:24, height=c(0.2,0.2,0.6,0.6), x=0.5, area=c(.1,.2,.1,.2))
    pp <- predict(fit, newdata=new, summary="samples")
    expect_equal(pp$magnitude[1,,], pp$magnitude[2,,])
    expect_equal(pp$magnitude[3,,], pp$magnitude[4,,])
    expect_equal(pp$proportion[1,,,], pp$proportion[2,,,])
    expect_equal(pp$proportion[3,,,], pp$proportion[4,,,])
    expect_equal(pp$count[2,,,], 2*pp$count[1,,,])
    obs <- predict(fit, newdata=new, target="observation", summary="samples")
    expect_true(any(obs$count[1,,,] != obs$count[2,,,]))
    expect_equal(as.numeric(apply(obs$count, c(1,3,4), sum)),
                  as.numeric(obs$magnitude * new$area))
  }
})

test_that("shared-node magnitude posterior agrees with numerical integration", {
  # A single latent value, several exposures and varying multipliers. With
  # IG(2,1) variance its marginal prior is t_4 with scale sqrt(1/2).
  # A tightly concentrated intercept prior makes its contribution negligible.
  plots <- data.frame(id=1:6, height=1, x=c(-1, .5, 1, 1.5, 2, 2),
                      area=c(1,1,2,2,1,3), count=c(0,1,5,7,4,13))
  data <- magCompData(plots, unit=id, tree_count=count, design=fixedArea(area))
  grid <- seq(-5,5,length.out=20001)
  log_density <- vapply(grid, function(f) sum(dnbinom(plots$count, size=10,
    mu=plots$area*exp(plots$x*f), log=TRUE)) + dt(f/sqrt(.5), df=4, log=TRUE), 0.0)
  weight <- exp(log_density-max(log_density))
  weight <- weight/sum(weight)
  mean_f <- sum(grid*weight)
  mean_variance <- sum((1 + grid^2/2)/1.5 * weight)
  for(method in c("gp", "nngp")){
    fit <- magComp(data, composition=NULL,
      magnitude=as.formula(paste0("~ x:",method,"(height, name='v')")),
      priors=magCompPriors(magnitude=list(coefficient_variance=1e-12)),
      starting=list(magnitude=list(terms=list(v=list(covariance=list(phi=fixed(1)))))),
      control=magCompControl(n_samples=6000,n_burn=1000,nb_size=10,chain_seeds=1,verbose=FALSE))
    draws <- samples(fit)$magnitude$terms$v
    expect_equal(nrow(draws$effects), 1L)
    expect_equal(mean(draws$effects), mean_f, tolerance=0.04)
    expect_equal(mean(draws$variance), mean_variance, tolerance=0.07)
  }
})

test_that("CV preserves shared-node starts when the representative observation is held out", {
  data <- shared_input_data()
  mag <- ~ x + x:nngp(height, name="v")
  comp <- magCompComposition(A ~ factor_gp(height, name="w", anchors="a"),
                            B | A ~ x:factor_nngp(height, name="u", anchors="a:s"))
  model <- magComp(data, magnitude=mag, composition=comp,
                   control=magCompControl(run_sampler=FALSE))
  starts <- magCompStarting(model)
  # These folds retain at least one observation at every height, while changing
  # which observation supplies a node's label. A second split omits whole sites.
  for(labels in list(rep(1:2,6), rep(c(1,1,2,2),3))){
    folds <- setNames(labels, data$unit_levels)
    cv <- magCompCV(data, magnitude=mag, composition=comp, starting=starts,
      fold_id=folds, keep_fits=TRUE,
      control=magCompControl(n_samples=6,n_burn=2,chain_seeds=1,verbose=FALSE))
    for(fit in cv$fits){
      term <- fit$magnitude$process_terms[[1]]
      old <- match(term$coords[,1], model$magnitude$process_terms[[1]]$coords[,1])
      expect_equal(unname(fit$starting$magnitude$terms$v$effects),
                    unname(starts$magnitude$terms$v$effects[old]))
      expect_equal(names(fit$starting$magnitude$terms$v$effects), term$location_levels)
      for(b in names(comp)){
        term <- fit$composition$allocation_components[[b]]$process_terms[[1]]
        original <- model$composition$allocation_components[[b]]$process_terms[[1]]
        old <- match(term$coords[,1], original$coords[,1])
        expect_equal(unname(fit$starting$composition$blocks[[b]]$terms[[1]]$scores),
          unname(starts$composition$blocks[[b]]$terms[[1]]$scores[old,,drop=FALSE]))
      }
    }
  }
})
