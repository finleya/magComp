sample_test_fit <- function(spatial = "gp", nested = FALSE, chains = 2L)
{
  rows <- expand.grid(plot = 1:6, species = c("a", "b", "c"), dbh = if(nested) 1:3 else 1)
  rows$x <- rep(c(-1, -0.6, -0.1, 0.3, 0.7, 1), length.out = nrow(rows))
  rows$y <- rep(c(0, 0.4, 0.1, 0.8, 0.5, 0.9), length.out = nrow(rows))
  rows$stand <- rep(c(1L, 1L, 1L, 2L, 2L, 3L), length.out = nrow(rows))
  rows$count <- 1L + (rows$plot + as.integer(rows$species) + rows$dbh) %% 5L
  rows$area <- 0.1
  dat <- do.call(magCompData, list(data = rows, unit = "plot",
    categories = if(nested) c("species", "dbh") else "species",
    tree_count = "count", design = fixedArea(area)))
  mag <- if(spatial == "gp") ~ x + iid(stand, name = "group") +
    x:iid(stand, name = "slope") + gp(x, y, name = "space") else
    ~ x + iid(stand, name = "group") + x:iid(stand, name = "slope") +
      nngp(x, y, name = "space", m = 3)
  comp <- if(spatial == "gp") ~ x + factor_iid(stand, K = 1, name = "group") +
    x:factor_iid(stand, K = 1, name = "slope") +
    factor_gp(x, y, K = 2, name = "space") else
    ~ x + factor_iid(stand, K = 1, name = "group") +
      x:factor_iid(stand, K = 1, name = "slope") +
      factor_nngp(x, y, K = 2, name = "space", m = 3)
  if(nested) comp <- magCompComposition(species ~ x,
    dbh | species ~ x + factor_iid(stand, K = 1, name = "group") +
      x:factor_iid(stand, K = 1, name = "slope") +
      factor_gp(x, y, K = 2, name = "space"))
  magComp(dat, magnitude = mag, composition = comp,
    control = magCompControl(n_samples = 12, n_burn = 4, n_chains = chains,
      chain_seeds = seq_len(chains), verbose = FALSE))
}

test_that("sample selection preserves pairing and pooled diagnostics retain chain identity", {
  fit <- sample_test_fit()
  raw <- samples(fit)
  selected <- samples(fit, start = 2, thin = 3, chains = c("chain2", "chain1"))
  expect_s3_class(selected, "magCompSamples")
  expect_equal(selected$magnitude$coefficients,
               raw$magnitude$coefficients[,c(2,5,8),c(2,1),drop = FALSE])
  expect_equal(selected$composition$blocks$species$terms$space$scores,
               raw$composition$blocks$species$terms$space$scores[,,c(2,5,8),c(2,1),drop = FALSE])
  meta <- attr(selected, "sample_info")
  expect_equal(meta$draw_ids, c(2,5,8))
  expect_equal(meta$iterations, c(6,9,12))
  expect_equal(meta$chain_ids, c("chain2", "chain1"))
  for(type in c("parameters", "effects", "combined")){
    args <- list(object = fit, type = type, start = 2, thin = 2, chains = 2:1)
    separate <- do.call(samples, args)
    pooled <- do.call(samples, c(args, list(combine_chains = TRUE)))
    expect_equal(summary(pooled), summary(separate))
    expect_equal(do.call(summary, args), summary(separate))
  }
  one <- samples(fit, start = 8, chains = 2)
  expect_equal(dim(one$magnitude$coefficients), c(2,1))
  expect_true(all(is.na(summary(one)$rhat)))
  expect_true(all(is.na(summary(one)$ess)))
  expect_equal(summary(fit, probs = 0.5)$q50, summary(raw)$q50)
  expect_setequal(unique(summary(fit, type = "effects")$parameter), c("effects", "scores"))
  expect_error(samples(fit, start = 9), "exceeds")
  expect_error(samples(fit, thin = 0), "positive integers")
  expect_error(samples(fit, chains = c(1,1)), "distinct existing")
  expect_error(samples(fit, chains = "missing"), "distinct existing")
  expect_error(samples(fit, group = "stand"), "only used")
  expect_error(samples(fit, type = "combined", group = "missing"), "fitted-unit column")
  expect_error(samples(fit, unknown = TRUE), "unknown")
  path <- tempfile(fileext = ".pdf")
  devices <- grDevices::dev.list()
  plot(samples(fit, combine_chains = TRUE), parameters = "variance", file = path)
  expect_identical(grDevices::dev.list(), devices)
  expect_gt(file.info(path)$size, 0)
  unlink(path)
})

test_that("combined coefficients reconstruct fitted linear predictors for GP and NNGP models", {
  for(spatial in c("gp", "nngp")){
    fit <- sample_test_fit(spatial)
    raw <- samples(fit)
    combined <- samples(fit, type = "combined")
    info <- attr(combined, "info")
    response <- fitted(fit, target = "expected", summary = "samples")
    x <- fit$data$unit_data$x
    get <- function(component, coefficient, unit, stick = "")
      combined[[which(info$component == component & info$coefficient == coefficient &
                        info$unit == unit & info$stick == stick)]]
    for(i in seq_along(fit$data$unit_levels)){
      unit <- fit$data$unit_levels[i]
      eta <- get("magnitude", "(Intercept)", unit) + x[i] * get("magnitude", "x", unit)
      expect_equal(unname(eta), unname(log(response$magnitude[i,,])), tolerance = 1e-10)
      remaining <- matrix(1, 8, 2)
      for(j in 1:2){
        stick <- rownames(raw$composition$blocks$species$coefficients)[j]
        eta <- get("composition", "(Intercept)", unit, stick) +
          x[i] * get("composition", "x", unit, stick)
        probability <- pnorm(eta)
        expect_equal(unname(remaining * probability), unname(response$proportion[i,j,,]),
                     tolerance = 1e-10)
        remaining <- remaining * (1 - probability)
      }
    }
    # Group means average matching unit draws before summarizing. Unequal group
    # sizes ensure the overall mean is over units, not over the three groups.
    for(group in c("stand", ".all")){
      grouped <- samples(fit, type = "combined", group = group)
      gi <- attr(grouped, "info")
      for(j in seq_len(nrow(gi))){
        keep_units <- if(group == ".all") fit$data$unit_levels else
          fit$data$unit_levels[fit$data$unit_data$stand == gi$group[j]]
        keep <- info$component == gi$component[j] & info$index == gi$index[j] &
          info$unit %in% keep_units
        oracle <- Reduce(`+`, combined[keep]) / sum(keep)
        expect_equal(grouped[[j]], oracle, tolerance = 1e-12)
      }
    }
  }
})

test_that("conditional blocks keep their parent-specific loadings and coefficients", {
  fit <- sample_test_fit(nested = TRUE, chains = 1L)
  raw <- samples(fit)
  combined <- samples(fit, type = "combined")
  info <- attr(combined, "info")
  block <- fit$data$allocation_blocks[["dbh|species"]]
  part <- raw$composition$blocks[["dbh|species"]]
  terms <- fit$composition$allocation_components[["dbh|species"]]$process_terms
  for(j in seq_len(nrow(part$coefficients))){
    parent <- (j - 1L) %/% 2L + 1L
    for(i in 1:6) for(coefficient in c("(Intercept)", "x")){
      block_row <- which(block$unit_index == i & block$parent_index == parent)
      oracle <- part$coefficients[j,coefficient,]
      for(term in terms){
        if((if(is.null(term$coef_name)) "(Intercept)" else term$coef_name) != coefficient) next
        t <- part$terms[[term$name]]
        g <- if(term$term_type == "factor_iid") term$group[block_row] else i
        for(k in seq_len(ncol(t$scores))) oracle <- oracle + t$loadings[j,k,] * t$scores[g,k,]
      }
      row <- which(info$block == "dbh|species" & info$stick == rownames(part$coefficients)[j] &
                     info$unit == fit$data$unit_levels[i] & info$coefficient == coefficient)
      expect_equal(as.numeric(combined[[row]]), as.numeric(oracle), tolerance = 1e-12)
      expect_equal(info$parent[row], block$parent_levels[parent])
    }
  }
  expect_equal(summary(fit, type = "combined", group = "stand"),
               summary(samples(fit, type = "combined", group = "stand", combine_chains = TRUE)))
})

test_that("prediction extraction preserves draws and reconstructs new-location predictors", {
  for(case in c("gp", "nngp", "child")){
    fit <- sample_test_fit(if(case == "nngp") "nngp" else "gp",
                           nested = case == "child", chains = if(case == "child") 1L else 2L)
    newdata <- fit$data$unit_data[c(1, 2, 1, 2), , drop = FALSE]
    newdata$plot <- 11:14
    newdata$x[3:4] <- c(1.4, 1.4)
    newdata$y[3:4] <- c(0.7, 0.7)
    newdata$stand[3:4] <- 4L
    set.seed(1)
    ordinary <- predict(fit, newdata, summary = "samples", n_draws = 4)
    set.seed(1)
    predicted <- predict(fit, newdata, summary = "samples", n_draws = 4,
                          keep_effects = TRUE)
    expect_identical(predicted[names(ordinary)], ordinary |> unclass())
    raw <- samples(predicted)
    combined <- samples(predicted, type = "combined")
    info <- attr(combined, "info")
    ns <- predicted$metadata$n_samples
    nc <- predicted$metadata$n_chains
    get <- function(component, coefficient, unit, block = "", stick = ""){
      index <- which(info$component == component & info$coefficient == coefficient &
                     info$unit == unit & info$block == block & info$stick == stick)
      if(length(index)) return(combined[[index]])
      if(component == "magnitude") return(matrix(raw$magnitude$coefficients[coefficient,,drop = FALSE], ns, nc))
      beta <- raw$composition$blocks[[block]]$coefficients
      matrix(if(nc == 1) beta[stick,coefficient,] else beta[stick,coefficient,,], ns, nc)
    }
    for(i in 1:4){
      unit <- predicted$metadata$unit_ids[i]
      eta <- get("magnitude", "(Intercept)", unit) + newdata$x[i] * get("magnitude", "x", unit)
      expect_equal(as.numeric(exp(eta)), as.numeric(matrix(predicted$magnitude, nrow = 4)[i,]), tolerance = 1e-10)
      probabilities <- list()
      for(b in names(raw$composition$blocks)){
        sticks <- dimnames(raw$composition$blocks[[b]]$coefficients)[[1L]]
        nstick <- ncol(fit$data$allocation_blocks[[b]]$counts) - 1L
        p <- list()
        for(parent in seq_len(length(sticks) / nstick)){
          remaining <- matrix(1, ns, nc)
          for(j in seq_len(nstick)){
            stick <- sticks[(parent - 1L) * nstick + j]
            eta <- get("composition", "(Intercept)", unit, b, stick) +
              newdata$x[i] * get("composition", "x", unit, b, stick)
            probability <- pnorm(eta)
            p[[length(p) + 1L]] <- remaining * probability
            remaining <- remaining * (1 - probability)
          }
          p[[length(p) + 1L]] <- remaining
        }
        probabilities[[b]] <- p
      }
      p <- probabilities[[1L]]
      if(length(probabilities) == 2L){
        child <- probabilities[[2L]]
        p <- lapply(seq_along(child), function(j) p[[(j - 1L) %/% 3L + 1L]] * child[[j]])
      }
      for(j in seq_along(p)){
        actual <- if(nc == 1) predicted$proportion[i,j,] else predicted$proportion[i,j,,]
        expect_equal(as.numeric(p[[j]]), as.numeric(actual), tolerance = 1e-10)
      }
    }
    expect_equal(summary(samples(predicted, type = "combined", group = "stand",
                                  combine_chains = TRUE)),
                 summary(samples(predicted, type = "combined", group = "stand")))
    selected <- samples(predicted, start = 2, thin = 2)
    expect_equal(attr(selected, "sample_info")$draw_ids, predicted$metadata$sample_ids[c(2,4)])
    expect_equal(attr(selected, "sample_info")$iterations,
                 fit$control$n_burn + predicted$metadata$sample_ids[c(2,4)] * fit$control$n_thin)
    expect_error(samples(ordinary), "keep_effects")
    expect_error(predict(fit, newdata, keep_effects = TRUE), "summary")
    expect_error(predict(fit, newdata, keep_effects = NA), "TRUE or FALSE")
    if(case == "gp"){
      set.seed(1)
      ordinary_counts <- predict(fit, newdata, target = "observation", summary = "samples")
      set.seed(1)
      retained_counts <- predict(fit, newdata, target = "observation", summary = "samples",
                                  keep_effects = TRUE)
      expect_identical(retained_counts[names(ordinary_counts)], unclass(ordinary_counts))
      expect_equal(summary(samples(retained_counts, type = "effects"))$mean,
                   summary(samples(retained_counts, type = "effects", combine_chains = TRUE))$mean)
    }
  }
})
