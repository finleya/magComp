test_that("positive diagonal priors reach leading and child iid/GP/NNGP terms", {
  rows <- expand.grid(plot = paste0("p", 1:8), species = c("a", "b", "c", "d"),
                       dbh = c("small", "large"))
  i <- as.integer(rows$plot)
  rows$count <- as.integer(i == 1L & rows$species == "a" & rows$dbh == "small")
  rows$area <- 1
  rows$group <- paste0("g", (i - 1L) %/% 2L + 1L)
  rows$x <- i / 8
  rows$y <- (i %% 3) / 3
  dat <- magCompData(rows, unit = plot, categories = c(species, dbh),
                      tree_count = count, design = fixedArea(area))
  priors <- magCompPriors(composition = list(loading_variance = 0.5,
    blocks = list(species = list(loading_variance = 1.5,
      terms = list(shape = list(loading_variance = 0.3))),
      "dbh|species" = list(loading_variance = 3))))
  for(kind in c("iid", "gp", "nngp")) for(link in c("probit", "logit")) {
    term <- paste0("factor_", kind, "(", if(kind == "iid") "group" else "x, y",
                   ", K = 2, name = 'shape', anchors = c('b:small', 'c:small')",
                   if(kind == "nngp") ", m = 3" else "",
                   ")")
    composition <- magCompComposition(
      species ~ factor_iid(group, K = 2, name = "shape", anchors = c("b", "c")),
      as.formula(paste("dbh | species ~", term)))
    fit <- magComp(dat, composition = composition, priors = priors,
      control = magCompControl(n_burn = 0, n_samples = 1600, nb_size = 20,
        chain_seeds = 619, composition_link = link, verbose = FALSE))
    samples <- samples(fit)$composition$blocks
    for(b in names(samples)) {
      loading <- samples[[b]]$terms$shape$loadings
      variance <- if(b == "species") 0.3 else 3
      for(k in 1:2) {
        d <- as.numeric(loading[k+1L,k,])
        expect_true(all(is.finite(d) & d > 0))
        # With zero trials the full conditional is exactly the half-normal
        # prior. Its CDF must turn these independent draws into uniforms.
        pit <- 2 * pnorm(d / sqrt(variance)) - 1
        expect_lt(abs(mean(pit) - 0.5), 6 / sqrt(12 * length(pit)))
        expect_lt(max(abs(sort(pit) - (seq_along(pit) - 0.5) / length(pit))), 0.06)
      }
      expect_true(all(loading[2,2,] == 0))
    }
    last <- magCompStarting(fit, use = "last")
    prepared <- magComp(dat, composition = composition, priors = priors, starting = last,
    control = magCompControl(run_sampler = FALSE))
    expect_equal(magCompStarting(prepared), last)
    expect_match(subset(magCompSettings(fit)$parameters,
                         parameter == "loadings")$prior[1], "restricted to")
  }
})
