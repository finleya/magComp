# Deliberately unequal densities and areas distinguish pooled composition from
# an area-weighted mean of local proportions. Two chains test draw alignment.
stand_prediction <- function(target = "observation")
{
  density <- array(c(10, 0, 3, 0, 40, 1,
                     0, 0, 0, 0, 0, 0,
                     20, 0, 4, 0, 20, 2,
                     30, 10, 2, 10, 30, 4), c(3, 2, 2, 2),
    dimnames = list(c("p2", "p1", "p3"), c("oak", "pine"),
                    c("draw7", "draw9"), c("chain2", "chain1")))
  structure(list(density = density, metadata = list(
    source = "predict", target = target, type = "density", summary = "samples",
    family = "fixed_area_stems", units = c(density = "stems per unit area"),
    unit_ids = dimnames(density)[[1]], exposure = c(p2 = 0.05, p1 = 0.1, p3 = 0.1),
    components = c("magnitude", "composition"), count_conditioning = "none",
    n_samples = 2L, n_chains = 2L, sample_ids = c(7L, 9L),
    chain_ids = c("chain2", "chain1"), n_draws = 4L)),
    class = c("magCompPrediction", "list"))
}

test_that("stand summaries preserve shared cell draws, weights and identities", {
  # p2 spans both stands. The rows also differ from prediction location order.
  areas <- data.frame(cell = c("p1", "p2", "p3", "p2"),
                      stand = c("B", "B", "A", "A"), ha = c(3, 0.4, 2, 0.6))
  for(target in c("expected", "observation")){
    pred <- stand_prediction(target)
    original <- pred
    total <- standSummary(pred, areas, unit = cell, group = stand,
                          represented_area = ha, type = "total", summary = "samples")
    # Independent arithmetic: the same p2 draw contributes to both stands.
    oracle <- pred$density[1:2,,,, drop = FALSE]
    dimnames(oracle)[[1]] <- c("B", "A")
    oracle[1,,,] <- 0.4 * pred$density[1,,,] + 3 * pred$density[2,,,]
    oracle[2,,,] <- 0.6 * pred$density[1,,,] + 2 * pred$density[3,,,]
    expect_equal(total$total, oracle)
    expect_identical(pred, original)
    expect_identical(total$metadata$target, target)
    expect_identical(total$metadata$chain_ids, pred$metadata$chain_ids)
    expect_identical(total$metadata$sample_ids, pred$metadata$sample_ids)
    expect_equal(total$metadata$represented_area, c(B = 3.4, A = 2.6))
    expect_equal(total$metadata$aggregation$exposure, pred$metadata$exposure)
    expect_equal(total$metadata$aggregation$data$unit, areas$cell)

    density <- standSummary(pred, areas, cell, stand, ha, summary = "samples")
    expected_density <- sweep(oracle, 1, c(3.4, 2.6), "/")
    expect_equal(density$density, expected_density)
    forest <- standSummary(pred, areas, cell, represented_area = ha, summary = "samples")
    expect_equal(as.numeric(forest$density),
                 as.numeric(apply(oracle, c(2, 3, 4), sum) / 6))
    expect_equal(standSummary(pred, areas, cell, stand, ha)$density,
                 apply(expected_density, c(1, 2), mean))
    q <- standSummary(pred, areas, cell, stand, ha, summary = "quantile", probs = 0.5)
    expect_equal(dim(q$density), c(2L, 2L, 1L))
    expect_equal(as.numeric(q$density),
                 as.numeric(apply(expected_density, c(1, 2), median)))

    reordered <- standSummary(pred, areas[c(2, 1, 4, 3), ], unit = "cell",
      group = "stand", represented_area = "ha", summary = "samples")
    expect_equal(reordered$density, density$density)
    subset <- standSummary(pred, areas[3, ], cell, stand, ha, summary = "samples")
    expect_equal(unname(subset$density), unname(pred$density[3,,,, drop = FALSE]))
    # Disconnected pieces of one cell in the same stand add their areas.
    pieces <- rbind(transform(areas, ha = ha / 2), transform(areas, ha = ha / 2))
    expect_equal(standSummary(pred, pieces, cell, stand, ha, summary = "samples")$density,
                 density$density)
  }
})

test_that("stand proportions pool abundance and condition on positive totals", {
  areas <- data.frame(cell = c("p2", "p1", "p3"), stand = c("B", "B", "A"),
                      ha = c(1, 3, 2))
  for(target in c("expected", "observation")){
    pred <- stand_prediction(target)
    prop <- standSummary(pred, areas, cell, stand, ha,
                         type = "proportion", summary = "samples")
    expect_equal(unname(prop$proportion[1,,1,1]), c(10, 120) / 130)
    expect_true(all(is.na(prop$proportion[,,2,1])))
    expect_equal(prop$metadata$zero_total_probability, c(B = 0.25, A = 0.25))
    expect_equal(prop$metadata$n_positive_draws, c(B = 3L, A = 3L))
    tab <- summary(prop)
    expect_true(all(tab$n_draws == 3))
    expect_equal(tab$mean, as.numeric(apply(prop$proportion, c(1, 2), mean, na.rm = TRUE)))
    expect_equal(tab$q50, as.numeric(apply(prop$proportion, c(1, 2), median, na.rm = TRUE)))
  }
  pred$density[] <- 0
  empty <- standSummary(pred, areas, cell, represented_area = ha,
                        type = "proportion", summary = "samples")
  expect_true(all(is.na(summary(empty)$mean)))
  expect_true(all(summary(empty)$n_draws == 0))
  expect_true(all(standSummary(pred, areas, cell, represented_area = ha)$density == 0))
  expect_true(all(standSummary(pred, areas, cell, represented_area = ha, type = "total")$total == 0))
})

test_that("stand summaries validate columns, IDs, areas and sampled input", {
  pred <- stand_prediction()
  areas <- data.frame(cell = pred$metadata$unit_ids, stand = "A", ha = 1)
  expect_error(standSummary(list()), "output from fitted.*predict")
  mean_only <- pred
  mean_only$metadata$summary <- "mean"
  expect_error(standSummary(mean_only, areas, cell, represented_area = ha), "requires retained draws")
  proportions <- pred
  proportions$density <- NULL
  proportions$proportion <- pred$density
  expect_error(standSummary(proportions, areas, cell, represented_area = ha), "proportions alone")
  expect_error(standSummary(pred), "data must be a nonempty data frame")
  expect_error(standSummary(pred, areas[FALSE, ], cell, represented_area = ha), "nonempty data frame")
  expect_error(standSummary(pred, areas), "unit is required")
  expect_error(standSummary(pred, areas, cell), "represented_area is required")
  expect_error(standSummary(pred, areas, absent, represented_area = ha), "unit column not found")
  expect_error(standSummary(pred, areas, cell, represented_area = "absent"), "represented_area column not found")
  expect_error(standSummary(pred, areas, c(cell, stand), represented_area = ha), "identify one column")
  expect_error(standSummary(pred, areas, cell, represented_area = 1), "variable name")
  expect_error(standSummary(pred, transform(areas, ha = "1"), cell, represented_area = ha), "must be numeric")
  for(area in list(0, -1, NA_real_, Inf))
    expect_error(standSummary(pred, transform(areas, ha = area), cell, represented_area = ha), "positive finite")
  expect_error(standSummary(pred, transform(areas, cell = "unknown"), cell, represented_area = ha), "must match prediction unit IDs")
  expect_error(standSummary(pred, transform(areas, cell = NA_character_), cell, represented_area = ha), "missing IDs")
  expect_error(standSummary(pred, transform(areas, cell = ""), cell, represented_area = ha), "empty IDs")
  expect_error(standSummary(pred, transform(areas, stand = NA_character_), cell, stand, ha), "missing or empty")
  expect_error(standSummary(pred, areas, cell, represented_area = ha, target = "expected"), "unused standSummary argument")
  expect_error(standSummary(pred, areas, cell, represented_area = ha, type = "count"), "arg.*should be")
})

test_that("single units and single-chain magnitude samples retain dimensions", {
  pred <- stand_prediction()
  single <- pred
  single$density <- pred$density[1,1,1,1, drop = FALSE]
  single$metadata$unit_ids <- "p2"
  single$metadata$n_samples <- single$metadata$n_chains <- single$metadata$n_draws <- 1L
  areas <- data.frame(cell = "p2", ha = 2)
  expect_equal(dim(standSummary(single, areas, cell, represented_area = ha,
                               summary = "samples")$density), c(1L, 1L, 1L, 1L))
  single$magnitude <- array(10, c(1, 1), dimnames = list("p2", "draw7"))
  single$density <- NULL
  expect_equal(standSummary(single, areas, cell, represented_area = ha, type = "total")$total,
               matrix(20, 1, 1, dimnames = list("all", "total")))
  expect_error(standSummary(single, areas, cell, represented_area = ha, type = "proportion"), "category density")
  # Numeric IDs are matched by value, and data-frame subclasses are accepted.
  single$metadata$unit_ids <- "2"
  dimnames(single$magnitude)[[1]] <- "2"
  areas <- structure(data.frame(cell = 2, ha = 2), class = c("tbl_df", "tbl", "data.frame"))
  expect_equal(as.numeric(standSummary(single, areas, cell, represented_area = ha)$density), 10)
})
