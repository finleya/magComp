test_that("unit-level model inputs cannot silently vary within units", {
  rows <- data.frame(plot = rep(c("a", "b"), each = 2),
    cls = rep(c("s", "l"), 2), z = c(1, 2, 3, 3),
    x = c(0, 0, 1, 1), y = c(0, 0, 2, 2), group = c("g", "h", "h", "h"))
  dat <- magCompData(rows, unit = plot, category = cls, design = fixedArea(1))
  expect_error(magComp(dat, magnitude = ~ log(z),
    control = magCompControl(run_sampler = FALSE)), "unit-level column 'z'")
  expect_error(magComp(dat, composition = ~ z,
    control = magCompControl(run_sampler = FALSE)), "unit-level column 'z'")
  expect_error(magComp(dat, magnitude = ~ iid(group),
    control = magCompControl(run_sampler = FALSE)), "unit-level column 'group'")
  expect_error(magComp(dat, composition = ~ factor_iid(group),
    control = magCompControl(run_sampler = FALSE)), "unit-level column 'group'")
  expect_s3_class(magComp(dat,
    control = magCompControl(run_sampler = FALSE)), "magComp") # Unused tree-level columns may vary.
  rows$x[2] <- 2
  expect_error(magCompData(rows, unit = plot, category = cls, coords = c(x, y)),
               "unit-level column 'x'")
  dat <- magCompData(rows, unit = plot, category = cls, design = fixedArea(1))
  expect_error(magComp(dat, magnitude = ~ gp(x, y),
    control = magCompControl(run_sampler = FALSE)), "unit-level column 'x'")
})

test_that("ordered categories require unambiguous levels and honor named flags", {
  rows <- expand.grid(plot = c("a", "b"),
    cls = c("10-20", "20-30", "100-110"), species = c("oak", "pine"),
    stringsAsFactors = FALSE)
  expect_error(magCompData(rows, unit = plot, category = cls, ordered = TRUE),
               "explicitly ordered levels")
  rows$cls <- factor(rows$cls, levels = c("10-20", "20-30", "100-110"))
  dat <- magCompData(rows, unit = plot, categories = c(species, cls),
                     ordered = c(cls = TRUE, species = FALSE))
  expect_identical(dat$category_levels$cls, levels(rows$cls))
  expect_identical(dat$category_ordered, c(species = FALSE, cls = TRUE))
  rows$cls <- c(100, 20, 10)[as.integer(rows$cls)]
  dat <- magCompData(rows, unit = plot, category = cls, ordered = TRUE)
  expect_identical(dat$category_levels$cls, c("10", "20", "100"))
})

test_that("count overflow is rejected before integer conversion", {
  rows <- data.frame(plot = c("a", "a"), cls = c("s", "l"), count = c(2^31, 1))
  expect_error(magCompData(rows, unit = plot, category = cls, tree_count = count),
               "tree_count exceeds")
  rows$count <- c(2^30, 2^30)
  expect_error(magCompData(rows, unit = plot, category = cls, tree_count = count),
               "counts per sampling unit exceed")
})

test_that("magCompData builds fixed-area tree-list summaries", {
  trees <- data.frame(
    plot = c("a", "a", "b", "b", "b"),
    cls = c("small", "large", "small", "small", "large"),
    x = c(0, 0, 1, 1, 1),
    y = c(0, 0, 1, 1, 1),
    area = c(0.1, 0.1, 0.2, 0.2, 0.2),
    tf = c(10, 10, 5, 5, 5),
    volume = c(1, 4, 2, 2, 5)
  )

  dat <- magCompData(
    trees,
    unit = plot,
    category = cls,
    coords = c(x, y),
    tree_factor = tf,
    design = fixedArea(area),
    attributes = volume
  )

  expect_s3_class(dat, "magCompData")
  expect_equal(dat$n_unit, 2L)
  expect_equal(dim(dat$counts), c(2L, 2L))
  expect_equal(unname(rowSums(dat$counts)), c(2, 3))
  expect_equal(unname(dat$area), c(0.1, 0.2))
  expect_equal(unname(rowSums(dat$expanded$volume)), c(50, 45))
  expect_equal(dat$category_names, "cls")
  expect_equal(length(dat$allocation_blocks), 1L)
  expect_equal(dat$allocation_blocks[[1L]]$type, "marginal")
  expect_equal(dat$allocation_blocks[[1L]]$counts, dat$counts)
})

test_that("magCompData builds one- and two-category allocation blocks", {
  trees <- data.frame(
    plot = c("p1", "p1", "p1", "p2", "p2", "p3", "p3", "p3"),
    spp = c("hw", "hw", "sw", "hw", "sw", "sw", "sw", "hw"),
    dbh_class = c("small", "large", "small", "large",
                  "large", "small", "large", "small"),
    area = 0.1
  )

  trees$dbh_class <- factor(trees$dbh_class, levels = c("large", "small"))
  dat_one <- magCompData(
    trees,
    unit = plot,
    categories = c(dbh_class),
    design = fixedArea(area)
  )
  expect_equal(dat_one$category_names, "dbh_class")
  expect_equal(dim(dat_one$counts), c(3L, 2L))
  expect_equal(length(dat_one$allocation_blocks), 1L)
  expect_equal(dat_one$allocation_blocks[[1L]]$levels,
               colnames(dat_one$counts))

  dat_two <- magCompData(
    trees,
    unit = plot,
    categories = c(spp, dbh_class),
    ordered = c(spp = FALSE, dbh_class = TRUE),
    design = fixedArea(area)
  )
  expect_equal(dat_two$category_names, c("spp", "dbh_class"))
  expect_equal(dat_two$category_ordered, c(spp = FALSE, dbh_class = TRUE))
  expect_equal(dat_two$n_category, 2L)
  expect_equal(length(dat_two$allocation_blocks), 2L)
  expect_equal(names(dat_two$allocation_blocks), c("spp", "dbh_class|spp"))
  expect_equal(dat_two$allocation_blocks[[1L]]$type, "marginal")
  expect_false(dat_two$allocation_blocks[[1L]]$ordered)
  expect_equal(dat_two$allocation_blocks[[2L]]$type, "conditional")
  expect_true(dat_two$allocation_blocks[[2L]]$ordered)
  expect_equal(dat_two$allocation_blocks[[2L]]$parent, "spp")
  expect_equal(dat_two$allocation_blocks[[2L]]$child, "dbh_class")
  expect_equal(nrow(dat_two$allocation_blocks[[2L]]$counts),
               dat_two$n_unit * length(dat_two$category_levels$spp))
  expect_equal(dat_two$allocation_blocks[[2L]]$coef_stick_count,
               length(dat_two$category_levels$spp) *
                 (length(dat_two$category_levels$dbh_class) - 1L))
  expect_equal(dat_two$stick_names,
               c(dat_two$category_levels$spp[-length(dat_two$category_levels$spp)],
                 paste(dat_two$category_levels$spp,
                       dat_two$category_levels$dbh_class[-length(dat_two$category_levels$dbh_class)],
                       sep = ":")))
  expect_equal(colnames(dat_two$counts),
               levels(interaction(dat_two$category_factors$spp,
                                  dat_two$category_factors$dbh_class,
                                  drop = FALSE,
                                  sep = ":",
                                  lex.order = TRUE)))
  expect_equal(sum(dat_two$counts),
               sum(dat_two$allocation_blocks[[1L]]$counts))
  expect_equal(sum(dat_two$counts),
               sum(dat_two$allocation_blocks[[2L]]$counts))

  counted <- trees
  counted$n <- c(2L, 0L, 1L, 3L, 1L, 0L, 4L, 1L)
  dat_counted <- magCompData(
    counted,
    unit = plot,
    categories = c(spp, dbh_class),
    tree_count = n,
    ordered = c(spp = FALSE, dbh_class = TRUE),
    design = fixedArea(area)
  )
  expect_equal(sum(dat_counted$counts), sum(counted$n))
  expect_equal(unname(rowSums(dat_counted$counts)),
               as.numeric(stats::xtabs(n ~ plot, counted)))
  expect_equal(sum(dat_counted$counts),
               sum(dat_counted$allocation_blocks[[1L]]$counts))
  expect_equal(sum(dat_counted$counts),
               sum(dat_counted$allocation_blocks[[2L]]$counts))

  expect_error(
    magCompData(
      trees,
      unit = plot,
      category = spp,
      categories = c(spp, dbh_class),
      design = fixedArea(area)
    ),
    "category or categories"
  )
})

test_that("magCompData accepts unit-level fixed-area exposures", {
  trees <- data.frame(
    plot = c("a", "a", "b", "b", "b"),
    cls = c("small", "large", "small", "small", "large"),
    area = c(0.08, 0.08, 0.10, 0.10, 0.10)
  )

  area <- c(a = 9, b = 9)
  dat_column <- magCompData(
    trees,
    unit = plot,
    category = cls,
    design = fixedArea(area)
  )

  expect_equal(unname(dat_column$area), c(0.08, 0.10))

  dat_named <- magCompData(
    trees,
    unit = plot,
    category = cls,
    design = fixedArea(c(a = 0.08, b = 0.10))
  )

  expect_equal(unname(dat_named$area), c(0.08, 0.10))

  dat_scalar <- magCompData(
    trees,
    unit = plot,
    category = cls,
    design = fixedArea(c(plot_area = 0.08))
  )

  expect_equal(unname(dat_scalar$area), c(0.08, 0.08))

  trees_no_area <- trees[,c("plot", "cls")]
  area_by_unit <- c(a = 0.08, b = 0.10)
  dat_external <- magCompData(
    trees_no_area,
    unit = plot,
    category = cls,
    design = fixedArea(area_by_unit)
  )

  expect_equal(unname(dat_external$area), c(0.08, 0.10))

  area_unnamed <- c(0.08, 0.10)
  expect_error(
    magCompData(
      trees_no_area,
      unit = plot,
      category = cls,
      design = fixedArea(area_unnamed)
    ),
    "named by sampling unit"
  )
})

test_that("magCompData rejects category-specific fixed-area exposures", {
  trees <- data.frame(
    plot = c("a", "a", "b", "b", "b"),
    cls = c("small", "large", "small", "small", "large"),
    area = c(0.02, 0.08, 0.02, 0.02, 0.08)
  )

  expect_error(
    magCompData(
      trees,
      unit = plot,
      category = cls,
      design = fixedArea(area)
    ),
    "constant within each sampling unit"
  )

  expect_error(
    magCompData(
      trees,
      unit = plot,
      category = cls,
      design = fixedArea(c(large = 0.08, small = 0.02))
    ),
    "sampling-unit names"
  )

  area_matrix <- matrix(c(0.07, 0.09, 0.03, 0.04), nrow = 2)
  expect_error(fixedArea(area_matrix),
               "area must be scalar or have one value per sampling unit")
})

test_that("magCompData links point-sampling BAF by column or named unit vector", {
  trees <- data.frame(
    point = c("a", "a", "b", "b"),
    cls = c("small", "large", "small", "large"),
    baf = c(10, 10, 20, 20)
  )

  baf <- c(a = 99, b = 99)
  dat_column <- magCompData(
    trees,
    unit = point,
    category = cls,
    design = pointSampling(baf)
  )

  expect_equal(unname(dat_column$baf), c(10, 20))

  trees_no_baf <- trees[,c("point", "cls")]
  baf_by_unit <- c(a = 10, b = 20)
  dat_external <- magCompData(
    trees_no_baf,
    unit = point,
    category = cls,
    design = pointSampling(baf_by_unit)
  )

  expect_equal(unname(dat_external$baf), c(10, 20))

  baf_unnamed <- c(10, 20)
  expect_error(
    magCompData(
      trees_no_baf,
      unit = point,
      category = cls,
      design = pointSampling(baf_unnamed)
    ),
    "named by sampling unit"
  )
})
