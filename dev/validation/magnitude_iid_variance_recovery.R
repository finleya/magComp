# Run with a freshly installed package, e.g. R_LIBS=/path/to/library Rscript <file>.
# A focused recovery experiment, not a substitute for repeated calibration.
dev_lib <- Sys.getenv("MAGCOMP_DEV_LIB", "")
if(nzchar(dev_lib)) .libPaths(c(dev_lib, .libPaths()))
library(magComp)
set.seed(1940)
n_stand <- 60L
n_block <- 30L
n <- n_stand * 6L
truth <- c(iid_1 = 0.16, iid_2 = 0.81)
stand <- rep(seq_len(n_stand), each = 6L)
block <- sample(rep(seq_len(n_block), length.out = n))
x <- rep(seq(-1, 1, length.out = 6L), n_stand)
intercept <- rnorm(n_stand, sd = sqrt(truth[1L]))
slope <- rnorm(n_block, sd = sqrt(truth[2L]))
mu <- exp(log(80) + 0.2 * x + intercept[stand] + slope[block] * x)
total <- rnbinom(n, mu = mu, size = 40)
first <- rbinom(n, total, 0.6)
rows <- data.frame(plot = rep(paste0("p", seq_len(n)), 2),
                    cls = factor(rep(c("a", "b"), each = n)),
                    count = c(first, total - first), area = 1,
                    stand = rep(paste0("s", stand), 2),
                    block = rep(paste0("b", block), 2), x = rep(x, 2))
dat <- magCompData(rows, unit = plot, category = cls, tree_count = count,
                    design = fixedArea(area))
fit <- magComp(dat, magnitude = ~ x + iid(stand) + x:iid(block),
                 control = magCompControl(n_samples = 3500, n_burn = 1500,
                                           n_chains = 2, chain_seeds = c(1941, 1942),
                                           nb_size = 40, verbose = FALSE))
samples <- samples(fit)
summary <- do.call(rbind, lapply(names(truth), function(term){
  draws <- as.numeric(samples$magnitude$terms[[term]]$variance)
  data.frame(term = term, truth = truth[[term]], mean = mean(draws),
               lower = unname(quantile(draws, 0.025)),
               upper = unname(quantile(draws, 0.975)))
}))
print(summary, row.names = FALSE)
diagnostics <- subset(summary(fit), parameter == "variance")
print(diagnostics)
stopifnot(all(summary$lower < summary$truth & summary$truth < summary$upper),
            summary$mean[2L] > 2 * summary$mean[1L],
            all(diagnostics$rhat < 1.05))
cat("Both distinct iid variances recovered; variance R-hat values below 1.05.\n")
