# Sampling lengths for the five models. All priors and starts use defaults.
model_settings <- tribble(
  ~model, ~n_samples, ~n_burn,
  "H1",      50000L,  20000L,
  "H2",      50000L,  20000L,
  "H3",      50000L,  20000L,
  "H4",      50000L,  20000L,
  "H5",      50000L,  20000L
)

# Optional longer runs for fits that do not pass the chain diagnostics.
# Leave empty to record these fits as unresolved for review.
extension_settings <- tibble(n_samples = integer(), n_burn = integer())
n_chains <- 3L
chain_seeds <- 1:3
n_thin <- 10L
prediction_draws <- 333L
diagnostic_draws <- 1000L
rhat_limit <- 1.05
fixed_size <- 1000
