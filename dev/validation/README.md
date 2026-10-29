# Numerical and simulation validation

These scripts supplement the routine regression tests in `tests/testthat/`.
Run them from the package root against an installed magComp package.
`MAGCOMP_DEV_LIB` can select an alternative R library. Generated results
are written under `dev/validation/outputs/`, which is excluded from Git.

| Script | Purpose |
|:--|:--|
| `covariance_updates.R` | Compare covariance updates with known uniform posterior targets, including transformation Jacobians and adaptation. |
| `covariance_family_recovery.R` | Simulate dense-GP populations and estimate Matérn, multiresolution, and Gneiting models using NNGPs. |
| `space_time_recovery.R` | Estimate a space-time model and predict held-out revisits. |
| `nngp_gp_prediction.R` | Compare implemented predictions with independent joint Gaussian calculations, including GP and NNGP cases. |
| `factor_term_recovery_battery.R` | Exercise iid, GP, and NNGP factor terms, conditional allocation, and alternative starts. |
| `magnitude_iid_variance_recovery.R` | Estimate separate variances for magnitude iid terms from simulated data. |
| `positive_diagonal_reference.R` | Compare the posterior of a positive loading with an independent quadrature calculation. |
| `positive_diagonal_recovery.R` | Estimate positive anchor loadings in simulated leading and conditional composition models. |
| `native_memory_cleanup.R` | Check native prediction error handling, RNG state, and memory cleanup. |

For example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript dev/validation/nngp_gp_prediction.R
```

Using one BLAS/OpenMP thread per worker avoids nested threading when a
script runs multiple chains or simulations concurrently. Simulation checks
can require considerably longer runs than the routine test suite.

The factor battery defaults to a short `smoke` profile for execution and
constraint checks. Set `MAGCOMP_FACTOR_BATTERY_PROFILE=standard` for its
larger simulations. Positive-diagonal experiments expose settings such as
`MAGCOMP_DIAG_BURN`, `MAGCOMP_DIAG_DRAWS`, `MAGCOMP_DIAG_REPS`, and
`MAGCOMP_DIAG_CASES`; see the script for their defaults. Short execution
checks do not establish convergence or interval calibration.

The memory-cleanup script accepts `MAGCOMP_MEMORY_MODE=normal` or `error`.
It can also be run under Valgrind:

```sh
MAGCOMP_CHECK_RSS=false R -d "valgrind --leak-check=full --error-exitcode=1" --vanilla -f dev/validation/native_memory_cleanup.R
```

These are focused checks of implementation and estimation behavior. Their
simulation results should be interpreted with their chain diagnostics and
the number of simulated populations.
