# Independent regression calculations for checking that model predictions
# include additional random/spatial contributions. These are test references,
# not a second prediction interface.
mc_regression_magnitude_samples <- function(object,
                                            newdata = object$data$unit_data)
{
  X <- model.matrix(object$magnitude$fixed_formula, newdata)
  beta <- object$backend$magnitude_coefficient_samples
  out <- array(exp(X %*% matrix(beta, nrow = ncol(X))),
               dim = c(nrow(X), dim(beta)[-1L]))
  dimnames(out) <- c(list(rownames(newdata)), dimnames(beta)[-1L])
  out
}

mc_regression_magnitude_mean <- function(object)
{
  matrix(apply(mc_regression_magnitude_samples(object), 1L, mean), ncol = 1L,
         dimnames = list(object$data$unit_levels, "total"))
}

mc_regression_proportion_mean <- function(object)
{
  stopifnot(length(object$data$allocation_blocks) == 1L)
  component <- object$composition$allocation_components[[1L]]
  X <- component$X
  beta <- object$backend$composition_block_coefficient_samples[[1L]]
  d <- dim(beta)
  n_chain <- if(length(d) == 4L) d[4L] else 1L
  dim(beta) <- c(d[1:3], n_chain)
  out <- matrix(0, nrow(X), d[1L] + 1L,
                dimnames = list(object$data$unit_levels,
                                colnames(object$data$counts)))
  inverse_link <- if(object$control$composition_link == "probit") pnorm else plogis
  for(ch in seq_len(n_chain)){
    for(s in seq_len(d[3L])){
      coef <- matrix(beta[,,s,ch], nrow = d[1L], ncol = d[2L])
      stick <- inverse_link(X %*% t(coef))
      remaining <- rep(1, nrow(X))
      for(j in seq_len(d[1L])){
        out[,j] <- out[,j] + remaining * stick[,j]
        remaining <- remaining * (1 - stick[,j])
      }
      out[,d[1L] + 1L] <- out[,d[1L] + 1L] + remaining
    }
  }
  out / (d[3L] * n_chain)
}
