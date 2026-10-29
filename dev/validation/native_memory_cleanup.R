# Small native error-path regression; no simulation fits or output files.
# Run in a fresh R process. MAGCOMP_DEV_LIB selects an isolated installation.
lib <- Sys.getenv("MAGCOMP_DEV_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(magComp)

mode <- Sys.getenv("MAGCOMP_MEMORY_MODE", "error")
stopifnot(mode %in% c("normal", "error"))
effect <- matrix(0, 1, 1)
n_new <- 131072L # A 1 MiB C++ workspace; the returned R matrix is one element.
invoke <- function(bad = FALSE) {
  .Call("magcomp_predict_iid_scalar_cpp", effect, 0L,
    if (bad) 0L else 1L, n_new, 1.0, PACKAGE = "magComp")
}
measure <- function(calls) {
  heap <- gc()
  rss <- if (file.exists("/proc/self/status")) {
    line <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
    as.numeric(sub("^VmRSS:\\s+(\\d+).*", "\\1", line)) / 1024
  } else NA_real_
  data.frame(mode = mode, calls = calls, rss_MiB = rss,
             R_heap_MiB = sum(heap[, 2]))
}

# Warm up R and the native allocator before measuring cumulative growth.
for (i in seq_len(20L)) invisible(invoke())
set.seed(57)
invisible(invoke())
expected_seed <- .Random.seed
set.seed(57)
snapshots <- list(measure(0L))
for (batch in seq_len(5L)) {
  for (i in seq_len(20L)) {
    set.seed(57)
    if (mode == "error") {
      message <- tryCatch({ invoke(TRUE); "no error" }, error = conditionMessage)
      stopifnot(identical(message, "invalid new iid scalar group index"))
    } else invisible(invoke())
    stopifnot(identical(.Random.seed, expected_seed))
  }
  snapshots[[batch + 1L]] <- measure(20L * batch)
}
measurements <- do.call(rbind, snapshots)
print(measurements, row.names = FALSE)

# Valgrind changes allocator/RSS behavior; use its leak report in that mode.
if (Sys.getenv("MAGCOMP_CHECK_RSS", "true") == "true" &&
    all(is.finite(measurements$rss_MiB))) {
  growth <- tail(measurements$rss_MiB, 1L) - measurements$rss_MiB[1L]
  stopifnot(growth < 32)
}
