# Shared fixtures used across test files.
# Functions defined here are auto-loaded by testthat before each test file.

# ── Small synthetic dataset ───────────────────────────────────────────────────
# 5 x 5 lateral grid, 5 depth layers → 125 rows; K = 3 classes.
# Z2 observed at 25 randomly chosen locations (20% coverage).
make_test_data <- function(seed = 1) {
  set.seed(seed)
  grid <- expand.grid(d = 1:5, x = 1:5, y = 1:5)
  grid$locID <- as.integer(interaction(grid$x, grid$y, drop = TRUE))
  grid$ID    <- seq_len(nrow(grid))
  grid$Z1    <- sample(1:3, nrow(grid), replace = TRUE)
  grid$Z2    <- NA_real_
  obs        <- sample(nrow(grid), 25)
  grid$Z2[obs] <- rnorm(25, mean = 10, sd = 2)
  grid[order(grid$locID, grid$d), ]
}

# Default variable mapping for the test dataset
test_variables <- list(
  loc   = "locID",
  xID   = "x", yID = "y", dID = "d",
  x     = "x", y   = "y", depth = "d",
  Z1    = "Z1", Z2  = "Z2"
)

# ── Mock MCMC sample matrix ───────────────────────────────────────────────────
# Produces a plausible posterior sample matrix (niter x params) without
# running MCMC, used to test extract_parameters() and run_mcmc_diagnostics().
make_mock_samples <- function(niter = 40, K = 3, N1 = 125, p = 2, seed = 42) {
  set.seed(seed)

  # column names matching NIMBLE monitor output
  sigma2_cols  <- paste0("sigma2[", 1:K, "]")
  lL_cols      <- paste0("lL[", 1:K, "]")
  lD_cols      <- paste0("lD[", 1:K, "]")
  alpha_cols   <- as.vector(outer(1:K, 1:p, function(k, j) sprintf("alpha[%d, %d]", k, j)))
  gamma_cols   <- as.vector(outer(1:K, 1:K, function(i, j) sprintf("gammaMat[%d, %d]", i, j)))
  Y1_cols      <- paste0("Y1[", 1:N1, "]")

  fixed_cols <- c("tau2", "sigma2_L", "sigma2_D", "h", "logProb_Z1", "logProb_Z2")

  all_cols <- c(sigma2_cols, lL_cols, lD_cols, alpha_cols, gamma_cols,
                fixed_cols, Y1_cols)

  mat <- matrix(0, nrow = niter, ncol = length(all_cols),
                dimnames = list(NULL, all_cols))

  # fill with vaguely realistic values
  for (k in seq_len(K)) {
    mat[, paste0("sigma2[", k, "]")] <- abs(rnorm(niter, c(5, 4, 7)[k], 0.5))
    mat[, paste0("lL[",    k, "]")] <- abs(rnorm(niter, c(1, 0.75, 1.5)[k], 0.1))
    mat[, paste0("lD[",    k, "]")] <- abs(rnorm(niter, c(2, 3, 4)[k], 0.2))
    for (j in seq_len(p)) {
      mat[, sprintf("alpha[%d, %d]", k, j)] <- rnorm(niter, c(10, 19, 7)[k] + j - 1, 0.3)
    }
    for (j in seq_len(K)) {
      mat[, sprintf("gammaMat[%d, %d]", k, j)] <- runif(niter, 0.1, 0.9)
    }
  }
  mat[, "tau2"]      <- abs(rnorm(niter, 1, 0.1))
  mat[, "sigma2_L"]  <- abs(rnorm(niter, 4, 0.3))
  mat[, "sigma2_D"]  <- abs(rnorm(niter, 4, 0.3))
  mat[, "h"]         <- abs(rnorm(niter, 1, 0.1))
  mat[, "logProb_Z1"] <- rnorm(niter, -300, 5)
  mat[, "logProb_Z2"] <- rnorm(niter, -200, 5)

  for (i in seq_len(N1)) {
    mat[, paste0("Y1[", i, "]")] <- sample(1:K, niter, replace = TRUE)
  }

  mat
}
