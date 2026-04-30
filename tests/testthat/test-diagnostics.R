K  <- 3
N1 <- 125
p  <- 2

make_params <- function(niter = 40, seed = 1) {
  mat <- make_mock_samples(niter = niter, K = K, N1 = N1, p = p, seed = seed)
  extract_parameters(mat)
}

# ── run_mcmc_diagnostics() ────────────────────────────────────────────────────

test_that("run_mcmc_diagnostics returns correct top-level names", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_named(diag, c("draws", "tables", "summaries", "plots"))
})

test_that("run_mcmc_diagnostics draws sublist has correct names", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_named(diag$draws, c("draws_array", "combined_matrix", "chain_matrices"))
})

test_that("run_mcmc_diagnostics tables sublist has correct names", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_named(diag$tables, c("diagnostics", "class_diagnostics",
                               "chain_diagnostics", "overall_diagnostics",
                               "worst_rhat", "worst_ess", "Y1count"))
})

test_that("run_mcmc_diagnostics diagnostics table has expected columns", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))
  tbl    <- diag$tables$diagnostics

  expect_true(all(c("param", "ess_bulk", "ess_tail", "rhat",
                    "sd", "mcse_mean", "class") %in% names(tbl)))
})

test_that("run_mcmc_diagnostics overall_diagnostics is a one-row tibble", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_equal(nrow(diag$tables$overall_diagnostics), 1L)
  expect_true("max_rhat" %in% names(diag$tables$overall_diagnostics))
  expect_true("min_bulk_ess" %in% names(diag$tables$overall_diagnostics))
})

test_that("run_mcmc_diagnostics worst_rhat has at most 10 rows", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_lte(nrow(diag$tables$worst_rhat), 10L)
})

test_that("run_mcmc_diagnostics plots sublist has expected names", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_true(all(c("trace_alpha_acf", "trace_cov_acf",
                    "trace_other_acf", "trace_gamma_acf",
                    "logProb", "Y1") %in% names(diag$plots)))
})

test_that("run_mcmc_diagnostics trace plots are ggplot objects", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_s3_class(diag$plots$logProb, "ggplot")
  expect_s3_class(diag$plots$Y1,      "ggplot")
  expect_s3_class(diag$plots$trace_cov_acf$trace, "ggplot")
  expect_s3_class(diag$plots$trace_cov_acf$acf,   "ggplot")
})

test_that("run_mcmc_diagnostics n_chains is reported correctly", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_equal(diag$summaries$n_chains, 1L)
})

test_that("run_mcmc_diagnostics handles multiple chains", {
  p1   <- make_params(seed = 1)
  p2   <- make_params(seed = 2)
  diag <- suppressMessages(run_mcmc_diagnostics(list(p1, p2)))

  expect_equal(diag$summaries$n_chains, 2L)
  expect_equal(dim(diag$draws$draws_array)[2], 2L)  # 2 chains
})

test_that("run_mcmc_diagnostics draws_array has correct dimensions", {
  niter  <- 40
  params <- make_params(niter = niter)
  diag   <- suppressMessages(run_mcmc_diagnostics(params))
  arr    <- diag$draws$draws_array

  # dims: iterations x chains x variables
  expect_equal(dim(arr)[1], niter)
  expect_equal(dim(arr)[2], 1L)
  expect_gt(dim(arr)[3], 0L)
})

test_that("run_mcmc_diagnostics key_stats_text is a character string", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  expect_type(diag$summaries$key_stats_text, "character")
  expect_true(nchar(diag$summaries$key_stats_text) > 0)
})

test_that("run_mcmc_diagnostics rhat values are finite and positive", {
  params <- make_params()
  diag   <- suppressMessages(run_mcmc_diagnostics(params))

  rhats <- diag$tables$diagnostics$rhat
  expect_true(all(is.finite(rhats)))
  expect_true(all(rhats > 0))
})
