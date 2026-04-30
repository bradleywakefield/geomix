K  <- 3
N1 <- 125
p  <- 2

# ── extract_parameters() ──────────────────────────────────────────────────────

test_that("extract_parameters returns correct top-level structure", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_named(out, c("params", "samples"))
  expect_type(out$params,  "list")
  expect_type(out$samples, "list")
})

test_that("extract_parameters params has expected entries", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)
  pr  <- out$params

  expect_true(all(c("sigma2", "tau2", "sigma2_L", "sigma2_D",
                    "lL", "lD", "a0", "a1", "gamma", "h",
                    "Y1", "Y1prob") %in% names(pr)))
})

test_that("extract_parameters sigma2 is a length-K vector of means", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_length(out$params$sigma2, K)
  expect_true(all(out$params$sigma2 > 0))
})

test_that("extract_parameters tau2 is a positive scalar", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_length(out$params$tau2, 1L)
  expect_gt(out$params$tau2, 0)
})

test_that("extract_parameters length scales are positive length-K vectors", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_length(out$params$lL, K)
  expect_length(out$params$lD, K)
  expect_true(all(out$params$lL > 0))
  expect_true(all(out$params$lD > 0))
})

test_that("extract_parameters Y1 has correct length and valid class labels", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_length(out$params$Y1, N1)
  expect_true(all(out$params$Y1 %in% seq_len(K)))
})

test_that("extract_parameters Y1prob is an N1 x K matrix summing to 1", {
  mat  <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out  <- extract_parameters(mat)
  prob <- out$params$Y1prob

  expect_equal(dim(prob), c(N1, K))
  expect_true(all(abs(rowSums(prob) - 1) < 1e-10))
  expect_true(all(prob >= 0 & prob <= 1))
})

test_that("extract_parameters Y1 agrees with Y1prob argmax", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expected <- max.col(out$params$Y1prob, ties.method = "first")
  expect_equal(out$params$Y1, expected)
})

test_that("extract_parameters samples alpha is iter x K x p array", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_equal(dim(out$samples$alpha), c(40L, K, p))
})

test_that("extract_parameters samples Y1 is iter x N1 matrix", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_equal(dim(out$samples$Y1), c(40L, N1))
  expect_true(all(out$samples$Y1 %in% seq_len(K)))
})

test_that("extract_parameters log probs are finite vectors", {
  mat <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out <- extract_parameters(mat)

  expect_length(out$samples$logProbZ1, 40L)
  expect_length(out$samples$logProbZ2, 40L)
  expect_true(all(is.finite(out$samples$logProb)))
})

test_that("extract_parameters applied to a list returns a list", {
  mat  <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  out  <- extract_parameters(list(chain1 = mat, chain2 = mat))

  expect_type(out, "list")
  expect_named(out, c("chain1", "chain2"))
  expect_named(out$chain1, c("params", "samples"))
})

test_that("extract_parameters posterior means are close to generating values", {
  # Tight samples → means should be near the generating centre values
  mat <- make_mock_samples(niter = 200, K = K, N1 = N1, p = p, seed = 7)
  out <- extract_parameters(mat)

  expect_equal(out$params$sigma2, colMeans(out$samples$sigma2), tolerance = 1e-10)
  expect_equal(out$params$tau2,   mean(out$samples$tau2),       tolerance = 1e-10)
})

# ── combine_chains() ──────────────────────────────────────────────────────────

test_that("combine_chains returns correct structure", {
  mat    <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  params <- extract_parameters(list(chain1 = mat, chain2 = mat))
  comb   <- combine_chains(params)

  expect_named(comb, c("params", "samples"))
})

test_that("combine_chains averages scalar params across chains", {
  mat    <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  p_list <- extract_parameters(list(chain1 = mat, chain2 = mat))

  # Two identical chains → combined == single-chain value
  comb <- combine_chains(p_list)
  expect_equal(comb$params$tau2,   p_list$chain1$params$tau2,   tolerance = 1e-10)
  expect_equal(comb$params$sigma2, p_list$chain1$params$sigma2, tolerance = 1e-10)
})

test_that("combine_chains stacks Y1 samples across chains", {
  mat    <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  p_list <- extract_parameters(list(chain1 = mat, chain2 = mat))
  comb   <- combine_chains(p_list)

  expect_equal(nrow(comb$samples$Y1), 80L)   # 40 + 40
  expect_equal(ncol(comb$samples$Y1), N1)
})

test_that("combine_chains stacks alpha array along iteration dimension", {
  mat    <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  p_list <- extract_parameters(list(chain1 = mat, chain2 = mat))
  comb   <- combine_chains(p_list)

  expect_equal(dim(comb$samples$alpha), c(80L, K, p))
})

test_that("combine_chains Y1 is derived from averaged Y1prob", {
  mat    <- make_mock_samples(niter = 40, K = K, N1 = N1, p = p)
  p_list <- extract_parameters(list(chain1 = mat, chain2 = mat))
  comb   <- combine_chains(p_list)

  expected <- apply(comb$params$Y1prob, 1, which.max)
  expect_equal(comb$params$Y1, as.integer(expected))
})
