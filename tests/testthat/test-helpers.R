# ── kcenter_select() (internal) ───────────────────────────────────────────────

test_that("kcenter_select returns correct structure", {
  set.seed(1)
  D   <- matrix(runif(6 * 10), nrow = 6, ncol = 10)
  out <- geomix:::kcenter_select(D, m = 3)

  expect_named(out, c("chosen_rows", "radius", "assignment"))
  expect_length(out$chosen_rows, 3L)
  expect_true(all(out$chosen_rows %in% seq_len(nrow(D))))
  expect_length(out$radius, 1L)
  expect_length(out$assignment, ncol(D))
})

test_that("kcenter_select chosen rows are unique", {
  set.seed(2)
  D   <- matrix(runif(20 * 15), nrow = 20, ncol = 15)
  out <- geomix:::kcenter_select(D, m = 5)

  expect_equal(length(out$chosen_rows), length(unique(out$chosen_rows)))
})

test_that("kcenter_select selects m = 1 correctly", {
  D   <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  out <- geomix:::kcenter_select(D, m = 1)

  expect_length(out$chosen_rows, 1L)
})

test_that("kcenter_select radius is non-negative", {
  set.seed(3)
  D   <- matrix(abs(rnorm(50)), nrow = 10, ncol = 5)
  out <- geomix:::kcenter_select(D, m = 4)

  expect_gte(out$radius, 0)
})

test_that("kcenter_select assignment indices are within chosen rows", {
  set.seed(4)
  D   <- matrix(runif(8 * 12), nrow = 8, ncol = 12)
  out <- geomix:::kcenter_select(D, m = 3)

  # assignment gives index into chosen_rows (1-based)
  expect_true(all(out$assignment %in% seq_along(out$chosen_rows)))
})

# ── box_cox / inv_box_cox (internal) ─────────────────────────────────────────

test_that("box_cox and inv_box_cox are inverses", {
  y      <- c(1, 2, 5, 10, 100)
  lambda <- 0.6
  expect_equal(geomix:::inv_box_cox(geomix:::box_cox(y, lambda), lambda),
               y, tolerance = 1e-10)
})

test_that("box_cox with lambda = 0.6 matches formula", {
  y      <- c(2, 4, 8)
  lambda <- 0.6
  expect_equal(geomix:::box_cox(y, lambda), (y^lambda - 1) / lambda,
               tolerance = 1e-12)
})

# ── offshore dataset ──────────────────────────────────────────────────────────

test_that("offshore dataset loads and has correct dimensions", {
  data(offshore)
  expect_equal(nrow(offshore), 8000L)
  expect_equal(ncol(offshore), 7L)
  expect_named(offshore, c("ID", "d", "x", "y", "locID", "Z1", "Z2"))
})

test_that("offshore dataset has correct variable ranges", {
  data(offshore)
  expect_equal(range(offshore$d),     c(1L, 20L))
  expect_equal(range(offshore$x),     c(1L, 20L))
  expect_equal(range(offshore$y),     c(1L, 20L))
  expect_equal(sort(unique(offshore$Z1)), c(1L, 2L, 3L))
  expect_equal(length(unique(offshore$locID)), 400L)
})

test_that("offshore dataset has correct Z2 coverage", {
  data(offshore)
  expect_equal(sum(!is.na(offshore$Z2)), 1600L)
  expect_equal(sum(is.na(offshore$Z2)),  6400L)
})

test_that("offshore Z2 values are positive", {
  data(offshore)
  z2 <- offshore$Z2[!is.na(offshore$Z2)]
  expect_true(all(z2 > 0))
})

# ── load_mcmc_samples() ───────────────────────────────────────────────────────

test_that("load_mcmc_samples reads and combines batch files correctly", {
  tmp <- tempdir()

  # Write two fake batch files for one chain
  chain_dir <- file.path(tmp, "GeoMix_1")
  dir.create(chain_dir, showWarnings = FALSE)

  mat1 <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5,
                 dimnames = list(NULL, paste0("p", 1:5)))
  mat2 <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5,
                 dimnames = list(NULL, paste0("p", 1:5)))

  saveRDS(mat1, file.path(chain_dir, "batch_1.rds"))
  saveRDS(mat2, file.path(chain_dir, "batch_2.rds"))

  samps <- load_mcmc_samples(tmp, name = "GeoMix")

  expect_type(samps, "list")
  expect_length(samps, 1L)
  expect_equal(nrow(samps[[1]]), 20L)
  expect_equal(ncol(samps[[1]]), 5L)
})

test_that("load_mcmc_samples index argument subsets batches", {
  tmp <- tempdir()

  chain_dir <- file.path(tmp, "GeoMix_sub")
  dir.create(chain_dir, showWarnings = FALSE)

  for (i in 1:4) {
    m <- matrix(rnorm(5), nrow = 1, ncol = 5)
    saveRDS(m, file.path(chain_dir, paste0("batch_", i, ".rds")))
  }

  samps <- load_mcmc_samples(tmp, name = "GeoMix_sub", index = 2:3)
  expect_equal(nrow(samps[[1]]), 2L)
})
