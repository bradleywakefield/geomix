test_that("setupVecchiaGeoMix returns correct structure", {
  skip_on_cran()

  set.seed(1)
  locs  <- as.data.frame(expand.grid(d = 1:5, x = 1:5))
  grps  <- seq_len(nrow(locs))

  V <- setupVecchiaGeoMix(locations = locs, m = 3, groups = grps)

  expect_named(V, c("constants", "groups", "groupNeighbours",
                    "groupMatrix", "groupLookup", "groupNum",
                    "groupDeps", "groupDepL"))

  expect_equal(V$constants$G, nrow(locs))
  expect_equal(V$constants$m, 3)
  expect_equal(nrow(V$groupNeighbours), nrow(locs))
  expect_equal(ncol(V$groupNeighbours), 3)
  expect_equal(length(V$groups), nrow(locs))
  expect_true(all(V$groups %in% seq_len(V$constants$G)))
})

test_that("setupVecchiaGeoMix depStructure options all return valid output", {
  skip_on_cran()

  set.seed(1)
  locs <- as.data.frame(expand.grid(d = 1:4, x = 1:4))
  grps <- seq_len(nrow(locs))

  for (dep in c("minimal", "previous", "full")) {
    V <- setupVecchiaGeoMix(locs, m = 3, groups = grps, depStructure = dep)
    expect_equal(V$constants$G, nrow(locs), info = dep)
    expect_equal(nrow(V$groupNeighbours), nrow(locs), info = dep)
  }
})

test_that("setupGeoMixModel returns correct structure on small dataset", {
  skip_on_cran()

  dat   <- make_test_data()
  setup <- setupGeoMixModel(
    data      = dat,
    K         = 3,
    dims      = c(5, 5, 5),
    variables = test_variables,
    aformula  = ~ d,
    m         = 3
  )

  expect_named(setup, c("constants", "constantsPotts", "data_list", "inits",
                         "df", "X", "lattice", "lattice_coords",
                         "controlGibbs", "controlHMC", "controlMCMC",
                         "grouping", "vecchia", "fix_lateral"))
})

test_that("setupGeoMixModel constants are correct", {
  skip_on_cran()

  dat   <- make_test_data()
  N_obs <- sum(!is.na(dat$Z2))

  setup <- setupGeoMixModel(
    data      = dat,
    K         = 3,
    dims      = c(5, 5, 5),
    variables = test_variables,
    aformula  = ~ d,
    m         = 3
  )

  expect_equal(setup$constants$N1, nrow(dat))
  expect_equal(setup$constants$N2, N_obs)
  expect_equal(setup$constants$K,  3L)
  expect_equal(setup$constants$p,  2L)   # intercept + d
})

test_that("setupGeoMixModel data_list has correct dimensions", {
  skip_on_cran()

  dat   <- make_test_data()
  N_obs <- sum(!is.na(dat$Z2))

  setup <- setupGeoMixModel(
    data      = dat,
    K         = 3,
    dims      = c(5, 5, 5),
    variables = test_variables,
    aformula  = ~ d,
    m         = 3
  )

  expect_equal(length(setup$data_list$Z2),  N_obs)
  expect_equal(length(setup$data_list$Z1),  nrow(dat))
  expect_equal(length(setup$data_list$dID), nrow(dat))
  expect_false(any(is.na(setup$data_list$Z2)))
})

test_that("setupGeoMixModel inits have correct lengths", {
  skip_on_cran()

  dat   <- make_test_data()
  setup <- setupGeoMixModel(
    data      = dat,
    K         = 3,
    dims      = c(5, 5, 5),
    variables = test_variables,
    aformula  = ~ d,
    m         = 3
  )

  expect_equal(length(setup$inits$Y1),     nrow(dat))
  expect_equal(length(setup$inits$sigma2), 3L)
  expect_equal(length(setup$inits$lL),     3L)
  expect_equal(length(setup$inits$lD),     3L)
  expect_true(all(setup$inits$Y1 %in% 1:3))
})

test_that("setupGeoMixModel custom hyperparams and inits are applied", {
  skip_on_cran()

  dat   <- make_test_data()
  setup <- setupGeoMixModel(
    data        = dat,
    K           = 3,
    dims        = c(5, 5, 5),
    variables   = test_variables,
    aformula    = ~ d,
    m           = 3,
    hyperparams = list(a_tau = 99, b_tau = 77),
    inits       = list(tau2 = 3.14)
  )

  expect_equal(setup$constants$a_tau, 99)
  expect_equal(setup$constants$b_tau, 77)
  expect_equal(setup$inits$tau2, 3.14)
})

test_that("setupGeoMixModel intercept-only formula gives p = 1", {
  skip_on_cran()

  dat   <- make_test_data()
  setup <- setupGeoMixModel(
    data      = dat,
    K         = 3,
    dims      = c(5, 5, 5),
    variables = test_variables,
    aformula  = ~ 1,
    m         = 3
  )

  expect_equal(setup$constants$p, 1L)
  expect_equal(ncol(setup$X), 1L)
})

test_that("setupGeoMixModel works on the offshore dataset", {
  skip_on_cran()

  data(offshore)
  setup <- setupGeoMixModel(
    data      = offshore,
    K         = 3,
    dims      = c(20, 20, 20),
    variables = list(
      loc = "locID", xID = "x", yID = "y", dID = "d",
      x = "x", y = "y", depth = "d", Z1 = "Z1", Z2 = "Z2"
    ),
    aformula  = ~ d,
    m         = 10
  )

  expect_equal(setup$constants$N1, 8000L)
  expect_equal(setup$constants$N2, 1600L)
  expect_equal(setup$constants$K,  3L)
})
