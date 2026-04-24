#' Estimate the Potts interaction parameter beta
#'
#' Fits a reduced Potts model to the observed class labels `Z1` in order to
#' estimate the spatial interaction parameter `beta`. The model is implemented
#' in \pkg{nimble} and sampled using a custom one-dimensional HMC sampler.
#'
#' The function constructs a Potts model with `Y1` fixed to the observed labels
#' `Z1`, runs MCMC for `beta`, and returns both the posterior samples and the
#' posterior mean estimate.
#'
#' @param geomix_setup A GeoMix setup object containing the constants, data,
#'   and control settings required to fit the Potts submodel. In particular,
#'   it must contain:
#'   \describe{
#'     \item{`controlGibbs`}{A list containing `N1`, `K`, `neighbourID`,
#'     `neighbourNum`, `weights`, and `penalty`.}
#'     \item{`data_list`}{A list containing the observed class labels `Z1`.}
#'     \item{`controlHMC`}{A list of control arguments passed to the custom
#'     HMC sampler `HMCSample1D_logFD`.}
#'   }
#'
#' @param niter Integer giving the total number of MCMC iterations per chain.
#'   Defaults to `2000`.
#'
#' @param nburnin Integer giving the number of burn-in iterations discarded
#'   from each chain. Defaults to `250`.
#'
#' @param thin Integer giving the thinning interval used when saving posterior
#'   samples. Defaults to `10`.
#'
#' @param nchains Integer giving the number of MCMC chains to run.
#'   Defaults to `4`.
#'
#' @return A named list with components:
#' \describe{
#'   \item{estimate}{Numeric scalar giving the posterior mean of the sampled
#'   `beta` values across all retained draws.}
#'   \item{samples}{Matrix or vector of retained posterior samples returned by
#'   `nimble::runMCMC()`. The exact structure depends on the value of
#'   `nchains` and the default return format of `runMCMC()`.}
#' }
#'
#' @details
#' The fitted model is
#' \deqn{
#' \beta \sim \mathrm{Uniform}(0.5, 5),
#' }
#' with latent field likelihood
#' \deqn{
#' Y_1 \sim \mathrm{Potts}(Z_1, K, \beta, \mathrm{weights},
#' \mathrm{neighbours}, \mathrm{neighbourNum}, \mathrm{penalty}),
#' }
#' where `Y1` is initialised at the observed labels `Z1`.
#'
#' The model is compiled with \pkg{nimble}, the default sampler for `beta`
#' is replaced by `HMCSample1D_logFD`, and posterior inference is based on
#' the retained MCMC samples after burn-in and thinning.
#'
#' This function is typically used to obtain a plug-in estimate of the Potts
#' interaction parameter for later use in the full GeoMix model.
#'
#' @examples
#' \dontrun{
#' beta_out <- estimate_beta(geomix_setup)
#'
#' beta_out$estimate
#' head(beta_out$samples)
#'
#' beta_out2 <- estimate_beta(
#'   geomix_setup,
#'   niter = 5000,
#'   nburnin = 500,
#'   thin = 5,
#'   nchains = 2
#' )
#' }
#'
#' @seealso [nimble::nimbleModel()], [nimble::runMCMC()]
#'
#' @export
estimate_beta <- function(geomix_setup, niter = 2000, nburnin = 250, thin = 10, nchains = 4){
  beta_code <- nimbleCode({
    beta ~ dunif(0.01,5)
    Y1[1:N1] ~ dPotts(Z1=Z1[1:N1], K = K, beta = beta, weights = weights[1:N1,1:M],
                      neighbours = neighbourID[1:N1,1:M], neighbourNum = neighbourNum[1:N1],
                      penalty = penalty[1:K,1:K])
  })

  beta_constants <- c(geomix_setup$controlGibbs[c("N1","K","neighbourID","neighbourNum","weights","penalty")],
                      geomix_setup$data_list["Z1"])
  beta_constants$M <- ncol(beta_constants$neighbourID)
  beta_data <- list(Y1 = geomix_setup$data_list$Z1)
  message("Building model for beta...")
  bmodel <- suppressMessages({
    nimbleModel(beta_code, constants = beta_constants, data = beta_data,
                        inits = list(beta=0.5), buildDerivs = F)
  })
  Cbmodel <- suppressMessages({compileNimble(bmodel)})
  message("Setting up samplers...")
  bconf <- configureMCMC(Cbmodel,print = FALSE)
  bconf$replaceSampler(
    target  = "beta",
    type    = HMCSample1D_logFD,
    control = geomix_setup$controlHMC
  )
  bmcmc <- buildMCMC(bconf)
  Cbmcmc <- suppressMessages({compileNimble(bmcmc)})
  message("Running samples...")
  beta_samples <- runMCMC(Cbmcmc,nburnin = nburnin, niter = niter, thin = thin, nchains = nchains)
  message("Done.")
  return(list(estimate = mean(beta_samples), samples = beta_samples))
}

