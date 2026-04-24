#' Estimate MAP covariance parameters for the grouped GP model
#'
#' Estimates covariance parameters for the GeoMix grouped Gaussian process
#' component by maximising the marginal posterior density of the observed
#' geotechnical response `Z2`, with regression coefficients integrated out.
#'
#' The optimisation is carried out on the log scale for positivity-constrained
#' parameters, using `optim()`. The latent class field `Y1` is fixed to the
#' observed labels `Z1`, and the function returns the resulting MAP estimates
#' for the covariance parameters.
#'
#' @param geomix_setup A GeoMix setup object containing model constants, data,
#'   and initial values. It must contain at least:
#'   \describe{
#'     \item{`constants`}{A list containing all constants required by the
#'     `dGPgroupvec_margAlpha` distribution, including `K`, `N1`, `N2`, `p`,
#'     `D`, `L`, `G`, `m`, `mG`, indexing objects, distance matrices, and
#'     prior hyperparameters.}
#'     \item{`data_list`}{A list containing `Z1`, `Z2`, `X`, `Z2_ind`,
#'     `dID`, `locID`, and any other required model data.}
#'     \item{`inits`}{A list containing initial values for `sigma2`, `tau2`,
#'     `sigma2_L`, `sigma2_D`, `lL`, and `lD`.}
#'   }
#'
#' @return A named list containing MAP estimates of the covariance parameters:
#' \describe{
#'   \item{lL}{Numeric vector of length `K` containing the estimated lateral
#'   length-scale parameters.}
#'   \item{lD}{Numeric vector of length `K` containing the estimated depth
#'   length-scale parameters.}
#'   \item{sigma2}{Numeric vector of length `K` containing the estimated
#'   class-specific process variance parameters.}
#'   \item{tau2}{Numeric scalar giving the estimated nugget or noise variance.}
#' }
#'
#' @details
#' The model fitted by this function includes priors for
#' \eqn{\tau^2}, \eqn{\sigma_L^2}, \eqn{\sigma_D^2}, and the class-specific
#' covariance parameters \eqn{l_{L,k}}, \eqn{l_{D,k}}, and \eqn{\sigma_k^2}.
#' The likelihood for `Z2` is defined through the custom
#' `dGPgroupvec_margAlpha` distribution, which integrates out the regression
#' coefficients `alpha`.
#'
#' Before optimisation:
#' \itemize{
#'   \item `Y1` is fixed to `Z1`
#'   \item only the covariance-related initial values are retained
#'   \item the model is compiled with \pkg{nimble}
#' }
#'
#' Optimisation is performed over the transformed parameter vector
#' \deqn{
#' (\log(lL_1 + 1), \ldots, \log(lL_K + 1),
#'   \log(lD_1 + 1), \ldots, \log(lD_K + 1),
#'   \log(\sigma^2_1), \ldots, \log(\sigma^2_K),
#'   \log(\tau^2)),
#' }
#' and the resulting optimiser output is back-transformed before being
#' returned.
#'
#' The objective function is the negative log-posterior contribution from the
#' `Z2` node only, evaluated using `Cmodel$calculate("Z2")`. If the evaluated
#' log-density is non-finite, a large penalty value is returned to stabilise
#' optimisation.
#'
#' Note that this function does not return the full `optim()` object, Hessian,
#' or convergence diagnostics. It only returns the transformed MAP estimates.
#'
#' @examples
#' \dontrun{
#' map_cov <- estimate_MAP_covariance(geomix_setup)
#'
#' map_cov$lL
#' map_cov$lD
#' map_cov$sigma2
#' map_cov$tau2
#' }
#'
#' @seealso [stats::optim()], [nimble::nimbleModel()],
#'   [nimble::compileNimble()]
#'
#' @export
estimate_MAP_covariance <- function(geomix_setup){
  code <- nimbleCode({
    # Priors
    tau2 ~ dinvgamma(a_tau,b_tau)
    sigma2_L ~ dinvgamma(a_L,b_L)
    sigma2_D ~ dinvgamma(a_D,b_D)

    #GP params
    for (k in 1:K) {
      lL[k] ~ T(dnorm(0, sd = sqrt(sigma2_L)), 0,)
      lD[k] ~ T(dnorm(0, sd = sqrt(sigma2_D)), 0,)
      sigma2[k] ~ dinvgamma(a_sigma,b_sigma)
    }
    Z2[1:N2] ~ dGPgroupvec_margAlpha(m_alpha = m_alpha[1:p],
                                     V_alpha =  Q_alpha[1:p, 1:p],
                                     sigma2 = sigma2[1:K], tau2 = tau2,
                                     lL = lL[1:K], lD = lD[1:K],
                                     LFlag = LFlag[1:K],
                                     Y1 = Y1[1:N1], X = X[1:N2, 1:p], K = K,
                                     Z2_ind = Z2_ind[1:N2],
                                     dID = dID[1:N1],locID = locID[1:N1],
                                     distD = distD[1:D, 1:D],distL = distL[1:L, 1:L],
                                     m = m, groupLookup = groupLookup[1:G,1:mG],
                                     groupNum = groupNum[1:G], groupNeighbours = groupNeighbours[1:G,1:m])

  })

  geomix_setup$data_list$Y1 <- geomix_setup$data_list$Z1
  geomix_setup$inits <- geomix_setup$inits[c('sigma2','tau2','sigma2_L','sigma2_D','lL','lD')]
  K <- geomix_setup$constants$K
  message('Constructing likelihood....')
  model <- suppressMessages({nimbleModel(
    code,
    constants   = geomix_setup$constants,
    data        = geomix_setup$data_list,
    inits       = geomix_setup$inits,
    buildDerivs = FALSE
  )})
  Cmodel <- suppressMessages({compileNimble(model)})

  compute_nll_logscale <- function(par, Cmodel, K) {
    lL     <- exp(par[1:K])
    lD     <- exp(par[(K + 1):(2 * K)])
    sigma2 <- exp(par[(2 * K + 1):(3 * K)])
    tau2   <- exp(par[3 * K + 1])

    Cmodel$lL[]     <- lL
    Cmodel$lD[]     <- lD
    Cmodel$sigma2[] <- sigma2
    Cmodel$tau2     <- tau2

    lp <- Cmodel$calculate("Z2")

    if(is.nan(lp) | !is.finite(lp)) return(1e12)
    -lp
  }

  par0 <- c(
    log(as.numeric(Cmodel$lL)+1),
    log(as.numeric(Cmodel$lD)+1),
    log(as.numeric(Cmodel$sigma2)),
    log(Cmodel$tau2)
  )
  message('Optimising likelihood....')
  fit <- optim(
    par = par0,
    fn = compute_nll_logscale,
    Cmodel = Cmodel,
    K = K,
    control = list(maxit = 10000, trace = 0)
  )
  par <- exp(fit$par)
  message('Done.')
  list(lL =  par[1:K],
       lD = par[(K + 1):(2 * K)],
       sigma2 = par[(2 * K + 1):(3 * K)],
       tau2 = par[3 * K + 1])
}

