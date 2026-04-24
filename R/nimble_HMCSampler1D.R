#' One-dimensional HMC sampler on the log scale using finite differences
#'
#' Internal NIMBLE sampler that updates a positive scalar parameter by applying
#' Hamiltonian Monte Carlo to its log-transformed value, with gradients
#' approximated by central finite differences.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Named list of sampler tuning parameters.
#'
#' @details
#' The sampler operates on `u = log(theta)` so that positivity constraints are
#' enforced automatically. The log-posterior on the transformed scale includes
#' the Jacobian adjustment, and gradients are approximated numerically using a
#' central finite-difference step.
#'
#' The `control` list is expected to contain:
#' \describe{
#'   \item{nLeap}{Number of leapfrog steps.}
#'   \item{eps}{Leapfrog step size.}
#'   \item{h}{Finite-difference step size on the log scale.}
#'   \item{mass}{Scalar mass parameter.}
#' }
#'
#' @return No return value. Updates the target node in `model` in place.
#'
#' @keywords internal
HMCSample1D_logFD <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,

  setup = function(model, mvSaved, target, control) {

    nLeap <- control$nLeap
    eps   <- control$eps
    h     <- control$h       # finite-diff step on u = log(theta)
    M     <- control$mass    # scalar mass

    targetExpanded <- model$expandNodeNames(target)
    targetNode <- targetExpanded[1]

    calcNodes <- model$getDependencies(targetExpanded)
  },

  run = function() {

    theta0 <- model[[targetNode]]
    u0 <- log(theta0)

    model[[targetNode]] <<- exp(u0)
    model$calculate(calcNodes)
    lp0 <- model$getLogProb(calcNodes) + u0

    r0 <- rnorm(1, 0, sqrt(M))
    r  <- r0
    u  <- u0

    ## grad U(u) via central diff on log-scale
    model[[targetNode]] <<- exp(u + h)
    model$calculate(calcNodes)
    lp_plus <- model$getLogProb(calcNodes) + (u + h)

    model[[targetNode]] <<- exp(u - h)
    model$calculate(calcNodes)
    lp_minus <- model$getLogProb(calcNodes) + (u - h)

    gLog <- (lp_plus - lp_minus) / (2.0 * h)
    gU   <- -gLog

    r <- r - 0.5 * eps * gU

    for(step in 1:nLeap) {

      u <- u + eps * (r / M)

      model[[targetNode]] <<- exp(u + h)
      model$calculate(calcNodes)
      lp_plus <- model$getLogProb(calcNodes) + (u + h)

      model[[targetNode]] <<- exp(u - h)
      model$calculate(calcNodes)
      lp_minus <- model$getLogProb(calcNodes) + (u - h)

      gLog <- (lp_plus - lp_minus) / (2.0 * h)
      gU   <- -gLog

      if(step < nLeap) {
        r <- r - eps * gU
      } else {
        r <- r - 0.5 * eps * gU
      }
    }

    model[[targetNode]] <<- exp(u)
    model$calculate(calcNodes)
    lp1 <- model$getLogProb(calcNodes) + u

    K0 <- 0.5 * (r0 * r0) / M
    K1 <- 0.5 * (r  * r ) / M

    logAlpha <- (lp1 - K1) - (lp0 - K0)

    if(log(runif(1, 0, 1)) < logAlpha) {
      ## accept: already set
    } else {
      model[[targetNode]] <<- exp(u0)
      model$calculate(calcNodes)
    }

    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    return()
  },

  methods = list(reset = function() {})
)
