#' One-dimensional HMC sampler for class-specific GP parameters on the log scale
#'
#' Internal NIMBLE sampler that updates a positive class-specific Gaussian
#' process parameter by applying Hamiltonian Monte Carlo to its log-transformed
#' value, with gradients approximated by central finite differences.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Named list of sampler tuning parameters and GP settings.
#'
#' @details
#' The sampler operates on `u = log(theta)` so that positivity constraints are
#' enforced automatically. The log-posterior on the transformed scale includes
#' the Jacobian adjustment and a class-specific grouped Gaussian process
#' likelihood contribution, evaluated using `evaluateGroupGPvec_kclass()`.
#' Gradients are approximated numerically using a central finite-difference
#' step on the log scale.
#'
#' The `control` list is expected to contain:
#' \describe{
#'   \item{nLeap}{Number of leapfrog steps.}
#'   \item{eps}{Leapfrog step size.}
#'   \item{h}{Finite-difference step size on the log scale.}
#'   \item{mass}{Scalar mass parameter.}
#'   \item{m}{Maximum conditioning set size used in grouped GP evaluation.}
#' }
#'
#' The class index is extracted from the target node name, for example from
#' targets such as `sigma2[k]`, `lL[k]`, or `lD[k]`.
#'
#' @return No return value. Updates the target node in `model` in place.
#'
#' @export
#' @keywords internal
HMCSample1D_kZ2 <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,

  setup = function(model, mvSaved, target, control) {

    nLeap <- control$nLeap
    eps   <- control$eps
    h     <- control$h       # finite-diff step on u = log(theta)
    M     <- control$mass    # scalar mass
    m <- control$m

    k <- as.numeric(sub(".*\\[(\\d+)\\].*", "\\1", target))
    targetExpanded <- model$expandNodeNames(target)
    targetNode <- targetExpanded[1]

    calcNodes <- model$getDependencies(targetExpanded)
  },

  run = function() {

    theta0 <- model[[targetNode]]
    u0 <- log(theta0)

    model[[targetNode]] <<- exp(u0)
    model$calculate(targetNode)
    tau2 <- model$tau2[1]
    lp0 <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                     alpha = model$alpha,
                                     sigma2 = model$sigma2, tau2 = tau2,
                                     lL = model$lL, lD = model$lD,
                                     LFlag = model$LFlag,
                                     Y1 = model$Y1, X = model$X,
                                     Z2_ind = model$Z2_ind,
                                     dID = model$dID,locID = model$locID,
                                     distD = model$distD, distL = model$distL,
                                     m = m, groupLookup = model$groupLookup,
                                     groupNum = model$groupNum,
                                     groupNeighbours = model$groupNeighbours) +
      model$getLogProb(targetNode) + u0

    r0 <- rnorm(1, 0, sqrt(M))
    r  <- r0
    u  <- u0

    ## grad U(u) via central diff on log-scale
    model[[targetNode]] <<- exp(u + h)
    model$calculate(targetNode)
    tau2 <- model$tau2[1]
    lp_plus <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                     alpha = model$alpha,
                                     sigma2 = model$sigma2, tau2 = tau2,
                                     lL = model$lL, lD = model$lD,
                                     LFlag = model$LFlag,
                                     Y1 = model$Y1, X = model$X,
                                     Z2_ind = model$Z2_ind,
                                     dID = model$dID,locID = model$locID,
                                     distD = model$distD, distL = model$distL,
                                     m = m, groupLookup = model$groupLookup,
                                     groupNum = model$groupNum,
                                     groupNeighbours = model$groupNeighbours) +
      model$getLogProb(targetNode) + (u + h)

    model[[targetNode]] <<- exp(u - h)
    model$calculate(targetNode)
    tau2 <- model$tau2[1]
    lp_minus <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                          alpha = model$alpha,
                                          sigma2 = model$sigma2, tau2 = tau2,
                                          lL = model$lL, lD = model$lD,
                                          LFlag = model$LFlag,
                                          Y1 = model$Y1, X = model$X,
                                          Z2_ind = model$Z2_ind,
                                          dID = model$dID,locID = model$locID,
                                          distD = model$distD, distL = model$distL,
                                          m = m, groupLookup = model$groupLookup,
                                          groupNum = model$groupNum,
                                          groupNeighbours = model$groupNeighbours) +
      model$getLogProb(targetNode) + (u - h)

    gLog <- (lp_plus - lp_minus) / (2.0 * h)
    gU   <- -gLog

    r <- r - 0.5 * eps * gU

    for(step in 1:nLeap) {

      u <- u + eps * (r / M)

      model[[targetNode]] <<- exp(u + h)
      model$calculate(targetNode)
      tau2 <- model$tau2[1]
      lp_plus <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                           alpha = model$alpha,
                                           sigma2 = model$sigma2, tau2 = tau2,
                                           lL = model$lL, lD = model$lD,
                                           LFlag = model$LFlag,
                                           Y1 = model$Y1, X = model$X,
                                           Z2_ind = model$Z2_ind,
                                           dID = model$dID,locID = model$locID,
                                           distD = model$distD, distL = model$distL,
                                           m = m, groupLookup = model$groupLookup,
                                           groupNum = model$groupNum,
                                           groupNeighbours = model$groupNeighbours) +
        model$getLogProb(targetNode) + (u + h)

      model[[targetNode]] <<- exp(u - h)
      model$calculate(targetNode)
      tau2 <- model$tau2[1]
      lp_minus <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                            alpha = model$alpha,
                                            sigma2 = model$sigma2, tau2 = tau2,
                                            lL = model$lL, lD = model$lD,
                                            LFlag = model$LFlag,
                                            Y1 = model$Y1, X = model$X,
                                            Z2_ind = model$Z2_ind,
                                            dID = model$dID,locID = model$locID,
                                            distD = model$distD, distL = model$distL,
                                            m = m, groupLookup = model$groupLookup,
                                            groupNum = model$groupNum,
                                            groupNeighbours = model$groupNeighbours) +
        model$getLogProb(targetNode) + (u - h)

      gLog <- (lp_plus - lp_minus) / (2.0 * h)
      gU   <- -gLog

      if(step < nLeap) {
        r <- r - eps * gU
      } else {
        r <- r - 0.5 * eps * gU
      }
    }

    model[[targetNode]] <<- exp(u)
    model$calculate(targetNode)
    tau2 <- model$tau2[1]
    lp1 <- evaluateGroupGPvec_kclass(k = k, Z2 = model$Z2,
                                         alpha = model$alpha,
                                         sigma2 = model$sigma2, tau2 = tau2,
                                         lL = model$lL, lD = model$lD,
                                         LFlag = model$LFlag,
                                         Y1 = model$Y1, X = model$X,
                                         Z2_ind = model$Z2_ind,
                                         dID = model$dID,locID = model$locID,
                                         distD = model$distD, distL = model$distL,
                                         m = m, groupLookup = model$groupLookup,
                                         groupNum = model$groupNum,
                                         groupNeighbours = model$groupNeighbours) +
      model$getLogProb(targetNode) + u

    K0 <- 0.5 * (r0 * r0) / M
    K1 <- 0.5 * (r  * r ) / M

    logAlpha <- (lp1 - K1) - (lp0 - K0)

    if(log(runif(1, 0, 1)) < logAlpha) {
      ## accept: already set
      model$calculate(calcNodes)
    } else {
      model[[targetNode]] <<- exp(u0)
      model$calculate(targetNode)
    }

    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    return()
  },

  methods = list(reset = function() {})
)
