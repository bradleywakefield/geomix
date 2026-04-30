#' Solve a linear assignment problem in R
#'
#' Internal helper that applies the Hungarian algorithm via
#' `clue::solve_LSAP()` to obtain a permutation maximising the total
#' assignment score.
#'
#' @param C Numeric cost or score matrix.
#'
#' @return A numeric vector giving the optimal assignment permutation.
#'
#' @export
#' @keywords internal
perm_R <- function(C) {
  as.numeric(clue::solve_LSAP(C, maximum = TRUE))
}
#' R interface for Hungarian assignment within NIMBLE
#'
#' Internal `nimbleRcall` wrapper exposing the Hungarian assignment solver
#' for use inside compiled NIMBLE functions.
#'
#' @param C Numeric cost or score matrix.
#'
#' @return A numeric vector giving the optimal assignment permutation.
#'
#' @export
#' @keywords internal
perm <- nimble::nimbleRcall(
  function(C = double(2)){},
  Rfun = "perm_R",returnType = double(1)
)
#' Hungarian assignment wrapper for NIMBLE
#'
#' Internal NIMBLE helper that calls the registered assignment solver to
#' compute an optimal label permutation.
#'
#' @param C Numeric cost or score matrix.
#'
#' @return A numeric vector giving the optimal assignment permutation.
#'
#' @export
#' @keywords internal
hungarian <- nimble::nimbleFunction(
  run = function(C = double(2)) {
    returnType(double(1))
    return(perm(C))
  }
)
#' Relabel latent classes using optimal permutation matching
#'
#' Internal helper used to post-process sampled latent labels by matching
#' current class labels to observed labels using the Hungarian algorithm.
#'
#' @param Y1 Numeric vector of latent class labels.
#' @param Z1 Numeric vector of observed class labels.
#' @param K Integer scalar giving the number of classes.
#' @param N1 Integer scalar giving the number of lattice sites.
#'
#' @return A numeric vector of relabelled latent class assignments.
#'
#' @export
#' @keywords internal
relabelY1 <- nimble::nimbleFunction(
  run = function(Y1 = double(1), Z1 = double(1),
                 K = integer(0), N1 = integer(0)) {
    returnType(double(1))
    Y1new <- numeric(N1)
    C <- matrix(0, K, K)
    for(i in 1:N1){
      C[Y1[i], Z1[i]] <-  C[Y1[i], Z1[i]] + 1
    }
    Y1perm <- perm(C)
    for(i in 1:N1) Y1new[i] <- Y1perm[Y1[i]]
    return(Y1new)
  })
#' Gibbs sampler for latent class labels without observation model on the ground model.
#'
#' Internal NIMBLE sampler that updates the latent class field `Y1` using
#' Gibbs steps combining Potts prior contributions and grouped Gaussian
#' process likelihood changes without the boundary-aware confusion matrix
#' observation model on the ground model information.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Named list of sampler constants and indexing structures.
#'
#' @details
#' Each site is updated conditionally on neighbouring labels and any
#' associated continuous observations. After sampling, labels are permuted
#' to a reference ordering using the Hungarian matching step.
#'
#' @return No return value. Updates `model$Y1` in place.
#'
#' @export
#' @keywords internal
GibbsSamplerLGFM <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Constants
    K <- control$K
    N1 <- control$N1
    N2 <- control$N2
    L <- control$L
    D <- control$D
    m <- control$m
    G <- control$G
    p <- control$p
    w <- control$w
    Z1 <- as.numeric(control$Z1)

    # Potts Constants
    beta <- control$beta
    neighbourID <- as.matrix(control$neighbourID)
    storage.mode(neighbourID) <- "double"
    neighbourNum <- as.numeric(control$neighbourNum)
    weights <- control$weights
    penalty <- control$penalty

    # GP Likelihood Constants
    Z2_ind <- as.numeric(control$Z2_ind)
    Z2_flag <- control$Z2_flag
    groupDepL <- control$groupDepL
    groupDeps <- control$groupDeps
    groups <- control$groups

    #Sampling Constants
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    tau2 <- model$tau2[1]
    Y1star <- model$Y1

    initLL <- evaluateGroupGPvec(model$Z2, alpha = model$alpha,
                                 sigma2 = model$sigma2, tau2 = tau2,
                                 lL = model$lL, lD = model$lD,
                                 LFlag = model$LFlag,
                                 Y1 = Y1star, X = model$X, K = K, Z2_ind = Z2_ind,
                                 dID = model$dID, locID = model$locID,
                                 distD = model$distD, distL = model$distL,
                                 m = m, groupLookup = model$groupLookup,
                                 groupNum = model$groupNum,
                                 groupNeighbours = model$groupNeighbours)

    for (s in 1:N1) {

      log_prior_terms <- rep(0,K)
      log_lik_terms1 <- rep(0,K)

      initialY1 <- Y1star[s]
      n_valid <- neighbourNum[s]
      ndi <- 0

      if(Z2_flag[s]==1){
        id <- which(Z2_ind==s)[1]
        ndi <- groupDepL[groups[id]]
        dpds <- numeric(ndi)
        oldLL <- matrix(0,nrow=ndi,ncol=K)
        if(ndi > 0){
          for(j in 1:ndi){
            dpds[j] <- groupDeps[groups[id],j]
            oldLL[j,] <- initLL[dpds[j],]
          }
        }
        diffLL <- diffGroupGPvec(i = id,
                                 upd_groups = dpds,
                                 orig_logdens = oldLL,
                                 Z2 =  model$Z2,
                                 alpha = model$alpha,
                                 sigma2 = model$sigma2, tau2 = tau2,
                                 lL = model$lL, lD = model$lD,
                                 LFlag = model$LFlag,
                                 Y1 = Y1star, X = model$X,
                                 K = K, Z2_ind = Z2_ind,
                                 dID = model$dID, locID = model$locID,
                                 distD = model$distD, distL = model$distL,
                                 m = m,  groupLookup = model$groupLookup,
                                 groupNum = model$groupNum,
                                 groupNeighbours = model$groupNeighbours)
      }

      for (k in 1:K) {
        # Potts prior
        penalty_term <- 0.0
        if(n_valid > 0){
          for (n in 1:n_valid) {
            neighbourVal <- Y1star[neighbourID[s, n]]
            penalty_term <- penalty_term + weights[s, n] * penalty[k, neighbourVal]
          }
        }
        log_prior_terms[k] <- - beta * penalty_term

        # GP likelihood change (only if observed in Z2 and a change in Y1 is proposed)
        if((Z2_flag[s] == 1) & (k != initialY1)){
          log_lik_terms1[k] <- sum(diffLL[,initialY1]) + sum(diffLL[,k])
        }
      }
      # Posterior
      log_post_terms <- log_prior_terms + log_lik_terms1
      Mtrick <- max(log_post_terms)
      exp_terms <- exp(log_post_terms - Mtrick)
      post_terms <- exp_terms / sum(exp_terms)
      Y1star[s] <- rcat(1, post_terms)
      #model$Y1[s] <<- Y1star[s]
      if((Z2_flag[s] == 1) & (initialY1 != Y1star[s]) & (ndi > 0)){
        for(g in 1:ndi){
          initLL[dpds[g],initialY1] <- oldLL[g,initialY1] + diffLL[g,initialY1]
          initLL[dpds[g],Y1star[s]] <- oldLL[g,Y1star[s]] + diffLL[g,Y1star[s]]
        }
      }
    }
    model$Y1 <<- relabelY1(Y1star,Z1,K,N1)
    #model$Y1 <<- Y1star

    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    returnType(void())
    return()
  },
  methods = list(reset = function() {})
)
