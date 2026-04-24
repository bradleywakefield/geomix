#' Gibbs sampler for latent class labels with confusion and GP updates
#'
#' Internal NIMBLE sampler that updates the latent class field `Y1` using
#' Gibbs steps combining Potts prior contributions, boundary-aware confusion
#' likelihood terms for `Z1`, and grouped Gaussian process likelihood changes
#' for `Z2`.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Named list of sampler constants and indexing structures.
#'
#' @details
#' For each site, the sampler constructs the local posterior over class
#' assignments using three components:
#' \itemize{
#'   \item the weighted Potts prior from neighbouring labels;
#'   \item the boundary-aware confusion contribution from observed labels `Z1`;
#'   \item the grouped Gaussian process likelihood contribution from `Z2`,
#'   when observed at the site.
#' }
#'
#' A local boundary window is determined from the current value of `h`, and
#' grouped Gaussian process log-density terms are updated incrementally for
#' affected groups only.
#'
#' @return No return value. Updates `model$Y1` in place.
#'
#' @keywords internal
GibbsSampler <- nimble::nimbleFunction(
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
    #h <- control$h
    w_tol <- abs(qnorm(control$w_tol))
    delta_depth <- control$delta_depth

    # Potts Constants
    beta <- control$beta
    neighbourID <- as.matrix(control$neighbourID)
    storage.mode(neighbourID) <- "double"
    neighbourNum <- as.numeric(control$neighbourNum)
    weights <- control$weights
    penalty <- control$penalty

    # Confusion Boundary Constants
    kappa <- control$kappa

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
    h <- model$h[1]
    Y1star <- model$Y1
    phiMat <- matrix(0.0, nrow = D-1, ncol = D)
    for (r in 1:(D-1)) {
      phiMat[r, ] <- pnorm(-0.5 * (sqrt(model$distD[r+1, ]) + sqrt(model$distD[r, ])) / h)
    }
    w <- max(ceiling(h*w_tol/delta_depth)-1,5)
    #cat("\n Window ",w,"\n")
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
      log_lik_terms2 <- rep(0,K)

      initialY1 <- Y1star[s]
      n_valid <- neighbourNum[s]
      ndi <- 0

      loc <- model$locID[s]
      loc0 <- max(model$locStack[loc, 1], s - w)
      loc1 <- min(model$locStack[loc, 2], s + w)
      nStack    <- loc1 - loc0 + 1
      Y1Stack   <- Y1star[loc0:loc1]
      Z1Stack   <- model$Z1[loc0:loc1]
      dIDStackC <- model$dID[loc0:loc1]
      dIDStackR <- model$dID[loc0:(loc1-1)]
      phiStack <- phiMat[dIDStackR, dIDStackC]
      log_lik_terms2 <- proposePiStack(s = s - loc0 + 1, loc0 = loc0, loc1 = loc1,
                                       Y1Stack = Y1Stack, gammaMat = model$gammaMat,
                                       Z1Stack = Z1Stack, phiStack = phiStack,
                                       K= K, kappa = kappa)

      if(Z2_flag[s]==1){
        id <- which(Z2_ind==s)[1]
        ndi <- groupDepL[groups[id]]
        dpds <- numeric(ndi, init = TRUE)
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
      log_post_terms <- log_prior_terms + log_lik_terms1 + log_lik_terms2
      Mtrick <- max(log_post_terms)
      exp_terms <- exp(log_post_terms - Mtrick)
      post_terms <- exp_terms / sum(exp_terms)
      Y1star[s] <- rcat(1, post_terms)
      model$Y1[s] <<- Y1star[s]
      if((Z2_flag[s] == 1) & (initialY1 != Y1star[s]) & (ndi > 0)){
        for(g in 1:ndi){
          initLL[dpds[g],initialY1] <- oldLL[g,initialY1] + diffLL[g,initialY1]
          initLL[dpds[g],Y1star[s]] <- oldLL[g,Y1star[s]] + diffLL[g,Y1star[s]]
        }
      }
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    returnType(void())
    return()
  },
  methods = list(reset = function() {})
)
