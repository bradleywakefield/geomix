#' Definition of nimbleList for grouped posterior terms
#'
#' Internal container storing posterior precision and linear terms.
#'
#' @keywords internal
alphaGroupPostTermsListDef <- nimble::nimbleList(
  Q = double(2),
  b = double(1)
)

#' Compute grouped posterior terms for regression coefficients
#'
#' Internal helper used by GeoMix samplers to calculate grouped Gaussian
#' posterior precision and linear terms for a class-specific regression
#' coefficient update.
#'
#' @param k Integer scalar giving the class index.
#' @param grp_inds Integer vector of indices for the current group.
#' @param nb_inds Numeric vector of neighbour indices used for conditioning.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param X Design matrix.
#' @param Z2 Numeric response vector.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Depth distance matrix.
#' @param distL Lateral distance matrix.
#'
#' @return A nimbleList containing:
#' \describe{
#'   \item{Q}{Posterior precision contribution matrix.}
#'   \item{b}{Posterior linear term contribution vector.}
#' }
#'
#' @keywords internal
alphaGroupPostTerms <- nimble::nimbleFunction(
  run = function(k = integer(0),
                 grp_inds = integer(1),
                 nb_inds = double(1),
                 sigma2 = double(1),
                 tau2 = double(0),
                 X = double(2), Z2 = double(1),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2)) {

    returnType(alphaGroupPostTermsListDef())

    nk <- length(nb_inds)
    ng <- length(grp_inds)
    p  <- dim(X)[2]

    K_gg <- matrix(nrow = ng, ncol = ng)

    Z2g <- matrix(Z2[grp_inds], ncol = 1)
    Xg  <- matrix(X[grp_inds, 1:p],nrow = ng)

    if(nk == 0){

      if(ng == 1){
        K_gg[1,1] <- sigma2[k] + tau2
      } else {
        if(LFlag[k] == 0){
          r2_gg <- distD[dID[grp_inds],  dID[grp_inds]]  / lD[k]^2
        }else{
          r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
            distD[dID[grp_inds],  dID[grp_inds]]  / lD[k]^2
        }
        r_gg <- sqrt(r2_gg)
        K_gg[1:ng,1:ng] <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
          tau2 * diag(ng)
      }

      LR  <- t(chol(K_gg))
      Xw   <- matrix(Xg, nrow = ng, ncol = p)
      Z2w <- matrix(Z2g, nrow = ng, ncol = 1)

    } else {

      Z2n <- matrix(Z2[nb_inds],nrow = nk, ncol = 1)
      Xn  <- matrix(X[nb_inds, 1:p],nrow = nk, ncol = p)

      K_nn <- matrix(nrow = nk, ncol = nk)
      K_ng <- matrix(nrow = nk, ncol = ng)

      if(LFlag[k] == 0){
        r2_ng <- distD[dID[nb_inds],  dID[grp_inds]]  / lD[k]^2
      }else{
        r2_ng <- distL[locID[nb_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[nb_inds],  dID[grp_inds]]  / lD[k]^2
      }
      r_ng <- sqrt(r2_ng)
      K_ng[1:nk,1:ng] <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

      if(LFlag[k] == 0){
        r2_gg <- distD[dID[grp_inds],  dID[grp_inds]]  / lD[k]^2
      }else{
        r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[grp_inds],  dID[grp_inds]]  / lD[k]^2
      }
      r_gg <- sqrt(r2_gg)
      K_gg[1:ng,1:ng] <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
        tau2 * diag(ng)

      if(LFlag[k] == 0){
        r2_nn <- distD[dID[nb_inds],  dID[nb_inds]]  / lD[k]^2
      }else{
        r2_nn <- distL[locID[nb_inds], locID[nb_inds]] / lL[k]^2 +
          distD[dID[nb_inds],  dID[nb_inds]]  / lD[k]^2
      }
      r_nn <- sqrt(r2_nn)
      K_nn[1:nk,1:nk] <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
        tau2 * diag(nk)

      Lnn <- t(chol(K_nn))

      W  <- forwardsolve(Lnn, K_ng)
      Bt <- backsolve(t(Lnn), W)
      B  <- t(Bt)

      Z2w <- Z2g - B %*% Z2n
      Xw   <- Xg  - B %*% Xn

      R  <- K_gg - t(W) %*% W
      LR <- t(chol(R))
    }

    Xpost <- forwardsolve(LR, Xw)
    Zpost <- forwardsolve(LR, Z2w)

    Q <- t(Xpost)%*%Xpost
    b <- (t(Xpost)%*%Zpost)[,1]

    out <- alphaGroupPostTermsListDef$new(Q = Q, b = b)
    return(out)
  }
)

#' Gibbs sampler for class-specific regression coefficients
#'
#' Internal NIMBLE sampler that updates the matrix of regression coefficients
#' `alpha` using conjugate Gibbs steps based on grouped Vecchia likelihood
#' contributions.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Optional named list of sampler control settings.
#'
#' @details
#' For each class, grouped posterior terms are accumulated across Vecchia
#' groups and combined with the prior to sample from the multivariate normal
#' full conditional distribution.
#'
#' @return No return value. Updates `model$alpha` in place.
#'
#' @keywords internal
alphaGibbsSampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    m_alpha       <- model$m_alpha
    Q_alpha       <- model$Q_alpha
    groupLookup   <- model$groupLookup
    groupNum      <- model$groupNum
    groupNeighbours <- model$groupNeighbours
    m <- dim(model$groupNeighbours)[2]
    Z2_ind        <- model$Z2_ind

    G <- length(groupNum)
    p <- length(m_alpha)
    K <- dim(model$alpha)[1]

    calcNodes <- model$getDependencies(target)
  },

  run = function() {
    returnType(void())
    tau2   <- model$tau2[1]
    Y1_valid <- model$Y1[Z2_ind]
    dIDv     <- model$dID[Z2_ind]
    locIDv   <- model$locID[Z2_ind]

    alpha <- matrix(nrow = K, ncol = p)
    Qpost <- array(dim = c(K,p,p))
    bpost <- matrix(nrow = K, ncol = p)

    b <- (Q_alpha %*% m_alpha)[,1]

    for(k in 1:K){
      Qpost[k,,] <- Q_alpha
      bpost[k,] <- b
    }

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j,1:grp_size]
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j,1:nb_size0]

      for(k in 1:K){
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          grp_indsk <- grp_inds[which(Y1_valid[grp_inds] == k)]   ### i -> jg -> i
          nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)] ### i -> jn -> i
          termList <- alphaGroupPostTerms(k = k,
                                          grp_inds = grp_indsk,
                                          nb_inds = nb_indsk,
                                          sigma2 = model$sigma2,
                                          tau2 = tau2,
                                          lL = model$lL, lD = model$lD,
                                          LFlag = model$LFlag,
                                          Z2 = model$Z2, X = model$X,
                                          dID = dIDv,
                                          locID = locIDv,
                                          distD = model$distD,
                                          distL = model$distL)
          Qpost[k,,] <- Qpost[k,,] + termList$Q
          bpost[k,] <- bpost[k,] + termList$b
        }
      }
    }

    for(k in 1:K){
      Lpost <- t(chol(Qpost[k,,]))
      y <- forwardsolve(Lpost, bpost[k,])
      mpost <- backsolve(t(Lpost), y)
      x <- rnorm(p,0,1)
      eps <- backsolve(t(Lpost),x)
      alpha[k,1:p] <- mpost + eps
    }

    model$alpha[1:K, 1:p] <<- alpha
    model$calculate(calcNodes)
    return()
  },  methods = list(reset = function() {})
)


