#' Evaluate a grouped Gaussian process log-density with alpha integrated out
#'
#' Internal NIMBLE helper that computes the conditional log-density
#' contribution for a class-specific observation group under the grouped
#' Gaussian process model after marginalising over regression coefficients.
#'
#' @param k Integer scalar giving the class index.
#' @param grp_inds Numeric vector of indices for the current group.
#' @param nb_inds Numeric vector of neighbour indices used for conditioning.
#' @param X Design matrix.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Z2 Numeric response vector.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#'
#' @return A numeric scalar giving the grouped log-density contribution.
#'
#' @keywords internal
groupGPvec_margAlpha <- nimble::nimbleFunction(
  run = function(k = integer(0),
                 grp_inds = double(1),
                 nb_inds = double(1),
                 X = double(2),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(1),
                 tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Z2 = double(1),
                 dID = double(1),
                 locID = double(1),
                 distD = double(2),
                 distL = double(2)) {
    returnType(double(0))
    nk <- length(nb_inds)
    ng <- length(grp_inds)
    p <- length(m_alpha)

    Xg <- X[grp_inds, 1:p]
    mg <- Xg %*% asCol(m_alpha)

    if(nk == 0){
      if(ng == 1){
        vg <- sigma2[k] + tau2 + inprod(Xg[1, 1:p], V_alpha[1:p, 1:p] %*% Xg[1, 1:p])
        if(vg <= 0) return(-Inf)
        logdens <- dnorm(Z2[grp_inds[1]], mg[1,1], sd = sqrt(vg), log = 1)
      }else{
        if(LFlag[k] == 0){
          r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }else{
          r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
            distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }
        r_gg <- sqrt(r2_gg)
        K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
          tau2 * diag(ng)
        K_gg <- K_gg + Xg %*% V_alpha[1:p, 1:p] %*% t(Xg)
        chol_i <- chol(K_gg)

        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]]
          mvecg[idx] <- mg[idx,1]
        }
        logdens <- dmnorm_chol(Z2g, mvecg, chol_i, prec_param = FALSE, log = 1)
      }
    } else {
      Xn <- X[nb_inds, 1:p]
      mn <- Xn %*% asCol(m_alpha)

      if(LFlag[k] == 0){
        r2_ng <- distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_ng <- distL[locID[nb_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }

      r_ng <- sqrt(r2_ng)
      K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

      if(nk == 1){
        K_ng <- matrix(K_ng, nrow = nk, ncol = ng)
      }
      K_ng <- K_ng + Xn %*% V_alpha[1:p, 1:p] %*% t(Xg)

      if(LFlag[k] == 0){
        r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }

      r_gg <- sqrt(r2_gg)
      K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
        tau2 * diag(ng)
      K_gg <- K_gg + Xg %*% V_alpha[1:p, 1:p] %*% t(Xg)

      if(LFlag[k] == 0){
        r2_nn <- distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }else{
        r2_nn <- distL[locID[nb_inds], locID[nb_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }

      r_nn <- sqrt(r2_nn)
      K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
        tau2 * diag(nk)
      K_nn <- K_nn + Xn %*% V_alpha[1:p, 1:p] %*% t(Xn)

      L <- t(chol(K_nn))
      v <- forwardsolve(L, Z2[nb_inds] - mn[,1])

      mu_i <- mg + t(K_ng) %*% backsolve(t(L), v)
      w <- forwardsolve(L, K_ng)
      cov_i <- K_gg - t(w) %*% w
      check <- logdet(cov_i)
      if(is.nan(check)) check <- -Inf
      if(check < -1e6) {
        cov_i <- K_gg - t(K_ng) %*% solve(K_nn, K_ng)
      }

      if(ng == 1){
        if(cov_i[1,1] <= 0) return(-Inf)
        logdens <- dnorm(Z2[grp_inds[1]], mean = mu_i[1,1], sd = sqrt(cov_i[1,1]), log = 1)
      } else {
        chol_i <- chol(cov_i)
        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]]
          mvecg[idx] <- mu_i[idx,1]
        }
        logdens <- dmnorm_chol(Z2g, mvecg, chol_i, prec_param = FALSE, log = 1)
      }
    }
    return(logdens)
  })
#' Evaluate grouped Gaussian process log-densities with alpha integrated out
#'
#' Internal NIMBLE helper that computes grouped Gaussian process log-density
#' contributions for each group and class combination after marginalising over
#' regression coefficients.
#'
#' @param Z2 Numeric response vector.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param K Integer scalar giving the number of classes.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param groupLookup Numeric matrix of group membership indices.
#' @param groupNum Numeric vector of group sizes.
#' @param groupNeighbours Numeric matrix of neighbour indices for each group.
#'
#' @return A numeric matrix of grouped log-density contributions.
#'
#' @keywords internal
evaluateGroupGPvec_margAlpha <- nimble::nimbleFunction(
  run = function(Z2 = double(1),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = double(2),
                 groupNum = double(1), groupNeighbours = double(2)) {
    returnType(double(2))
    G <- length(groupNum)
    logdens <- matrix(0, nrow = G, ncol = K)
    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

      for(k in 1:K){
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          grp_indsk <- grp_inds[which(Y1_valid[grp_inds] == k)]
          nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)]

          logdens[j,k] <- groupGPvec_margAlpha(k = k,
                                               grp_inds = grp_indsk,
                                               nb_inds = nb_indsk,
                                               X = X,
                                               m_alpha = m_alpha,
                                               V_alpha = V_alpha,
                                               sigma2 = sigma2,
                                               tau2 = tau2,
                                               lL = lL, lD = lD,
                                               LFlag = LFlag,
                                               Z2 = Z2,
                                               dID = dIDv,
                                               locID = locIDv,
                                               distD = distD,
                                               distL = distL)
        }
      }
    }
    return(logdens)
  }
)


#' Density for the grouped Gaussian process model with alpha integrated out
#'
#' Internal NIMBLE distribution for the grouped Gaussian process likelihood
#' after marginalising over regression coefficients.
#'
#' @param x Numeric response vector.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param K Integer scalar giving the number of classes.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param groupLookup Numeric matrix of group membership indices.
#' @param groupNum Numeric vector of group sizes.
#' @param groupNeighbours Numeric matrix of neighbour indices for each group.
#' @param log Integer scalar; if `1`, return the log-density.
#'
#' @return A numeric scalar giving the density or log-density.
#'
#' @keywords internal
dGPgroupvec_margAlpha <- nimble::nimbleFunction(
  run = function(x = double(1),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = double(2),
                 groupNum = double(1), groupNeighbours = double(2),
                 log = integer(0)) {
    returnType(double(0))
    log_terms <- evaluateGroupGPvec_margAlpha(Z2 = x,
                                              m_alpha = m_alpha,
                                              V_alpha = V_alpha,
                                              sigma2 = sigma2, tau2 = tau2,
                                              lL = lL, LFlag = LFlag, lD = lD,
                                              Y1 = Y1, X = X, K = K, Z2_ind = Z2_ind,
                                              dID = dID, locID = locID,
                                              distD = distD, distL = distL,
                                              m = m, groupLookup = groupLookup,
                                              groupNum = groupNum,
                                              groupNeighbours = groupNeighbours)
    logDens <- sum(log_terms)
    if(log) return(logDens) else return(exp(logDens))
  })

#' Random generation from the grouped Gaussian process model with alpha integrated out
#'
#' Internal NIMBLE random-generation function for simulating from the grouped
#' Gaussian process model after marginalising over regression coefficients.
#'
#' @param n Integer scalar giving the number of draws.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param K Integer scalar giving the number of classes.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param groupLookup Numeric matrix of group membership indices.
#' @param groupNum Numeric vector of group sizes.
#' @param groupNeighbours Numeric matrix of neighbour indices for each group.
#'
#' @return A numeric vector of simulated responses.
#'
#' @keywords internal
rGPgroupvec_margAlpha <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = double(2),
                 groupNum = double(1), groupNeighbours = double(2)) {

    returnType(double(1))

    if(n != 1) print("rGPgroupvec_margAlpha only supports n = 1")

    N2 <- length(Z2_ind)
    Z2 <- numeric(N2, init = FALSE)

    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    G <- length(groupNum)
    p <- length(m_alpha)

    for(j in 1:G){

      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]

      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

      for(k in 1:K){

        kcheck <- sum(Y1_valid[grp_inds] == k)

        if(kcheck > 0){

          grp_indsk <- grp_inds[which(Y1_valid[grp_inds] == k)]
          nb_indsk  <- nb_inds[which(Y1_valid[nb_inds] == k)]

          ng <- length(grp_indsk)
          nk <- length(nb_indsk)

          Xg <- X[grp_indsk, 1:p]
          mg <- Xg %*% asCol(m_alpha)

          if(nk == 0){

            # unconditional draw

            if(LFlag[k] == 0){
              r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            } else {
              r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }

            r_gg <- sqrt(r2_gg)
            K_gg <- sigma2[k] * (1 + sqrt(3)*r_gg) * exp(-sqrt(3)*r_gg) +
              tau2 * diag(ng)

            K_gg <- K_gg + Xg %*% V_alpha[1:p,1:p] %*% t(Xg)

            chol_i <- chol(K_gg)
            z <- rnorm(ng)
            draw <- mg + chol_i %*% z

          } else {

            # conditional draw

            Xn <- X[nb_indsk, 1:p]
            mn <- Xn %*% asCol(m_alpha)

            # K_nn
            if(LFlag[k] == 0){
              r2_nn <- distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            } else {
              r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }

            r_nn <- sqrt(r2_nn)
            K_nn <- sigma2[k] * (1 + sqrt(3)*r_nn) * exp(-sqrt(3)*r_nn) +
              tau2 * diag(nk)

            K_nn <- K_nn + Xn %*% V_alpha[1:p,1:p] %*% t(Xn)

            # K_ng
            if(LFlag[k] == 0){
              r2_ng <- distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            } else {
              r2_ng <- distL[locIDv[nb_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            }

            r_ng <- sqrt(r2_ng)
            K_ng <- sigma2[k] * (1 + sqrt(3)*r_ng) * exp(-sqrt(3)*r_ng)

            if(nk == 1) K_ng <- matrix(K_ng, nrow = nk, ncol = ng)

            K_ng <- K_ng + Xn %*% V_alpha[1:p,1:p] %*% t(Xg)

            # K_gg
            if(LFlag[k] == 0){
              r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            } else {
              r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }

            r_gg <- sqrt(r2_gg)
            K_gg <- sigma2[k] * (1 + sqrt(3)*r_gg) * exp(-sqrt(3)*r_gg) +
              tau2 * diag(ng)

            K_gg <- K_gg + Xg %*% V_alpha[1:p,1:p] %*% t(Xg)

            # conditional mean + cov
            L <- t(chol(K_nn))
            v <- forwardsolve(L, Z2[nb_indsk] - mn[,1])

            mu_i <- mg + t(K_ng) %*% backsolve(t(L), v)
            w <- forwardsolve(L, K_ng)
            cov_i <- K_gg - t(w) %*% w

            # jitter for safety
            for(ii in 1:ng) cov_i[ii,ii] <- cov_i[ii,ii] + 1e-8

            chol_i <- chol(cov_i)
            z <- rnorm(ng)
            draw <- mu_i + chol_i %*% z
          }

          # assign
          for(ii in 1:ng){
            Z2[grp_indsk[ii]] <- draw[ii,1]
          }
        }
      }
    }

    return(Z2)
  }
)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Evaluate a grouped Gaussian process log-density for one class with alpha integrated out
#'
#' Internal NIMBLE helper that computes the conditional log-density
#' contribution for a single class-specific observation group under the
#' grouped Gaussian process model after marginalising over regression
#' coefficients.
#'
#' @param k Integer scalar giving the class index.
#' @param grp_inds Numeric vector of indices for the current group.
#' @param nb_inds Numeric vector of neighbour indices used for conditioning.
#' @param X Design matrix.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric scalar process variance for class `k`.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric scalar lateral range parameter.
#' @param lD Numeric scalar depth range parameter.
#' @param LFlag Numeric scalar indicating lateral correlation usage.
#' @param Z2 Numeric response vector.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#'
#' @return A numeric scalar giving the grouped log-density contribution.
#'
#' @keywords internal
groupGPvec_margAlphak <- nimble::nimbleFunction(
  run = function(k = integer(0),
                 grp_inds = double(1),
                 nb_inds = double(1),
                 X = double(2),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(0),
                 tau2 = double(0),
                 lL = double(0), lD = double(0),
                 LFlag = double(0),
                 Z2 = double(1),
                 dID = double(1),
                 locID = double(1),
                 distD = double(2),
                 distL = double(2)) {
    returnType(double(0))
    nk <- length(nb_inds)
    ng <- length(grp_inds)
    p <- length(m_alpha)

    Xg <- X[grp_inds, 1:p]
    mg <- Xg %*% asCol(m_alpha)

    if(nk == 0){
      if(ng == 1){
        vg <- sigma2 + tau2 + inprod(Xg[1, 1:p], V_alpha[1:p, 1:p] %*% Xg[1, 1:p])
        if(vg <= 0) return(-Inf)
        logdens <- dnorm(Z2[grp_inds[1]], mg[1,1], sd = sqrt(vg), log = 1)
      }else{
        if(LFlag == 0){
          r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD^2
        }else{
          r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL^2 +
            distD[dID[grp_inds], dID[grp_inds]] / lD^2
        }
        r_gg <- sqrt(r2_gg)
        K_gg <- sigma2 * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
          tau2 * diag(ng)
        K_gg <- K_gg + Xg %*% V_alpha[1:p, 1:p] %*% t(Xg)
        chol_i <- chol(K_gg)

        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]]
          mvecg[idx] <- mg[idx,1]
        }
        logdens <- dmnorm_chol(Z2g, mvecg, chol_i, prec_param = FALSE, log = 1)
      }
    } else {
      Xn <- X[nb_inds, 1:p]
      mn <- Xn %*% asCol(m_alpha)

      if(LFlag == 0){
        r2_ng <- distD[dID[nb_inds], dID[grp_inds]] / lD^2
      }else{
        r2_ng <- distL[locID[nb_inds], locID[grp_inds]] / lL^2 +
          distD[dID[nb_inds], dID[grp_inds]] / lD^2
      }

      r_ng <- sqrt(r2_ng)
      K_ng <- sigma2 * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

      if(nk == 1){
        K_ng <- matrix(K_ng, nrow = nk, ncol = ng)
      }
      K_ng <- K_ng + Xn %*% V_alpha[1:p, 1:p] %*% t(Xg)

      if(LFlag == 0){
        r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD^2
      }else{
        r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL^2 +
          distD[dID[grp_inds], dID[grp_inds]] / lD^2
      }

      r_gg <- sqrt(r2_gg)
      K_gg <- sigma2 * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
        tau2 * diag(ng)
      K_gg <- K_gg + Xg %*% V_alpha[1:p, 1:p] %*% t(Xg)

      if(LFlag == 0){
        r2_nn <- distD[dID[nb_inds], dID[nb_inds]] / lD^2
      }else{
        r2_nn <- distL[locID[nb_inds], locID[nb_inds]] / lL^2 +
          distD[dID[nb_inds], dID[nb_inds]] / lD^2
      }

      r_nn <- sqrt(r2_nn)
      K_nn <- sigma2 * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
        tau2 * diag(nk)
      K_nn <- K_nn + Xn %*% V_alpha[1:p, 1:p] %*% t(Xn)

      L <- t(chol(K_nn))
      v <- forwardsolve(L, Z2[nb_inds] - mn[,1])

      mu_i <- mg + t(K_ng) %*% backsolve(t(L), v)
      w <- forwardsolve(L, K_ng)
      cov_i <- K_gg - t(w) %*% w
      check <- logdet(cov_i)
      if(is.nan(check)) check <- -Inf
      if(check < -1e6) {
        cov_i <- K_gg - t(K_ng) %*% solve(K_nn, K_ng)
      }

      if(ng == 1){
        if(cov_i[1,1] <= 0) return(-Inf)
        logdens <- dnorm(Z2[grp_inds[1]], mean = mu_i[1,1], sd = sqrt(cov_i[1,1]), log = 1)
      } else {
        chol_i <- chol(cov_i)
        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]]
          mvecg[idx] <- mu_i[idx,1]
        }
        logdens <- dmnorm_chol(Z2g, mvecg, chol_i, prec_param = FALSE, log = 1)
      }
    }
    return(logdens)
  })
#' Density for one-class grouped Gaussian process model with alpha integrated out
#'
#' Internal NIMBLE distribution for the grouped Gaussian process likelihood
#' for a single class after marginalising over regression coefficients.
#'
#' @param x Numeric response vector.
#' @param k Integer scalar giving the class index.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric scalar process variance for class `k`.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric scalar lateral range parameter.
#' @param lD Numeric scalar depth range parameter.
#' @param LFlag Numeric scalar indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param groupLookup Numeric matrix of group membership indices.
#' @param groupNum Numeric vector of group sizes.
#' @param groupNeighbours Numeric matrix of neighbour indices for each group.
#' @param log Integer scalar; if `1`, return the log-density.
#'
#' @return A numeric scalar giving the density or log-density.
#'
#' @keywords internal
dGPgroupvec_margAlphak <- nimble::nimbleFunction(
  run = function(x = double(1),
                 k = integer(0),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(0), tau2 = double(0),
                 lL = double(0), lD = double(0),
                 LFlag = double(0),
                 Y1 = double(1), X = double(2),
                 Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = double(2),
                 groupNum = double(1), groupNeighbours = double(2),
                 log = integer(0)) {
    returnType(double(0))
    G <- length(groupNum)
    logDens <- 0
    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

        logDens <- logDens + groupGPvec_margAlphak(k = k,
                                            grp_inds = grp_inds,
                                            nb_inds = nb_inds,
                                            X = X,
                                            m_alpha = m_alpha,
                                            V_alpha = V_alpha,
                                            sigma2 = sigma2,
                                            tau2 = tau2,
                                            lL = lL, lD = lD,
                                            LFlag = LFlag,
                                            Z2 = Z2,
                                            dID = dIDv,
                                            locID = locIDv,
                                            distD = distD,
                                            distL = distL)
      }
    if(log) return(logDens) else return(exp(logDens))
  })
#' Random generation from one-class grouped Gaussian process model with alpha integrated out
#'
#' Internal NIMBLE random-generation function for simulating from the grouped
#' Gaussian process model for a single class after marginalising over
#' regression coefficients.
#'
#' @param n Integer scalar giving the number of draws.
#' @param k Integer scalar giving the class index.
#' @param m_alpha Numeric prior mean vector for regression coefficients.
#' @param V_alpha Numeric prior covariance matrix for regression coefficients.
#' @param sigma2 Numeric scalar process variance for class `k`.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric scalar lateral range parameter.
#' @param lD Numeric scalar depth range parameter.
#' @param LFlag Numeric scalar indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param groupLookup Numeric matrix of group membership indices.
#' @param groupNum Numeric vector of group sizes.
#' @param groupNeighbours Numeric matrix of neighbour indices for each group.
#'
#' @return A numeric vector of simulated responses.
#'
#' @keywords internal
rGPgroupvec_margAlphak <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 k = integer(0),
                 m_alpha = double(1),
                 V_alpha = double(2),
                 sigma2 = double(0), tau2 = double(0),
                 lL = double(0), lD = double(0),
                 LFlag = double(0),
                 Y1 = double(1), X = double(2),
                 Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = double(2),
                 groupNum = double(1), groupNeighbours = double(2)) {

    returnType(double(1))

    if(n != 1) print("rGPgroupvec_margAlpha only supports n = 1")

    N2 <- length(Z2_ind)
    Z2 <- numeric(N2, init = FALSE)

    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    G <- length(groupNum)
    p <- length(m_alpha)

    for(j in 1:G){

      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]

      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

      kcheck <- sum(Y1_valid[grp_inds] == k)

      if(kcheck > 0){

        grp_indsk <- grp_inds[which(Y1_valid[grp_inds] == k)]
        nb_indsk  <- nb_inds[which(Y1_valid[nb_inds] == k)]

        ng <- length(grp_indsk)
        nk <- length(nb_indsk)

        Xg <- X[grp_indsk, 1:p]
        mg <- Xg %*% asCol(m_alpha)

        if(nk == 0){

          # unconditional draw

          if(LFlag == 0){
            r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD^2
          } else {
            r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL^2 +
              distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD^2
          }

          r_gg <- sqrt(r2_gg)
          K_gg <- sigma2 * (1 + sqrt(3)*r_gg) * exp(-sqrt(3)*r_gg) +
            tau2 * diag(ng)

          K_gg <- K_gg + Xg %*% V_alpha[1:p,1:p] %*% t(Xg)

          chol_i <- chol(K_gg)
          z <- rnorm(ng)
          draw <- mg + chol_i %*% z

        } else {

          # conditional draw

          Xn <- X[nb_indsk, 1:p]
          mn <- Xn %*% asCol(m_alpha)

          # K_nn
          if(LFlag == 0){
            r2_nn <- distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD^2
          } else {
            r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / lL^2 +
              distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD^2
          }

          r_nn <- sqrt(r2_nn)
          K_nn <- sigma2 * (1 + sqrt(3)*r_nn) * exp(-sqrt(3)*r_nn) +
            tau2 * diag(nk)

          K_nn <- K_nn + Xn %*% V_alpha[1:p,1:p] %*% t(Xn)

          # K_ng
          if(LFlag == 0){
            r2_ng <- distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD^2
          } else {
            r2_ng <- distL[locIDv[nb_indsk], locIDv[grp_indsk]] / lL^2 +
              distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD^2
          }

          r_ng <- sqrt(r2_ng)
          K_ng <- sigma2 * (1 + sqrt(3)*r_ng) * exp(-sqrt(3)*r_ng)

          if(nk == 1) K_ng <- matrix(K_ng, nrow = nk, ncol = ng)

          K_ng <- K_ng + Xn %*% V_alpha[1:p,1:p] %*% t(Xg)

          # K_gg
          if(LFlag == 0){
            r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD^2
          } else {
            r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL^2 +
              distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD^2
          }

          r_gg <- sqrt(r2_gg)
          K_gg <- sigma2 * (1 + sqrt(3)*r_gg) * exp(-sqrt(3)*r_gg) +
            tau2 * diag(ng)

          K_gg <- K_gg + Xg %*% V_alpha[1:p,1:p] %*% t(Xg)

          # conditional mean + cov
          L <- t(chol(K_nn))
          v <- forwardsolve(L, Z2[nb_indsk] - mn[,1])

          mu_i <- mg + t(K_ng) %*% backsolve(t(L), v)
          w <- forwardsolve(L, K_ng)
          cov_i <- K_gg - t(w) %*% w

          # jitter for safety
          for(ii in 1:ng) cov_i[ii,ii] <- cov_i[ii,ii] + 1e-8

          chol_i <- chol(cov_i)
          z <- rnorm(ng)
          draw <- mu_i + chol_i %*% z
        }

        # assign
        for(ii in 1:ng){
          Z2[grp_indsk[ii]] <- draw[ii,1]
        }
      }
    }

    return(Z2)
  }
)

