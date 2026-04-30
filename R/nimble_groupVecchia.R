#' Evaluate a grouped Gaussian process log-density contribution
#'
#' Internal NIMBLE helper that computes the conditional log-density
#' contribution for a class-specific observation group under the grouped
#' Gaussian process model.
#'
#' @param k Integer scalar giving the class index.
#' @param grp_inds Numeric vector of indices for the current group.
#' @param nb_inds Numeric vector of neighbour indices used for conditioning.
#' @param mvec Numeric mean vector.
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
#' @export
#' @keywords internal
groupGPvec <- nimble::nimbleFunction(
  run = function(k = integer(0),
                 grp_inds = integer(1),
                 nb_inds = double(1),
                 mvec = double(1),
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

    if(nk == 0){
      if(ng == 1){
        logdens <- dnorm(Z2[grp_inds[1]], mvec[grp_inds[1]], sd = sqrt(sigma2[k] + tau2), log = 1)
      }else{
        if(LFlag[k] == 0){
          r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }else{
          r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
            distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }
        r_gg <- sqrt(r2_gg)
        K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
          tau2*diag(ng)
        chol_i <- chol(K_gg)

        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]]
          mvecg[idx] <- mvec[grp_inds[idx]]
        }
        logdens <- dmnorm_chol(Z2g, mvecg, chol_i, prec_param = FALSE, log = 1)
      }
    } else {
      if(LFlag[k] == 0){
        r2_ng <- distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_ng <- distL[locID[nb_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }

      r_ng <- sqrt(r2_ng)
      K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

      if(nk==1){
        K_ng <- matrix(K_ng,nrow = nk, ncol = ng)
      }

      if(LFlag[k] == 0){
        r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }

      r_gg <- sqrt(r2_gg)
      K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
        tau2*diag(ng)

      if(LFlag[k] == 0){
        r2_nn <- distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }else{
        r2_nn <- distL[locID[nb_inds], locID[nb_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }

      r_nn <- sqrt(r2_nn)
      K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
        tau2*diag(nk)

      L <- t(chol(K_nn))
      v <- forwardsolve(L, Z2[nb_inds] - mvec[nb_inds])

      mu_i <- t(K_ng) %*% backsolve(t(L), v)
      w <- forwardsolve(L, K_ng)
      cov_i <- K_gg - t(w) %*% w
      check <- logdet(cov_i)
      if(is.nan(check)) check <- -Inf
      if(check < -1e6) {
        cov_i <- K_gg - t(K_ng) %*% solve(K_nn, K_ng)
      }

      if(ng == 1){
        if(cov_i[1,1] <= 0) return(-Inf)
        logdens <- dnorm(Z2[grp_inds[1]] - mvec[grp_inds[1]], mean = mu_i[1,1] , sd = sqrt(cov_i[1,1]), log = 1)
      } else {
        chol_i <- chol(cov_i)
        Z2g <- numeric(ng)
        mvecg <- numeric(ng)
        for(idx in 1:ng) {
          Z2g[idx] <- Z2[grp_inds[idx]] - mvec[grp_inds[idx]]
        }
        logdens <- dmnorm_chol(Z2g, mu_i[,1], chol_i, prec_param = FALSE, log = 1)
      }
    }
    return(logdens)
  })
#' Evaluate grouped Gaussian process log-densities by group and class
#'
#' Internal NIMBLE helper that computes grouped Gaussian process log-density
#' contributions for each group and class combination.
#'
#' @param Z2 Numeric response vector.
#' @param alpha Numeric matrix of class-specific regression coefficients.
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
#' @export
#' @keywords internal
evaluateGroupGPvec <- nimble::nimbleFunction(
  run = function(Z2 = double(1),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(2))
    G <- length(groupNum)
    logdens <- matrix(0,nrow = G, ncol = K)
    Y1_valid <- Y1[Z2_ind]  ### s -> i
    dIDv <- dID[Z2_ind]   ### s -> i
    locIDv <- locID[Z2_ind]  ### s -> i

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j,1:grp_size]          ### j -> i
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j,1:nb_size0]    ### j -> i

      for(k in 1:K){
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          grp_which_k <- which(Y1_valid[grp_inds] == k)
          grp_indsk   <- nimInteger(length(grp_which_k))
          for(ii in 1:length(grp_which_k)) { grp_indsk[ii] <- grp_inds[grp_which_k[ii]] }
          nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)] ### i -> jn -> i
          mvec <- X %*% asCol(alpha[k,]) ### i

          logdens[j,k] <-groupGPvec(k = k,
                                    grp_inds = grp_indsk,
                                    nb_inds = nb_indsk,
                                    mvec = mvec[,1],
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
#' Evaluate grouped Gaussian process log-density for one class
#'
#' Internal NIMBLE helper that computes the total grouped Gaussian process
#' log-density contribution for a single class.
#'
#' @param k Integer scalar giving the class index.
#' @param Z2 Numeric response vector.
#' @param alpha Numeric matrix of class-specific regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
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
#' @return A numeric scalar giving the total grouped log-density for class `k`.
#'
#' @export
#' @keywords internal
evaluateGroupGPvec_kclass <- nimble::nimbleFunction(
  run = function(k = integer(0), Z2 = double(1),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(0))
    G <- length(groupNum)
    logdens <- 0
    Y1_valid <- Y1[Z2_ind]  ### s -> i
    dIDv <- dID[Z2_ind]   ### s -> i
    locIDv <- locID[Z2_ind]  ### s -> i
    mvec <- X %*% asCol(alpha[k,]) ### i

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j,1:grp_size]          ### j -> i
      nb_size0 <- min(m, cum_grp_size)

      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j,1:nb_size0]    ### j -> i

        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          grp_which_k <- which(Y1_valid[grp_inds] == k)
          grp_indsk   <- nimInteger(length(grp_which_k))
          for(ii in 1:length(grp_which_k)) { grp_indsk[ii] <- grp_inds[grp_which_k[ii]] }
          nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)] ### i -> jn -> i

          logdens <- logdens + groupGPvec(k = k,
                                    grp_inds = grp_indsk,
                                    nb_inds = nb_indsk,
                                    mvec = mvec[,1],
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
    return(logdens)
  }
)

#' Compute grouped Gaussian process log-density differences
#'
#' Internal NIMBLE helper used during latent label updates to compute the
#' change in grouped Gaussian process log-density for affected groups under
#' alternative class assignments.
#'
#' @param i Integer scalar giving the observed-index position being updated.
#' @param upd_groups Numeric vector of groups affected by the update.
#' @param orig_logdens Numeric matrix of original grouped log-densities.
#' @param Z2 Numeric response vector.
#' @param alpha Numeric matrix of class-specific regression coefficients.
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
#' @return A numeric matrix of grouped log-density differences.
#'
#' @export
#' @keywords internal
diffGroupGPvec <- nimble::nimbleFunction(
  run = function(i = integer(0),
                 upd_groups = double(1),
                 orig_logdens = double(2),
                 Z2 = double(1),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(2))

    nG <- length(upd_groups)

    Y1_valid <- Y1[Z2_ind]  ### s -> i
    dIDv <- dID[Z2_ind]   ### s -> i
    locIDv <- locID[Z2_ind]  ### s -> i

    k_old <- Y1_valid[i]

    upd_logdens <- matrix(0,nrow=nG,ncol =K)

    for(uj in 1:nG){
      j <- upd_groups[uj]
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j,1:grp_size]          ### j -> i
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j,1:nb_size0]       ### j -> i

      for(k in 1:K){
        if(k != k_old) Y1_valid[i] <- k else Y1_valid[i] <- 0
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          grp_which_k <- which(Y1_valid[grp_inds] == k)
          grp_indsk   <- nimInteger(length(grp_which_k))
          for(ii in 1:length(grp_which_k)) { grp_indsk[ii] <- grp_inds[grp_which_k[ii]] }
          nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)] ### i -> jn -> i
          mvec <- X %*% asCol(alpha[k,]) ### i

          upd_logdens[uj,k] <- groupGPvec(k = k,
                                     grp_inds = grp_indsk,
                                     nb_inds = nb_indsk,
                                     mvec = mvec[,1],
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
    return(upd_logdens - orig_logdens) # Returns change in log dens
  }
)
#' Density for the grouped Gaussian process model
#'
#' Internal NIMBLE distribution for the grouped Gaussian process likelihood
#' used for continuous observations in GeoMix.
#'
#' @param x Numeric response vector.
#' @param alpha Numeric matrix of class-specific regression coefficients.
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
#' @export
#' @keywords internal
dGPgroupvec <- nimble::nimbleFunction(
  run = function(x = double(1),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2),
                 log = integer(0)){
    returnType(double(0))
    log_terms <- evaluateGroupGPvec(Z2 = x,alpha = alpha,
                                    sigma2 = sigma2, tau2 = tau2, lL = lL, LFlag = LFlag, lD = lD,
                                    Y1 = Y1, X = X, K = K, Z2_ind = Z2_ind,
                                    dID = dID, locID = locID,distD = distD,distL = distL,
                                    m = m,  groupLookup = groupLookup,
                                    groupNum = groupNum, groupNeighbours = groupNeighbours)
    logDens <- sum(log_terms)
    if(log) return(logDens) else return(exp(logDens))
  })
#' Random generation from the grouped Gaussian process model
#'
#' Internal NIMBLE random-generation function for simulating from the grouped
#' Gaussian process model one sample at a time.
#'
#' @param n Integer scalar giving the number of draws.
#' @param alpha Numeric matrix of class-specific regression coefficients.
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
#' @export
#' @keywords internal
rGPgroupvec <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(1))

    if(n != 1) stop("rGPgroupvec only generates one sample at a time")

    N <- length(Z2_ind)
    Z2_out <- numeric(N)

    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    G <- length(groupNum)

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

      for(k in 1:K){
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          # filter grp_inds
          grp_indsk_tmp <- numeric(length(grp_inds))
          count <- 0
          for(idx in 1:length(grp_inds)) {
            if(Y1_valid[grp_inds[idx]] == k) {
              count <- count + 1
              grp_indsk_tmp[count] <- grp_inds[idx]
            }
          }
          grp_indsk <- grp_indsk_tmp[1:count]

          # filter nb_inds
          nb_indsk_tmp <- numeric(length(nb_inds))
          count_nb <- 0
          for(idx in 1:length(nb_inds)) {
            if(Y1_valid[nb_inds[idx]] == k) {
              count_nb <- count_nb + 1
              nb_indsk_tmp[count_nb] <- nb_inds[idx]
            }
          }
          nb_indsk <- nb_indsk_tmp[1:count_nb]

          mvec <- (X %*% t(alpha[k,]))[,1]

          nk <- length(nb_indsk)
          ng <- length(grp_indsk)

          if(nk == 0){
            if(ng == 1){
              Z2_out[grp_indsk[1]] <- rnorm(1,
                                            mean = mvec[grp_indsk[1]],
                                            sd = sqrt(sigma2[k] + tau2))
            } else {
              if(LFlag[k]==0){
                r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
              }else{
                r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                  distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
              }
              r_gg <- sqrt(r2_gg)
              K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
                tau2*diag(ng)
              chol_i <- chol(K_gg)
              draws <- rmnorm_chol(1, mvec[grp_indsk], chol_i, prec_param = FALSE)
              for(idx in 1:ng) {
                Z2_out[grp_indsk[idx]] <- draws[idx]
              }
            }
          } else {
            # conditional case
            if(LFlag[k]==0){
              r2_ng <- distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            }else{
              r2_ng <- distL[locIDv[nb_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            }

            r_ng <- sqrt(r2_ng)
            K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

            if(LFlag[k]==0){
              r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }else{
              r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }
            r_gg <- sqrt(r2_gg)
            K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
              tau2*diag(ng)

            if(LFlag[k]==0){
              r2_nn <- distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }else{
              r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }

            r_nn <- sqrt(r2_nn)
            K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
              tau2*diag(nk)

            L <- t(chol(K_nn))
            v <- forwardsolve(L, Z2_out[nb_indsk]) # conditional on already-sampled nb
            tmp_mu <- t(K_ng) %*% backsolve(t(L), v)

            # flatten tmp_mu into a vector
            mu_vec <- numeric(ng)
            for(ii in 1:ng) mu_vec[ii] <- tmp_mu[ii, 1]

            w <- forwardsolve(L, K_ng)
            cov_i <- K_gg - t(w) %*% w

            if(ng == 1){
              if(cov_i[1,1] <= 0) return(Z2_out)
              Z2_out[grp_indsk[1]] <- rnorm(1,
                                            mean = mu_vec[1] + mvec[grp_indsk[1]],
                                            sd = sqrt(cov_i[1,1]))
            } else {
              chol_i <- chol(cov_i)
              draws <- rmnorm_chol(1,
                                   mu_vec + mvec[grp_indsk],
                                   chol_i,
                                   prec_param = FALSE)
              for(idx in 1:ng) {
                Z2_out[grp_indsk[idx]] <- draws[idx]
              }
            }
          }
        }
      }
    }
    return(Z2_out)
  })

#' Density for the grouped Gaussian process model when p = 1
#'
#' Internal NIMBLE distribution for the grouped Gaussian process likelihood
#' used for continuous observations in GeoMix.
#'
#' @param x Numeric response vector.
#' @param alpha Numeric vector of class-specific regression coefficients.
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
#' @export
#' @keywords internal
dGPgroupvecP1 <- nimble::nimbleFunction(
  run = function(x = double(1),
                 alpha = double(1),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(1),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2),
                 log = integer(0)){
    returnType(double(0))
    alphaMat <- matrix(alpha,nrow = K, ncol = 1)
    XMat <- matrix(X, nrow = length(X), ncol = 1)
    log_terms <- evaluateGroupGPvec(Z2 = x,alpha = alphaMat,
                                    sigma2 = sigma2, tau2 = tau2, lL = lL, LFlag = LFlag, lD = lD,
                                    Y1 = Y1, X = XMat, K = K, Z2_ind = Z2_ind,
                                    dID = dID, locID = locID,distD = distD,distL = distL,
                                    m = m,  groupLookup = groupLookup,
                                    groupNum = groupNum, groupNeighbours = groupNeighbours)
    logDens <- sum(log_terms)
    if(log) return(logDens) else return(exp(logDens))
  })
#' Random generation from the grouped Gaussian process model when p =1
#'
#' Internal NIMBLE random-generation function for simulating from the grouped
#' Gaussian process model one sample at a time.
#'
#' @param n Integer scalar giving the number of draws.
#' @param alpha Numeric matrix of class-specific regression coefficients.
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
#' @export
#' @keywords internal
rGPgroupvecP1 <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 alpha = double(1),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(1),
                 K = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0), groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(1))

    if(n != 1) stop("rGPgroupvec only generates one sample at a time")
    alphaM <- matrix(alpha, nrow = K, ncol = 1)
    XM <- matrix(X, nrow = length(X), ncol = 1)
    N <- length(Z2_ind)
    Z2_out <- numeric(N)

    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    G <- length(groupNum)

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j, 1:grp_size]
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j, 1:nb_size0]

      for(k in 1:K){
        kcheck <- sum(Y1_valid[grp_inds] == k)
        if(kcheck > 0){
          # filter grp_inds
          grp_indsk_tmp <- numeric(length(grp_inds))
          count <- 0
          for(idx in 1:length(grp_inds)) {
            if(Y1_valid[grp_inds[idx]] == k) {
              count <- count + 1
              grp_indsk_tmp[count] <- grp_inds[idx]
            }
          }
          grp_indsk <- grp_indsk_tmp[1:count]

          # filter nb_inds
          nb_indsk_tmp <- numeric(length(nb_inds))
          count_nb <- 0
          for(idx in 1:length(nb_inds)) {
            if(Y1_valid[nb_inds[idx]] == k) {
              count_nb <- count_nb + 1
              nb_indsk_tmp[count_nb] <- nb_inds[idx]
            }
          }
          nb_indsk <- nb_indsk_tmp[1:count_nb]

          mvec <- (XM %*% t(alphaM[k,]))[,1]

          nk <- length(nb_indsk)
          ng <- length(grp_indsk)

          if(nk == 0){
            if(ng == 1){
              Z2_out[grp_indsk[1]] <- rnorm(1,
                                            mean = mvec[grp_indsk[1]],
                                            sd = sqrt(sigma2[k] + tau2))
            } else {
              if(LFlag[k]==0){
                r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
              }else{
                r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                  distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
              }
              r_gg <- sqrt(r2_gg)
              K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
                tau2*diag(ng)
              chol_i <- chol(K_gg)
              draws <- rmnorm_chol(1, mvec[grp_indsk], chol_i, prec_param = FALSE)
              for(idx in 1:ng) {
                Z2_out[grp_indsk[idx]] <- draws[idx]
              }
            }
          } else {
            # conditional case
            if(LFlag[k]==0){
              r2_ng <- distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            }else{
              r2_ng <- distL[locIDv[nb_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[grp_indsk]] / lD[k]^2
            }

            r_ng <- sqrt(r2_ng)
            K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)

            if(LFlag[k]==0){
              r2_gg <- distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }else{
              r2_gg <- distL[locIDv[grp_indsk], locIDv[grp_indsk]] / lL[k]^2 +
                distD[dIDv[grp_indsk], dIDv[grp_indsk]] / lD[k]^2
            }
            r_gg <- sqrt(r2_gg)
            K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
              tau2*diag(ng)

            if(LFlag[k]==0){
              r2_nn <- distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }else{
              r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }

            r_nn <- sqrt(r2_nn)
            K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
              tau2*diag(nk)

            L <- t(chol(K_nn))
            v <- forwardsolve(L, Z2_out[nb_indsk]) # conditional on already-sampled nb
            tmp_mu <- t(K_ng) %*% backsolve(t(L), v)

            # flatten tmp_mu into a vector
            mu_vec <- numeric(ng)
            for(ii in 1:ng) mu_vec[ii] <- tmp_mu[ii, 1]

            w <- forwardsolve(L, K_ng)
            cov_i <- K_gg - t(w) %*% w

            if(ng == 1){
              if(cov_i[1,1] <= 0) return(Z2_out)
              Z2_out[grp_indsk[1]] <- rnorm(1,
                                            mean = mu_vec[1] + mvec[grp_indsk[1]],
                                            sd = sqrt(cov_i[1,1]))
            } else {
              chol_i <- chol(cov_i)
              draws <- rmnorm_chol(1,
                                   mu_vec + mvec[grp_indsk],
                                   chol_i,
                                   prec_param = FALSE)
              for(idx in 1:ng) {
                Z2_out[grp_indsk[idx]] <- draws[idx]
              }
            }
          }
        }
      }
    }
    return(Z2_out)
  })
nimble::registerDistributions(list(
  dGPgroupvec = list(
    BUGSdist = "dGPgroupvec(alpha, sigma2, tau2, lL, lD, LFlag, Y1, X, K, Z2_ind, dID, locID, distD, distL, m, groupLookup, groupNum, groupNeighbours)",
    types = c("value = double(1)",
              "alpha = double(2)",
              "sigma2 = double(1)", "tau2 = double(0)",
              "lL = double(1)", "lD = double(1)",
              "LFlag = double(1)",
              "Y1 = double(1)", "X = double(2)",
              "K = integer(0)", "Z2_ind = double(1)",
              "dID = double(1)", "locID = double(1)",
              "distD = double(2)", "distL = double(2)",
              "m = integer(0)", "groupLookup = integer(2)",
              "groupNum = integer(1)", "groupNeighbours = integer(2)"),
    discrete = FALSE
  ),
  dGPgroupvecP1 = list(
    BUGSdist = "dGPgroupvecP1(alpha, sigma2, tau2, lL, lD, LFlag, Y1, X, K, Z2_ind, dID, locID, distD, distL, m, groupLookup, groupNum, groupNeighbours)",
    types = c("value = double(1)",
              "alpha = double(1)",
              "sigma2 = double(1)", "tau2 = double(0)",
              "lL = double(1)", "lD = double(1)",
              "LFlag = double(1)",
              "Y1 = double(1)", "X = double(1)",
              "K = integer(0)", "Z2_ind = double(1)",
              "dID = double(1)", "locID = double(1)",
              "distD = double(2)", "distL = double(2)",
              "m = integer(0)", "groupLookup = integer(2)",
              "groupNum = integer(1)", "groupNeighbours = integer(2)"),
    discrete = FALSE
  )
))
#' Predict missing values under the grouped Gaussian process model
#'
#' Internal NIMBLE helper that computes predictive means and variances for
#' missing or unobserved locations under the grouped Gaussian process model.
#'
#' @param alpha Numeric matrix of class-specific regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param K Integer scalar giving the number of classes.
#' @param G Integer scalar giving the number of prediction groups.
#' @param Z2 Numeric response vector.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param pred_inds Numeric vector of prediction indices.
#' @param pred_groups Numeric vector assigning prediction indices to groups.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param predNeighbours Numeric matrix of neighbour indices for prediction groups.
#' @param predL Numeric vector giving neighbour counts for prediction groups.
#'
#' @return A numeric matrix with predictive means and variances.
#'
#' @export
#' @keywords internal
predictMissingGPvec <- nimble::nimbleFunction(
  run = function(alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), G = integer(0),
                 Z2 = double(1), Z2_ind = double(1),
                 pred_inds = double(1),       # indices for prediction
                 pred_groups = double(1),       # corresponding group for each pred_ind
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0),
                 predNeighbours = double(2),
                 predL = double(1)) {
    returnType(double(2))
    # Return matrix: rows = length(pred_inds), cols = [mean, variance]
    N1 <- length(Y1)
    np <- length(pred_inds)
    out <- matrix(0.0, nrow = N1, ncol = 2)
    out[Z2_ind,1] <- Z2
    # restrict Z2 etc. to observed subset
    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    for(g in 1:G){
      for(k in 1:K){
        grp_indsk <- pred_inds[which(pred_groups == g & Y1[pred_inds] == k)]

        if(length(grp_indsk) > 0){
          nb_size0 <- predL[g]
          if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- predNeighbours[g,1:nb_size0]
          nb_indsk  <- nb_inds[which(Y1_valid[nb_inds] == k)]
          nb_size0 <- min(length(nb_indsk),m)
          nb_indsk  <- nb_indsk[1:nb_size0]

          ng <- length(grp_indsk)
          nk <- length(nb_indsk)
          mvec <- (X %*% asCol(alpha[k,]))[,1]

          # restrict mean vector to observed sites
          mvec_valid <- mvec[Z2_ind]

          # If no neighbours
          if(nk == 0){
            out[grp_indsk,1] <- mvec[grp_indsk]
            out[grp_indsk,2] <- sigma2[k] + tau2
          }else{
            # --- Build covariance matrices ---
            if(LFlag[k]==0){
              r2_gg <- distD[dID[grp_indsk], dID[grp_indsk]] / lD[k]^2
            }else{
              r2_gg <- distL[locID[grp_indsk], locID[grp_indsk]] / lL[k]^2 +
                distD[dID[grp_indsk], dID[grp_indsk]] / lD[k]^2
            }
            r_gg <- sqrt(r2_gg)
            K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
              tau2*diag(ng)
            K_gg <- matrix(K_gg, nrow = ng, ncol = ng)

            if(LFlag[k]==0){
              r2_ng <- distD[dIDv[nb_indsk], dID[grp_indsk]] / lD[k]^2
            }else{
              r2_ng <- distL[locIDv[nb_indsk], locID[grp_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dID[grp_indsk]] / lD[k]^2
            }

            r_ng <- sqrt(r2_ng)
            K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)
            K_ng <- matrix(K_ng, nrow = nk, ncol = ng)

            if(LFlag[k]==0){
              r2_nn <- distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }else{
              r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / lL[k]^2 +
                distD[dIDv[nb_indsk], dIDv[nb_indsk]] / lD[k]^2
            }
            r_nn <- sqrt(r2_nn)
            K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
              tau2*diag(nk)
            K_nn <- matrix(K_nn, nrow = nk, ncol = nk)

            # --- Conditioning ---
            L <- t(chol(K_nn))
            v <- forwardsolve(L, Z2[nb_indsk] - mvec_valid[nb_indsk])
            alpha <- backsolve(t(L), v)
            mu_pred <- mvec[grp_indsk] + t(K_ng) %*% alpha
            w <- forwardsolve(L, K_ng)
            cov_pred <- K_gg - t(w) %*% w
            var_pred <- diag(cov_pred)
            out[grp_indsk,1] <- mu_pred
            out[grp_indsk,2] <- var_pred
          }
        }
      }
    }
    return(out)
  })
#' Compute gradient contributions for grouped Gaussian process regression terms
#'
#' Internal NIMBLE helper that computes the class-specific gradient of the
#' grouped Gaussian process log-posterior with respect to regression
#' coefficients.
#'
#' @param k Integer scalar giving the class index.
#' @param grp_inds Numeric vector of indices for the current group.
#' @param nb_inds Numeric vector of neighbour indices used for conditioning.
#' @param alpha Numeric matrix of class-specific regression coefficients.
#' @param X Design matrix.
#' @param prior_mu Numeric prior mean vector.
#' @param prior_Q Numeric prior precision matrix.
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
#' @return A numeric vector giving the gradient with respect to `alpha[k, ]`.
#'
#' @export
#' @keywords internal
groupGPvecGrad <- nimble::nimbleFunction(
  run = function(k = integer(0),
                 grp_inds = integer(1),
                 nb_inds = double(1),
                 alpha = double(2),
                 X = double(2),
                 prior_mu = double(1),
                 prior_Q = double(2),
                 sigma2 = double(1),
                 tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Z2 = double(1),
                 dID = double(1),
                 locID = double(1),
                 distD = double(2),
                 distL = double(2)) {
    returnType(double(1))
    nk <- length(nb_inds)
    ng <- length(grp_inds)
    p <- dim(X)[2]
    grad <- numeric(p)

    mvec <- X %*% asCol(alpha[k,])
    Z2r <- Z2 - mvec[,1]

    if(nk == 0){
      if(ng == 1){
        ind <- grp_inds[1]
        grad[1:p] <- X[ind,1:p]*Z2r[ind]/(sigma2[k] + tau2)
      }else{
        if(LFlag[k]==0){
          r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }else{
          r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
            distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
        }
        r_gg <- sqrt(r2_gg)
        K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
          tau2*diag(ng)
        L <- t(chol(K_gg))
        Z2g <- Z2r[grp_inds]
        Xg <- X[grp_inds,]
        w <- forwardsolve(L,Z2g)
        v <- backsolve(t(L),w)
        grad_mat <- t(Xg) %*% asCol(v)
        grad[1:p] <- grad_mat[1:p,1]
      }
    } else {
      if(LFlag[k]==0){
        r2_ng <- distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_ng <- distL[locID[nb_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[grp_inds]] / lD[k]^2
      }
      r_ng <- sqrt(r2_ng)
      K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)
      K_ng <- matrix(K_ng,nrow = nk, ncol = ng)

      if(LFlag[k]==0){
        r2_gg <- distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }else{
        r2_gg <- distL[locID[grp_inds], locID[grp_inds]] / lL[k]^2 +
          distD[dID[grp_inds], dID[grp_inds]] / lD[k]^2
      }

      r_gg <- sqrt(r2_gg)
      K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg) +
        tau2*diag(ng)

      if(LFlag[k]==0){
        r2_nn <- distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }else{
        r2_nn <- distL[locID[nb_inds], locID[nb_inds]] / lL[k]^2 +
          distD[dID[nb_inds], dID[nb_inds]] / lD[k]^2
      }
      r_nn <- sqrt(r2_nn)
      K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) +
        tau2*diag(nk)

      L <- t(chol(K_nn))
      vr <- forwardsolve(L, Z2r[nb_inds])
      Br <- t(K_ng) %*% backsolve(t(L), vr)

      eG <- Z2r[grp_inds] - Br

      vX <- forwardsolve(L, (X[nb_inds,]))
      BX <- t(K_ng) %*% backsolve(t(L), vX)

      M <- (X[grp_inds,]) - BX
      wL <- forwardsolve(L, K_ng)

      cov_i <- K_gg - t(wL) %*% wL
      check <- logdet(cov_i)
      if(is.nan(check)) check <- -Inf
      if(check < -1e6) {
        cov_i <- K_gg - t(K_ng) %*% solve(K_nn, K_ng)
      }

      Lc <- t(chol(cov_i))
      vC <- forwardsolve(Lc,eG)
      wc <- backsolve(t(Lc),vC)

      grad_mat <- t(M) %*% wc
      grad[1:p] <- grad_mat[1:p,1]
    }
    prior_grad <- prior_Q%*%(prior_mu - alpha[k,])
    post_grad <- grad + prior_grad[1:p,1]
    return(post_grad)
  })
#' Evaluate grouped Gaussian process gradient for one class
#'
#' Internal NIMBLE helper that accumulates grouped gradient contributions for
#' a single class under the grouped Gaussian process model.
#'
#' @param Z2 Numeric response vector.
#' @param alpha Numeric matrix of class-specific regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param prior_mu Numeric prior mean vector.
#' @param prior_Q Numeric prior precision matrix.
#' @param k Integer scalar giving the class index.
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
#' @return A numeric vector giving the accumulated gradient for class `k`.
#'
#' @export
#' @keywords internal
evalGPvecGrad_k <- nimble::nimbleFunction(
  run = function(Z2 = double(1),
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 prior_mu = double(1), prior_Q = double(2),
                 k = integer(0), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2),distL = double(2),
                 m = integer(0),  groupLookup = integer(2),
                 groupNum = integer(1), groupNeighbours = integer(2)) {
    returnType(double(1))
    G <- length(groupNum)
    p <- dim(X)[2]
    loggrad <- numeric(p)
    Y1_valid <- Y1[Z2_ind]  ### s -> i
    dIDv <- dID[Z2_ind]   ### s -> i
    locIDv <- locID[Z2_ind]  ### s -> i

    for(j in 1:G){
      grp_size <- groupNum[j]
      if(j > 1) cum_grp_size <- sum(groupNum[1:(j-1)]) else cum_grp_size <- 0
      grp_inds <- groupLookup[j,1:grp_size]          ### j -> i
      nb_size0 <- min(m, cum_grp_size)
      if(nb_size0 == 0) nb_inds <- numeric(0) else nb_inds <- groupNeighbours[j,1:nb_size0]    ### j -> i

      kcheck <- sum(Y1_valid[grp_inds] == k)
      if(kcheck > 0){
        grp_which_k <- which(Y1_valid[grp_inds] == k)
        grp_indsk   <- nimInteger(length(grp_which_k))
        for(ii in 1:length(grp_which_k)) { grp_indsk[ii] <- grp_inds[grp_which_k[ii]] }
        nb_indsk <- nb_inds[which(Y1_valid[nb_inds] == k)] ### i -> jn -> i

        loggrad <- groupGPvecGrad(k = k,
                                  grp_inds = grp_indsk,
                                  nb_inds = nb_indsk,
                                  X = X,
                                  alpha = alpha,
                                  sigma2 = sigma2,
                                  tau2 = tau2,
                                  prior_mu = prior_mu,
                                  prior_Q = prior_Q,
                                  lL = lL, lD = lD,
                                  LFlag = LFlag,
                                  Z2 = Z2,
                                  dID = dIDv,
                                  locID = locIDv,
                                  distD = distD,
                                  distL = distL) + loggrad
      }
    }
    return(loggrad)
  }
)
#' Definition of nimbleList for Gaussian process prediction outputs
#'
#' Internal container storing predictive means, standard deviations, and
#' covariance matrices for grouped Gaussian process predictions.
#'
#' @export
#' @keywords internal
predictOutputListDef <- nimble::nimbleList(
  mean = double(1),
  sd   = double(1),
  cov  = double(2)
)
#' Predict grouped Gaussian process responses
#'
#' Internal NIMBLE helper that computes predictive means, standard deviations,
#' and covariance matrices for selected prediction locations under the grouped
#' Gaussian process model.
#'
#' @param pred_inds Numeric vector of prediction indices.
#' @param pred_groups Numeric vector assigning prediction indices to groups.
#' @param alpha Numeric matrix of class-specific regression coefficients.
#' @param sigma2 Numeric vector of class process variances.
#' @param tau2 Numeric scalar noise variance.
#' @param lL Numeric vector of lateral range parameters.
#' @param lD Numeric vector of depth range parameters.
#' @param LFlag Numeric vector indicating lateral correlation usage.
#' @param Y1 Numeric vector of latent class labels.
#' @param X Design matrix.
#' @param K Integer scalar giving the number of classes.
#' @param G Integer scalar giving the number of prediction groups.
#' @param Z2 Numeric response vector.
#' @param Z2_ind Numeric vector of indices for observed `Z2` values.
#' @param dID Numeric vector of depth identifiers.
#' @param locID Numeric vector of location identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param distL Numeric matrix of lateral distances.
#' @param m Integer scalar giving the maximum conditioning set size.
#' @param predNeighbours Numeric matrix of neighbour indices for prediction groups.
#' @param predL Numeric vector giving neighbour counts for prediction groups.
#' @param nugget Integer scalar indicating whether to include the nugget term.
#'
#' @return A nimbleList containing:
#' \describe{
#'   \item{mean}{Predictive mean vector.}
#'   \item{sd}{Predictive standard deviation vector.}
#'   \item{cov}{Predictive covariance matrix.}
#' }
#'
#' @export
#' @keywords internal
predictGPvec <- nimble::nimbleFunction(
  run = function(pred_inds = double(1),       # indices for prediction
                 pred_groups = double(1),     # corresponding group for each pred_ind
                 alpha = double(2),
                 sigma2 = double(1), tau2 = double(0),
                 lL = double(1), lD = double(1),
                 LFlag = double(1),
                 Y1 = double(1), X = double(2),
                 K = integer(0), G = integer(0),
                 Z2 = double(1), Z2_ind = double(1),
                 dID = double(1), locID = double(1),
                 distD = double(2), distL = double(2),
                 m = integer(0),
                 predNeighbours = double(2),
                 predL = double(1),
                 nugget = integer(0)) {

    returnType(predictOutputListDef())  # declare type with ()

    Np <- length(pred_inds)
    out_mean <- numeric(Np)
    out_sd   <- numeric(Np)
    out_cov  <- matrix(0.0, Np, Np)     # initialise to zeros

    # restrict to observed subset
    Y1_valid <- Y1[Z2_ind]
    dIDv <- dID[Z2_ind]
    locIDv <- locID[Z2_ind]

    for(g in 1:G){
      for(k in 1:K){
        pred_grp_indsk <- which((pred_groups == g) & (Y1[pred_inds] == k))
        ng <- length(pred_grp_indsk)
        if(ng > 0){
          grp_indsk <- pred_inds[pred_grp_indsk]
          nb_size0 <- predL[g]
          if(nb_size0 == 0) {
            nb_inds <- numeric(0)
          } else {
            nb_inds <- predNeighbours[g, 1:nb_size0]
          }
          nb_indsk_vec <- which(Y1_valid[nb_inds] == k)
          nb_check <- length(nb_indsk_vec)
          if(nb_check > 0){
            nb_indsk  <- nb_inds[nb_indsk_vec]
            nb_size0  <- min(nb_check, m)
            nb_indsk  <- nb_indsk[1:nb_size0]
            nk <- length(nb_indsk)
          }else{
            nb_indsk <- numeric(0)
            nk <- 0
          }
          ng <- length(grp_indsk)

          mvec <- (X %*% asCol(alpha[k,]))[,1]
          mvec_valid <- mvec[Z2_ind]

          # base K_gg
          if(LFlag[k]==0){
            r2_gg <- distD[dID[grp_indsk],   dID[grp_indsk]]   / (lD[k]^2)
          }else{
            r2_gg <- distL[locID[grp_indsk], locID[grp_indsk]] / (lL[k]^2) +
              distD[dID[grp_indsk],   dID[grp_indsk]]   / (lD[k]^2)
          }
          r_gg <- sqrt(r2_gg)
          K_gg <- sigma2[k] * (1 + sqrt(3) * r_gg) * exp(-sqrt(3) * r_gg)
          K_gg <- matrix(K_gg, ng, ng)
          if(nugget != 0) K_gg <- K_gg + tau2 * diag(ng)

          if(nk == 0){
            mu_pred <- mvec[grp_indsk]
            if(nugget != 0){
              sd_pred <- rep(sqrt(sigma2[k] + tau2), ng)
            } else {
              sd_pred <- rep(sqrt(sigma2[k]), ng)
            }
            out_mean[pred_grp_indsk] <- mu_pred
            out_sd[pred_grp_indsk]   <- sd_pred
            out_cov[pred_grp_indsk, pred_grp_indsk] <- K_gg

          } else {
            # cross-covs
            if(LFlag[k]==0){
              r2_ng <- distD[dIDv[nb_indsk],   dID[grp_indsk]]   / (lD[k]^2)
            }else{
              r2_ng <- distL[locIDv[nb_indsk], locID[grp_indsk]] / (lL[k]^2) +
                distD[dIDv[nb_indsk],   dID[grp_indsk]]   / (lD[k]^2)
            }
            r_ng <- sqrt(r2_ng)
            K_ng <- sigma2[k] * (1 + sqrt(3) * r_ng) * exp(-sqrt(3) * r_ng)
            K_ng <- matrix(K_ng, nk, ng)

            if(LFlag[k]==0){
              r2_nn <- distD[dIDv[nb_indsk],   dIDv[nb_indsk]]   / (lD[k]^2)
            }else{
              r2_nn <- distL[locIDv[nb_indsk], locIDv[nb_indsk]] / (lL[k]^2) +
                distD[dIDv[nb_indsk],   dIDv[nb_indsk]]   / (lD[k]^2)
            }
            r_nn <- sqrt(r2_nn)
            K_nn <- sigma2[k] * (1 + sqrt(3) * r_nn) * exp(-sqrt(3) * r_nn) + tau2 * diag(nk)
            K_nn <- matrix(K_nn, nk, nk)

            # conditioning
            U <- chol(K_nn)  # upper-tri
            y <- Z2[nb_indsk] - mvec_valid[nb_indsk]
            z     <- forwardsolve(t(U), y)
            A <- backsolve(U, z)

            mu_pred_col <- asCol(mvec[grp_indsk]) + t(K_ng) %*% asCol(A)
            mu_pred     <- mu_pred_col[, 1]

            Q <- forwardsolve(t(U), K_ng)
            P <- backsolve(U, Q)
            cov_pred <- K_gg - t(K_ng) %*% P

            sd_pred <- numeric(ng)
            for(i in 1:ng) sd_pred[i] <- sqrt(cov_pred[i,i])

            out_mean[pred_grp_indsk] <- mu_pred
            out_sd[pred_grp_indsk]   <- sd_pred
            out_cov[pred_grp_indsk,  pred_grp_indsk] <- cov_pred
          }
        }
      }
    }

    # instantiate nimbleList in compiled mode
    out <- predictOutputListDef$new(mean = out_mean, sd = out_sd, cov = out_cov)
    return(out)
  }
)
#' Sample from grouped Gaussian process predictions
#'
#' Internal NIMBLE helper that generates predictive samples from grouped
#' Gaussian process mean and covariance outputs.
#'
#' @param m Numeric vector of predictive means.
#' @param cov Numeric covariance matrix.
#' @param Y1 Numeric vector of latent class labels.
#' @param groups Numeric vector of group assignments.
#' @param K Integer scalar giving the number of classes.
#' @param G Integer scalar giving the number of groups.
#'
#' @return A numeric vector of sampled predictive values.
#'
#' @export
#' @keywords internal
sampleGPvec <- nimble::nimbleFunction(
  run = function(m = double(1),
                 cov = double(2),
                 Y1 = double(1),
                 groups = double(1),
                 K = integer(0),
                 G = integer(0)) {
    returnType(double(1))
    N <- length(m)
    out <- numeric(N)
    for(k in 1:K){
      for(g in 1:G){
        inds <- which(groups == g & Y1 == k)
        n <- length(inds)
        if(n>0){
          U <- chol(cov[inds,inds])
          temp <-  rnorm(length(inds))%*%U
          out[inds] <- m[inds] + temp[1,]
        }
      }
    }
    return(out)
  })



