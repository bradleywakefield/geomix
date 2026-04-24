#' Density for a weighted Potts model using pseduo-likelihood approximation
#'
#' Internal NIMBLE distribution for a weighted Potts model over a discrete
#' latent field with user-supplied neighbourhood structure and penalty matrix.
#'
#' @param x Numeric vector of class labels at all sites.
#' @param Z1 Numeric vector providing the reference field used in the local
#' conditional calculations.
#' @param K Integer scalar giving the number of classes.
#' @param beta Numeric scalar controlling Potts interaction strength.
#' @param weights Numeric matrix of neighbour-specific weights.
#' @param neighbours Numeric matrix of neighbour indices for each site.
#' @param neighbourNum Numeric vector giving the number of neighbours per site.
#' @param penalty Numeric matrix of pairwise class penalties.
#' @param log Integer scalar; if `1`, return the log-density.
#'
#' @return A numeric scalar giving the density or log-density.
#'
#' @keywords internal
dPotts <- nimble::nimbleFunction(
  run = function(x = double(1), Z1 = double(1), K = integer(0), beta = double(0),
                 weights = double(2),  neighbours = double(2),
                 neighbourNum = double(1), penalty = double(2), log = integer(0)) {
    returnType(double(0))
    N1 <- length(Z1)
    logProb <- 0.0
    for(s in 1:N1){
      M <- neighbourNum[s]
      neighID <- neighbours[s,1:M]
      penalty_terms <- numeric(K)
      ks <- x[s]
      for(k in 1:K){
        for(j in 1:M){
          sn <- neighID[j]
          kn <- Z1[sn]
          penalty_terms[k] <- penalty_terms[k] - beta * weights[s,j] * penalty[k,kn]
        }
      }
      logProb <- logProb + penalty_terms[ks] - max(penalty_terms) -
        log(sum(exp(penalty_terms - max(penalty_terms))))
    }
    if(log == 1) return(logProb) else return(exp(logProb))
  })

#' Random generation from a weighted Potts model using pseudo-likelihood approximation
#'
#' Internal NIMBLE sampler for generating a discrete field from a weighted
#' Potts model using sequential conditional updates.
#'
#' @param n Integer scalar giving the number of draws.
#' @param Z1 Numeric vector of initial or reference class labels.
#' @param K Integer scalar giving the number of classes.
#' @param beta Numeric scalar controlling Potts interaction strength.
#' @param weights Numeric matrix of neighbour-specific weights.
#' @param neighbours Numeric matrix of neighbour indices for each site.
#' @param neighbourNum Numeric vector giving the number of neighbours per site.
#' @param penalty Numeric matrix of pairwise class penalties.
#'
#' @return A numeric vector of simulated class labels.
#'
#' @keywords internal
rPotts <- nimble::nimbleFunction(
  run = function(n = integer(0), Z1 = double(1), K = integer(0), beta = double(0),
                 weights = double(2), neighbours = double(2),
                 neighbourNum = double(1), penalty = double(2)) {
    returnType(double(1))
    N1 <- length(Z1)
    Z1new <- Z1

    for(s in 1:N1){
      if(Z1[s] > 0){
        logProb <- rep(0, K)
        M <- neighbourNum[s]
        neighID <- neighbours[s,1:M]
        for(k in 1:K){
          for(j in 1:M){
            sn <- neighID[j]
            kn <- Z1new[sn]
            logProb[k] <- logProb[k] - beta * weights[s,j] * penalty[k,kn]
          }
        }
        Z1new[s] <- rcat(n = 1, prob = exp(logProb - max(logProb)))  # numerically stable
      }
    }
    return(Z1new)
  })

nimble::registerDistributions(list(
  dPotts = list(
    BUGSdist = "dPotts(Z1, K, beta, weights, neighbours, neighbourNum, penalty)",
    types = c("value = double(1)",
              "Z1 = double(1)",
              "K = integer(0)",
              "beta = double(0)",
              "weights = double(2)",
              "neighbours = double(2)",
              "neighbourNum = double(1)",
              "penalty = double(2)"),
    discrete = TRUE
  )
))
