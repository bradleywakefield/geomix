#' Compute confusion probabilities for a vertical stack
#'
#' Internal helper used to construct the stack-specific misclassification
#' probability matrix for observed labels given the latent class field.
#'
#' @param loc0 Integer scalar giving the first index of the stack.
#' @param loc1 Integer scalar giving the last index of the stack.
#' @param Y1Stack Numeric vector of latent class labels for the stack.
#' @param gammaMat Numeric matrix of class-specific misclassification weights.
#' @param phiStack Numeric matrix of boundary proximity weights for the stack.
#' @param K Integer scalar giving the number of classes.
#' @param kappa Numeric scalar controlling the strength of the confusion model.
#'
#' @return A numeric matrix of stack-specific misclassification probabilities.
#'
#' @keywords internal
computePiStack <- nimble::nimbleFunction(
  run = function(loc0 = integer(0), loc1 = integer(0), Y1Stack = double(1),
                 gammaMat = double(2), phiStack = double(2), K= integer(0),
                 kappa = double(0)) {
    returnType(double(2))

    n <- loc1 - loc0 + 1
    Imat <- matrix(0, nrow = n, ncol = K)
    gmat <- matrix(0, nrow = n, ncol = K)

    for(j in 1:n) {
      k <- Y1Stack[j]
      Imat[j, k] <- 1
      gmat[j, ] <- gammaMat[k, ]
    }

    J <- abs(Imat[1:(n-1), ] - Imat[2:n, ])
    phi <- t(J) %*% phiStack
    phi_check <- matrix(phi <= 1, nrow = K, ncol = n)
    phi <- phi*phi_check + (1-phi_check)

    misterms <- gmat * (1 - Imat) * t(phi) * kappa + (1-kappa)/(K-1)
    corterms <- 1 - misterms %*% matrix(1, nrow=K,ncol=K)

    piMat <- misterms + Imat * corterms

    return(piMat)
  }
)
#' Propose stack-wise confusion probabilities for a site update
#'
#' Internal helper used when updating a single latent class label within a
#' vertical stack. Computes normalised log proposal probabilities over all
#' possible class assignments at the selected site.
#'
#' @param s Integer scalar giving the stack position being updated.
#' @param loc0 Integer scalar giving the first index of the stack.
#' @param loc1 Integer scalar giving the last index of the stack.
#' @param Y1Stack Numeric vector of latent class labels for the stack.
#' @param gammaMat Numeric matrix of class-specific misclassification weights.
#' @param Z1Stack Numeric vector of observed class labels for the stack.
#' @param phiStack Numeric matrix of boundary proximity weights for the stack.
#' @param K Integer scalar giving the number of classes.
#' @param kappa Numeric scalar controlling the strength of the confusion model.
#'
#' @return A numeric vector of normalised log proposal probabilities.
#'
#' @keywords internal
proposePiStack <- nimble::nimbleFunction(
  run = function(s = integer(0), loc0 = integer(0), loc1 = integer(0),
                 Y1Stack = double(1), gammaMat = double(2), Z1Stack = double(1),
                 phiStack = double(2), K= integer(0), kappa = double(0)) {
    returnType(double(1))

    n <- loc1 - loc0 + 1
    Imat <- matrix(0, nrow = n, ncol = K)
    gmat <- matrix(0, nrow = n, ncol = K)
    Z1mat <- matrix(0, nrow = n, ncol = K)
    for(j in 1:n) {
      k <- Y1Stack[j]
      l <- Z1Stack[j]
      Imat[j, k] <- 1
      gmat[j, ] <- gammaMat[k, ]
      Z1mat[j,l] <- 1
    }
    logPropProb <- numeric(K)
    for(k in 1:K){
      Imat[s, ] <- (1:K == k)*1
      J <- abs(Imat[1:(n-1), ] - Imat[2:n, ])
      phi <- t(J) %*% phiStack
      phi_check <- matrix(phi <= 1, nrow = K, ncol = n)
      phi <- phi*phi_check + (1-phi_check)

      misterms <- gmat * (1 - Imat) * t(phi) * kappa + (1-kappa)/(K-1)
      corterms <- 1 - misterms %*% matrix(1, nrow=K,ncol=K)

      piMat <- misterms + Imat * corterms

      logPropProb[k] <- sum(log((piMat*Z1mat)%*%matrix(1,nrow = K)))
    }
    M <- max(logPropProb)
    return(logPropProb-M - log(sum(exp(logPropProb-M))))
  }
)
#' Density for the boundary-aware confusion model
#'
#' Internal NIMBLE distribution for observed class labels under the
#' boundary-aware confusion model used in GeoMix.
#'
#' @param x Numeric vector of observed class labels.
#' @param Y1 Numeric vector of latent class labels.
#' @param kappa Numeric scalar controlling the strength of the confusion model.
#' @param gammaMat Numeric matrix of class-specific misclassification weights.
#' @param h Numeric scalar controlling the boundary bandwidth.
#' @param dID Numeric vector of depth identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param locStack Numeric matrix giving start and end indices for each stack.
#' @param log Integer scalar; if `1`, return the log-density.
#'
#' @return A numeric scalar giving the density or log-density.
#'
#' @export
#' @keywords internal
dConfuse <- nimble::nimbleFunction(
  run = function(x = double(1),               # Z1[1:N1]
                 Y1 = double(1),               # Y1[1:N1]
                 kappa = double(0),
                 gammaMat = double(2),
                 h = double(0),
                 dID = double(1),
                 distD = double(2),
                 locStack = double(2),
                 log = integer(0)) {
    returnType(double(0))
    L <- dim(locStack)[1]
    K <- dim(gammaMat)[1]
    D <- dim(distD)[1]
    logProb <- 0.0
    #phiMat <- pnorm(-0.5*(distD[2:D,] + distD[1:(D-1),])/h)
    phiMat <- matrix(0.0, nrow = D-1, ncol = D)
    for (r in 1:(D-1)) {
      phiMat[r, ] <- pnorm(-0.5 * (sqrt(distD[r+1, ]) + sqrt(distD[r, ])) / h)
    }
    for (i in 1:L) {
      loc0 <- locStack[i, 1]
      loc1 <- locStack[i, 2]
      n    <- loc1 - loc0 + 1
      Y1Stack   <- Y1[loc0:loc1]
      Z1Stack   <- x[loc0:loc1]
      dIDStackR <- dID[loc0:(loc1-1)]
      dIDStackC <- dID[loc0:loc1]
      phiStack <- phiMat[dIDStackR, dIDStackC]
      piMat <- computePiStack(loc0 = loc0, loc1 = loc1,
                              Y1Stack = Y1Stack,
                              gammaMat  = gammaMat,
                              phiStack = phiStack,
                              K = K, kappa = kappa)
      for(j in 1:n){
        l <- Z1Stack[j]
        logProb <- logProb + log(piMat[j,l])
      }
    }
    if (log == 1) return(logProb) else return(exp(logProb))
  }
)
#' Random generation from the boundary-aware confusion model
#'
#' Internal NIMBLE sampler for generating observed class labels from the
#' boundary-aware confusion model conditional on the latent class field.
#'
#' @param n Integer scalar giving the number of draws.
#' @param Y1 Numeric vector of latent class labels.
#' @param kappa Numeric scalar controlling the strength of the confusion model.
#' @param gammaMat Numeric matrix of class-specific misclassification weights.
#' @param h Numeric scalar controlling the boundary bandwidth.
#' @param dID Numeric vector of depth identifiers.
#' @param distD Numeric matrix of depth distances.
#' @param locStack Numeric matrix giving start and end indices for each stack.
#'
#' @return A numeric vector of simulated observed class labels.
#'
#' @export
#' @keywords internal
rConfuse <- nimble::nimbleFunction(
  run = function(n = integer(0),               # Z1[1:N1]
                 Y1 = double(1),               # Y1[1:N1]
                 kappa = double(0),
                 gammaMat = double(2),
                 h = double(0),
                 dID = double(1),
                 distD = double(2),
                 locStack = double(2)) {
    returnType(double(1))
    L <- dim(locStack)[1]
    N1 <- length(Y1)
    K <- dim(gammaMat)[1]
    D <- dim(distD)[1]
    phiMat <- matrix(0.0, nrow = D-1, ncol = D)
    for (r in 1:(D-1)) {
      phiMat[r, ] <- pnorm(-0.5 * (sqrt(distD[r+1, ]) + sqrt(distD[r, ])) / h)
    }
    x <- numeric(N1)
    for (i in 1:L) {
      loc0 <- locStack[i, 1]
      loc1 <- locStack[i, 2]
      N    <- loc1 - loc0 + 1
      Y1Stack   <- Y1[loc0:loc1]
      dIDStackR <- dID[loc0:(loc1-1)]
      dIDStackC <- dID[loc0:loc1]
      phiStack <- phiMat[dIDStackR, dIDStackC]
      piMat <- computePiStack(loc0 = loc0, loc1 = loc1,
                              Y1Stack = Y1Stack,
                              gammaMat  = gammaMat,
                              phiStack = phiStack,
                              K = K, kappa = kappa)
      for(j in 1:N) x[loc0+j-1] <- rcat(1,piMat[j,])
    }
    return(x)
  }
)

nimble::registerDistributions(list(
  dConfuse = list(
    BUGSdist = "dConfuse(Y1, kappa, gammaMat, h, dID, distD, locStack)",
    types = c(
      "value = double(1)",
      "Y1 = double(1)",
      "kappa = double(0)",
      "gammaMat = double(2)",
      "h = double(0)",
      "dID = double(1)",
      "distD = double(2)",
      "locStack = double(2)"
    ),
    discrete = TRUE
  )
))
