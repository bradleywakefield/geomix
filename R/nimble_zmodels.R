#' NIMBLE model code for the GeoMix model
#'
#' Internal `nimbleCode` object defining the full GeoMix hierarchical model,
#' including priors, latent class structure, Gaussian process component,
#' and boundary-aware confusion model for observed class labels.
#'
#' @details
#' The model includes:
#' \itemize{
#'   \item priors for Gaussian process variance, range, and noise parameters;
#'   \item Dirichlet priors for the confusion matrix rows;
#'   \item class-specific regression and covariance parameters;
#'   \item the grouped Gaussian process likelihood for `Z2`;
#'   \item the latent class field `Y1`, updated through a custom Gibbs sampler;
#'   \item the confusion-model likelihood for observed labels `Z1`.
#' }
#'
#' This object is used internally when constructing the GeoMix NIMBLE model.
#'
#' @noRd
code <- nimble::nimbleCode({
  # Priors
  tau2 ~ dinvgamma(a_tau,b_tau)
  sigma2_L ~ dinvgamma(a_L,b_L)
  sigma2_D ~ dinvgamma(a_D,b_D)
  h ~ T(dnorm(0, sd = sigma_h), 0,)

  #Confusion Matrix
  for(k in 1:K){
    gammaMat[k,1:K] ~ ddirch(gamma0[1:K])
  }
  #GP params
  for (k in 1:K) {
    lL[k] ~ T(dnorm(0, sd = sqrt(sigma2_L)), 0,)
    lD[k] ~ T(dnorm(0, sd = sqrt(sigma2_D)), 0,)
    alpha[k,1:p] ~ dmnorm(mean = m_alpha[1:p],prec = Q_alpha[1:p, 1:p])
    sigma2[k] ~ dinvgamma(a_sigma,b_sigma)
  }
  Z2[1:N2] ~ dGPgroupvec(alpha = alpha[1:K, 1:p],
                         sigma2 = sigma2[1:K], tau2 = tau2,
                         lL = lL[1:K], lD = lD[1:K],
                         LFlag = LFlag[1:K],
                         Y1 = Y1[1:N1], X = X[1:N2, 1:p], K = K,
                         Z2_ind = Z2_ind[1:N2],
                         dID = dID[1:N1],locID = locID[1:N1],
                         distD = distD[1:D, 1:D],distL = distL[1:L, 1:L],
                         m = m, groupLookup = groupLookup[1:G,1:mG],
                         groupNum = groupNum[1:G], groupNeighbours = groupNeighbours[1:G,1:m])
  Y1[1:N1] ~ dDummyCat(N = N1, logProb = Y1LogProb)  ## Y1LogProb updated in Gibbs Sampler
  Z1[1:N1] ~ dConfuse(Y1 = Y1[1:N1], kappa = kappa,
                      gammaMat = gammaMat[1:K,1:K], h = h, dID = dID[1:N1],
                      distD = distD[1:D, 1:D],locStack = locStack[1:L,1:2])
})

#' NIMBLE model code for the GeoMix model for p = 1
#'
#' Internal `nimbleCode` object defining the full GeoMix hierarchical model,
#' including priors, latent class structure, Gaussian process component,
#' and boundary-aware confusion model for observed class labels in an intercept only model (p=1).
#'
#' @details
#' The model includes:
#' \itemize{
#'   \item priors for Gaussian process variance, range, and noise parameters;
#'   \item Dirichlet priors for the confusion matrix rows;
#'   \item class-specific regression and covariance parameters;
#'   \item the grouped Gaussian process likelihood for `Z2`;
#'   \item the latent class field `Y1`, updated through a custom Gibbs sampler;
#'   \item the confusion-model likelihood for observed labels `Z1`.
#' }
#'
#' This object is used internally when constructing the GeoMix NIMBLE model.
#'
#' @noRd
codep1 <- nimble::nimbleCode({
  # Priors
  tau2 ~ dinvgamma(a_tau,b_tau)
  sigma2_L ~ dinvgamma(a_L,b_L)
  sigma2_D ~ dinvgamma(a_D,b_D)
  h ~ T(dnorm(0, sd = sigma_h), 0,)

  #Confusion Matrix
  for(k in 1:K){
    gammaMat[k,1:K] ~ ddirch(gamma0[1:K])
  }
  #GP params
  for (k in 1:K) {
    lL[k] ~ T(dnorm(0, sd = sqrt(sigma2_L)), 0,)
    lD[k] ~ T(dnorm(0, sd = sqrt(sigma2_D)), 0,)
    alpha[k,1] ~ dnorm(mean = m_alpha,sd = sqrt(1/Q_alpha[1, 1]))
    sigma2[k] ~ dinvgamma(a_sigma,b_sigma)
  }
  Z2[1:N2] ~ dGPgroupvecP1(alpha = alpha[1:K, 1],
                         sigma2 = sigma2[1:K], tau2 = tau2,
                         lL = lL[1:K], lD = lD[1:K],
                         LFlag = LFlag[1:K],
                         Y1 = Y1[1:N1], X = X[1:N2, 1], K = K,
                         Z2_ind = Z2_ind[1:N2],
                         dID = dID[1:N1],locID = locID[1:N1],
                         distD = distD[1:D, 1:D],distL = distL[1:L, 1:L],
                         m = m, groupLookup = groupLookup[1:G,1:mG],
                         groupNum = groupNum[1:G], groupNeighbours = groupNeighbours[1:G,1:m])
  Y1[1:N1] ~ dDummyCat(N = N1, logProb = Y1LogProb)  ## Y1LogProb updated in Gibbs Sampler
  Z1[1:N1] ~ dConfuse(Y1 = Y1[1:N1], kappa = kappa,
                      gammaMat = gammaMat[1:K,1:K], h = h, dID = dID[1:N1],
                      distD = distD[1:D, 1:D],locStack = locStack[1:L,1:2])
})
#' NIMBLE model code for GeoMix without the boundary-aware confusion mechanism (LGFM)
#'
#' Internal `nimbleCode` object defining the LGFM-style model used in GeoMix,
#' consisting of the latent class field and grouped Gaussian process component
#' without the observed-label confusion model for `Z1`.
#'
#' @details
#' The model includes:
#' \itemize{
#'   \item priors for Gaussian process variance, range, and noise parameters;
#'   \item class-specific regression and covariance parameters;
#'   \item the grouped Gaussian process likelihood for `Z2`;
#'   \item the latent class field `Y1`, updated through a custom Gibbs sampler.
#' }
#'
#' This object is used internally when fitting the LGFM variant of the model.
#'
#' @noRd
codeLGFM <- nimble::nimbleCode({
  # Priors
  tau2 ~ dinvgamma(a_tau,b_tau)
  sigma2_L ~ dinvgamma(a_L,b_L)
  sigma2_D ~ dinvgamma(a_D,b_D)

  #GP params
  for (k in 1:K) {
    lL[k] ~ T(dnorm(0, sd = sqrt(sigma2_L)), 0,)
    lD[k] ~ T(dnorm(0, sd = sqrt(sigma2_D)), 0,)
    alpha[k,1:p] ~ dmnorm(mean = m_alpha[1:p],prec = Q_alpha[1:p, 1:p])
    sigma2[k] ~ dinvgamma(a_sigma,b_sigma)
  }
  Z2[1:N2] ~ dGPgroupvec(alpha = alpha[1:K, 1:p],
                         sigma2 = sigma2[1:K], tau2 = tau2,
                         lL = lL[1:K], lD = lD[1:K],
                         LFlag = LFlag[1:K],
                         Y1 = Y1[1:N1], X = X[1:N2, 1:p], K = K,
                         Z2_ind = Z2_ind[1:N2],
                         dID = dID[1:N1],locID = locID[1:N1],
                         distD = distD[1:D, 1:D],distL = distL[1:L, 1:L],
                         m = m, groupLookup = groupLookup[1:G,1:mG],
                         groupNum = groupNum[1:G], groupNeighbours = groupNeighbours[1:G,1:m])
  Y1[1:N1] ~ dDummyCat(N = N1, logProb = Y1LogProb)  ## Y1LogProb updated in Gibbs Sampler
})

#' NIMBLE model code for intercept-only GeoMix without the boundary-aware confusion mechanism (LGFM)
#'
#' Internal `nimbleCode` object defining the LGFM-style model used in GeoMix,
#' consisting of the latent class field and grouped Gaussian process component
#' without the observed-label confusion model for `Z1` when p=1.
#'
#' @details
#' The model includes:
#' \itemize{
#'   \item priors for Gaussian process variance, range, and noise parameters;
#'   \item class-specific regression and covariance parameters;
#'   \item the grouped Gaussian process likelihood for `Z2`;
#'   \item the latent class field `Y1`, updated through a custom Gibbs sampler.
#' }
#'
#' This object is used internally when fitting the LGFM variant of the model.
#'
#' @noRd
codeLGFMp1 <- nimble::nimbleCode({
  # Priors
  tau2 ~ dinvgamma(a_tau,b_tau)
  sigma2_L ~ dinvgamma(a_L,b_L)
  sigma2_D ~ dinvgamma(a_D,b_D)

  #GP params
  for (k in 1:K) {
    lL[k] ~ T(dnorm(0, sd = sqrt(sigma2_L)), 0,)
    lD[k] ~ T(dnorm(0, sd = sqrt(sigma2_D)), 0,)
    alpha[k,1] ~ dnorm(mean = m_alpha,sd = sqrt(1/Q_alpha[1, 1]))
    sigma2[k] ~ dinvgamma(a_sigma,b_sigma)
  }
  Z2[1:N2] ~ dGPgroupvecP1(alpha = alpha[1:K, 1],
                         sigma2 = sigma2[1:K], tau2 = tau2,
                         lL = lL[1:K], lD = lD[1:K],
                         LFlag = LFlag[1:K],
                         Y1 = Y1[1:N1], X = X[1:N2, 1], K = K,
                         Z2_ind = Z2_ind[1:N2],
                         dID = dID[1:N1],locID = locID[1:N1],
                         distD = distD[1:D, 1:D],distL = distL[1:L, 1:L],
                         m = m, groupLookup = groupLookup[1:G,1:mG],
                         groupNum = groupNum[1:G], groupNeighbours = groupNeighbours[1:G,1:m])
  Y1[1:N1] ~ dDummyCat(N = N1, logProb = Y1LogProb)  ## Y1LogProb updated in Gibbs Sampler
})
