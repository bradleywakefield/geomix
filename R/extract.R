#' Extract posterior parameter summaries from MCMC samples
#'
#' Extracts parameter samples and posterior mean summaries from an MCMC sample
#' matrix, returning them in a structured list suitable for downstream GeoMix
#' workflows. The function recognises common GeoMix parameter names, reshapes
#' them where needed, and computes simple posterior summaries.
#'
#' Supported parameters include class-specific covariance parameters, regression
#' coefficients, latent labels, confusion matrix basis function weights, variance
#' hyperparameters, and stored log-probability terms.
#'
#' If a list of sample objects is supplied, the function is applied
#' recursively to each element.
#'
#' @param samples Either:
#'   \describe{
#'     \item{matrix-like object}{An MCMC sample object coercible to a matrix,
#'     with rows corresponding to iterations and columns named according to the
#'     monitored parameter names.}
#'     \item{list}{A list of such sample objects, for example one per chain.}
#'   }
#'
#' @return If `samples` is a list, a list of outputs obtained by applying
#'   `extract_parameters()` to each element.
#'
#'   Otherwise, a named list with components:
#'   \describe{
#'     \item{params}{Posterior summaries, typically posterior means or modal
#'     class allocations.}
#'     \item{samples}{Extracted posterior sample objects in structured form.}
#'   }
#'
#'   The `params` component may contain:
#'   \describe{
#'     \item{a0, a1, ...}{Posterior means of regression coefficient matrices,
#'     stored as vectors of length `K` for each covariate.}
#'     \item{sigma2}{Posterior means of the class-specific process variances.}
#'     \item{tau2}{Posterior mean of the nugget or noise variance.}
#'     \item{sigma2_L}{Posterior mean of the hypervariance controlling `lL`.}
#'     \item{sigma2_D}{Posterior mean of the hypervariance controlling `lD`.}
#'     \item{lL}{Posterior means of the class-specific lateral length scales.}
#'     \item{lD}{Posterior means of the class-specific depth length scales.}
#'     \item{gamma}{Posterior mean confusion matrix basis function weights, 
#'     returned as a `K x K` numeric matrix when `gammaMat[...]` samples are 
#'     present.}
#'     \item{h}{Posterior mean of the boundary bandwidth parameter.}
#'     \item{Y1}{Sitewise posterior modal class labels, obtained from the
#'     highest posterior class probability at each site.}
#'     \item{Y1prob}{Matrix of posterior class probabilities for each site,
#'     with rows indexing sites and columns indexing classes.}
#'   }
#'
#'   The `samples` component may contain:
#'   \describe{
#'     \item{alpha}{Three-dimensional array of regression coefficient samples
#'     with dimensions iteration × class × coefficient.}
#'     \item{a0, a1, ...}{Per-coefficient regression sample matrices with rows
#'     corresponding to iterations and columns to classes.}
#'     \item{sigma2}{Matrix of class-specific variance samples.}
#'     \item{tau2}{Vector of sampled `tau2` values.}
#'     \item{sigma2_L}{Vector of sampled `sigma2_L` values.}
#'     \item{sigma2_D}{Vector of sampled `sigma2_D` values.}
#'     \item{lL}{Matrix of class-specific lateral length-scale samples.}
#'     \item{lD}{Matrix of class-specific depth length-scale samples.}
#'     \item{gamma}{Matrix of sampled `gammaMat[...]` entries.}
#'     \item{h}{Vector of sampled `h` values.}
#'     \item{logProbZ1}{Vector of stored `Z1` log-probabilities, if present.}
#'     \item{logProbZ2}{Vector of stored `Z2` log-probabilities, if present.}
#'     \item{logProb}{Vector of total log-probabilities formed as
#'     `logProbZ1 + logProbZ2` when available.}
#'     \item{Y1}{Matrix of sampled latent class labels with rows corresponding
#'     to iterations and columns to sites.}
#'   }
#'
#' @details
#' Parameter extraction is based on regular expression matching of column names.
#' The function currently recognises the following naming conventions:
#' \itemize{
#'   \item `sigma2[k]`
#'   \item `alpha[k, j]`
#'   \item `Y1[i]`
#'   \item `lD[k]`
#'   \item `lL[k]`
#'   \item `gammaMat[i, j]`
#'   \item `tau2`
#'   \item `sigma2_L`
#'   \item `sigma2_D`
#'   \item `h`
#'   \item names containing `logProb_Z1` or `logProb_Z2`
#' }
#'
#' For regression coefficients `alpha[k, j]`, the function additionally
#' constructs coefficient-specific matrices `a0`, `a1`, and so on, where the
#' suffix corresponds to `j - 1`.
#'
#' For latent labels `Y1`, posterior class probabilities are computed sitewise
#' by averaging class indicators over MCMC iterations. The returned `Y1`
#' summary in `params` is the sitewise modal class based on these estimated
#' posterior probabilities.
#'
#'
#' @examples
#' \dontrun{
#' data(offshore)
#'
#' setup <- setupGeoMixModel(
#'   data      = offshore,
#'   K         = 3,
#'   dims      = c(20, 20, 20),
#'   variables = list(
#'     loc = "locID", xID = "x", yID = "y", dID = "d",
#'     x = "x", y = "y", depth = "d", Z1 = "Z1", Z2 = "Z2"
#'   ),
#'   aformula  = ~ d,
#'   m         = 10
#' )
#'
#' fit <- run_chains(setup, nchains = 2, seed = 42)
#'
#' # Extract from a single chain's sample matrix
#' params1 <- extract_parameters(fit$samples$GeoMix_1)
#'
#' params1$params$sigma2   # posterior mean class variances
#' params1$params$tau2     # posterior mean noise variance
#' params1$params$lL       # posterior mean lateral length scales
#' params1$params$lD       # posterior mean depth length scales
#' params1$params$Y1[1:10] # posterior modal class at first 10 sites
#' params1$params$Y1prob[1:5, ]  # class probabilities at first 5 sites
#'
#' # Extract from all chains at once
#' params_list <- extract_parameters(fit$samples)
#' params_list$GeoMix_1$params$sigma2
#' params_list$GeoMix_2$params$sigma2
#' }
#'
#' @seealso [run_chains()], [combine_chains()], [run_mcmc_diagnostics()]
#'
#' @export
extract_parameters <- function(samples) {
  if(is.list(samples)){
    return(lapply(samples,extract_parameters))
  }else{
    
    # Coerce once
    if (!is.matrix(samples)) {
      samples <- as.matrix(samples)
    }
    
    params <- colnames(samples)
    niter  <- nrow(samples)
    
    # ---- column locations ----
    sigma2_idx <- grep("^sigma2\\[", params)
    alpha_idx  <- grep("^alpha\\[", params)
    Y1_idx     <- grep("^Y1\\[", params)
    lD_idx     <- grep("^lD\\[", params)
    lL_idx     <- grep("^lL\\[", params)
    gamma_idx  <- grep("^gammaMat\\[", params)
    
    tau2_idx      <- match("tau2", params)
    s2L_idx      <- match("sigma2_L", params)
    s2D_idx      <- match("sigma2_D", params)
    h_idx         <- match("h", params)
    logProbZ2_idx <- grep("logProb_Z2", params)
    logProbZ1_idx <- grep("logProb_Z1", params)
    
    K  <- length(sigma2_idx)
    N1 <- length(Y1_idx)
    
    # ---- output containers ----
    output_samples <- list()
    output_params  <- list()
    
    # ---- alpha ----
    if (length(alpha_idx) > 0) {
      alpha_names <- params[alpha_idx]
      
      # Parse indices from names like alpha[1, 2]
      alpha_k <- as.integer(sub("^alpha\\[(\\d+),\\s*(\\d+)\\]$", "\\1", alpha_names))
      alpha_p <- as.integer(sub("^alpha\\[(\\d+),\\s*(\\d+)\\]$", "\\2", alpha_names))
      
      p <- max(alpha_p)
      alpha_mat <- samples[, alpha_idx, drop = FALSE]
      
      alpha_arr <- array(NA_real_, dim = c(niter, K, p))
      
      for (j in seq_len(p)) {
        cols_j <- which(alpha_p == j)
        ord_j  <- order(alpha_k[cols_j])
        mat_j  <- alpha_mat[, cols_j[ord_j], drop = FALSE]
        colnames(mat_j) <- as.character(seq_len(K))
        
        output_samples[[paste0("a", j - 1)]] <- mat_j
        output_params[[paste0("a", j - 1)]]  <- colMeans(mat_j)
        alpha_arr[, , j] <- mat_j
      }
      
      output_samples[["alpha"]] <- alpha_arr
    }
    
    # ---- sigma2 ----
    if (length(sigma2_idx) > 0) {
      sigma2_names <- params[sigma2_idx]
      sigma2_k <- as.integer(sub("^sigma2\\[(\\d+)\\]$", "\\1", sigma2_names))
      ord <- order(sigma2_k)
      
      sigma2_mat <- samples[, sigma2_idx[ord], drop = FALSE]
      colnames(sigma2_mat) <- as.character(seq_len(K))
      
      output_samples[["sigma2"]] <- sigma2_mat
      output_params[["sigma2"]]  <- colMeans(sigma2_mat)
    }
    
    # ---- tau2 ----
    if (!is.na(tau2_idx)) {
      tau2_vec <- samples[, tau2_idx]
      output_samples[["tau2"]] <- tau2_vec
      output_params[["tau2"]]  <- mean(tau2_vec)
    }
    
    # ---- sigma2_L ----
    if (!is.na(s2L_idx)) {
      s2L_vec <- samples[, s2L_idx]
      output_samples[["sigma2_L"]] <- s2L_vec
      output_params[["sigma2_L"]]  <- mean(s2L_vec)
    }
    
    # ---- sigma2_D ----
    if (!is.na(s2D_idx)) {
      s2D_vec <- samples[, s2D_idx]
      output_samples[["sigma2_D"]] <- s2D_vec
      output_params[["sigma2_D"]]  <- mean(s2D_vec)
    }
    
    # ---- lD ----
    if (length(lD_idx) > 0) {
      lD_names <- params[lD_idx]
      lD_k <- as.integer(sub("^lD\\[(\\d+)\\]$", "\\1", lD_names))
      ord <- order(lD_k)
      
      lD_mat <- samples[, lD_idx[ord], drop = FALSE]
      colnames(lD_mat) <- as.character(seq_len(ncol(lD_mat)))
      
      output_samples[["lD"]] <- lD_mat
      output_params[["lD"]]  <- colMeans(lD_mat)
    }
    
    # ---- lL ----
    if (length(lL_idx) > 0) {
      lL_names <- params[lL_idx]
      lL_k <- as.integer(sub("^lL\\[(\\d+)\\]$", "\\1", lL_names))
      ord <- order(lL_k)
      
      lL_mat <- samples[, lL_idx[ord], drop = FALSE]
      colnames(lL_mat) <- as.character(seq_len(ncol(lL_mat)))
      
      output_samples[["lL"]] <- lL_mat
      output_params[["lL"]]  <- colMeans(lL_mat)
    }
    
    # ---- gamma ----
    if (length(gamma_idx) > 0) {
      gamma_names <- params[gamma_idx]
      gamma_i <- as.integer(sub("^gammaMat\\[(\\d+),\\s*(\\d+)\\]$", "\\1", gamma_names))
      gamma_j <- as.integer(sub("^gammaMat\\[(\\d+),\\s*(\\d+)\\]$", "\\2", gamma_names))
      ord <- order(gamma_i, gamma_j)
      
      gamma_mat_samples <- samples[, gamma_idx[ord], drop = FALSE]
      output_samples[["gamma"]] <- gamma_mat_samples
      
      gamma_mean <- colMeans(gamma_mat_samples)
      output_params[["gamma"]] <- matrix(gamma_mean, nrow = K, ncol = K, byrow = TRUE)
      
      if (!is.na(h_idx)) {
        h_vec <- samples[, h_idx]
        output_samples[["h"]] <- h_vec
        output_params[["h"]]  <- mean(h_vec)
      }
    }
    
    # ---- log probabilities ----
    if (length(logProbZ2_idx) > 0) {
      output_samples[["logProbZ2"]] <- samples[, logProbZ2_idx[1]]
    }
    if (length(logProbZ1_idx) > 0) {
      output_samples[["logProbZ1"]] <- samples[, logProbZ1_idx[1]]
    }
    if (length(logProbZ1_idx) > 0 || length(logProbZ2_idx) > 0) {
      lp1 <- if (length(logProbZ1_idx) > 0) samples[, logProbZ1_idx[1]] else 0
      lp2 <- if (length(logProbZ2_idx) > 0) samples[, logProbZ2_idx[1]] else 0
      output_samples[["logProb"]] <- lp1 + lp2
    }
    
    # ---- Y1 ----
    if (length(Y1_idx) > 0) {
      Y1_names <- params[Y1_idx]
      Y1_site  <- as.integer(sub("^Y1\\[(\\d+)\\]$", "\\1", Y1_names))
      ord <- order(Y1_site)
      
      Y1_mat <- samples[, Y1_idx[ord], drop = FALSE]
      colnames(Y1_mat) <- paste0("Y1[", seq_len(ncol(Y1_mat)), "]")
      output_samples[["Y1"]] <- Y1_mat
      
      # assumes labels are positive integers, typically 1:K
      classes <- seq_len(K)
      Y1_prob_mat <- sapply(classes, function(k) colMeans(Y1_mat == k))
      
      if (is.vector(Y1_prob_mat)) {
        Y1_prob_mat <- matrix(Y1_prob_mat, ncol = 1)
      }
      
      colnames(Y1_prob_mat) <- as.character(classes)
      rownames(Y1_prob_mat) <- colnames(Y1_mat)
      
      output_params[["Y1"]] <- unname(max.col(Y1_prob_mat, ties.method = "first"))
      output_params[["Y1prob"]] <- Y1_prob_mat
    }
    return(list(params = output_params, samples = output_samples))
  }
}