#' Produce posterior predictions from GeoMix MCMC samples
#'
#' Generates posterior predictive summaries for unobserved `Z2` locations using
#' saved GeoMix MCMC samples. Predictions are computed conditional on each
#' retained posterior draw using compiled prediction routines, then aggregated
#' across iterations and optionally across chains.
#'
#' The function can return full posterior predictive samples or only posterior
#' mean and standard deviation summaries.
#'
#' @param samples Either:
#'   \describe{
#'     \item{matrix-like object}{Posterior sample matrix for a single chain,
#'     with rows corresponding to iterations and columns to monitored
#'     parameters.}
#'     \item{list}{A list of such sample objects, typically one per chain.}
#'   }
#'
#' @param geomix_setup A GeoMix setup object containing the model data,
#'   constants, design matrix, grouping structure, and distance objects needed
#'   for prediction.
#'
#' @param nugget Logical indicating whether predictive variances should include
#'   the nugget or noise variance component. Defaults to `TRUE`.
#'
#' @param include_samples Logical indicating whether to generate full posterior
#'   predictive samples. If `FALSE`, only predictive
#'   means and standard deviations are returned using moment identities.
#'   Defaults to `TRUE`.
#'
#' @param run_parallel Logical indicating whether posterior draws should be
#'   processed in parallel using `parallel::mclapply()`. Defaults to `FALSE`.
#'
#' @param mc.cores Integer giving the number of CPU cores to use when
#'   `run_parallel = TRUE`. Defaults to `parallel::detectCores()`.
#'
#' @param predict_index Optional integer vector giving indices of locations to
#'   predict. Defaults to `NULL`, in which case all rows with missing `Z2`
#'   values in `geomix_setup$df$Z2` are predicted.
#'
#' @return A named list with components:
#' \describe{
#'   \item{mean}{Numeric vector of posterior predictive means at each predicted
#'   location.}
#'   \item{sd}{Numeric vector of posterior predictive standard deviations at
#'   each predicted location.}
#'   \item{samples}{Matrix of posterior predictive samples with rows
#'   corresponding to prediction locations and columns to retained posterior
#'   draws, or `NULL` when `include_samples = FALSE`.}
#' }
#'
#' @details
#' For each retained posterior draw, the function:
#' \enumerate{
#'   \item extracts model parameters
#'   \item constructs prediction inputs for the requested locations
#'   \item evaluates predictive moments
#'   \item optionally draws predictive samples
#' }
#'
#' When `include_samples = TRUE`, posterior predictive summaries are computed
#' directly from the simulated predictive draws.
#'
#' When `include_samples = FALSE`, predictive uncertainty is approximated using
#' the law of total variance:
#' \deqn{
#' \mathrm{Var}(Y \mid data)
#' =
#' E[\mathrm{Var}(Y \mid \theta, data)]
#' +
#' \mathrm{Var}(E[Y \mid \theta, data]).
#' }
#'
#' If multiple chains are supplied:
#' \itemize{
#'   \item predictive samples are concatenated across chains when
#'   `include_samples = TRUE`
#'   \item posterior means and standard deviations are averaged across chains
#'   when `include_samples = FALSE`
#' }
#'
#' The function includes internal validation checks for prediction dimensions,
#' finite outputs, and compiled routine failures.
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
#' # Predict all 6400 unobserved Z2 locations
#' pred <- produce_prediction(fit$samples, setup)
#' length(pred$mean)  # 6400
#' range(pred$mean)
#' range(pred$sd)
#'
#' # Attach predictions back to the data frame
#' pred_df <- setup$df |>
#'   dplyr::filter(is.na(Z2)) |>
#'   dplyr::mutate(pred_mean = pred$mean, pred_sd = pred$sd)
#'
#' # Fast moment-only prediction (no predictive draws)
#' pred_fast <- produce_prediction(
#'   fit$samples,
#'   setup,
#'   include_samples = FALSE
#' )
#'
#' # Predict a subset of locations
#' pred_sub <- produce_prediction(
#'   fit$samples,
#'   setup,
#'   predict_index = which(is.na(offshore$Z2))[1:100]
#' )
#' }
#'
#' @seealso [run_chains()], [extract_parameters()], [parallel::mclapply()]
#'
#' @export
produce_prediction <- function(samples,
                               geomix_setup,
                               nugget = TRUE,
                               include_samples = TRUE,
                               run_parallel = F,
                               mc.cores = NULL,
                               predict_index = NULL) {

  is_sample_list <- is.list(samples)
  if(run_parallel & is.null(mc.cores)) mc.cores <- parallel::detectCores()
  build_prediction_inputs <- function(params, geomix_setup, predict_index, nugget = TRUE, include_samples = TRUE) {
    if (is.null(predict_index)) {
      predict_index <- which(is.na(geomix_setup$df$Z2))
    }

    # Precompute once
    predNeighbours <- cbind(
      geomix_setup$data_list$groupLookup,
      geomix_setup$data_list$groupNeighbours
    )
    predNeighbours <- t(apply(predNeighbours, 1, function(x) c(x[x > 0], x[x == 0])))
    predL <- rowSums(predNeighbours > 0)

    pred_group <- geomix_setup$df[[geomix_setup$grouping]][predict_index]

    nsamples <- dim(params$samples$alpha)[1]

    if(is.null(geomix_setup$data_list$LFlag)){
      LFlag <- rep(1,geomix_setup$constants$K)
    }else{
      LFlag <- geomix_setup$data_list$LFlag
    }
    # Pre-extract static pieces (so we don’t repeatedly index geomix_setup)
    static <- list(
      pred_inds = predict_index,
      pred_groups = pred_group,
      K = geomix_setup$constants$K,
      G = geomix_setup$constants$G,
      m = geomix_setup$constants$m,
      Z2 = geomix_setup$data_list$Z2,
      X = geomix_setup$X,
      Z2_ind = geomix_setup$constants$Z2_ind,
      dID = geomix_setup$data_list$dID,
      locID = geomix_setup$data_list$locID,
      LFlag = LFlag,
      distD = geomix_setup$data_list$distD,
      distL = geomix_setup$data_list$distL,
      predNeighbours = predNeighbours,
      predL = predL,
      nugget = as.integer(nugget),
      include_samples = include_samples
    )

    # Build per-iteration inputs
    inputs <- vector("list", nsamples)

    for (j in seq_len(nsamples)) {
      inputs[[j]] <- c(static, list(
        alpha  = params$samples$alpha[j, , ],
        sigma2 = params$samples$sigma2[j, ],
        tau2   = params$samples$tau2[j],
        lL     = params$samples$lL[j, ],
        lD     = params$samples$lD[j, ],
        Y1     = params$samples$Y1[j, ]
      ))
    }

    return(inputs)
  }

  predict_one_chain <- function(samples_chain) {
    if (!is.matrix(samples_chain)) {
      samples_chain <- as.matrix(samples_chain)
    }

    compiled <- .get_compiled_predictors()
    CpredictGPvec <- compiled$CpredictGPvec
    CsampleGPvec  <- compiled$CsampleGPvec

    cat("\nProcessing samples...\n")
    params <- extract_parameters(samples_chain)
    nsamples <- nrow(samples_chain)

    inputs <- build_prediction_inputs(params,geomix_setup,predict_index,nugget)
    predict_sample_iter <- function(x, include_samples = TRUE, iter = NA_integer_) {
      samps <- NULL

      # basic checks
      stopifnot(is.list(x))
      stopifnot(length(x$pred_inds) == length(x$pred_groups))
      stopifnot(all(is.finite(x$pred_inds)))
      stopifnot(all(x$pred_inds >= 1))
      stopifnot(all(x$pred_inds <= length(x$Y1)))

      pars <- tryCatch(
        CpredictGPvec(
          pred_inds = x$pred_inds,
          pred_groups = x$pred_groups,
          alpha = x$alpha,
          sigma2 = x$sigma2,
          tau2 = x$tau2,
          lL = x$lL,
          lD = x$lD,
          LFlag = x$LFlag,
          Y1 = x$Y1,
          K = x$K,
          G = x$G,
          m = x$m,
          Z2 = x$Z2,
          X = x$X,
          Z2_ind = x$Z2_ind,
          dID = x$dID,
          locID = x$locID,
          distD = x$distD,
          distL = x$distL,
          predNeighbours = x$predNeighbours,
          predL = x$predL,
          nugget = as.integer(x$nugget)
        ),
        error = function(e) {
          stop(sprintf("CpredictGPvec failed at iter %s: %s", iter, e$message))
        }
      )

      # validate predictor output
      if (is.null(pars$mean) || is.null(pars$sd) || is.null(pars$cov)) {
        stop(sprintf("Missing output from CpredictGPvec at iter %s", iter))
      }

      if (any(!is.finite(pars$mean))) {
        stop(sprintf("Non-finite pars$mean at iter %s", iter))
      }

      if (any(!is.finite(pars$sd))) {
        stop(sprintf("Non-finite pars$sd at iter %s", iter))
      }

      if (any(!is.finite(pars$cov))) {
        stop(sprintf("Non-finite pars$cov at iter %s", iter))
      }

      n_pred <- length(x$pred_inds)

      if (length(pars$mean) != n_pred) {
        stop(sprintf(
          "Length mismatch at iter %s: length(mean)=%d, n_pred=%d",
          iter, length(pars$mean), n_pred
        ))
      }

      if (!is.matrix(pars$cov) || nrow(pars$cov) != n_pred || ncol(pars$cov) != n_pred) {
        stop(sprintf(
          "Covariance dimension mismatch at iter %s: got %d x %d, expected %d x %d",
          iter, nrow(pars$cov), ncol(pars$cov), n_pred, n_pred
        ))
      }

      if (include_samples) {
        samps <- tryCatch(
          CsampleGPvec(
            m = pars$mean,
            cov = pars$cov,
            Y1 = x$Y1[x$pred_inds],
            groups = x$pred_groups,
            K = x$K,
            G = x$G
          ),
          error = function(e) {
            stop(sprintf("CsampleGPvec failed at iter %s: %s", iter, e$message))
          }
        )
      }

      list(
        mean = pars$mean,
        sd = pars$sd,
        samples = samps
      )
    }
    cat("\nComputing prediction terms and running samples...\n")
    if(run_parallel){
      post_params <- parallel::mclapply(inputs,predict_sample_iter,mc.cores = mc.cores)
    }else{
      post_params <- lapply(inputs,predict_sample_iter)
    }
    out <- list()

    if (include_samples) {
      post_samples <- post_mean <- do.call(cbind, lapply(post_params, `[[`, "samples"))

      out$samples <- post_samples
      out$mean <- rowMeans(post_samples)
      out$sd <- apply(post_samples, 1, stats::sd)
    } else {
      post_mean <- do.call(cbind, lapply(post_params, `[[`, "mean"))
      post_sd   <- do.call(cbind, lapply(post_params, `[[`, "sd"))

      mean_cond_var <- rowMeans(post_sd^2)
      var_cond_mean <- apply(post_mean, 1, stats::var)

      out$samples <- NULL
      out$mean <- rowMeans(post_mean)
      out$sd <- sqrt(mean_cond_var + var_cond_mean)
    }

    cat("\nDone\n")
    out
  }

  if (is_sample_list) {
    pred_list <- vector("list", length(samples))

    for (ch in seq_along(samples)) {
      cat(paste0("\nChain ", ch, "\n"))
      pred_list[[ch]] <- predict_one_chain(samples[[ch]])
    }

    if (include_samples) {
      all_samples <- do.call(cbind, lapply(pred_list, `[[`, "samples"))

      return(list(
        mean = rowMeans(all_samples),
        sd = apply(all_samples, 1, stats::sd),
        samples = all_samples
      ))
    } else {
      means_mat <- do.call(cbind, lapply(pred_list, `[[`, "mean"))
      sds_mat   <- do.call(cbind, lapply(pred_list, `[[`, "sd"))

      return(list(
        mean = rowMeans(means_mat),
        sd = rowMeans(sds_mat),
        samples = NULL
      ))
    }
  }

  predict_one_chain(samples)
}
