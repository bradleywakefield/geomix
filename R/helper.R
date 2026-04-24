#' Apply the Box-Cox transformation
#'
#' Internal helper for transforming positive responses.
#'
#' @param y Numeric vector.
#' @param lambda Transformation parameter. Defaults to `0.6`.
#'
#' @return Numeric vector of transformed values.
#' @noRd
box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda

#' Apply the inverse Box-Cox transformation
#'
#' Internal helper for back-transforming Box-Cox transformed values.
#'
#' @param y Numeric vector on the transformed scale.
#' @param lambda Transformation parameter. Defaults to `0.6`.
#'
#' @return Numeric vector on the original scale.
#' @noRd
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)

#' Compute the modal value
#'
#' Internal helper returning the most frequent value in a vector.
#'
#' @param x A vector.
#'
#' @return A length-one character vector giving the modal value.
#' @noRd
get_mode <- function(x) names(sort(table(x), decreasing = TRUE))[1]

#' Compute the empirical mode proportion
#'
#' Internal helper returning the proportion of observations equal to the mode.
#'
#' @param x A vector.
#'
#' @return Numeric scalar giving the relative frequency of the modal value.
#' @noRd
get_mode_prob <- function(x) max(table(x))/length(x)

#' Retrieve cached compiled NIMBLE prediction routines
#'
#' Internal helper that lazily compiles and caches NIMBLE prediction
#' functions used by GeoMix posterior prediction workflows. Compilation is
#' performed only on the first call in an R session, with compiled objects
#' reused on subsequent calls.
#'
#' @details
#' This helper is intended to avoid repeated calls to
#' `nimble::compileNimble()` for prediction routines, which can be expensive.
#' The compiled objects are session-specific and are recreated automatically
#' after restarting R.
#'
#' It compiles and returns cached versions of:
#' \describe{
#'   \item{CpredictGPvec}{Compiled prediction routine based on `predictGPvec`.}
#'   \item{CsampleGPvec}{Compiled sampling routine based on `sampleGPvec`.}
#' }
#'
#' @return A named list containing compiled NIMBLE functions.
#'
#' @keywords internal
#' @noRd
.get_compiled_predictors <- local({
  CpredictGPvec <- NULL
  CsampleGPvec <- NULL

  function() {
    if (is.null(CpredictGPvec)) {
      CpredictGPvec <<- nimble::compileNimble(predictGPvec)
    }
    if (is.null(CsampleGPvec)) {
      CsampleGPvec <<- nimble::compileNimble(sampleGPvec)
    }
    list(
      CpredictGPvec = CpredictGPvec,
      CsampleGPvec = CsampleGPvec
    )
  }
})
