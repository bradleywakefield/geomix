#' Synthetic offshore CPT dataset
#'
#' A synthetic dataset representing a 20 x 20 x 20 subsurface lattice designed
#' to illustrate the GeoMix workflow for offshore site characterisation. The
#' domain contains 400 CPT locations on a regular lateral grid, each measured
#' at 20 depth increments, giving 8 000 rows in total.
#'
#' Ground model units (`Z1`) are observed at all locations. The geotechnical
#' property (`Z2`) is observed at 1 600 locations (20% of sites), with the
#' remaining 6 400 treated as prediction targets.
#'
#' @format A data frame with 8 000 rows and 7 columns:
#' \describe{
#'   \item{ID}{Integer row identifier (1–8000).}
#'   \item{d}{Depth layer index (1–20).}
#'   \item{x}{Lateral x grid index (1–20).}
#'   \item{y}{Lateral y grid index (1–20).}
#'   \item{locID}{CPT location identifier (1–400).}
#'   \item{Z1}{Observed ground model unit, an integer in \{1, 2, 3\}.}
#'   \item{Z2}{Observed geotechnical property (cone tip resistance). Numeric,
#'     with \code{NA} at unobserved locations.}
#' }
#'
#' @examples
#' data(offshore)
#' head(offshore)
#' table(offshore$Z1)
#' sum(!is.na(offshore$Z2))  # number of observed CPT measurements
#'
#' @source Simulated from a three-class GeoMix model with known parameters.
#'   Full simulation details and true parameter values are available in
#'   \code{system.file("extdata", "offshoreFull.rds", package = "geomix")}.
"offshore"
