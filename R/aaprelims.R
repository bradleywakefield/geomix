#' @keywords internal
"_PACKAGE"
#' @import nimble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import stringr
#' @import purrr
#' @import future
#' @import parallel
#' @import posterior
#' @importFrom bayesplot mcmc_trace mcmc_acf mcmc_rank_overlay
#' @importFrom tibble tibble
#' @importFrom GpGp order_maxmin
#' @importFrom FNN get.knnx
#' @importFrom proxy dist
#' @importFrom clue solve_LSAP
#' @importFrom abind abind
#' @importFrom methods is
#' @importFrom stats dnorm model.matrix median optim pnorm qnorm rnorm
#' @importFrom utils capture.output head modifyList tail
#' @importFrom magrittr %>%

# Suppress R CMD check notes for NSE variables used in dplyr/ggplot2 pipelines
utils::globalVariables(c(
  "latticeID", "rowID", "loc", "x", "y", "minD", "maxD",
  "ind", "iter", "chain", "Y1count", "param", "mcse_over_sd",
  "Z2", "logProb", "name", "value", "class", "count"
))

NULL

# confirm_run <- function(x) {
#   message = paste0("\n Welcome to the GeoMix Workflow.",
#                    "\n Note the runtime on the MCMC is substantial.",
#                    "\n To generate results you can preload previous MCMC chain runs.",
#                    "\n For ",x," \n Would you like to re-run the MCMC chains again?\n[y/n]:")
#   response <- tolower(trimws(readline(prompt = message)))
#
#   if (response %in% c("y", "yes")) {
#     out <- T
#   }else{
#     out <- F
#   }
#   out
# }
#
# # pkgs <- c("posterior", "bayesplot", "ggplot2", "coda", "dplyr", "tidyr")
# # to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
# # if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
#
#
