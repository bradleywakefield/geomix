#' Load saved MCMC sample batches from one or more chain directories
#'
#' Loads saved MCMC sample batches from subdirectories within a given path and
#' combines the selected batches within each directory by row-binding them.
#' Subdirectories are filtered by a common prefix, typically corresponding to a
#' model name such as `"GeoMix"`.
#'
#' This function is useful when MCMC output has been saved in separate batch
#' files across multiple directories, for example one directory per chain.
#'
#' @param path Character string giving the parent directory containing the
#'   saved MCMC output directories.
#'
#' @param name Character string giving the prefix used to identify relevant
#'   subdirectories. Only directories whose names begin with this prefix are
#'   loaded. Defaults to `"GeoMix"`.
#'
#' @param index Optional integer vector specifying which batch files to load
#'   from each matching directory, after batch files have been ordered by their
#'   numeric batch index. Defaults to `NULL`, in which case all available
#'   batches are loaded.
#'
#' @param thin Optional positive integer. If supplied, every `thin`-th row is
#'   retained from the combined sample matrix (i.e. rows `thin`, `2*thin`,
#'   `3*thin`, ...). Defaults to `NULL` (no thinning).
#'
#' @return A list with one element per matching subdirectory. Each element is
#'   the result of row-binding the selected batch objects read via `readRDS()`.
#'   In typical use, each element will be a posterior sample matrix or similar
#'   object with rows corresponding to retained MCMC iterations.
#'
#' @details
#' The function:
#' \itemize{
#'   \item lists all subdirectories in `path`
#'   \item keeps only those whose names begin with `name`
#'   \item finds files matching `"batch"` within each selected directory
#'   \item orders batch files by the numeric suffix extracted from names such as
#'   `batch_1`, `batch_2`, and so on
#'   \item reads the selected files using `readRDS()`
#'   \item combines them using `rbind()`
#' }
#'
#' The return value is not explicitly named by directory, so if directory names
#' are needed downstream it may be helpful to assign them manually.
#'
#' Note that `index` is interpreted relative to the ordered batch list within
#' each directory.
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
#' # Run chains with batch saving enabled
#' run_chains(
#'   setup,
#'   nchains     = 2,
#'   path        = "results/mcmc",
#'   controlMCMC = list(
#'     niter        = 10000,
#'     thin         = 5,
#'     nbatches     = 10,
#'     save_batches = TRUE,
#'     retain_draws = FALSE
#'   ),
#'   seed = 42
#' )
#'
#' # Reload all saved batches
#' samples <- load_mcmc_samples("results/mcmc", name = "GeoMix")
#'
#' # Reload only the last five batches (useful after resuming a run)
#' samples_sub <- load_mcmc_samples(
#'   path  = "results/mcmc",
#'   name  = "GeoMix",
#'   index = 6:10
#' )
#'
#' # Reload all batches but keep every 5th iteration
#' samples_thin <- load_mcmc_samples(
#'   path = "results/mcmc",
#'   name = "GeoMix",
#'   thin = 5
#' )
#' }
#'
#' @seealso [run_chains()], [extract_parameters()]
#'
#' @export
load_mcmc_samples <- function(path, name = "GeoMix", index = NULL, thin = NULL){
  dirs <- list.dirs(path, full.names = FALSE)
  dirs <- dirs[grepl(paste0("^",name), dirs)]
  samples <- lapply(dirs, function(dir){
    batches <- list.files(file.path(path, dir), pattern = "batch", full.names = TRUE)
    batches <- batches[order(as.numeric(str_remove(str_extract(batches, "batch\\_\\d+"), "batch\\_")))]
    nbatches <- length(batches)
    if(is.null(index)) index <- 1:nbatches
    mat <- do.call(rbind, lapply(batches[index], readRDS))
    if(!is.null(thin)){
      thin <- as.integer(thin)
      keep <- seq(thin, nrow(mat), by = thin)
      mat <- mat[keep, , drop = FALSE]
    }
    mat
  })
}

#' Combine extracted MCMC outputs across chains
#'
#' Combines a list of extracted MCMC outputs, typically one per chain, into a
#' single object containing averaged posterior summaries and concatenated
#' posterior samples.
#'
#' The function assumes each element of `x` has the same structure as the
#' output from `extract_parameters()`, with top-level components `params` and
#' `samples`.
#'
#' Posterior summaries in `params` are averaged across chains, while posterior
#' samples in `samples` are combined by stacking iterations across chains.
#'
#' @param x A list of extracted parameter objects, usually one per chain. Each
#'   element must contain:
#'   \describe{
#'     \item{params}{A named list of posterior summaries.}
#'     \item{samples}{A named list of posterior sample objects.}
#'   }
#'
#' @return A named list with components:
#' \describe{
#'   \item{params}{A list of averaged posterior summaries across chains.}
#'   \item{samples}{A list of posterior sample objects combined across chains.}
#' }
#'
#' The returned `params` component contains the elementwise average of each
#' summary in `x[[i]]$params`, excluding `Y1`, which is reconstructed from the
#' averaged posterior class probabilities `Y1prob`.
#'
#' Specifically:
#' \describe{
#'   \item{Y1prob}{Averaged across chains elementwise.}
#'   \item{Y1}{Recomputed as the sitewise most probable class using the rows of
#'   `Y1prob`.}
#' }
#'
#' The returned `samples` component contains chain-combined posterior draws for
#' each named sample object:
#' \describe{
#'   \item{vectors and matrices}{Combined using `rbind()`.}
#'   \item{arrays}{Combined along the first dimension using
#'   `abind::abind(..., along = 1)`.}
#' }
#'
#' @details
#' This function is intended for post-processing after per-chain outputs have
#' already been converted into a structured form. It assumes that all chains
#' have the same parameter names, sample names, and compatible dimensions.
#'
#' Averaging of `params` is performed elementwise using `purrr::map2()`. This
#' works naturally for numeric scalars, vectors, matrices, and arrays so long
#' as the corresponding objects have the same dimensions in every chain.
#'
#' The latent class summary `Y1` is not averaged directly. Instead, it is
#' reconstructed from the averaged posterior probability matrix `Y1prob` by
#' choosing the class with highest posterior probability at each site.
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
#' # Extract per-chain posteriors then combine
#' params_list <- extract_parameters(fit$samples)
#' combined    <- combine_chains(params_list)
#'
#' combined$params$sigma2   # averaged posterior mean variances
#' combined$params$Y1[1:10] # pooled modal class at first 10 sites
#' table(combined$params$Y1)  # class counts across the full lattice
#'
#' dim(combined$samples$alpha)  # (total iterations) x K x p
#' }
#'
#' @seealso [extract_parameters()], [run_mcmc_diagnostics()], [produce_prediction()]
#'
#' @export
combine_chains <- function(x) {
  
  n_chains <- length(x)
  
  # ---------------------------
  # 1. Average params
  # ---------------------------
  
  param_names <- setdiff(names(x[[1]]$params), "Y1")
  
  avg_params <- x |>
    map(~ .x$params[param_names]) |>
    reduce(~ map2(.x, .y, `+`)) |>
    map(~ .x / n_chains)
  
  avg_params$Y1 <- unname(apply(avg_params$Y1prob, 1, \(z) as.integer(which.max(z))))
  
  # ---------------------------
  # 2. Combine samples
  # ---------------------------
  
  sample_names <- names(x[[1]]$samples)
  
  combined_samples <- sample_names |>
    set_names() |>
    map(function(nm) {
      objs <- map(x, ~ .x$samples[[nm]])
      
      # handle arrays vs matrices/vectors
      if (length(dim(objs[[1]])) <= 2) {
        do.call(rbind, objs)
      } else {
        # for arrays (e.g. alpha: iter x K x p)
        abind::abind(objs, along = 1)
      }
    })
  
  # ---------------------------
  # Return
  # ---------------------------
  
  list(
    params  = avg_params,
    samples = combined_samples
  )
}