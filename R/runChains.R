#' Run One or More GeoMix MCMC Chains
#'
#' Runs one or more Markov chain Monte Carlo (MCMC) chains for a fitted
#' GeoMix model object. Chains may be run sequentially or in parallel,
#' with optional batch-wise saving to disk and
#' optional continuation from a previously saved state.
#'
#' The function constructs the NIMBLE model, configures custom samplers,
#' compiles the model and MCMC objects, and executes the requested number
#' of iterations split across batches.
#'
#' @param geomix_setup A list containing all model inputs required to build and
#' run the model. Typically includes elements such as `constants`,
#' `data_list`, `inits`, `controlGibbs`, `controlHMC`, and optionally
#' `controlMCMC`.
#'
#' @param nchains Integer. Number of independent chains to run.
#' Default is `1`.
#'
#' @param path Optional character string giving a directory where batch outputs
#' should be saved when `save_batches = TRUE`.
#' Default is `NULL`.
#'
#' @param controlMCMC Optional named list of MCMC controls. Values supplied here
#' override defaults and any values stored in `geomix_setup$controlMCMC`.
#' Supported fields are:
#' \describe{
#'   \item{niter}{Total number of MCMC iterations per chain.}
#'   \item{thin}{Thinning interval passed to `runMCMC()`.}
#'   \item{nbatches}{Number of batches used to split iterations.}
#'   \item{save_batches}{Logical; save each batch to disk.}
#'   \item{retain_draws}{Logical; retain posterior draws in memory.}
#' }
#'
#' @param LGFM Logical. If `TRUE`, runs the LGFM model (GeoMix without Z1 data model)
#' instead of GeoMix. Default is `FALSE`.
#'
#' @param run_parallel Logical. If `TRUE`, chains are run in parallel using
#' \pkg{future} with multisession workers. Default is `FALSE`.
#'
#' @param load_previous_state Logical. If `TRUE`, attempts to resume each chain
#' from previously saved batch files in `path`. Requires
#' `save_batches = TRUE`. Default is `FALSE`.
#'
#' @param mc.cores Optional integer giving the number of parallel workers when
#' `run_parallel = TRUE`. Defaults to the minimum of `nchains` and available
#' cores.
#'
#' @param seed Integer random seed used to initialise chains. Chain `j` uses
#' seed `seed + j`.
#'
#' @details
#' Iterations are divided into approximately equal batch sizes using
#' `nbatches`. After each batch:
#' \itemize{
#'   \item posterior draws may be saved to disk;
#'   \item draws may be retained in memory;
#'   \item in-memory sample storage is cleared to reduce RAM usage.
#' }
#'
#' When `load_previous_state = TRUE`, the function loads the most recent saved
#' batch for each chain, reconstructs parameter initial values, and continues
#' sampling from that state.
#'
#' Parallel execution uses `future::plan(multisession)`.
#'
#' @return
#' If `retain_draws = TRUE`, returns a list with components:
#' \describe{
#'   \item{chains}{Per-chain metadata and results.}
#'   \item{samples}{Named list of posterior draw matrices, one per chain.}
#' }
#'
#' Otherwise returns only the per-chain result list.
#'
#' Each chain result contains:
#' \describe{
#'   \item{chain}{Chain index.}
#'   \item{name}{Chain name `GeoMix_j` (or `LGFM_j` when `LGFM == TRUE`).}
#'   \item{path}{Output directory used for saved batches.}
#'   \item{batches_run}{Batch indices completed.}
#'   \item{batch_sizes}{Iterations per batch.}
#'   \item{controlMCMC}{Resolved MCMC control list.}
#'   \item{samples}{Posterior draws matrix, if retained.}
#' }
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
#' # Single chain
#' fit <- run_chains(setup, nchains = 1, seed = 42)
#'
#' # Two chains in parallel
#' fit <- run_chains(setup, nchains = 2, run_parallel = TRUE, seed = 42)
#'
#' # Access posterior samples
#' dim(fit$samples$GeoMix_1)   # iterations x parameters
#'
#' # Save intermediate batches to disk (useful for long runs)
#' fit <- run_chains(
#'   setup,
#'   nchains     = 2,
#'   path        = "results/mcmc",
#'   controlMCMC = list(
#'     niter        = 10000,
#'     thin         = 5,
#'     nbatches     = 10,
#'     save_batches = TRUE
#'   ),
#'   seed = 42
#' )
#'
#' # Resume a run from the last saved batch
#' fit2 <- run_chains(
#'   setup,
#'   nchains              = 2,
#'   path                 = "results/mcmc",
#'   load_previous_state  = TRUE,
#'   controlMCMC          = list(save_batches = TRUE)
#' )
#' }
#'
#' @seealso [setupGeoMixModel()], [extract_parameters()], [load_mcmc_samples()],
#'   \code{\link[nimble:nimble-package]{nimble}}
#'
#' @export
run_chains <- function(geomix_setup,
                       nchains = 1,
                       path = NULL,
                       controlMCMC = NULL,
                       LGFM = FALSE,
                       run_parallel = FALSE,
                       load_previous_state = FALSE,
                       mc.cores = NULL,
                       seed = 16) {
  message('Processing data...',appendLF = F)
  `%||%` <- function(a,b) if(!is.null(a)) a else b
  make_batch_sizes <- function(niter, nbatches) {
    q <- niter %/% nbatches
    r <- niter %% nbatches
    c(rep(q + 1, r), rep(q, nbatches - r))
  }
  default_controlMCMC <- list(
    niter = 2000,
    thin = 1,
    nbatches = 1,
    save_batches = FALSE,
    retain_draws = TRUE
  )

  control_arg <- controlMCMC

  controlMCMC <- modifyList(
    default_controlMCMC,
    geomix_setup$controlMCMC %||% list()
  )

  controlMCMC <- modifyList(
    controlMCMC,
    control_arg %||% list()
  )

  if (run_parallel) {
    mc.cores <- mc.cores %||% min(nchains, parallel::detectCores())
    mc.cores <- min(as.integer(mc.cores), nchains)
  }

  required_fields <- c("niter", "thin", "nbatches", "save_batches", "retain_draws" )
  missing_fields <- setdiff(required_fields, names(controlMCMC))
  if (length(missing_fields) > 0) {
    stop("controlMCMC is missing: ", paste(missing_fields, collapse = ", "))
  }

  if (!controlMCMC$save_batches & !controlMCMC$retain_draws) {
    stop("At least one of controlMCMC$save_batches or controlMCMC$retain_draws must be TRUE.")
  }

  if (controlMCMC$save_batches & is.null(path)) {
    stop("path must be supplied when controlMCMC$save_batches = TRUE.")
  }

  if (controlMCMC$niter <= 0 | controlMCMC$thin <= 0 | controlMCMC$nbatches <= 0) {
    stop("controlMCMC$niter, controlMCMC$thin, and controlMCMC$nbatches must all be positive.")
  }

  if (nchains <= 0) {
    stop("nchains must be positive.")
  }

  if (controlMCMC$nbatches > controlMCMC$niter) {
    warning("controlMCMC$nbatches > controlMCMC$niter, reducing nbatches to niter.")
    controlMCMC$nbatches <- controlMCMC$niter
  }

  batch_sizes <- make_batch_sizes(
    niter = controlMCMC$niter,
    nbatches = controlMCMC$nbatches
  )
  inputs <- list(
    seed = seed,
    LGFM = LGFM,
    path = path,
    load_previous_state = load_previous_state,
    geomix_setup = geomix_setup,
    controlMCMC = controlMCMC,
    batch_sizes = batch_sizes
  )

  run_one_chain <- function(chain, inputs) {
    geomix_local <- inputs$geomix_setup
    batch_start <- 1
    set.seed(inputs$seed + chain)

    if (inputs$LGFM) {
      name <- paste0("LGFM_", chain)
    } else {
      name <- paste0("GeoMix_", chain)
    }

    if (inputs$controlMCMC$save_batches) {
      dir <- file.path(inputs$path, name)
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    } else {
      dir <- NULL
    }
    if (inputs$LGFM) {
      local_GibbsSampler <- GibbsSamplerLGFM
      if(geomix_local$constants$p == 1){
        local_code <- codeLGFMp1
      }else{
        local_code <- codeLGFM
      }
      geomix_local$controlGibbs$Z1 <- geomix_local$data_list$Z1
    }else{
      local_GibbsSampler <- GibbsSampler
      if(geomix_local$constants$p == 1){
        local_code <- codep1
      }else{
        local_code <- code
      }
    }

    if (!inputs$load_previous_state) {
      geomix_local$inits$Y1 <- sample(
        1:geomix_local$constants$K,
        size = geomix_local$constants$N1,
        replace = TRUE
      )
      if (inputs$controlMCMC$save_batches) {
        saveRDS(geomix_local, file.path(dir, "geomix_setup.rds"))
      }
    } else {
      if (is.null(dir)) {
        stop("load_previous_state = TRUE requires controlMCMC$save_batches = TRUE and a valid path.")
      }

      geomix_local <- readRDS(file.path(dir, "geomix_setup.rds"))

      batch_files <- list.files(dir, pattern = "^batch_\\d+\\.rds$", full.names = TRUE)
      if (length(batch_files) == 0) {
        stop("No previous batch files found in ", dir)
      }

      batch_numbers <- as.numeric(stringr::str_extract(basename(batch_files), "\\d+"))
      samps <- readRDS(batch_files[which.max(batch_numbers)])
      last <- samps[nrow(samps), ]

      grab_node <- function(x, var) {
        x[grep(paste0("^", var, "\\["), names(x))]
      }

      batch_start <- max(batch_numbers, na.rm = TRUE) + 1
      K <- geomix_local$constants$K

      geomix_local$inits <- list(
        tau2     = unname(last["tau2"]),
        h        = unname(last["h"]),
        sigma2_L = unname(last["sigma2_L"]),
        sigma2_D = unname(last["sigma2_D"]),
        alpha    = matrix(
          unname(grab_node(last, "alpha")),
          nrow = K,
          ncol = geomix_local$constants$p,
          byrow = TRUE
        ),
        sigma2   = unname(grab_node(last, "sigma2")),
        lL       = unname(grab_node(last, "lL")),
        lD       = unname(grab_node(last, "lD")),
        gammaMat = matrix(
          unname(grab_node(last, "gammaMat")),
          nrow = K,
          ncol = K,
          byrow = TRUE
        ),
        Y1LogProb = 1,
        Y1 = unname(grab_node(last, "Y1"))
      )

      if (inputs$controlMCMC$save_batches) {
        saveRDS(geomix_local, file.path(dir, "geomix_setup_reloaded.rds"))
      }
    }
    message('Done.')
    message('Building model...', appendLF = F)
    model <- suppressMessages({nimble::nimbleModel(
      local_code,
      constants   = geomix_local$constants,
      data        = geomix_local$data_list,
      inits       = geomix_local$inits,
      buildDerivs = FALSE
    )})

    Cmodel <- suppressMessages({nimble::compileNimble(model)})
    message('Done.')
    message('Setting up samplers...',appendLF = F)
    conf <- nimble::configureMCMC(model, print = FALSE)
    conf$replaceSampler(
      target  = "Y1",
      type    = local_GibbsSampler,
      control = geomix_local$controlGibbs
    )
    conf$replaceSampler(
      target  = "alpha",
      type    = alphaGibbsSampler,
      control = list(m_alpha = geomix_local$constants$m_alpha,
                     Q_alpha = geomix_local$constants$Q_alpha)
    )

    scalarTargets <- c("tau2", "sigma2_L", "sigma2_D", "h")
    scalarTargetsExpanded <- model$expandNodeNames(scalarTargets)
    for (tar in scalarTargetsExpanded) {
      conf$replaceSampler(
        target  = tar,
        type    = HMCSample1D_logFD,
        control = geomix_local$controlHMC
      )
    }

    Z2scalarTargets <- c("sigma2", "lL", "lD")
    Z2scalarTargetsExpanded <- model$expandNodeNames(Z2scalarTargets)
    noUpdate <- character(0)
    if (!is.null(geomix_local$fix_lateral)) {
      if (geomix_local$fix_lateral) {
        noUpdate <- paste0("lL[", 1:geomix_local$constants$K, "]")
      } else if (!is.null(geomix_local$data_list$LFlag)) {
        noUpdate <- paste0("lL[", which(geomix_local$data_list$LFlag == 0), "]")
      }
    }

    for (tar in Z2scalarTargetsExpanded) {
      if (tar %in% noUpdate) {
        conf$replaceSampler(
          target  = tar,
          type    = DummySampler,
          control = NULL
        )
      } else {
        conf$replaceSampler(
          target  = tar,
          type    = HMCSample1D_kZ2,
          control = geomix_local$controlHMC
        )
      }
    }

    conf$addMonitors(c("Y1", "lL", "lD"), print = FALSE)
    if (inputs$LGFM) {
      conf$addMonitors(c("logProb_Z2"), print = FALSE)
    } else {
      conf$addMonitors(c("logProb_Z1", "logProb_Z2"), print = FALSE)
    }
    mcmc  <- nimble::buildMCMC(conf)
    Cmcmc <- suppressMessages(nimble::compileNimble(mcmc, project = model))
    message('Done.')
    message('Running samples...')
    bseq <- seq_len(inputs$controlMCMC$nbatches) + batch_start - 1
    retained_batches <- if (inputs$controlMCMC$retain_draws) vector("list", length(bseq)) else NULL

    cat("Chain", chain, ": batch sizes =", paste(inputs$batch_sizes, collapse = ", "), "\n")

    for (i in seq_along(bseq)) {
      b <- bseq[i]
      Cmcmc$run(inputs$batch_sizes[i], thin = inputs$controlMCMC$thin)

      batch <- as.matrix(Cmcmc$mvSamples)

      if (inputs$controlMCMC$save_batches) {
        saveRDS(batch, file = file.path(dir, paste0("batch_", b, ".rds")))
      }

      if (inputs$controlMCMC$retain_draws) {
        retained_batches[[i]] <- batch
      }

      Cmcmc$mvSamples$resize(0)

      cat("Chain", chain, ": completed batch", b, "(", inputs$batch_sizes[i], "iterations)\n")
    }
    message('Done.')
    retained_draws_mat <- if (inputs$controlMCMC$retain_draws) do.call(rbind, retained_batches) else NULL

    list(
      chain = chain,
      name = name,
      path = dir,
      batches_run = bseq,
      batch_sizes = inputs$batch_sizes,
      controlMCMC = inputs$controlMCMC,
      samples = retained_draws_mat
    )
  }

  chains <- seq_len(nchains)

  if (nchains == 1 | !isTRUE(run_parallel)) {
    out <- lapply(chains, run_one_chain, inputs = inputs)
  } else {
    future::plan(future::multisession, workers = mc.cores)
    on.exit(future::plan(future::sequential), add = TRUE)

    out <- future.apply::future_lapply(
      X = chains,
      FUN = run_one_chain,
      inputs = inputs,
      future.seed = TRUE,
      future.packages = "geomix"
    )
  }

  names(out) <- paste0(if (LGFM) "LGFM_" else "GeoMix_", chains)

  if (controlMCMC$retain_draws) {
    draw_list <- lapply(out, `[[`, "samples")
    names(draw_list) <- names(out)

    return(list(
      chains = out,
      samples = draw_list
    ))
  }

  out
}
