#' Run MCMC diagnostics for GeoMix model outputs
#'
#' Computes and summarises convergence diagnostics, effective sample sizes,
#' Monte Carlo standard errors, and representative plots for one or more
#' MCMC chains produced by GeoMix model fitting routines. The function is
#' designed for posterior samples stored in `params` objects containing a
#' `$samples` element.
#'
#' It supports multiple chains, combines draws into `posterior` formats,
#' produces diagnostic tables, rank plots, trace plots, autocorrelation plots,
#' log-posterior traces, and latent class count summaries.
#'
#' By default, high-dimensional latent fields such as `Y1` are excluded from
#' standard convergence summaries, but selected summaries of sampled class
#' counts are provided.
#'
#' @param params_list A single fitted GeoMix parameter object with a `$samples`
#'   component, or a list of such objects representing multiple MCMC chains.
#'
#' @param name Character string giving the model name used in plot titles.
#'   Defaults to `"GeoMix"`.
#'
#' @param exclude Optional character vector of parameter names to exclude from
#'   diagnostics after chain matrices are constructed. Defaults to `NULL`.
#'
#' @param plot_variables Character vector specifying parameter groups to include
#'   in representative diagnostic plots. Defaults to
#'   `c("alpha0","alpha1","sigma2","lL","lD","tau2","sigma2_L","sigma2_D","gamma","h","logProb")`.
#'
#' @param plot_classes Integer vector of class indices used when selecting
#'   class-specific parameters (such as `sigma2[k]`, `lL[k]`, `lD[k]`) for
#'   plotting. Defaults to `c(1,2)`.
#'
#' @param plot_gamma Integer matrix with two columns specifying off-diagonal
#'   gamma parameters to include in plots. Each row corresponds to one
#'   `gammaMat[i,j]` pair. Defaults to
#'   `rbind(c(1,2), c(1,3), c(2,1), c(2,3))`.
#'
#' @param Y1index Optional integer vector giving selected lattice indices of
#'   latent labels `Y1` to use when producing class-count trace plots.
#'   Defaults to all columns of sampled `Y1`.
#'
#' @return A named list containing MCMC draws, diagnostic tables,
#' summaries, and plots.
#'
#' \describe{
#'   \item{draws}{Posterior draws in multiple formats.}
#'   \describe{
#'     \item{draws_array}{A `posterior::draws_array` object with dimensions
#'     iterations × chains × parameters.}
#'     \item{combined_matrix}{A stacked
#'     `posterior::draws_matrix` containing all chains combined by row.}
#'     \item{chain_matrices}{A list of numeric matrices, one per chain,
#'     with rows corresponding to iterations and columns to parameters.}
#'   }
#'
#'   \item{tables}{Diagnostic summary tables.}
#'   \describe{
#'     \item{diagnostics}{Parameter-level diagnostics including `rhat`,
#'     bulk/tail ESS, posterior SD, MCSE, quantiles, and parameter class.}
#'     \item{class_diagnostics}{Grouped summaries by parameter class
#'     (e.g. `sigma2`, `lL`, `gamma`).}
#'     \item{chain_diagnostics}{Per-chain summaries for each parameter,
#'     including ESS, MCSE, and posterior quantiles.}
#'     \item{overall_diagnostics}{One-row summary of global convergence
#'     statistics across all monitored parameters.}
#'     \item{worst_rhat}{Top 10 parameters ranked by largest `rhat`.}
#'     \item{worst_ess}{Top 10 parameters ranked by lowest bulk ESS.}
#'   }
#'
#'   \item{summaries}{General metadata and printed summary text.}
#'   \describe{
#'     \item{key_stats_text}{Character string containing the printed
#'     diagnostic report.}
#'     \item{n_chains}{Number of MCMC chains analysed.}
#'   }
#'
#'   \item{plots}{Diagnostic `ggplot2` / `bayesplot` objects.}
#'   \describe{
#'     \item{trace_acf}{List containing representative trace plots,
#'     autocorrelation plots, and selected parameter names.}
#'     \item{rank}{Rank histogram overlay plot for selected parameters.}
#'     \item{logProb}{Trace plots of log-posterior components by chain.}
#'     \item{Y1}{Trace plots of sampled latent class counts over iterations.}
#'   }
#' }
#'
#' @details
#' Diagnostics are based on functions from the \pkg{posterior} package,
#' including:
#' \itemize{
#'   \item rank-normalised split-\eqn{\hat{R}}
#'   \item bulk effective sample size
#'   \item tail effective sample size
#'   \item Monte Carlo standard error
#' }
#'
#' For non-LGFM models, off-diagonal gamma parameters and boundary parameter
#' `h` are automatically included where available.
#'
#' A printed summary is emitted via `message()` when the function runs.
#'
#' @examples
#' \dontrun{
#' # Single chain
#' out <- run_mcmc_diagnostics(params)
#'
#' # Multiple chains
#' out <- run_mcmc_diagnostics(
#'   params_list = list(chain1, chain2, chain3),
#'   name = "GeoMix",
#'   plot_classes = c(1, 4)
#' )
#'
#' # View worst mixing parameters
#' out$tables$worst_rhat
#'
#' # Print trace plots
#' out$plots$trace_acf$trace
#' }
#'
#' @seealso [posterior::rhat()], [posterior::ess_bulk()],
#'   [bayesplot::mcmc_trace()]
#'
#' @export
#'
run_mcmc_diagnostics <- function(params_list, name = "GeoMix",
                                 exclude = NULL,
                                 plot_variables = c("alpha0","alpha1","sigma2","lL","lD","tau2","sigma2_L","sigma2_D","gamma","h","logProb"),
                                 plot_classes = c(1,2),
                                 plot_gamma = rbind(c(1,2),c(1,3),c(2,1),c(2,3)),
                                 Y1index = NULL) {

  ## ============================================================
  ## GeoMix MCMC diagnostics (supplementary)
  ## - Handles multiple chains stored as a list of params objects
  ## - Key checks for continuous / low-dim parameters
  ## - Ignores Y1 (too large) by design
  ## - Returns tables, plots, summaries, and combined multi-chain draws
  ## ============================================================

  default_theme <-   theme(legend.position = "bottom",
                           legend.key.size = unit(0.5,"cm"),
                           axis.text = element_text(size=8),
                           axis.title = element_text(size=10),
                           legend.title = element_text(size=10),
                           plot.subtitle = element_text(size=8),
                           legend.title.align = 0.5,
                           legend.spacing.x  = unit(0.02, "cm"),
                           legend.margin = margin(0,0,0,0),
                           plot.margin = margin(0, 0, 0, 0),
                           strip.background = element_blank(),
                           strip.placement = "inside",   # optional
                           strip.switch.pad.wrap = unit(0.0, "cm"),
                           strip.text = element_text(size = 8,margin = margin(0,0,0,0)),
                           legend.text = element_text(size=9))

  ## ------------------------------------------------------------
  ## Helpers
  ## ------------------------------------------------------------

  ## Convert vector / matrix / array into plain iteration x parameter matrix
  to_iter_matrix <- function(x, prefix) {
    if (is.null(dim(x))) {
      dm <- matrix(x, ncol = 1)
      colnames(dm) <- prefix
      return(dm)
    }

    if (length(dim(x)) == 2) {
      ## matrix: [iter, j]
      dm <- x
      colnames(dm) <- paste0(prefix, "[", seq_len(ncol(dm)), "]")
      return(dm)
    }

    if (length(dim(x)) == 3) {
      ## array: [iter, j, q]
      it <- dim(x)[1]
      j  <- dim(x)[2]
      q  <- dim(x)[3]
      dm <- matrix(NA_real_, nrow = it, ncol = j * q)
      nm <- character(j * q)
      idx <- 1L
      for (jj in seq_len(j)) {
        for (qq in seq_len(q)) {
          dm[, idx] <- x[, jj, qq]
          nm[idx] <- paste0(prefix, "[", jj, ",", qq, "]")
          idx <- idx + 1L
        }
      }
      colnames(dm) <- nm
      return(dm)
    }

    stop("Unsupported dims for ", prefix)
  }

  ## Convert one chain's samples into a named iteration x parameter matrix
  build_chain_matrix <- function(params, name) {
    s <- params$samples

    mats <- list(
      tau2   = to_iter_matrix(s$tau2,   "tau2"),
      sigma2_L   = to_iter_matrix(s$sigma2_L,   "sigma2_L"),
      sigma2_D   = to_iter_matrix(s$sigma2_D,   "sigma2_D"),
      sigma2 = to_iter_matrix(s$sigma2, "sigma2"),
      lL     = to_iter_matrix(s$lL,     "lL"),
      lD     = to_iter_matrix(s$lD,     "lD"),
      alpha0 = to_iter_matrix(s$a0,     "alpha0"),
      alpha1 = to_iter_matrix(s$a1,     "alpha1")
    )

    if (name != "LGFM") {
      gamma_mat <- as.matrix(s$gamma)
      h_mat     <- to_iter_matrix(s$h, "h")

      ## Drop diagonal gamma terms if present as gammaMat[i, i]
      if (!is.null(colnames(gamma_mat))) {
        drop_cols <- grepl("^gammaMat\\[(\\d+), \\1\\]$", colnames(gamma_mat))
        if (any(drop_cols)) {
          gamma_mat <- gamma_mat[, !drop_cols, drop = FALSE]
        }
      }

      mats$gamma <- gamma_mat
      mats$h     <- h_mat
    }

    out <- do.call(cbind, mats)
    out
  }

  ## Parameter classes
  param_class <- function(x) {
    cls <- rep(NA_character_, length(x))

    cls[grepl("^(alpha0|alpha1|a0|a1)\\[", x)] <- "alpha"
    cls[grepl("^sigma2\\[", x)] <- "sigma2"
    cls[grepl("^sigma2\\_L", x)] <- "Other variances"
    cls[grepl("^sigma2\\_D", x)] <- "Other variances"
    cls[grepl("^lL\\[", x)]     <- "lL"
    cls[grepl("^lD\\[", x)]     <- "lD"
    cls[grepl("^gamma(Mat)?\\[|^gamma\\[", x)] <- "gamma"
    cls[grepl("^tau2(\\[|$)", x)] <- "Other variances"
    cls[grepl("^h(\\[|$)", x)] <- "h"

    cls
  }

  ## Build posterior draws array: iterations x chains x variables
  build_draws_array <- function(chain_mats) {
    pnames <- colnames(chain_mats[[1]])
    n_iter <- nrow(chain_mats[[1]])
    n_chain <- length(chain_mats)
    n_var <- length(pnames)

    arr <- array(NA_real_, dim = c(n_iter, n_chain, n_var),
                 dimnames = list(
                   iteration = NULL,
                   chain = paste0("chain", seq_len(n_chain)),
                   variable = pnames
                 ))

    for (ch in seq_len(n_chain)) {
      if (!identical(colnames(chain_mats[[ch]]), pnames)) {
        stop("Column names differ across chains.")
      }
      if (nrow(chain_mats[[ch]]) != n_iter) {
        stop("Number of iterations differs across chains.")
      }
      arr[, ch, ] <- as.matrix(chain_mats[[ch]])
    }

    posterior::as_draws_array(arr)
  }

  ## Combined per-parameter diagnostics using all chains
  diag_table <- function(draws_arr) {
    var_names <- posterior::variables(draws_arr)

    ess_b <- apply(draws_arr,3,posterior::ess_bulk)
    ess_t <- apply(draws_arr,3,posterior::ess_tail)
    rhat  <- apply(draws_arr,3,posterior::rhat)
    mcse_m <- apply(draws_arr,3,posterior::mcse_mean)
    mcse_s <- apply(draws_arr,3,posterior::mcse_sd)
    sds   <- posterior::summarise_draws(draws_arr, sd = sd)$sd

    qs <- posterior::summarise_draws(
      draws_arr,
      q025 = ~quantile(.x, probs = 0.025),
      q50  = ~quantile(.x, probs = 0.50),
      q975 = ~quantile(.x, probs = 0.975)
    )

    tibble::tibble(
      param = var_names,
      ess_bulk = as.numeric(ess_b),
      ess_tail = as.numeric(ess_t),
      rhat = as.numeric(rhat),
      sd = as.numeric(sds),
      mcse_mean = as.numeric(mcse_m),
      mcse_sd = as.numeric(mcse_s),
      q025 = as.numeric(qs$`2.5%`),
      q50  = as.numeric(qs$`50%`),
      q975 = as.numeric(qs$`97.5%`)
    ) %>%
      mutate(
        mcse_over_sd = mcse_mean / sd,
        class = param_class(param)
      ) %>%
      arrange(desc(rhat), desc(mcse_over_sd))
  }

  ## Per-chain summaries for each parameter
  chain_diag_table <- function(draws_arr) {
    n_chain <- posterior::nchains(draws_arr)
    var_names <- posterior::variables(draws_arr)

    out <- vector("list", n_chain)

    for (ch in seq_len(n_chain)) {
      darr_ch <- posterior::subset_draws(draws_arr, chain = ch)
      ess_b <- apply(darr_ch,3,posterior::ess_bulk)
      ess_t <- apply(darr_ch,3,posterior::ess_tail)
      mcse_m <- apply(darr_ch,3,posterior::mcse_mean)
      mcse_s <- apply(darr_ch,3,posterior::mcse_sd)
      sds <- posterior::summarise_draws(darr_ch, sd = sd)$sd

      qs <- posterior::summarise_draws(
        darr_ch,
        q025 = ~quantile(.x, probs = 0.025),
        q50  = ~quantile(.x, probs = 0.50),
        q975 = ~quantile(.x, probs = 0.975)
      )

      out[[ch]] <- tibble::tibble(
        chain = ch,
        param = var_names,
        ess_bulk = as.numeric(ess_b),
        ess_tail = as.numeric(ess_t),
        sd = as.numeric(sds),
        mcse_mean = as.numeric(mcse_m),
        q025 = as.numeric(qs$`2.5%`),
        q50  = as.numeric(qs$`50%`),
        q975 = as.numeric(qs$`97.5%`),
        class = param_class(var_names)
      ) %>%
        mutate(mcse_over_sd = mcse_mean / sd)
    }

    dplyr::bind_rows(out)
  }

  ## Class-level summaries
  class_summary_table <- function(diag_all) {
    diag_all %>%
      group_by(class) %>%
      summarise(
        count = n(),
        med_rhat = median(rhat, na.rm = TRUE),
        max_rhat = max(rhat, na.rm = TRUE),
        med_ess = median(ess_bulk, na.rm = TRUE),
        min_ess = min(ess_bulk, na.rm = TRUE),
        med_tail = median(ess_tail, na.rm = TRUE),
        min_tail = min(ess_tail, na.rm = TRUE),
        med_mcse_over_sd = median(mcse_over_sd, na.rm = TRUE),
        max_mcse_over_sd = max(mcse_over_sd, na.rm = TRUE),
        .groups = "drop"
      )
  }

  ## Overall multi-chain summary
  overall_summary <- function(diag_all) {
    tibble::tibble(
      n_param = nrow(diag_all),
      n_rhat_gt_1_01 = sum(diag_all$rhat > 1.01, na.rm = TRUE),
      n_rhat_gt_1_05 = sum(diag_all$rhat > 1.05, na.rm = TRUE),
      max_rhat = max(diag_all$rhat, na.rm = TRUE),
      median_rhat = median(diag_all$rhat, na.rm = TRUE),
      min_bulk_ess = min(diag_all$ess_bulk, na.rm = TRUE),
      median_bulk_ess = median(diag_all$ess_bulk, na.rm = TRUE),
      min_tail_ess = min(diag_all$ess_tail, na.rm = TRUE),
      median_tail_ess = median(diag_all$ess_tail, na.rm = TRUE),
      max_mcse_over_sd = max(diag_all$mcse_over_sd, na.rm = TRUE),
      median_mcse_over_sd = median(diag_all$mcse_over_sd, na.rm = TRUE)
    )
  }

  ## Plot representative trace and ACF across multiple chains
  plot_trace_acf <- function(draws_arr, pars, title_prefix = "") {
    keep <- intersect(pars, posterior::variables(draws_arr))
    if (length(keep) == 0) return(NULL)
    dsel <- posterior::subset_draws(draws_arr, variable = keep)

    p1 <- bayesplot::mcmc_trace(dsel) +
      ggtitle(paste0(title_prefix, "Trace plots"))+
      geom_line(alpha = 0.2)

    p2 <- bayesplot::mcmc_acf(dsel, lags = 50) +
      ggtitle(paste0(title_prefix, "Autocorrelation"))

    list(
      trace = p1,
      acf = p2,
      selected_pars = keep
    )
  }

  ## Rank histograms are useful for checking chain mixing
  plot_rank_overlay <- function(draws_arr, pars, title_prefix = "") {
    keep <- intersect(pars, posterior::variables(draws_arr))
    if (length(keep) == 0) return(NULL)

    dsel <- posterior::subset_draws(draws_arr, variable = keep)
    bayesplot::mcmc_rank_overlay(dsel) +
      ggtitle(paste0(title_prefix, "Rank plots"))
  }

  ## ------------------------------------------------------------
  ## 1) Standardise input: allow either one params object or list
  ## ------------------------------------------------------------
  if (!is.list(params_list) || is.null(params_list$samples)) {
    ## Could still be a list of params objects
    is_params_list <- is.list(params_list) &&
      all(vapply(params_list, function(x) is.list(x) && !is.null(x$samples), logical(1)))

    if (!is_params_list) {
      stop("params_list must be either a single params object with $samples or a list of such objects.")
    }
  }

  if (!is.null(params_list$samples)) {
    params_list <- list(params_list)
  }

  n_chains <- length(params_list)

  ## ------------------------------------------------------------
  ## 2) Build one matrix per chain and combine
  ## ------------------------------------------------------------
  chain_mats <- lapply(params_list, build_chain_matrix, name = name)
  if(!is.null(exclude)) chain_mats <- map(chain_mats,~.x[,setdiff(colnames(.x),exclude)])
  draws_arr <- build_draws_array(chain_mats)

  ## Optional combined matrix if you want stacked iterations across chains
  draws_combined_mat <- do.call(rbind, chain_mats)
  draws_combined_mat <- posterior::as_draws_matrix(draws_combined_mat)

  ## ------------------------------------------------------------
  ## 3) Multi-chain numeric checks
  ## ------------------------------------------------------------
  diag_all <- diag_table(draws_arr)
  class_diags <- class_summary_table(diag_all)
  overall_diags <- overall_summary(diag_all)
  chain_diags <- chain_diag_table(draws_arr)

  worst_rhat <- diag_all %>%
    arrange(desc(rhat), desc(mcse_over_sd)) %>%
    select(param, class, rhat, ess_bulk, ess_tail, mcse_over_sd) %>%
    head(10)

  worst_ess <- diag_all %>%
    arrange(ess_bulk, ess_tail) %>%
    select(param, class, rhat, ess_bulk, ess_tail, mcse_over_sd) %>%
    head(10)

  key_stats_text <- paste(
    capture.output({
      cat("\nMulti-chain diagnostic summary\n")
      cat("----------------------------------------\n")
      cat("Number of chains: ", n_chains, "\n", sep = "")
      cat("Iterations per chain: ", posterior::ndraws(draws_arr) / n_chains, "\n", sep = "")
      cat("Number of monitored parameters: ", nrow(diag_all), "\n", sep = "")

      cat("\nOverall summary:\n")
      print(overall_diags)

      cat("\nClass-specific summary:\n")
      print(class_diags)

      cat("\nWorst 10 parameters by R-hat:\n")
      print(worst_rhat)

      cat("\nWorst 10 parameters by bulk ESS:\n")
      print(worst_ess)
    }),
    collapse = "\n"
  )

  ## ------------------------------------------------------------
  ## 4) Representative multi-chain plots
  ## ------------------------------------------------------------
  class_vars <- c("alpha0","alpha1","sigma2","lL","lD")
  other_vars <- c("tau2","sigma2_L","sigma2_D","h")
  class_vars <- plot_variables[(plot_variables %in% class_vars)]
  other_vars <- plot_variables[(plot_variables %in% other_vars)]

  rep_pars <- c(paste0(rep(paste0(class_vars,"["),each=length(plot_classes)),
                       rep(plot_classes,length(class_vars)),"]"),other_vars)
  if("gamma" %in% plot_variables){
  gamma_names <- posterior::variables(draws_arr)
  gamma_str <- paste0("^gammaMat\\[",plot_gamma[,1],", ",plot_gamma[,2],"\\]$")
  gamma_keep <- sapply(gamma_str,function(x) grep(x, gamma_names, value = TRUE))
  rep_pars <- unique(c(rep_pars, gamma_keep))
  }

  plots <- list(
    trace_acf = plot_trace_acf(draws_arr, rep_pars, title_prefix = paste0(name, ": ")),
    rank = plot_rank_overlay(draws_arr, rep_pars, title_prefix = paste0(name, ": "))
  )

  ## ------------------------------------------------------------
  ## 5) Additional plots
  ## ------------------------------------------------------------

  if(is.null(Y1index)) Y1index <- 1:ncol(params_list[[1]]$samples$Y1)

  logProb <- imap(params_list, ~ data.frame(
    iter      = seq_along(.x$samples$logProb),
    logProb   = .x$samples$logProb,
    logProbZ1 = .x$samples$logProbZ1,
    logProbZ2 = .x$samples$logProbZ2,
    chain     = .y
  )) %>% bind_rows()

  Y1_count <- lapply(seq_along(params_list), function(i) {
    Y1 <- params_list[[i]]$samples$Y1[,Y1index]
    K  <- max(Y1)

    counts <- t(vapply(
      seq_len(nrow(Y1)),
      function(r) tabulate(Y1[r, ], nbins = K),
      integer(K)
    ))

    as.data.frame(counts) |>
      dplyr::mutate(iter = seq_len(nrow(Y1)), chain = i)
  })
  Y1_count <- dplyr::bind_rows(Y1_count) %>%
    pivot_longer(starts_with("V"),names_prefix = "V",
                 names_to = "class", values_to = "count")

  logProbplot <- logProb %>%
    pivot_longer(starts_with("logProb"),names_prefix = "logProb") %>%
    mutate(name = if_else(name=="","All",name)) %>%
    ggplot()+
    geom_line(aes(x=iter,y=value, col = factor(chain)), alpha = 0.6, linewidth=0.2)+
    facet_wrap(vars(name),scales = "free_y")+
    theme_bw()+default_theme +
    theme(legend.position = "bottom")+
    scale_x_continuous(breaks = c(0,250,500),
                       labels = \(x) if_else(x==500,
                                             paste0(format(x*10, big.mark = ","),'  '),
                                             format(x*10, big.mark = ",") ))+
    labs(y="Count",x="Iteration", col = "Chain")

  Y1countplot <-
    ggplot(Y1_count)+
    geom_line(aes(x=iter,y=count, col = factor(chain)), alpha = 0.6, linewidth=0.2)+
    facet_wrap(vars(paste0("Y1 = ",class)),scales = "free_y")+
    theme_bw()+default_theme +
    labs(y="Count",x="Iteration", col = "Chain")

  plots <- c(plots, list(logProb = logProbplot, Y1 = Y1countplot))
  ## ------------------------------------------------------------
  ## 6) Return everything
  ## ------------------------------------------------------------
  out <- list(
    draws = list(
      draws_array = draws_arr,
      combined_matrix = draws_combined_mat,
      chain_matrices = chain_mats
    ),
    tables = list(
      diagnostics = diag_all,
      class_diagnostics = class_diags,
      chain_diagnostics = chain_diags,
      overall_diagnostics = overall_diags,
      worst_rhat = worst_rhat,
      worst_ess = worst_ess
    ),
    summaries = list(
      key_stats_text = key_stats_text,
      n_chains = n_chains
    ),
    plots = plots
  )
  message(out$summaries$key_stats_text)
  return(out)
}
