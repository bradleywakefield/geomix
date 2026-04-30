#' Set up data, constants, and initial values for a GeoMix model
#'
#' Prepares the full collection of objects required to fit a GeoMix model,
#' including model constants, Potts prior constants, processed data inputs,
#' initial values, lattice information, grouped Vecchia structures, and MCMC
#' control settings.
#'
#' This function constructs the lattice representation of the domain, builds the
#' design matrix from a model formula, organises observed and missing responses,
#' prepares distance matrices for the Gaussian process component, sets up the
#' Potts neighbourhood structure, and optionally creates grouped Vecchia
#' approximations if these are not supplied directly.
#'
#' @param data A data frame containing the observed data and coordinate/index
#'   variables (see variables for required columns) required by the model.
#' @param K Integer giving the number of latent classes.
#' @param dims Integer vector or list giving the lattice dimensions used to
#'   define neighbourhood structure.
#' @param variables Named list specifying the names of variables in
#'   `data`. Defaults to list names. Supported entries include:
#'   \describe{
#'     \item{loc}{Location/profile identifier. }
#'     \item{xID, yID, dID}{Integer lattice index variables for x, y, and depth.}
#'     \item{x, y, depth}{Observed coordinate variables.}
#'     \item{Z1}{Observed categorical label variable.}
#'     \item{Z2}{Observed continuous response variable.}
#'     \item{groups}{Optional grouping variable for Vecchia ordering/groups.}
#'   }
#' @param aformula A model formula used to construct the design matrix for the
#'   mean structure of the continuous response. Defaults to `~1`.
#' @param beta Numeric Potts interaction parameter. Defaults to `1.29`.
#' @param diagonals Logical or numeric indicator passed to the neighbourhood
#'   construction to determine whether diagonal neighbours are included.
#' @param kappa Numeric mixing parameter controlling the contribution of the
#'   structured confusion component.
#' @param w_tol Numeric tolerance used in the confusion matrix numeric approximation.
#' @param vecchia Optional precomputed Vecchia object. If `NULL`, this function
#'   constructs one using `setupVecchiaGeoMix()`.
#' @param hyperparams Optional named list of hyperparameters. Missing values are
#'   filled using defaults. Supported entries and defaults are:
#'   \describe{
#'     \item{a_tau}{Shape parameter for the inverse-gamma prior on `tau2`.
#'       Default: `4`.}
#'     \item{b_tau}{Scale parameter for the inverse-gamma prior on `tau2`.
#'       Default: `2`.}
#'     \item{a_sigma}{Shape parameter for the inverse-gamma prior on each
#'       class-specific `sigma2[k]`. Default: `3`.}
#'     \item{b_sigma}{Scale parameter for the inverse-gamma prior on each
#'       class-specific `sigma2[k]`. Default: `4`.}
#'     \item{a_L}{Shape parameter for the inverse-gamma prior on `sigma2_L`.
#'       Default: `3`.}
#'     \item{b_L}{Scale parameter for the inverse-gamma prior on `sigma2_L`.
#'       Default: `5`.}
#'     \item{a_D}{Shape parameter for the inverse-gamma prior on `sigma2_D`.
#'       Default: `3`.}
#'     \item{b_D}{Scale parameter for the inverse-gamma prior on `sigma2_D`.
#'       Default: `8`.}
#'     \item{sigma_h}{Scale parameter controlling the prior for `h`.
#'       Default: `1`.}
#'     \item{gamma0}{Prior vector for each row of `gammaMat`.
#'       Default: `rep(1, K)`.}
#'     \item{m_alpha}{Prior mean vector for regression coefficients `alpha[k, ]`.
#'       Default: `rep(0, p)`.}
#'     \item{Q_alpha}{Prior precision matrix for regression coefficients
#'       `alpha[k, ]`. Default: `0.01 * diag(p)`.}
#'   }
#' @param fix_lateral Logical indicating whether lateral covariance parameters
#'   should be treated as fixed. If `TRUE`, entries of `inits$lL` equal to `0`
#'   are used to disable lateral length scales for the corresponding classes.
#' @param inits Optional named list of initial values for model parameters.
#'   Missing values are filled using defaults. Supported entries and defaults are:
#'   \describe{
#'     \item{Y1}{Initial latent class allocations. Default:
#'       random draw `sample(1:K, N1, replace = TRUE)`.}
#'     \item{sigma2}{Initial class-specific process variances.
#'       Default: `rep(5, K)`.}
#'     \item{tau2}{Initial measurement noise variance.
#'       Default: `2`.}
#'     \item{Y1LogProb}{Initial stored log-probability value for `Y1`.
#'       Default: `1`.}
#'     \item{gammaMat}{Initial class confusion / transition matrix.
#'       Default: `matrix(1 / K, nrow = K, ncol = K)`.}
#'     \item{sigma2_L}{Initial global variance hyperparameter for lateral
#'       range effects. Default: `4`.}
#'     \item{sigma2_D}{Initial global variance hyperparameter for depth
#'       range effects. Default: `4`.}
#'     \item{lL}{Initial class-specific lateral range parameters.
#'       Default: `rep(1, K)`.}
#'     \item{lD}{Initial class-specific depth range parameters.
#'       Default: `rep(1, K)`.}
#'     \item{alpha}{Initial regression coefficient matrix with one row per
#'       class. Default: `matrix(0, nrow = K, ncol = p)`.}
#'     \item{h}{Initial boundary bandwidth / smoothing parameter.
#'       Default: `1`.}
#'   }
#' @param penalty Optional `K x K` penalty matrix used in the Potts/confusion
#'   setup. Defaults to `1 - diag(K)`.
#' @param weights Optional matrix of neighbour weights. If `NULL`, equal weights
#'   are used.
#' @param m Integer controlling the Vecchia conditioning size and also passed
#'   into the default HMC control list. Defaults to 100.
#' @param HMC_control Optional named list of HMC tuning parameters. Missing
#'   entries are filled using defaults. Supported entries and defaults are:
#'   \describe{
#'     \item{nLeap}{Number of leapfrog steps per HMC proposal.
#'       Default: `10`.}
#'     \item{eps}{Leapfrog step size.
#'       Default: `0.015`.}
#'     \item{h}{Finite-difference step size used for numerical gradient
#'       approximation (if applicable). Default: `1e-3`.}
#'     \item{mass}{Mass matrix scalar (identity-scaled momentum variance).
#'       Default: `1.0`.}
#'     \item{m}{Neighbourhood / approximation size passed through from the
#'       main function argument `m`. Default: value of `m`.}
#'   }
#' @param mcmc_control Optional named list of overall MCMC control parameters.
#'   Missing entries are filled using defaults. Supported entries and defaults are:
#'   \describe{
#'     \item{niter}{Total number of MCMC iterations.
#'       Default: `2750`.}
#'     \item{nburnin}{Number of initial burn-in iterations discarded.
#'       Default: `250`.}
#'     \item{nthin}{Thinning interval for retained samples.
#'       Default: `10`.}
#'   }
#' @details
#' The function performs the following main tasks:
#' \enumerate{
#'   \item Defines default variable names and hyperparameters.
#'   \item Builds the design matrix from `aformula`.
#'   \item Creates a lattice and aligns the supplied data to that lattice.
#'   \item Identifies observed `Z2` indices and constructs distance matrices.
#'   \item Builds or accepts a grouped Vecchia approximation.
#'   \item Constructs Potts neighbourhood information and associated constants.
#'   \item Creates data lists, constants, and initial values for downstream
#'   model fitting.
#' }
#'
#' The returned object is intended to be passed to subsequent model-building,
#' MCMC, or sampler configuration functions in the GeoMix workflow.
#'
#' @return A named list with components:
#' \describe{
#'   \item{constants}{Model constants including dimensions and hyperparameters.}
#'   \item{constantsPotts}{Constants specific to the Potts prior and Gibbs
#'   updates.}
#'   \item{data_list}{Processed data required for model fitting.}
#'   \item{inits}{Initial values for the model parameters.}
#'   \item{df}{Processed version of the input data.}
#'   \item{X}{Design matrix built from `aformula`.}
#'   \item{lattice}{Full lattice data frame.}
#'   \item{lattice_coords}{Ranked coordinate/index subset used in ordering and Vecchia
#'   construction.}
#'   \item{controlGibbs}{Combined constants used for Gibbs updates.}
#'   \item{controlHMC}{HMC control settings.}
#'   \item{controlMCMC}{MCMC control settings.}
#'   \item{grouping}{Grouping variable used for Vecchia construction.}
#'   \item{vecchia}{Vecchia setup object.}
#'   \item{fix_lateral}{Logical flag passed through from the input.}
#' }
#'
#' @examples
#' \dontrun{
#' data(offshore)
#'
#' # Basic setup with depth as a linear mean covariate
#' setup <- setupGeoMixModel(
#'   data      = offshore,
#'   K         = 3,
#'   dims      = c(20, 20, 20),
#'   variables = list(
#'     loc   = "locID",
#'     xID   = "x",
#'     yID   = "y",
#'     dID   = "d",
#'     x     = "x",
#'     y     = "y",
#'     depth = "d",
#'     Z1    = "Z1",
#'     Z2    = "Z2"
#'   ),
#'   aformula = ~ d,
#'   m        = 10
#' )
#'
#' # Inspect key constants
#' setup$constants$N1   # total lattice size (8000)
#' setup$constants$N2   # observed Z2 locations (1600)
#' setup$constants$K    # number of classes (3)
#'
#' # Override hyperparameters and initial values
#' setup2 <- setupGeoMixModel(
#'   data        = offshore,
#'   K           = 3,
#'   dims        = c(20, 20, 20),
#'   variables   = list(
#'     loc = "locID", xID = "x", yID = "y", dID = "d",
#'     x = "x", y = "y", depth = "d", Z1 = "Z1", Z2 = "Z2"
#'   ),
#'   aformula    = ~ d,
#'   m           = 10,
#'   hyperparams = list(a_tau = 3, b_tau = 1),
#'   inits       = list(tau2 = 1, sigma2 = c(4, 4, 7)),
#'   mcmc_control = list(niter = 5000, thin = 5)
#' )
#' }
#'
#' @seealso [setupVecchiaGeoMix()], [run_chains()]
#'
#' @export
setupGeoMixModel <- function(data, K, dims, variables = NULL,
                             aformula = ~1,
                             beta = 1.29, diagonals = 0,
                             kappa = 0.90, w_tol = 0.001,
                             vecchia = NULL, hyperparams = NULL,
                             fix_lateral = F, inits = NULL,
                             penalty = NULL, weights = NULL, m =100,
                             HMC_control = NULL, mcmc_control = NULL) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  # === Defaults ===
  default_variables <- list(
    loc = "loc", xID = "xID", yID = "yID", dID = "dID",
    x = "x", y = "y", depth = "d",
    Z1 = "Z1", Z2 = "Z2",
    groups = NULL
  )
  variables <- modifyList(default_variables, variables %||% list())

  #if(identical(aformula,~1)) aformula <- ~ zeros
  X <- model.matrix(aformula,data)
  p <- ncol(X)

  default_hyperparams <- list(
    a_tau = 4, b_tau = 2,
    a_sigma = 3, b_sigma = 4,
    a_L = 3, b_L = 5,
    a_D = 3, b_D = 8,
    sigma_h = 1,
    gamma0 = rep(1, K),
    m_alpha = rep(0,p),  Q_alpha = 0.01*diag(p)
  )

  hyperparams <- modifyList(default_hyperparams, hyperparams %||% list())

  penalty <- penalty %||% (1 - diag(K))

  default_mcmc <- list(
    niter = 2500,
    thin = 10,
    nbatches = 1,
    save_batches = FALSE,
    retain_draws = TRUE
  )
  mcmc_control <- modifyList(default_mcmc, mcmc_control %||% list())

  default_HMC <- list(nLeap = 10, eps = 0.015, h = 1e-3, mass = 1.0, m = m)
  HMC_control <- modifyList(default_HMC, HMC_control %||% list())

  lattice <- mutate(expand.grid(lapply(dims,seq,from = 1)),latticeID = 1:n(),.before = everything())
  colnames(lattice) <- c("latticeID",variables$dID,variables$xID,variables$yID)
  lattice <- suppressMessages(left_join(lattice,select(data,c(variables$dID,variables$xID,variables$yID,variables$loc))))

  data <- data %>%
    left_join(select(lattice,variables$dID,variables$xID,variables$yID,latticeID),
              by = c(variables$dID,variables$xID,variables$yID)) %>%
    relocate(latticeID) %>% ungroup() %>% arrange(latticeID) %>%
    mutate(rowID = 1:n(),.before = everything())

  N1 <- nrow(data)
  Z2ind <- match(data$rowID[which(!is.na(data[[variables$Z2]]))],1:N1)
  N2 <- length(Z2ind)

  depth_vec <- sort(unique(data[[variables$depth]]))

  rank_dist <- select(data,variables$dID,variables$xID,variables$yID)

  grouping <- variables$groups
  if(is.null(vecchia)){
    if(!is.null(variables$groups)){
      groups_full <- as.numeric(factor(data[[variables$groups]]))
      groups <- groups_full[Z2ind]
      Z2ind <- Z2ind[order(groups)]
      groups <- groups[order(groups)]
    }else{
      Z2ind <- Z2ind[GpGp::order_maxmin(rank_dist[Z2ind,])]
      groups <- 1:N2
      notZ2ind <- (1:N1)[-Z2ind]
      notZ2ind <- notZ2ind[GpGp::order_maxmin(rank_dist[notZ2ind,])]
      groups_full <-c()
      groups_full[Z2ind] <- 1:N2
      groups_full[notZ2ind] <- (N2+1):N1
      variables$groups <- "VecchiaGroups"
      data$VecchiaGroups <- groups_full
      grouping <- variables$groups
    }
    vecchia <- setupVecchiaGeoMix(rank_dist[Z2ind,], m = m, groups = groups)
  }else{
    HMC_control$m  <- vecchia$constants$m
  }

  #Confusion Setup
  loc_df <- data %>%
    select(rowID,loc=variables$loc,x=variables$x, y=variables$y) %>%
    group_by(loc,x,y) %>%
    summarise(minD = min(rowID),
              maxD = max(rowID),.groups="drop") %>%
    arrange(loc)

  locStack <- as.matrix(select(loc_df,minD,maxD))
  udistD <- as.matrix(dist(depth_vec))
  udistL <- as.matrix(dist(loc_df[,c("x","y")]))
  dID <- match(data[[variables$depth]],depth_vec)
  locID <- match(data[[variables$loc]],loc_df$loc)

  #Potts Prior
  neighbours0 <- getNeighbours(dims = dims, diagonals = diagonals)
  neigh_df <- neighbours0[data$latticeID,]
  lattice_to_sample <- c(0,replace_na(match(lattice$latticeID, data$latticeID),0))
  neigh_df <- apply(neigh_df,2,function(i)lattice_to_sample[i+1])
  neighbours <- t(apply(neigh_df,1,function(x) c(x[x>0], x[x==0])))
  neighbourNum <- apply(neighbours>0,1,sum)
  weights <- weights %||% matrix(1, nrow = nrow(data), ncol = ncol(neighbours))

  constantsPotts <- list(
    w_tol = w_tol,
    delta_depth = min(diff(depth_vec)),
    neighbourID = neighbours,
    neighbourNum = neighbourNum,
    weights = weights,
    penalty = penalty,
    beta = beta,
    Z2_flag = (!is.na(data[[variables$Z2]])),
    groups = vecchia$groups,
    groupDeps =  vecchia$groupDeps,
    groupDepL =   vecchia$groupDepL
  )

  constants <- list(
    N2 = N2,
    N1 = N1,
    K = K,
    p = p,
    D = nrow(udistD),
    L = nrow(udistL),
    Z2_ind = Z2ind,
    kappa = kappa
  )
  constants <- c(constants,vecchia$constants)

  constantsM <- c(constants, hyperparams)

  initialY1 <- sample(1:K,N1,replace = T)

  default_inits <- list(
    Y1 = initialY1,
    sigma2 = rep(5,K),
    tau2 = 2,
    Y1LogProb = 1,
    gammaMat = matrix(1/K,nrow=K,ncol=K),
    sigma2_L = 4, sigma2_D = 4,
    lL = rep(1,K), lD = rep(1,K),
    alpha = matrix(0,ncol=p,nrow=K),
    h = 1
  )
  inits <- modifyList(default_inits, inits %||% list())

  if(fix_lateral){
    FlagL <- as.integer(inits$lL != 0)
  }else{
    FlagL <- rep(1L,K)
  }
  
  data_list <- list(
    Z2 = data[[variables$Z2]][Z2ind],
    Z1 = data[[variables$Z1]],
    dID = dID,
    locID = locID,
    distL = udistL^2,
    distD = udistD^2,
    LFlag = FlagL,
    locStack = locStack,
    X=array(X[Z2ind,],dim = c(N2,p)),
    groupLookup = matrix(as.integer(vecchia$groupLookup), nrow = nrow(vecchia$groupLookup)),
    groupNum = as.integer(vecchia$groupNum),
    groupNeighbours = matrix(as.integer(vecchia$groupNeighbours), nrow = nrow(vecchia$groupNeighbours))
  )

  return(list(
    constants = constantsM,
    constantsPotts = constantsPotts,
    data_list = data_list,
    inits = inits,
    df = data,
    X = X,
    lattice = lattice,
    lattice_coords = rank_dist,
    controlGibbs = c(constants,constantsPotts),
    controlHMC = HMC_control,
    controlMCMC = mcmc_control,
    grouping = grouping,
    vecchia = vecchia,
    fix_lateral = fix_lateral
  ))
}
