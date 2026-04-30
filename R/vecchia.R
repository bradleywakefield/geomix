#' Construct grouped Vecchia structures for GeoMix
#'
#' Builds the grouping, dependency, lookup, and neighbour structures used by the
#' grouped Vecchia approximation within GeoMix. The function converts an
#' observation-level grouping vector into group-level indexing objects and
#' determines which previous groups each group may condition on.
#'
#' Depending on the chosen dependency structure, the function creates either a
#' user-supplied dependency graph or one of several built-in sequential group
#' dependency schemes. It then selects up to `m` candidate neighbour
#' observations for each group from its allowable predecessor groups.
#'
#' @param locations Numeric matrix or data frame of observation coordinates,
#'   with one row per observation and columns representing spatial coordinates.
#'   These coordinates are used for nearest-neighbour searches.
#'
#' @param m Integer giving the maximum number of conditioning neighbours used
#'   for each group in the Vecchia approximation.
#'
#' @param groups Vector of length `N` giving the group membership for each
#'   observation. Values need not be consecutive; they are internally converted
#'   to sequential integer group labels.
#'
#' @param groupMatrix Optional user-supplied group dependency matrix. Each row
#'   corresponds to a group, and non-zero entries indicate predecessor groups
#'   that the row group may condition on. If supplied, this overrides
#'   `depStructure`.
#'
#' @param depStructure Character string specifying the default dependency
#'   structure when `groupMatrix = NULL`. Supported options are:
#'   \describe{
#'     \item{"minimal"}{Chooses the smallest number of previous groups whose
#'      combined size reaches at least `m` observations where possible.}
#'     \item{"previous"}{Each group depends only on the immediately previous
#'       group.}
#'     \item{"full"}{Each group depends on all previous groups.}
#'   }
#'   Default is `"minimal"`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Re-indexes supplied groups into consecutive integers.
#'   \item Computes group sizes and the total number of groups.
#'   \item Constructs a group-level dependency matrix.
#'   \item Builds a lookup table mapping groups to observation indices.
#'   \item For each group, selects up to `m` candidate neighbour observations
#'   from allowable predecessor groups using nearest-neighbour search.
#'   \item Computes recursive dependency closures used in downstream likelihood
#'   updates.
#' }
#'
#' Neighbour candidates are selected using Euclidean distances in `locations`.
#' When too many candidates are available, a reduced representative subset is
#' chosen using `kcenter_select()`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{constants}{List with structural constants:
#'     `m` (neighbour count), `G` (number of groups), and
#'     `mG` (maximum group size).}
#'   \item{groups}{Re-indexed integer group labels for each observation.}
#'   \item{groupNeighbours}{`G x m` matrix of selected neighbour observation
#'     indices for each group.}
#'   \item{groupMatrix}{Group dependency matrix used in the construction.}
#'   \item{groupLookup}{Lookup matrix mapping groups to observation indices.}
#'   \item{groupNum}{Vector of group sizes.}
#'   \item{groupDeps}{Expanded recursive dependency matrix.}
#'   \item{groupDepL}{Number of non-zero dependencies for each group.}
#' }
#'
#' @examples
#' \dontrun{
#' data(offshore)
#'
#' # Extract observed Z2 locations and their index coordinates from offshore
#' obs <- subset(offshore, !is.na(Z2))
#' locs <- obs[, c("d", "x", "y")]   # coordinate columns used for ordering
#' grps <- seq_len(nrow(obs))         # one observation per group (no grouping)
#'
#' V <- setupVecchiaGeoMix(
#'   locations    = locs,
#'   m            = 10,
#'   groups       = grps,
#'   depStructure = "minimal"
#' )
#'
#' # Key output components
#' V$constants$G    # number of groups
#' V$constants$m    # conditioning set size
#' dim(V$groupNeighbours)  # G x m neighbour index matrix
#'
#' # Pass a precomputed Vecchia object to setupGeoMixModel() to skip
#' # recomputation when experimenting with different MCMC settings
#' setup <- setupGeoMixModel(
#'   data      = offshore,
#'   K         = 3,
#'   dims      = c(20, 20, 20),
#'   variables = list(
#'     loc = "locID", xID = "x", yID = "y", dID = "d",
#'     x = "x", y = "y", depth = "d", Z1 = "Z1", Z2 = "Z2"
#'   ),
#'   aformula  = ~ d,
#'   vecchia   = V
#' )
#' }
#'
#' @seealso [setupGeoMixModel()]
#'
#' @export
setupVecchiaGeoMix <- function(locations,
                               m,
                               groups,
                               groupMatrix = NULL,
                               depStructure = "minimal"
){
  N <- length(groups)
  groups <- as.numeric(factor(groups))
  groupNum <- as.numeric(table(sort(groups)))
  G <- length(groupNum)

  if(!is.null(groupMatrix)){
    maxDeps <- ncol(groupMatrix)
    nDeps <- apply(groupMatrix>0,1,sum)
  }else if(depStructure == "full"){
    maxDeps <- G-1
    nDeps <- 0:(G-1)
    groupMatrix <- matrix(t(sapply(1:G,function(g) c(g:1,rep(0,G-g))))[,-1],ncol = maxDeps)
  }else if(depStructure == "previous"){
    maxDeps <- 1
    nDeps <- c(0,rep(1,G-1))
    groupMatrix <- matrix(1:G-1,ncol=maxDeps)
  }else if(depStructure == "minimal"){
    groupMask <- rbind(NA,t(sapply(1:(G-1),function(g)c(cumsum(groupNum[g:1]),rep(NA,G-g-1)))))
    find_minimum_group <- function(x){
      y <- x
      y[y < m] <- NA
      y <- which.min(y)
      if(length(y) == 0) y <- which.max(x[1:G])
      if(length(y) == 0) y <- 0
      return(y)
    }
    nDeps <- apply(groupMask,1,find_minimum_group)
    end <- 1:G-nDeps
    groupMatrix <- t(sapply(1:G,function(g) c(g:end[g],rep(0,max(nDeps)-nDeps[g]))))
    maxDeps <- ncol(groupMatrix) - 1
    groupMatrix <- matrix(groupMatrix[,-1],ncol = maxDeps)
  }

  groupLookup <- data.frame(id=1:N,groups) %>%
    group_by(groups) %>%
    mutate(ind = 1:n()) %>%
    pivot_wider(names_from = ind, values_from = id, values_fill = 0) %>%
    arrange(groups) %>% ungroup() %>% select(!groups) %>% as.matrix()
  mG <- ncol(groupLookup)

  groupNeighbours <- matrix(0,nrow = G, ncol = m)
  for(g in 1:G){
    ng <- groupNum[g]
    group_inds <- groupLookup[g,1:ng]
    if(nDeps[g]>0){
      query_inds <- unlist(lapply(groupMatrix[g,1:nDeps[g]],function(j) groupLookup[j,1:groupNum[j]]))
      ncands <- length(query_inds)
      if(ncands < m){
        candidates <- query_inds
      }else{
        NNlist <- FNN::get.knnx(
          query = locations[group_inds, , drop = FALSE],
          data  = locations[query_inds, , drop = FALSE],
          k     = m
        )
        NNDF <- data.frame(
          node     = rep(seq_len(ng), times = m),
          node_id  = rep(group_inds,  times = m),
          order    = rep(seq_len(m), each = ng),
          cand_id  = query_inds[c(NNlist$nn.ind)],
          dist     = c(NNlist$nn.dist)
        )
        max_order <- ceiling(m/ng) + 1
        NNDF_tmp <- filter(NNDF,order <=  max_order)
        candidates <-  unique(NNDF_tmp$cand_id)
        ncands_upd <- length(candidates)
        while(ncands_upd < m){
          max_order <- max_order + 1
          NNDF_tmp <- filter(NNDF,order <=  max_order)
          candidates <-  sort(unique(NNDF_tmp$cand_id))
          ncands_upd <- length(candidates)
        }
        if(ncands_upd > m){
          D <- proxy::dist(
            locations[candidates, , drop = FALSE],
            locations[group_inds, , drop = FALSE]
          )
          candidates <- candidates[kcenter_select(D,m)$chosen_rows]
        }
      }
      ncands <- length(candidates)
      groupNeighbours[g,1:ncands] <- candidates
    }
  }
  if(maxDeps > 5){
    CgetDependencies <- suppressMessages(compileNimble(getDependencies))
    groupDeps <- CgetDependencies(cbind(1:G,groupMatrix),G,maxDeps+1)
  }else{
    groupDeps <- getDependencies(cbind(1:G,groupMatrix),G,maxDeps+1)
  }
  groupDepL <- apply(groupDeps != 0,1,sum)
  if(ncol(groupLookup) == 1){
    groupLookup <- cbind(groupLookup, 0)
    mG <- mG + 1
  }
  out <- list(constants = list(m = m, G=G, mG = mG), groups = groups,
              groupNeighbours = groupNeighbours, groupMatrix = groupMatrix,
              groupLookup = groupLookup, groupNum = groupNum, groupDeps = groupDeps, groupDepL = groupDepL)
  return(out)
}
