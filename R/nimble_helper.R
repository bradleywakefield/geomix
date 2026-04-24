#' Construct neighbour indices for a lattice
#'
#' Internal NIMBLE helper that generates a neighbour matrix for a regular
#' lattice with optional inclusion of diagonal neighbours.
#'
#' @param dims Numeric vector giving the lattice dimensions.
#' @param diagonals Integer scalar; if `1`, include diagonal neighbours,
#' otherwise include only axis-aligned neighbours.
#'
#' @return A numeric matrix whose rows give neighbour indices for each lattice
#' site.
#'
#' @keywords internal
getNeighbours <- nimble::nimbleFunction(
  run = function(dims = double(1),
                 diagonals = integer(0)) {
    returnType(double(2))

    N <- prod(dims)
    d <- length(dims)

    x <- array(1:N, dim = dims, nDim = d)
    y <- array(0, dim = dims+2, nDim = d)

    yindices <- matrix(0, nrow = N, ncol = d)
    for (i in 1:d) {
      if (i == 1) block_size <- 1 else block_size <- prod((dims[1:(i-1)]))
      if (i == d) repeat_size <- 1 else repeat_size <- prod((dims[(i+1):d]))
      pattern <- rep(rep(1:dims[i], each = block_size), times = repeat_size)
      yindices[,i] <- pattern
    }
    yvec <- rep(1,N)
    dims0 <-  c(1, dims[1:(d-1)]+2)
    for(j in 1:d) yvec <- yvec + yindices[,j]*prod(dims0[1:j])
    y[yvec] <- x

    n_offsets <- 3^d
    offsets <- matrix(0, nrow=n_offsets, ncol=d)
    for (i in 1:d) {
      block_size <- 3^(i-1)
      repeat_size <- 3^(d-i)
      pattern <- rep(rep(-1:1, each=block_size), times=repeat_size)
      offsets[,i] <- pattern
    }
    if(diagonals==1){
      valid_offset <- numeric(0)
      for(j in 1:n_offsets) if(sum(abs(offsets)[j,])!=0) valid_offset <- c(valid_offset,j)
      n_offsets <- n_offsets - 1
      f_offsets <- offsets[valid_offset,]
    }else{
      valid_offset <- numeric(0)
      for(j in 1:n_offsets) if(sum(abs(offsets)[j,])==1) valid_offset <- c(valid_offset,j)
      f_offsets <- offsets[valid_offset,]
      n_offsets <- length(valid_offset)
    }

    valid_pts <- which(y>0)
    neighbours <- matrix(0,nrow = N, ncol = n_offsets)
    numNeighbours <- numeric(N)
    dims2 <- c(1,dims+2)
    for(i in 1:N){
      yindex <- rep(valid_pts[i],n_offsets)
      for(j in 1:d){
        yindex <- yindex+f_offsets[,j]*prod(dims2[1:j])
      }
      ypos <- which(y[yindex]>0)
      yzero <- which(y[yindex]==0)
      yord1 <- c(ypos,yzero)
      yord <- yindex[yord1]
      neighbours[i,] <- y[yord]
    }

    return(neighbours)
  })
#' Identify dependency relationships from neighbour sets
#'
#' Internal NIMBLE helper that computes reverse dependencies from a neighbour
#' matrix, recording which rows depend on each index.
#'
#' @param vneighbours Numeric matrix of neighbour indices.
#' @param N2 Integer scalar giving the number of indexed units.
#' @param m Integer scalar giving the maximum number of neighbours per unit.
#'
#' @return A numeric matrix of dependency indices.
#'
#' @keywords internal
getDependencies <- nimble::nimbleFunction(
  run = function(vneighbours = double(2),
                 N2 = integer(0),
                 m = integer(0)) {
    returnType(double(2))
    dependencies <- matrix(0, nrow = N2, ncol = N2)
    maxDep <- 0
    for (i in 1:N2) {
      count <- 0
      for (k in 1:m) {
        for (j in 1:N2) {
          if (vneighbours[j,k] == i) {
            count <- count + 1
            dependencies[i,count] <- j
          }
        }
      }
      if (count > maxDep) {
        maxDep <- count
      }
      if (i %% ceiling(N2*0.01) == 0 | i == N2) {
        cat("-")
      }
    }
    # Only keep the first maxDep columns
    result <- dependencies[ , 1:maxDep]
    return(result)
  }
)
#' Dummy categorical density
#'
#' Internal NIMBLE distribution used as a placeholder likelihood for latent
#' categorical variables whose log-probability is supplied externally.
#'
#' @param x Numeric vector of latent class labels.
#' @param N Integer scalar giving the length of `x`.
#' @param logProb Numeric scalar giving the externally computed log-density.
#' @param log Integer scalar; if `1`, return the log-density.
#'
#' @return A numeric scalar giving the density or log-density.
#'
#' @keywords internal
dDummyCat <- nimble::nimbleFunction(
  run = function(x = double(1),
                 N = integer(0),
                 logProb = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(log == 1) return(logProb) else return(exp(logProb))
  }
)
#' Dummy categorical random generator
#'
#' Internal NIMBLE random-generation function paired with `dDummyCat`. This
#' acts as a placeholder sampler and returns a zero vector of the requested
#' length.
#'
#' @param n Integer scalar giving the number of draws.
#' @param N Integer scalar giving the length of the output vector.
#' @param logProb Numeric scalar giving the externally computed log-density.
#'
#' @return A numeric vector of length `N`.
#'
#' @keywords internal
rDummyCat <- nimble::nimbleFunction(
  run = function(n = integer(0), N = integer(0), logProb = double(0)) {
    returnType(double(1))
    x <- rep(0,N)
    return(x[1:N])
  }
)
#' Compute pseudo-likelihood Potts log-probability approximation
#'
#' Internal NIMBLE helper that evaluates the weighted Potts interaction term
#' for a discrete latent field under a supplied neighbourhood structure.
#'
#' @param beta Numeric scalar controlling Potts interaction strength.
#' @param Y1 Numeric vector of class labels.
#' @param weights Numeric matrix of neighbour-specific weights.
#' @param neighbours Numeric matrix of neighbour indices for each site.
#' @param neighbourNum Numeric vector giving the number of neighbours per site.
#' @param penalty Numeric matrix of pairwise class penalties.
#'
#' @return A numeric scalar giving the Potts log-probability contribution.
#'
#' @keywords internal
computePotts <- nimble::nimbleFunction(
  run = function(beta = double(0), Y1 = double(1), weights = double(2),
                 neighbours = double(2), neighbourNum = double(1), penalty = double(2)) {
    returnType(double(0))
    N1 <- length(Y1)
    logProb <- 0.0
    for(s in 1:N1){
      M <- neighbourNum[s]
      neighID <- neighbours[s,1:M]
      k <- Y1[s]
      for(j in 1:M){
        sn <- neighID[j]
        kn <- Y1[sn]
        logProb <- logProb - beta*weights[s,j]*penalty[k,kn]

      }
    }
    return(logProb)
  })
#' Sample from a weighted Potts model using the psedo-likelihood approximation
#'
#' Internal NIMBLE helper that updates a discrete latent field using
#' sequential conditional draws under a weighted Potts model.
#'
#' @param beta Numeric scalar controlling Potts interaction strength.
#' @param Y1 Numeric vector of current class labels.
#' @param K Integer scalar giving the number of classes.
#' @param weights Numeric matrix of neighbour-specific weights.
#' @param neighbours Numeric matrix of neighbour indices for each site.
#' @param neighbourNum Numeric vector giving the number of neighbours per site.
#' @param penalty Numeric matrix of pairwise class penalties.
#'
#' @return A numeric vector of updated class labels.
#'
#' @keywords internal
samplePotts <- nimble::nimbleFunction(
  run = function(beta = double(0), Y1 = double(1), K = integer(0), weights = double(2),
                 neighbours = double(2), neighbourNum = double(1), penalty = double(2)) {
    returnType(double(1))
    N1 <- length(Y1)
    Y1star <- Y1
    for(s in 1:N1){
      if(Y1[s] > 0){
        logProb <- rep(0,K)
        M <- neighbourNum[s]
        neighID <- neighbours[s,1:M]
        for(k in 1:K){
          for(j in 1:M){
            sn <- neighID[j]
            kn <- Y1star[sn]
            logProb[k] <- logProb[k] - beta*weights[s,j]*penalty[k,kn]
          }
        }
        Y1star[s] <- rcat(n=1, prob = exp(logProb - max(logProb)))
      }
    }
    return(Y1star)
  })
