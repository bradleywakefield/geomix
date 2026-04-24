#' Select representative rows using a k-center criterion
#'
#' Internal helper that selects `m` rows of a distance matrix to provide
#' representative coverage of the columns. At each step, the algorithm chooses
#' the row that minimises the maximum distance from each column to its nearest
#' selected row.
#'
#' @param D Numeric distance matrix with candidate rows in rows and target
#'   points in columns.
#' @param m Integer number of rows to select.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{chosen_rows}{Indices of the selected rows.}
#'     \item{radius}{Maximum distance from any column to its nearest selected
#'       row after selection.}
#'     \item{assignment}{For each column, the index of its nearest selected row
#'       among `chosen_rows`.}
#'   }
#' @noRd
kcenter_select <- function(D, m) {
  
  nc <- nrow(D)
  nn <- ncol(D)
  
  chosen <- integer(0)
  
  # best distance from each node to chosen candidates
  best <- rep(Inf, nn)
  
  for (t in seq_len(m)) {
    score <- rep(Inf, nc)
    
    for (i in seq_len(nc)) {
      if (i %in% chosen) next
      new_best <- pmin(best, D[i, ])
      score[i] <- max(new_best)
    }
    
    i_star <- which.min(score)
    
    chosen <- c(chosen, i_star)
    best <- pmin(best, D[i_star, ])
  }
  
  list(
    chosen_rows = chosen,
    radius = max(best),
    assignment = apply(D[chosen, , drop = FALSE], 2, which.min)
  )
}