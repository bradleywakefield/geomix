#' Placeholder sampler that leaves a target unchanged
#'
#' Internal NIMBLE sampler used when a model node should remain fixed during
#' MCMC. The current state and associated log-probabilities are copied to the
#' saved modelValues object without proposing an update.
#'
#' @param model Compiled or uncompiled NIMBLE model object.
#' @param mvSaved NIMBLE modelValues object.
#' @param target Character string giving the sampler target node.
#' @param control Optional named list of sampler control settings.
#'
#' @details
#' Useful for disabling updates of selected parameters while preserving the
#' expected sampler interface within a configured MCMC.
#'
#' @return No return value. Copies the current model state to `mvSaved`.
#'
#' @keywords internal
DummySampler <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    returnType(void())
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    return()
  },
  methods = list(reset = function() {})
)
