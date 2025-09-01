
#' Blend two survival curves on the cumulative-hazard scale
#'
#' @param S_obs Observed-side survival.
#' @param S_ext External-side survival.
#' @param w Weight in [0,1] on the time grid.
#' @return Blended survival \eqn{S_{\mathrm{obs}}^{1-w} S_{\mathrm{ext}}^w}.
#' @export
blend_curves <- function(S_obs, S_ext, w) {
  stopifnot(length(S_obs) == length(S_ext), length(w) == length(S_obs))
  S_obs^(1 - w) * S_ext^w
}
