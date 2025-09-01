
#' Extrapolation summary table (annual survival and RMST)
#'
#' @description
#' Produce a tidy table over the **extrapolation window** (from `t_obs`
#' to the end of the grid), reporting the blended survival at each whole
#' year and the RMST of the blended curve from `t_obs` up to that year.
#'
#' @param fit A list returned by [survblendr_extrapolate()].
#' @param by Step size for reporting (time units). Default 1 (yearly).
#' @return A `data.frame` with columns:
#'   - `time`: reporting time points (>= `t_obs`).
#'   - `S_blend`: blended survival at `time`.
#'   - `S_obs`, `S_ext`: component survivals (for reference).
#'   - `RMST_blend`: integral of blended survival on `[t_obs, time]`.
#' @details
#' Survival values are taken on the nearest grid point of `fit$time`.
#' @export
survblendr_extrap_table <- function(fit, by = 1) {
  stopifnot(is.list(fit), !is.null(fit$time), !is.null(fit$S_blend))
  time <- as.numeric(fit$time)
  t0   <- if (!is.null(fit$t_obs)) fit$t_obs else min(time)

  t_end <- seq(ceiling(t0), max(time), by = by)
  if (!length(t_end)) t_end <- max(time)
  
  pick_idx <- vapply(t_end, function(tt)
    which.min(abs(time - tt)), integer(1))
  
  getv <- function(v) if (is.null(v)) rep(NA_real_, length(pick_idx)) else as.numeric(v)[pick_idx]
  
  S_bl <- getv(fit$S_blend)
  S_ob <- getv(fit$S_obs)
  S_ex <- getv(fit$S_ext)
  
  rmst_vec <- vapply(t_end, function(tt)
    rmst_numeric(time, as.numeric(fit$S_blend), t0 = t0, t_star = tt),
    numeric(1))
  
  out <- data.frame(
    time = t_end,
    S_blend = S_bl,
    S_obs   = S_ob,
    S_ext   = S_ex,
    RMST_blend = rmst_vec,
    check.names = FALSE
  )
  rownames(out) <- NULL
  out
}
