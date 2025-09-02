
#' @title Summarize survival and extrapolation RMST at target times
#' @name survblendr_summary_table
#'
#' @description
#' Returns `S_obs`, `S_ext`, `S_blend` and **extrapolation RMST** for each
#' requested time in `years`. Extrapolation RMST is defined as
#' \eqn{\int_{t_{\mathrm{obs}}}^{t} S(u)\,du}. For `years` below `t_obs`,
#' RMST is 0 by definition.
#'
#' RMST is computed on the **reporting grid** (`years`) using the trapezoidal
#' rule with the same survival values you see in the table, so increments
#' match the displayed S at the year endpoints and remain monotone non-decreasing.
#'
#' @param fit A list returned by [survblendr_extrapolate()].
#' @param years Numeric vector of target times (e.g., `10:25`).
#' @return A data.frame with columns:
#'   `time, S_obs, S_ext, S_blend, RMST_obs, RMST_ext, RMST_blend`.
#' @export
survblendr_summary_table <- function(fit, years) {
  stopifnot(all(c("time","S_obs","S_ext","S_blend") %in% names(fit)))
  
  times <- as.numeric(fit$time)
  Sobs  <- pmin(pmax(as.numeric(fit$S_obs),   0), 1)
  Sext  <- pmin(pmax(as.numeric(fit$S_ext),   0), 1)
  Sbl   <- pmin(pmax(as.numeric(fit$S_blend), 0), 1)
  
  t0 <- if (!is.null(fit$t_obs)) fit$t_obs else min(times)
  
  s_at <- function(S, tt) {
    vals <- stats::approx(times, S, xout = tt, rule = 2)$y
    pmin(pmax(vals, 0), 1)
  }
  
  yrs_in  <- as.numeric(years)
  yrs_ord <- order(yrs_in)
  yrs     <- yrs_in[yrs_ord]
  grid    <- sort(unique(c(t0, yrs)))
  
  Sobs_g <- s_at(Sobs, grid)
  Sext_g <- s_at(Sext, grid)
  Sbl_g  <- s_at(Sbl,  grid)
  
  cumtrapz_from_t0 <- function(x, y) {
    i0 <- match(t0, x)
    n  <- length(x)
    area <- numeric(n)
    if (n >= 2) {
      for (i in seq.int(i0 + 1L, n)) {
        area[i] <- area[i - 1L] + (x[i] - x[i - 1L]) * (y[i] + y[i - 1L]) / 2
      }
    }
    area
  }
  
  A_obs <- cumtrapz_from_t0(grid, Sobs_g)
  A_ext <- cumtrapz_from_t0(grid, Sext_g)
  A_bl  <- cumtrapz_from_t0(grid, Sbl_g)
  
  idx <- match(yrs, grid)
  RMST_obs   <- A_obs[idx]
  RMST_ext   <- A_ext[idx]
  RMST_blend <- A_bl[idx]
  
  RMST_obs[order(yrs)]   <- cummax(RMST_obs[order(yrs)])
  RMST_ext[order(yrs)]   <- cummax(RMST_ext[order(yrs)])
  RMST_blend[order(yrs)] <- cummax(RMST_blend[order(yrs)])
  
  S_obs_t <- s_at(Sobs, yrs)
  S_ext_t <- s_at(Sext, yrs)
  S_bl_t  <- s_at(Sbl,  yrs)
  
  out <- data.frame(
    time       = yrs_in,
    S_obs      = S_obs_t[match(yrs_in, yrs)],
    S_ext      = S_ext_t[match(yrs_in, yrs)],
    S_blend    = S_bl_t [match(yrs_in, yrs)],
    RMST_obs   = RMST_obs  [match(yrs_in, yrs)],
    RMST_ext   = RMST_ext  [match(yrs_in, yrs)],
    RMST_blend = RMST_blend[match(yrs_in, yrs)],
    row.names  = NULL,
    check.names = FALSE
  )
  
  out
}
