
#' Summarize survival and RMST over calendar times
#'
#' @description
#' Create a table of survival at selected times and the Restricted Mean
#' Survival Time (RMST) **from 0 up to each time** for the observed-side
#' PEM curve, the external Gompertz tail, and the blended curve.
#'
#' @param fit A list returned by [aswb_extrapolate()].
#' @param years Numeric vector of target times (e.g., `10:25`).
#'
#' @return A data.frame with columns:
#' `time, S_obs, S_ext, S_blend, RMST_obs, RMST_ext, RMST_blend`.
#' @export
aswb_summary_table <- function(fit, years) {
  stopifnot(all(c("time","S_obs","S_ext","S_blend") %in% names(fit)))
  times <- as.numeric(fit$time)
  Sobs  <- as.numeric(fit$S_obs)
  Sext  <- as.numeric(fit$S_ext)
  Sbl   <- as.numeric(fit$S_blend)
  
  # --- helpers ---------------------------------------------------------------
  # trapezoid rule
  trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, 1)) / 2)
  
  # linear on log-S tail extrapolation (stable for survival tails)
  s_at <- function(x, y, x0) {
    eps <- 1e-12
    if (x0 <= max(x) && x0 >= min(x)) {
      return(stats::approx(x, y, xout = x0, rule = 2)$y)
    }
    # extrapolate in log-S using last two points
    if (x0 > max(x)) {
      x1 <- x[length(x) - 1]; x2 <- x[length(x)]
      y1 <- pmax(y[length(y) - 1], eps); y2 <- pmax(y[length(y)], eps)
    } else { # x0 < min(x), use first two points (rare in practice)
      x1 <- x[1]; x2 <- x[2]
      y1 <- pmax(y[1], eps); y2 <- pmax(y[2], eps)
    }
    b <- (log(y2) - log(y1)) / (x2 - x1)
    a <- log(y2) - b * x2
    val <- exp(a + b * x0)
    pmin(pmax(val, 0), 1)
  }
  
  # RMST from 0 to t (allow t between/ beyond grid with interpolation/extrapolation)
  rmst_to <- function(t, x, y) {
    if (t <= min(x)) {
      return(t * s_at(x, y, t))  # degenerate short span
    }
    idx <- which(x <= t)
    x2  <- c(x[idx], t)
    y2  <- c(y[idx], s_at(x, y, t))
    trapz(x2, y2)
  }
  # ---------------------------------------------------------------------------
  
  S_obs_t <- vapply(years, function(tt) s_at(times, Sobs, tt), numeric(1))
  S_ext_t <- vapply(years, function(tt) s_at(times, Sext, tt), numeric(1))
  S_bl_t  <- vapply(years, function(tt) s_at(times, Sbl,  tt), numeric(1))
  
  RMST_obs <- vapply(years, function(tt) rmst_to(tt, times, Sobs), numeric(1))
  RMST_ext <- vapply(years, function(tt) rmst_to(tt, times, Sext), numeric(1))
  RMST_bl  <- vapply(years, function(tt) rmst_to(tt, times, Sbl ), numeric(1))
  
  out <- data.frame(
    time = years,
    S_obs = S_obs_t,
    S_ext = S_ext_t,
    S_blend = S_bl_t,
    RMST_obs = RMST_obs,
    RMST_ext = RMST_ext,
    RMST_blend = RMST_bl,
    row.names = NULL
  )
  out
}
