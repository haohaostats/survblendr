
#' Numerical RMST on [t0, t*] with endpoint interpolation
#'
#' @description
#' Computes RMST over [t0, t_star] by linearly interpolating S at the
#' interval endpoints and applying the trapezoidal rule. With S >= 0,
#' the result is non-decreasing in t_star (up to tiny FP noise).
#'
#' @param times Strictly increasing time grid.
#' @param S Survival values on the same grid (in [0, 1]).
#' @param t0 Lower limit (e.g., t_obs for extrapolation RMST).
#' @param t_star Upper limit.
#' @return Scalar RMST on [t0, t_star].
#' @export
rmst_numeric <- function(times, S, t0, t_star) {
  if (t_star <= t0) return(0)
  stopifnot(is.numeric(times), length(times) == length(S))
  stopifnot(all(diff(times) > 0))
  
  # Clamp survival for numerical robustness
  S <- pmin(pmax(S, 0), 1)
  
  s_at <- function(x, y, x0) stats::approx(x, y, xout = x0, rule = 2)$y
  
  inside <- times > t0 & times < t_star
  x2 <- c(t0, times[inside], t_star)
  y2 <- c(s_at(times, S, t0), S[inside], s_at(times, S, t_star))
  
  # Final clamp in case interpolation produced tiny excursions
  y2 <- pmin(pmax(y2, 0), 1)
  
  # Trapezoidal rule
  sum(diff(x2) * (head(y2, -1) + tail(y2, 1)) / 2)
}

