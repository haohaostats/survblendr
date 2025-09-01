
#' Monotone spline weight with temperature calibration (SCAM + logistic)
#'
#' Fits a monotone increasing P-spline (via \pkg{scam}, \code{bs="mpi"}) to
#' z(t) = log H_obs(t) - log H_ext(t), then maps to w(t) = logit^{-1}{eta(t)/tau}.
#' The temperature tau is increased multiplicatively until the blended hazard
#' h_blend = (1-w) h_obs + w h_ext + (H_ext - H_obs) * w'(t) is nonnegative on the grid.
#'
#' @param times Numeric grid (strictly increasing).
#' @param S_obs Numeric mean observed-side survival on `times`.
#' @param S_ext Numeric mean external-side survival on `times`.
#' @param J Basis dimension for `s(..., bs="mpi")`. Default 10.
#' @param tau_step Multiplicative step for temperature. Default 1.5.
#' @param tau_max Max temperature. Default 1e6.
#'
#' @return List with `pi` (weights), `dpi`, `tau`, `eta`,
#'   `hmin` (min blended hazard), and `ok`.
#' @importFrom scam scam
#' @importFrom mgcv s
#' @importFrom stats plogis predict
#' @export
weight_fit_scam_tau <- function(times, S_obs, S_ext, J = 10, tau_step = 1.5, tau_max = 1e6) {
  if (!requireNamespace("scam", quietly = TRUE))
    stop("Package 'scam' is required. Please install.packages('scam').")
  
  eps <- .Machine$double.eps
  H_obs <- -log(pmax(S_obs, eps))
  H_ext <- -log(pmax(S_ext, eps))
  z <- log(pmax(H_obs, eps)) - log(pmax(H_ext, eps))
  dat <- data.frame(t = times, z = z)
  fit <- scam::scam(z ~ s(t, k = J, bs = "mpi"), data = dat, family = stats::gaussian())
  
  eta  <- as.numeric(stats::predict(fit, newdata = dat, type = "link"))
  deta <- num_deriv(times, eta)
  
  h_obs  <- hazard_from_S(times, S_obs)
  h_ext  <- hazard_from_S(times, S_ext)
  H_diff <- H_ext - H_obs
  denom  <- pmax(abs(H_diff), 1e-10)
  
  plogis_tau <- function(eta, tau) 1 / (1 + exp(-eta / tau))
  
  tau <- 1; ok <- FALSE
  while (!ok && tau <= tau_max) {
    pi_tau  <- plogis_tau(eta, tau)
    dpi_tau <- (1 / tau) * pi_tau * (1 - pi_tau) * deta
    bound   <- ((1 - pi_tau) * h_obs + pi_tau * h_ext) / denom
    ok <- all(dpi_tau <= bound + 1e-12)
    if (!ok) tau <- tau * tau_step
  }
  
  pi_tau  <- plogis_tau(eta, tau)
  dpi_tau <- (1 / tau) * pi_tau * (1 - pi_tau) * deta
  h_blend <- (1 - pi_tau) * h_obs + pi_tau * h_ext + H_diff * dpi_tau
  
  list(
    pi = pmin(pmax(pi_tau, 0), 1),
    dpi = dpi_tau,
    tau = tau,
    eta = eta,
    ok  = all(h_blend >= -1e-10),
    hmin = min(h_blend, na.rm = TRUE),
    scam = fit
  )
}
