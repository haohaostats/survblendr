
#' One-shot extrapolation: PEM(INLA) + anchored Gompertz + monotone spline weight
#' @export
aswb_extrapolate <- function(
    df, t_obs, t_max, interval = 0.5,
    nsim_inla = 1000, nsim_ext = 1000,
    inla_threads = 1,
    # --- anchor controls exposed to users ---
    anchor_t = t_max,
    anchor_mean_Sa = 0.10, anchor_ess = 50,
    anchor_type = c("beta","fixed"),
    # Gompertz gamma prior
    loggamma_mean = log(0.08), loggamma_sd = 0.5,
    # weight (SCAM) controls
    scam_k = 15, tau_init = 1, tau_max = 100, tau_step = 1.2, tol = -1e-10,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  anchor_type <- match.arg(anchor_type)
  stopifnot(anchor_t <= t_max)
  
  # short-term observation-side
  pem <- pem_fit_inla(
    df, t_obs = t_obs, t_max = t_max, interval = interval,
    nsim = nsim_inla, inla_threads = inla_threads
  )
  times      <- as.numeric(pem$time)
  S_obs_mean <- as.numeric(colMeans(pem$S_obs_draws))
  
  # external tail with chosen anchor
  ext <- draw_ext_gompertz(
    times,
    t_a = anchor_t, S_a_mean = anchor_mean_Sa, ess_Sa = anchor_ess,
    loggamma_mean = loggamma_mean, loggamma_sd = loggamma_sd,
    nsim = nsim_ext, seed = seed, anchor_type = anchor_type
  )
  S_ext_mean <- as.numeric(colMeans(ext$S_ext_draws))
  
  # monotone spline weight with temperature calibration
  wfit <- weight_fit_scam_tau(
    times, S_obs = S_obs_mean, S_ext = S_ext_mean,
    J = scam_k, tau_step = tau_step, tau_max = tau_max
  )
  w <- wfit$pi
  S_blend <- blend_curves(S_obs_mean, S_ext_mean, w)
  
  structure(list(
    time    = times,
    S_obs   = S_obs_mean,
    S_ext   = S_ext_mean,
    S_blend = S_blend,
    w       = w,
    weight  = wfit,
    train   = pem$train,
    t_obs   = t_obs,
    anchor  = list(t = anchor_t, S_mean = anchor_mean_Sa, ess = anchor_ess, type = anchor_type)
  ), class = "aswb_fit")
}
