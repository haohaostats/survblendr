
#' One-shot extrapolation: PEM(INLA) + anchored Gompertz + monotone spline weight
#' @export
aswb_extrapolate <- function(
    df, t_obs, t_max, interval = 0.5,
    nsim_inla = 1000, nsim_ext = 1000,
    inla_threads = 1,
    anchor_t = t_max,
    anchor_mean_Sa = 0.10, anchor_ess = 50,
    anchor_type = c("beta","fixed"),
    loggamma_mean = log(0.08), loggamma_sd = 0.5,
    scam_k = 15, tau_init = 1, tau_max = 100, tau_step = 1.2, tol = -1e-10,
    seed = NULL
) {
  anchor_type <- match.arg(anchor_type)
  
  # 🔒 1) 顶层一次性设置“局部”随机种子（仅在本函数生命周期内生效）
  if (!is.null(seed)) withr::local_seed(as.integer(seed))
  # 若用户没给 seed，R 会自动初始化；不再强行 set.seed()
  
  # 🔒 2) 锁定数值库线程，避免并行规约的非确定性
  withr::local_envvar(c(
    OMP_NUM_THREADS      = "1",
    MKL_NUM_THREADS      = "1",
    OPENBLAS_NUM_THREADS = "1"
  ))
  
  stopifnot(anchor_t <= t_max)
  
  # 观测侧
  pem <- pem_fit_inla(
    df, t_obs = t_obs, t_max = t_max, interval = interval,
    nsim = nsim_inla, inla_threads = inla_threads,
    seed_draws = if (is.null(seed)) NULL else as.integer(seed) + 100L
  )
  
  times      <- as.numeric(pem$time)
  S_obs_mean <- as.numeric(colMeans(pem$S_obs_draws))
  
  # 外部尾部（⚠️ 此处不再单独设 seed）
  ext <- draw_ext_gompertz(
    times,
    t_a = anchor_t, S_a_mean = anchor_mean_Sa, ess_Sa = anchor_ess,
    loggamma_mean = loggamma_mean, loggamma_sd = loggamma_sd,
    nsim = nsim_ext, anchor_type = anchor_type
  )
  S_ext_mean <- as.numeric(colMeans(ext$S_ext_draws))
  
  # 权重
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
