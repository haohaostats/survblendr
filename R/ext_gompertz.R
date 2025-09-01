
#' Draw external Gompertz survival curves given an anchor constraint
#' @param times grid
#' @param t_a anchor time
#' @param S_a_mean anchor survival mean S(t_a); if anchor_type="fixed" this is exact S(t_a)
#' @param ess_Sa equivalent sample size for Beta prior on S(t_a) (ignored if fixed)
#' @param anchor_type "beta" (default) or "fixed"
#' @param loggamma_mean,loggamma_sd log-normal prior for Gompertz gamma
#' @export
draw_ext_gompertz <- function(times, t_a, S_a_mean, ess_Sa = 50,
                              loggamma_mean, loggamma_sd,
                              nsim = 2000,
                              anchor_type = c("beta", "fixed")) {
  anchor_type <- match.arg(anchor_type)
  
  # 直接抽样；随机性已由 aswb_extrapolate() 的 local_seed 统一控制
  if (anchor_type == "fixed" || is.infinite(ess_Sa)) {
    Sa <- rep(S_a_mean, nsim)
  } else {
    alpha <- S_a_mean * ess_Sa
    beta  <- (1 - S_a_mean) * ess_Sa
    Sa    <- stats::rbeta(nsim, alpha, beta)
  }
  gamma <- stats::rlnorm(nsim, meanlog = loggamma_mean, sdlog = loggamma_sd)
  gamma <- pmax(gamma, 1e-6)
  Sa    <- clip01(Sa)
  
  rho <- - gamma * log(Sa) / (exp(gamma * t_a) - 1)
  rho <- pmax(rho, 1e-10)
  
  tt <- matrix(times, nrow = nsim, ncol = length(times), byrow = TRUE)
  g  <- matrix(gamma, nrow = nsim, ncol = length(times))
  r  <- matrix(rho,   nrow = nsim, ncol = length(times))
  S  <- exp(-(r / g) * (exp(g * tt) - 1))
  
  list(S_ext_draws = S, gamma = gamma, rho = rho, Sa = Sa, t_a = t_a)
}
