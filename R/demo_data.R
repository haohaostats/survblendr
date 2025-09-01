
#' Simulate Weibull data with target administrative + random censoring
#'
#' \code{S(20)=S20} calibration on Weibull(shape, scale).
#'
#' @param n Sample size.
#' @param shape Weibull shape.
#' @param S20 Target survival at t=20.
#' @param censor_rate Target overall censoring rate in [0,1].
#' @param t_max Administrative censoring time (default 20).
#' @param seed Optional RNG seed.
#'
#' @return Data frame with \code{id, time, status, true_time}.
#' @export
demo_sim_weibull <- function(n = 2000, shape = 1.6, S20 = 0.20,
                             censor_rate = 0.30, t_max = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  scale <- t_max / (-log(S20))^(1/shape)
  Ttrue <- rweibull(n, shape = shape, scale = scale)
  time  <- pmin(Ttrue, t_max)
  status <- as.integer(Ttrue <= t_max)
  
  n_target <- round(n * censor_rate)
  n_admin  <- sum(status == 0L)
  n_need   <- max(0L, n_target - n_admin)
  
  idx_elig <- which(status == 1L & time < t_max)
  if (length(idx_elig) && n_need > 0L) {
    pick <- sample(idx_elig, min(n_need, length(idx_elig)))
    ctim <- runif(length(pick), min = tiny_time, max = time[pick])
    earlier <- ctim < time[pick]
    time[pick[earlier]]   <- ctim[earlier]
    status[pick[earlier]] <- 0L
  }
  time <- pmax(time, tiny_time)
  data.frame(id = seq_len(n), time = time, status = status, true_time = Ttrue)
}
