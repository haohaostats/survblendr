
#' Numerical RMST on [t0, t*] by trapezoidal rule
#' @param times grid
#' @param S survival on grid
#' @param t0 lower limit
#' @param t_star upper limit
#' @return scalar
#' @export
rmst_numeric <- function(times, S, t0, t_star) {
  idx <- which(times >= t0 & times <= t_star)
  trapz(times[idx], S[idx])
}

#' Error metrics for extrapolated survival curves
#' @param times grid
#' @param S_true true survival
#' @param S_model model survival
#' @param t0 observation horizon
#' @param eval_times vector of evaluation times (default 15, 20)
#' @return list of metrics
#' @export
compute_metrics <- function(times, S_true, S_model, t0, eval_times = c(15, 20)) {
  surv_err <- sapply(eval_times, function(tstar) {
    i <- which.min(abs(times - tstar))
    abs(S_true[i] - S_model[i])
  })
  names(surv_err) <- paste0("AE_Surv_", eval_times, "y")
  
  rmst_true  <- sapply(eval_times, function(tstar) rmst_numeric(times, S_true,  t0, tstar))
  rmst_model <- sapply(eval_times, function(tstar) rmst_numeric(times, S_model, t0, tstar))
  rmst_err   <- abs(rmst_true - rmst_model)
  names(rmst_err) <- paste0("AE_RMST_", eval_times, "y")
  
  idx_ex <- which(times > t0)
  mae <- mean(abs(S_true[idx_ex] - S_model[idx_ex]))
  fd  <- abs(rmst_err[2] - rmst_err[1])
  
  list(surv_err = surv_err, rmst_err = rmst_err, MAE = mae, FD = fd,
       RMST_true = rmst_true, RMST_model = rmst_model)
}


