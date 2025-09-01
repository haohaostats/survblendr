
#' Example survival dataset (Weibull, admin follow-up to 10)
#'
#' @description
#' A simulated cohort used as the package example. Time is administratively
#' truncated at 10 (same unit as your analysis; e.g., years). Overall censoring
#' is ~40%. Generated via [demo_sim_weibull()] with parameters
#' `shape = 1.6`, `S20 = 0.20`, `t_max = 10`, `n = 1500`, `seed = 1`.
#'
#' @format A data frame with 1,500 rows and 3 columns:
#' \describe{
#'   \item{id}{Integer subject id (1..n).}
#'   \item{time}{Nonnegative event/censoring time (administratively truncated at 10).}
#'   \item{status}{Event indicator (1 = event, 0 = censored).}
#' }
#'
#' @usage data(aswb_demo)
#' @keywords datasets
#' @examples
#' data(aswb_demo)
#' head(aswb_demo)
#' fit <- aswb_extrapolate(aswb_demo, t_obs = 10, t_max = 20, interval = 0.5,
#'                         nsim_inla = 200, nsim_ext = 200, inla_threads = 1)
#' p <- plot_curves(fit); print(p)
"aswb_demo"
