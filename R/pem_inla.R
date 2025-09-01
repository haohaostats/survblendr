
#' Piecewise-Exponential (PEM) short-term fit via INLA
#'
#' Fits a baseline hazard with \code{inla.coxph(..., control.hazard = list(model="rw1"))}
#' on administratively truncated data (at \code{t_obs}). Returns posterior draws
#' of the observed-side survival curve on a regular grid.
#'
#' @param df Data frame with columns \code{time}, \code{status}.
#' @param t_obs Administrative cut at observation horizon.
#' @param t_max Maximum time grid for posterior survival.
#' @param interval Grid step (must divide \code{t_max}).
#' @param nsim Number of posterior draws.
#' @param inla_threads Integer threads for INLA (passed to \code{INLA::inla.setOption}).
#'
#' @return List with \code{time}, \code{S_obs_draws}, \code{train}, \code{inla}.
#' @importFrom INLA inla.coxph inla inla.posterior.sample inla.setOption
#' @export
pem_fit_inla <- function(df, t_obs, t_max, interval = 1, nsim = 2000, inla_threads = 1) {
  # truncate to t_obs
  train <- df
  train$status[train$time > t_obs] <- 0L
  train$time[train$time > t_obs]   <- t_obs
  train$time[!is.finite(train$time)] <- t_obs
  train$time <- pmax(train$time, tiny_time)
  
  if (sum(train$status == 1L) == 0L)
    stop("No events before t_obs; cannot estimate short-term baseline hazard.")
  
  stopifnot(t_max > 0, interval > 0, (t_max %% interval) == 0)
  cutpoints <- seq(0, t_max, by = interval)
  
  # hazard with RW1 prior
  p <- INLA::inla.coxph(
    INLA::inla.surv(time, status) ~ -1,
    data = train,
    control.hazard = list(
      constr = FALSE,
      cutpoints = cutpoints,
      model = "rw1",
      hyper = list(prec = list(prior = "loggamma", param = c(0.01, 0.01))),
      scale.model = TRUE
    )
  )
  
  old_threads <- tryCatch(INLA::inla.setOption(num.threads = inla_threads), error = function(e) NULL)
  on.exit(try(INLA::inla.setOption(num.threads = old_threads), silent = TRUE), add = TRUE)
  
  m <- INLA::inla(
    p$formula, family = p$family,
    data = c(as.list(p$data), p$data.list),
    control.compute = list(config = TRUE, dic = TRUE),
    control.fixed = list(mean = 0, prec = 1/1000)
  )
  
  sbh <- m$summary.random$baseline.hazard
  if (is.null(sbh)) stop("INLA: baseline.hazard summary missing; check control.hazard.")
  
  K_expected <- length(cutpoints) - 1L
  K_got <- nrow(sbh)
  if (K_got == K_expected + 1L) sbh <- sbh[seq_len(K_expected), , drop = FALSE]
  if (nrow(sbh) != K_expected) stop("Baseline segments mismatch.")
  
  # joint posterior sample -> baseline hazard draws
  use_marginal_fallback <- FALSE
  jp <- try(INLA::inla.posterior.sample(n = nsim, result = m), silent = TRUE)
  if (inherits(jp, "try-error") || length(jp) == 0L) use_marginal_fallback <- TRUE
  
  h0_draws <- NULL
  if (!use_marginal_fallback) {
    nm <- names(jp[[1]]$latent)
    idx_base <- grep("^baseline\\.hazard", nm, ignore.case = TRUE)
    if (!length(idx_base)) {
      contents <- m$misc$configs$contents
      if (!is.null(contents)) {
        row_bh <- which(tolower(contents$tag) %in% c("baseline.hazard","baseline hazard","baselinehazard"))
        if (length(row_bh)) {
          row_bh <- row_bh[1]
          start  <- contents$start[row_bh]
          len    <- min(contents$length[row_bh], K_expected)
          idx_base <- start:(start + len - 1L)
        }
      }
    }
    if (length(idx_base) != K_expected) {
      use_marginal_fallback <- TRUE
    } else {
      h0_draws <- t(vapply(jp, function(s) exp(as.numeric(s$latent[idx_base])),
                           FUN.VALUE = numeric(K_expected)))
    }
  }
  
  if (use_marginal_fallback) {
    mr <- m$marginals.random$baseline.hazard
    if (is.null(mr) || length(mr) < K_expected)
      stop("Marginal fallback unavailable or too short.")
    if (length(mr) > K_expected) mr <- mr[seq_len(K_expected)]
    h0_draws <- matrix(NA_real_, nrow = nsim, ncol = K_expected)
    for (k in seq_len(K_expected)) h0_draws[, k] <- exp(INLA::inla.rmarginal(nsim, mr[[k]]))
  }
  
  H0_mat <- t(apply(h0_draws, 1, cumsum)) * interval
  S_obs_draws <- cbind(1, exp(-H0_mat))  # include S(0)=1
  tt <- seq(0, by = interval, length.out = ncol(S_obs_draws))
  list(time = tt, S_obs_draws = S_obs_draws, train = train,
       inla = list(model = m, cutpoints = cutpoints))
}
