

#' Plot blended vs components with extrapolation marker (+ KM step; seamless at t0)
#'
#' Shows Kaplan–Meier (step), Cox(PEM), External(Gompertz) and Blended(Spline).
#' For Cox/External, draw solid up to t0 and dashed from t0 onward, with no visual gap.
#'
#' @param fit object from [survblendr_extrapolate()], must contain $time, $S_obs, $S_ext, $S_blend;
#'   if `fit$train` has `time`/`status`, a KM curve is added automatically.
#' @param add_true optional numeric vector (same length as `fit$time`) for simulations.
#' @param t0 optional; if NULL, inferred from `fit$t_obs` or max(fit$train$time).
#' @param shade_extrap logical; shade t >= t0.
#' @export
plot_curves <- function(
    fit, add_true = NULL, t0 = NULL, shade_extrap = TRUE,
    t0_color = "#B22222", t0_linetype = "longdash", t0_linewidth = 1.8,
    t0_label = "Extrapolation starts", show_t0_label = TRUE
) {
  stopifnot(is.list(fit), !is.null(fit$time))
  time <- as.numeric(fit$time)
  
  t0_line <- if (!is.null(t0)) {
    as.numeric(t0)
  } else if (!is.null(fit$t_obs)) {
    as.numeric(fit$t_obs)
  } else if (!is.null(fit$train) && !is.null(fit$train$time)) {
    suppressWarnings(max(as.numeric(fit$train$time), na.rm = TRUE))
  } else NA_real_
  
  cols <- c(
    "True"                = "#000000",
    "Kaplan–Meier"        = "#56B4E9",
    "Cox (PEM)"           = "#4D4D4D",
    "External (Gompertz)" = "#9467BD",
    "Blended (Spline)"    = "#009E73"
  )
  legend_levels <- c("True","Kaplan–Meier","Cox (PEM)","External (Gompertz)","Blended (Spline)")
  
  ser <- list(
    "Cox (PEM)"           = fit$S_obs,
    "External (Gompertz)" = fit$S_ext,
    "Blended (Spline)"    = fit$S_blend
  )
  if (!is.null(add_true)) ser[["True"]] <- add_true
  
  long <- do.call(
    rbind,
    lapply(names(ser), function(nm) {
      y <- ser[[nm]]; if (is.null(y)) return(NULL)
      data.frame(time = time, S = as.numeric(y), series = nm, stringsAsFactors = FALSE)
    })
  )
  if (is.null(long) || !nrow(long)) stop("Nothing to plot.")
  long <- long[is.finite(long$S), , drop = FALSE]
  long$series <- factor(long$series, levels = legend_levels)
  
  km_df <- NULL
  if (!is.null(fit$train) && all(c("time","status") %in% names(fit$train))) {
    sf <- survival::survfit(survival::Surv(time, status) ~ 1, data = fit$train)
    if (is.finite(t0_line)) {
      keep <- sf$time <= t0_line
      km_time <- c(0, sf$time[keep])
      km_surv <- c(1, sf$surv[keep])
    } else {
      km_time <- c(0, sf$time)
      km_surv <- c(1, sf$surv)
    }
    km_df <- data.frame(time = km_time, S = km_surv, series = "Kaplan–Meier")
  }
  
  split_series <- function(df, nm) {
    d <- df[df$series == nm, , drop = FALSE]
    if (!nrow(d)) return(list(pre = NULL, post = NULL))
    if (!is.finite(t0_line)) return(list(pre = d, post = NULL))
    if (!any(abs(d$time - t0_line) < 1e-12)) {
      i <- max(which(d$time < t0_line))
      if (is.finite(i) && i < nrow(d)) {
        t1 <- d$time[i]; t2 <- d$time[i + 1]
        s1 <- d$S[i];    s2 <- d$S[i + 1]
        s0 <- s1 + (s2 - s1) * (t0_line - t1) / (t2 - t1)
        d <- rbind(d, data.frame(time = t0_line, S = s0, series = nm))
        d <- d[order(d$time), , drop = FALSE]
      }
    }
    list(pre  = d[d$time <= t0_line, , drop = FALSE],
         post = d[d$time >= t0_line, , drop = FALSE])
  }
  
  cx <- split_series(long, "Cox (PEM)")
  ex <- split_series(long, "External (Gompertz)")
  bl <- long[long$series == "Blended (Spline)", , drop = FALSE]
  tr <- long[long$series == "True",              , drop = FALSE]
  
  p <- ggplot2::ggplot()
  
  if (isTRUE(shade_extrap) && is.finite(t0_line)) {
    p <- p + ggplot2::annotate("rect",
                               xmin = t0_line, xmax = max(long$time, na.rm = TRUE),
                               ymin = -Inf,    ymax = Inf, alpha = 0.08)
  }
  
  if (!is.null(km_df)) {
    p <- p + ggplot2::geom_step(
      data = km_df,
      ggplot2::aes(x = time, y = S, color = "Kaplan–Meier"),
      direction = "hv", linewidth = 1.4, na.rm = TRUE
    )
  }
  
  # True
  if (nrow(tr)) p <- p + ggplot2::geom_line(
    data = tr,
    ggplot2::aes(x = time, y = S, color = "True"),
    linewidth = 1.4, na.rm = TRUE
  )
  
  if (nrow(bl)) p <- p + ggplot2::geom_line(
    data = bl,
    ggplot2::aes(x = time, y = S, color = "Blended (Spline)"),
    linewidth = 1.4, na.rm = TRUE
  )
  
  if (nrow(cx$pre))  p <- p + ggplot2::geom_line(
    data = cx$pre,
    ggplot2::aes(x = time, y = S, color = "Cox (PEM)"),
    linewidth = 1.4, linetype = "solid", lineend = "butt", na.rm = TRUE
  )
  if (nrow(cx$post)) p <- p + ggplot2::geom_line(
    data = cx$post,
    ggplot2::aes(x = time, y = S, color = "Cox (PEM)"),
    linewidth = 1.4, linetype = "dashed", lineend = "butt", na.rm = TRUE
  )
  
  if (nrow(ex$pre))  p <- p + ggplot2::geom_line(
    data = ex$pre,
    ggplot2::aes(x = time, y = S, color = "External (Gompertz)"),
    linewidth = 1.4, linetype = "solid", lineend = "butt", na.rm = TRUE
  )
  if (nrow(ex$post)) p <- p + ggplot2::geom_line(
    data = ex$post,
    ggplot2::aes(x = time, y = S, color = "External (Gompertz)"),
    linewidth = 1.4, linetype = "dashed", lineend = "butt", na.rm = TRUE
  )
  
  p <- p +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Time", y = "Survival", color = NULL) +
    ggplot2::scale_color_manual(values = cols, breaks = legend_levels, drop = FALSE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.box = "horizontal") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, override.aes = list(linewidth = 1.8)))
  
  if (is.finite(t0_line)) {
    p <- p + ggplot2::geom_vline(xintercept = t0_line,
                                 linetype = t0_linetype,
                                 linewidth = t0_linewidth,
                                 color = t0_color)
    if (isTRUE(show_t0_label)) {
      p <- p + ggplot2::annotate("text", x = t0_line, y = 0.04,
                                 label = paste0(t0_label, " \u2192"),
                                 hjust = 0, vjust = 0, size = 3.4, color = t0_color)
    }
  }
  
  p
}
