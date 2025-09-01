
#' Package utilities (numeric helpers)
#' @keywords internal

#' Clamp to [0, 1]
#' @param x numeric
#' @return numeric
#' @noRd
clip01 <- function(x) pmin(pmax(x, 0), 1)

#' @keywords internal
#' @noRd
tiny_time <- 1e-6

#' Trapezoidal rule integral
#' @param x numeric increasing grid
#' @param y numeric values on x
#' @return scalar integral
#' @noRd
trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, 1)) / 2)

#' Approximate hazard from survival on a regular grid
#' @param times numeric increasing
#' @param S survival in [0,1]
#' @return numeric hazard approximation
#' @noRd
hazard_from_S <- function(times, S) {
  stopifnot(length(times) == length(S))
  eps <- .Machine$double.eps
  H <- -log(pmax(S, eps))
  G <- length(times)
  h <- numeric(G)
  if (G == 1) return(h)
  h[1] <- (H[2] - H[1]) / (times[2] - times[1])
  if (G > 2) {
    for (g in 2:(G - 1)) {
      h[g] <- (H[g + 1] - H[g - 1]) / (times[g + 1] - times[g - 1])
    }
  }
  h[G] <- (H[G] - H[G - 1]) / (times[G] - times[G - 1])
  h
}

#' @param x numeric increasing
#' @param y numeric values
#' @return numeric derivative dy/dx
#' @noRd
num_deriv <- function(x, y) {
  G <- length(x); dy <- numeric(G)
  if (G == 1) return(dy)
  dy[1] <- (y[2] - y[1]) / (x[2] - x[1])
  if (G > 2) {
    for (g in 2:(G - 1)) {
      dy[g] <- (y[g + 1] - y[g - 1]) / (x[g + 1] - x[g - 1])
    }
  }
  dy[G] <- (y[G] - y[G - 1]) / (x[G] - x[G - 1])
  dy
}
