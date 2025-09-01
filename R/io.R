
#' Path to the packaged example CSV
#' @return A file path to `inst/extdata/survblendr_demo.csv`.
#' @export
survblendr_example_data_path <- function() {
  system.file("extdata", "survblendr_demo.csv",
              package = "survblendr", mustWork = TRUE)
}


#' Read user-provided survival data
#'
#' @description
#' Reads a CSV and returns a standardized data.frame with columns
#' `id` (optional), `time`, `status`. `status` is coerced to 0/1 with
#' `event_value` indicating the event label in the file.
#'
#' @param file CSV path.
#' @param time_col Column name for time.
#' @param status_col Column name for status.
#' @param id_col Optional id column name; if NULL, a sequence is created.
#' @param event_value Value in `status_col` that means "event" (default 1).
#' @return A data.frame with columns `id`, `time`, `status`.
#' @export
survblendr_read_csv <- function(file,
                          time_col = "time",
                          status_col = "status",
                          id_col = NULL,
                          event_value = 1) {
  df <- utils::read.csv(file, check.names = FALSE)
  if (!time_col %in% names(df) || !status_col %in% names(df)) {
    stop("Columns `", time_col, "` and `", status_col, "` must exist in the CSV.")
  }
  id <- if (!is.null(id_col) && id_col %in% names(df)) df[[id_col]] else seq_len(nrow(df))
  time <- suppressWarnings(as.numeric(df[[time_col]]))
  status_raw <- df[[status_col]]
  status <- as.integer(status_raw == event_value)
  
  if (any(!is.finite(time) | time < 0)) stop("Non-finite or negative times found.")
  if (any(!status %in% c(0L, 1L))) stop("Status could not be coerced to 0/1.")
  
  out <- data.frame(id = as.integer(id), time = time, status = status)
  out$time <- pmax(out$time, .Machine$double.eps) # avoid exact zeros
  out
}
