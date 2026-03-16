#' Launch the interactive simulation explorer
#'
#' Opens a Shiny app that lets you tweak all simulation parameters with
#' sliders and immediately visualise the resulting repeat array. Requires
#' the \pkg{shiny} package.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' launch_explorer()
#' }
launch_explorer <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("Install the 'shiny' package to use the explorer app", call. = FALSE)

  app_dir <- system.file("shiny", "explorer", package = "reCYCLing")
  if (app_dir == "")
    stop("Could not find the explorer app. Try re-installing reCYCLing.",
         call. = FALSE)

  shiny::runApp(app_dir, ...)
}

#' Launch the ABC parameter inference dashboard
#'
#' Opens a Shiny app that runs SMC-ABC parameter inference in real time.
#' You can watch the convergence live, inspect posterior distributions of
#' parameters as they narrow, and explore the best-fitting simulations
#' interactively. Requires the \pkg{shiny} package.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' launch_abc()
#' }
launch_abc <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("Install the 'shiny' package to use the ABC app", call. = FALSE)

  app_dir <- system.file("shiny", "abc", package = "reCYCLing")
  if (app_dir == "")
    stop("Could not find the ABC app. Try re-installing reCYCLing.",
         call. = FALSE)

  shiny::runApp(app_dir, ...)
}
