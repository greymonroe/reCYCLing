#' Integer axis-break generator with fixed step
#'
#' Returns a function that, given a numeric range, produces integer break
#' positions at a fixed interval. Useful as the \code{breaks} argument in
#' \pkg{ggplot2} scale functions.
#'
#' @param step Numeric. Step size between breaks. Default 1.
#'
#' @return A function \code{function(x)} that returns a numeric vector of
#'   break positions.
#' @export
#'
#' @examples
#' br <- int_breaks(step = 5)
#' br(c(0, 23))  # 0 5 10 15 20
int_breaks <- function(step = 1) {
  function(x) seq(floor(min(x, na.rm = TRUE)),
                  ceiling(max(x, na.rm = TRUE)), by = step)
}

#' Pretty integer axis-break generator
#'
#' Returns a function that uses \code{\link[base]{pretty}} to choose
#' approximately \code{n} integer break positions within a range. Useful as
#' the \code{breaks} argument in \pkg{ggplot2} scale functions.
#'
#' @param n Integer. Approximate number of breaks. Default 5.
#'
#' @return A function \code{function(x)} that returns a numeric vector of
#'   integer break positions.
#' @export
#'
#' @examples
#' br <- int_breaks_pretty(n = 4)
#' br(c(0, 100))
int_breaks_pretty <- function(n = 5) {
  function(x) {
    b <- pretty(range(x, na.rm = TRUE), n)
    b[b == round(b)]
  }
}
