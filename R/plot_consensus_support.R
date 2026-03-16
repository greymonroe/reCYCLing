#' Plot per-position consensus support across the array
#'
#' Draws a histogram of the proportion of monomers that agree with the
#' consensus base at each alignment position. A distribution concentrated near
#' 1 indicates strong conservation; a spread distribution indicates
#' heterogeneity—which can result from a high mutation rate, multiple
#' co-existing sequence variants, or both.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param counts Optional. Pre-computed output from
#'   \code{\link{counts_long_nogap}}. If \code{NULL} (default), it is computed
#'   internally. Pass a pre-computed value when calling multiple plot functions
#'   to avoid redundant computation.
#' @param bins Integer. Number of histogram bins. Default 80.
#' @param base_size Numeric. Base font size. Default 11.
#' @param fill Character. Fill colour for histogram bars.
#'   Default \code{"seagreen"}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_consensus_support(sim)
#' }
plot_consensus_support <- function(sim, counts = NULL, bins = 80,
                                   base_size = 11, fill = "seagreen") {
  monos <- as.data.table(sim$monomers)

  if (is.null(counts))
    counts <- counts_long_nogap(monos$bponly)
  counts_dt  <- as.data.table(counts)
  cons_dt    <- counts_dt[symbol == consensus][order(pos)]

  ggplot(cons_dt, aes(x = prop)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = fill, colour = "white", linewidth = 0.1,
                   alpha = 0.7) +
    geom_density(linewidth = 0.6, colour = "black") +
    theme_classic(base_size = base_size) +
    xlab("Consensus base proportion") +
    ylab("Density")
}
