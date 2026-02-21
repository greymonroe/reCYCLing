#' Plot per-position consensus support across the array
#'
#' Draws a histogram of the proportion of monomers that agree with the
#' consensus base at each alignment position. A distribution concentrated near
#' 1 indicates strong conservation; a spread distribution indicates
#' heterogeneity—which can result from a high mutation rate, multiple
#' co-existing sequence variants, or both.
#'
#' @param ps_results Output from \code{\link{run_sim_ps}}.
#' @param i Integer. Replicate index to plot.
#' @param counts Optional. Pre-computed output from
#'   \code{\link{counts_long_nogap}} for this replicate. If \code{NULL}
#'   (default), it is computed internally. Pass a pre-computed value when
#'   calling multiple plot functions on the same replicate to avoid redundant
#'   computation.
#' @param bins Integer. Number of histogram bins. Default 50.
#' @param base_size Numeric. Base font size. Default 6.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' plot_consensus_support(sim, i = 1)
#' }
plot_consensus_support <- function(ps_results, i, counts = NULL,
                                   bins = 50, base_size = 6) {
  monos <- as.data.table(ps_results$monomers_list[[i]])

  if (is.null(counts))
    counts <- counts_long_nogap(monos$bponly)
  counts_dt  <- as.data.table(counts)
  cons_dt    <- counts_dt[symbol == consensus][order(pos)]

  ggplot(cons_dt, aes(x = prop)) +
    geom_histogram(bins = bins) +
    theme_classic(base_size = base_size) +
    xlab("Consensus base proportion") +
    ylab("Count")
}
