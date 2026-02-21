#' Combined multi-panel summary plot for one simulation replicate
#'
#' Assembles five diagnostic panels into a single figure using
#' \pkg{cowplot}. The panels are produced by the individual plot functions
#' and can also be generated independently:
#' \itemize{
#'   \item \strong{Fingerprint} (\code{\link{plot_repeat_fingerprint}}): dot-plot
#'     of identical monomer pairs.
#'   \item \strong{Pair distances} (\code{\link{plot_pair_distances}}): histogram
#'     of \eqn{|num_1 - num_2|} on a log\eqn{_{10}} scale.
#'   \item \strong{Copy numbers} (\code{\link{plot_copy_numbers}}): histogram of
#'     copies per unique monomer sequence.
#'   \item \strong{Mutation load} (\code{\link{plot_mutation_load}}): histogram
#'     of per-monomer mismatches to the consensus.
#'   \item \strong{Consensus support} (\code{\link{plot_consensus_support}}):
#'     histogram of per-position consensus proportions.
#' }
#'
#' @param ps_results Output from \code{\link{run_sim_ps}}.
#' @param i Integer. Replicate index to plot.
#' @param title Character or \code{NULL}. Optional title placed above the
#'   fingerprint panel. Default \code{NULL}.
#' @param rel_widths_bottom Numeric vector of length 4. Relative widths of the
#'   four bottom-row panels (pair distances, copy numbers, mutation load,
#'   consensus support). Default \code{c(1, 1, 1, 1)}.
#' @param rel_heights Numeric vector of length 2. Relative heights of the
#'   fingerprint row and the bottom row. Default \code{c(2, 1)}.
#'
#' @return A \pkg{cowplot} grid object (a \code{ggplot2} grob).
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' plot_ps_summary(sim, i = 1, title = "mu = 3e-5")
#' }
plot_ps_summary <- function(ps_results, i,
                             title = NULL,
                             rel_widths_bottom = c(1, 1, 1, 1),
                             rel_heights = c(2, 1)) {

  # Pre-compute consensus once; pass to the two panels that need it to
  # avoid running counts_long_nogap() twice.
  monos  <- as.data.table(ps_results$monomers_list[[i]])
  counts <- counts_long_nogap(monos$bponly)

  p_dist <- plot_pair_distances(ps_results, i)
  p_fp   <- plot_repeat_fingerprint(ps_results, i,
                                    rotate = FALSE, ptsize = 0.01,
                                    gridalpha = 0)
  p_cn   <- plot_copy_numbers(ps_results, i)
  p_ml   <- plot_mutation_load(ps_results, i, counts = counts)
  p_cs   <- plot_consensus_support(ps_results, i, counts = counts)

  if (!is.null(title))
    p_fp <- p_fp + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  bottom_row <- cowplot::plot_grid(p_dist, p_cn, p_ml, p_cs,
                                   nrow = 1,
                                   rel_widths = rel_widths_bottom)
  cowplot::plot_grid(p_fp, bottom_row,
                     ncol = 1,
                     rel_heights = rel_heights)
}
