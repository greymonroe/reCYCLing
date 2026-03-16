#' Combined multi-panel summary plot
#'
#' Assembles five diagnostic panels into a single figure using
#' \pkg{cowplot}: the repeat fingerprint on top (full width), and four
#' histogram panels in a 2x2 grid below.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param title Character or \code{NULL}. Optional title placed above the
#'   fingerprint panel. Default \code{NULL}.
#' @param rel_heights Numeric vector of length 2. Relative heights of the
#'   fingerprint row and the histogram grid. Default \code{c(2, 1.5)}.
#'
#' @return A \pkg{cowplot} grid object (a \code{ggplot2} grob).
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_ps_summary(sim, title = "mu = 3e-5")
#' }
plot_ps_summary <- function(sim, title = NULL, rel_heights = c(2, 1.5)) {

  monos  <- as.data.table(sim$monomers)
  counts <- counts_long_nogap(monos$bponly)

  # Use smaller base_size for the multi-panel layout
  p_fp   <- plot_repeat_fingerprint(sim, rotate = FALSE, ptsize = 0.01,
                                    gridalpha = 0)
  p_dist <- plot_pair_distances(sim, base_size = 10)
  p_cn   <- plot_copy_numbers(sim, base_size = 10)
  p_ml   <- plot_mutation_load(sim, counts = counts, base_size = 10)
  p_cs   <- plot_consensus_support(sim, counts = counts, base_size = 10)

  if (!is.null(title))
    p_fp <- p_fp + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  # 2x2 grid for the histograms
  hist_grid <- cowplot::plot_grid(p_dist, p_cn, p_ml, p_cs,
                                  nrow = 2, ncol = 2)

  cowplot::plot_grid(p_fp, hist_grid,
                     ncol = 1,
                     rel_heights = rel_heights)
}
