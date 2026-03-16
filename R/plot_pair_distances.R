#' Plot distance distribution between identical repeat pairs
#'
#' Draws a histogram of the positional distance \eqn{|num_1 - num_2|} between
#' every pair of identical monomers in the array, using a log\eqn{_{10}}
#' x-axis. This distribution reflects the duplication history: a peak near
#' distance 1 indicates predominantly local (tandem) duplications, while a
#' broad or bimodal distribution suggests distal duplications as well.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param bins Integer. Number of histogram bins. Default 80.
#' @param base_size Numeric. Base font size passed to
#'   \code{theme_classic}. Default 11.
#' @param fill Character. Fill colour for histogram bars.
#'   Default \code{"steelblue"}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_pair_distances(sim)
#' }
plot_pair_distances <- function(sim, bins = 80, base_size = 11,
                                fill = "steelblue") {
  ps <- as.data.table(sim$ps)
  # Compute pairs on the fly if sim was run with compute_pairs = FALSE
  if (nrow(ps) == 0 && nrow(as.data.table(sim$monomers)) > 1)
    ps <- pairs_identical(as.data.table(sim$monomers))
  ps[, dist := abs(num1 - num2)]
  ps <- ps[dist > 0]

  ggplot(ps, aes(x = dist)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = fill, colour = "white", linewidth = 0.1,
                   alpha = 0.7) +
    geom_density(linewidth = 0.6, colour = "black") +
    scale_x_log10() +
    theme_classic(base_size = base_size) +
    xlab("Distance between identical pairs") +
    ylab("Density")
}
