#' Plot distance distribution between identical repeat pairs
#'
#' Draws a histogram of the positional distance \eqn{|num_1 - num_2|} between
#' every pair of identical monomers in the array, using a log\eqn{_{10}}
#' x-axis. This distribution reflects the duplication history: a peak near
#' distance 1 indicates predominantly local (tandem) duplications, while a
#' broad or bimodal distribution suggests distal duplications as well.
#'
#' @param ps_results Output from \code{\link{run_sim_ps}}.
#' @param i Integer. Replicate index to plot.
#' @param bins Integer. Number of histogram bins. Default 50.
#' @param base_size Numeric. Base font size passed to
#'   \code{theme_classic}. Default 6.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' plot_pair_distances(sim, i = 1)
#' }
plot_pair_distances <- function(ps_results, i, bins = 50, base_size = 6) {
  ps <- as.data.table(ps_results$ps_list[[i]])

  ggplot(ps, aes(x = abs(num1 - num2))) +
    geom_histogram(bins = bins) +
    scale_x_log10() +
    theme_classic(base_size = base_size) +
    xlab("|distance between identical pairs|") +
    ylab("Count")
}
