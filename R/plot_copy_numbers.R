#' Plot copy-number distribution of unique monomers
#'
#' Draws a histogram of the number of copies of each distinct monomer sequence,
#' using a log\eqn{_{10}} x-axis. A highly right-skewed distribution (many
#' singletons, few high-copy sequences) suggests limited homogenisation;
#' a distribution shifted toward higher copy numbers indicates strong
#' concerted evolution.
#'
#' @param ps_results Output from \code{\link{run_sim_ps}}.
#' @param i Integer. Replicate index to plot.
#' @param bins Integer. Number of histogram bins. Default 50.
#' @param base_size Numeric. Base font size. Default 6.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' plot_copy_numbers(sim, i = 1)
#' }
plot_copy_numbers <- function(ps_results, i, bins = 50, base_size = 6) {
  monos   <- as.data.table(ps_results$monomers_list[[i]])
  copy_dt <- data.table(N = as.integer(table(monos$bponly)))

  ggplot(copy_dt, aes(x = N)) +
    geom_histogram(bins = bins) +
    scale_x_log10() +
    theme_classic(base_size = base_size) +
    xlab("Monomer copy number") +
    ylab("Count")
}
