#' Plot copy-number distribution of unique monomers
#'
#' Draws a histogram of the number of copies of each distinct monomer sequence,
#' using a log\eqn{_{10}} x-axis. A highly right-skewed distribution (many
#' singletons, few high-copy sequences) suggests limited homogenisation;
#' a distribution shifted toward higher copy numbers indicates strong
#' concerted evolution.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param bins Integer. Number of histogram bins. Default 80.
#' @param base_size Numeric. Base font size. Default 11.
#' @param fill Character. Fill colour for histogram bars.
#'   Default \code{"darkorange"}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_copy_numbers(sim)
#' }
plot_copy_numbers <- function(sim, bins = 80, base_size = 11,
                              fill = "darkorange") {
  monos   <- as.data.table(sim$monomers)
  copy_dt <- data.table(N = as.integer(table(monos$bponly)))

  ggplot(copy_dt, aes(x = N)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = fill, colour = "white", linewidth = 0.1,
                   alpha = 0.7) +
    geom_density(linewidth = 0.6, colour = "black") +
    scale_x_log10() +
    theme_classic(base_size = base_size) +
    xlab("Monomer copy number") +
    ylab("Density")
}
