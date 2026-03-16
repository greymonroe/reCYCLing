#' Plot per-monomer mutation load relative to the array consensus
#'
#' For each monomer in the array, counts the number of positions that differ
#' from the consensus sequence (computed via \code{\link{counts_long_nogap}})
#' and draws a histogram of these mismatch counts. Higher average load
#' indicates greater sequence divergence within the array; the distribution
#' shape reflects the interplay between mutation rate and the homogenising
#' effect of duplication.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param counts Optional. Pre-computed output from
#'   \code{\link{counts_long_nogap}}. If \code{NULL} (default), it is computed
#'   internally. Pass a pre-computed value when calling multiple plot functions
#'   to avoid redundant computation.
#' @param bins Integer. Number of histogram bins. Default 80.
#' @param base_size Numeric. Base font size. Default 11.
#' @param fill Character. Fill colour for histogram bars.
#'   Default \code{"firebrick"}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_mutation_load(sim)
#'
#' # Pre-compute counts to share with plot_consensus_support()
#' ct <- counts_long_nogap(sim$monomers$bponly)
#' plot_mutation_load(sim, counts = ct)
#' }
plot_mutation_load <- function(sim, counts = NULL, bins = 80, base_size = 11,
                               fill = "firebrick") {
  monos <- as.data.table(sim$monomers)

  if (is.null(counts))
    counts <- counts_long_nogap(monos$bponly)
  counts_dt  <- as.data.table(counts)
  cons_dt    <- counts_dt[symbol == consensus][order(pos)]

  load_vec <- mu_load_from_consensus(monos$bponly, cons_dt)
  load_dt  <- data.table(load = load_vec)

  ggplot(load_dt, aes(x = load)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = fill, colour = "white", linewidth = 0.1,
                   alpha = 0.7) +
    geom_density(linewidth = 0.6, colour = "black") +
    theme_classic(base_size = base_size) +
    xlab("Mismatches to consensus") +
    ylab("Density")
}
