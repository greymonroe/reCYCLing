#' Plot per-monomer mutation load relative to the array consensus
#'
#' For each monomer in the array, counts the number of positions that differ
#' from the consensus sequence (computed via \code{\link{counts_long_nogap}})
#' and draws a histogram of these mismatch counts. Higher average load
#' indicates greater sequence divergence within the array; the distribution
#' shape reflects the interplay between mutation rate and the homogenising
#' effect of duplication.
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
#' plot_mutation_load(sim, i = 1)
#'
#' # Pre-compute counts to share with plot_consensus_support()
#' mono <- sim$monomers_list[[1]]
#' ct   <- counts_long_nogap(mono$bponly)
#' plot_mutation_load(sim, i = 1, counts = ct)
#' }
plot_mutation_load <- function(ps_results, i, counts = NULL,
                               bins = 50, base_size = 6) {
  monos <- as.data.table(ps_results$monomers_list[[i]])

  if (is.null(counts))
    counts <- counts_long_nogap(monos$bponly)
  counts_dt  <- as.data.table(counts)
  cons_dt    <- counts_dt[symbol == consensus][order(pos)]

  load_vec <- mu_load_from_consensus(monos$bponly, cons_dt)
  load_dt  <- data.table(load = load_vec)

  ggplot(load_dt, aes(x = load)) +
    geom_histogram(bins = bins) +
    theme_classic(base_size = base_size) +
    xlab("Mismatches relative to consensus") +
    ylab("Count")
}
