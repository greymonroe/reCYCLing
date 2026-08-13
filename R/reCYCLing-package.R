#' reCYCLing: Repeat Element Cyclical Evolutionary Simulations in Genomes
#'
#' Simulates the evolution of tandem repeat arrays—such as satellite
#' DNA-like monomers—through local duplication, distal duplication, chunk
#' deletion, and per-base mutation. Companion visualization functions let
#' you explore the resulting arrays from multiple angles.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item Run one or more replicates with \code{\link{run_sim_ps}}.
#'   \item Inspect the output object (see \code{vignette("reCYCLing-intro")}).
#'   \item Plot individual panels: \code{\link{plot_pair_distances}},
#'         \code{\link{plot_repeat_fingerprint}},
#'         \code{\link{plot_copy_numbers}},
#'         \code{\link{plot_mutation_load}},
#'         \code{\link{plot_consensus_support}}.
#'   \item Produce a combined summary with \code{\link{plot_ps_summary}}.
#' }
#'
#' @docType package
#' @name reCYCLing-package
#' @aliases reCYCLing
#'
#' @import data.table
#' @import ggplot2
#' @importFrom parallel mclapply
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rgamma rgeom rnorm rpois runif quantile
#' @useDynLib reCYCLing, .registration = TRUE
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table column-name variables used
# throughout the package via non-standard evaluation.
utils::globalVariables(c(
  ".", ".I", ".N", ".SD",
  "angle", "bponly", "chrom", "chrom1", "chrom2",
  "consensus", "count", "dir", "factorlevel", "fill",
  "hap", "hjust", "id", "id1", "id2",
  "lab", "label", "load", "meanN", "N",
  "num", "num1", "num2",
  "pos", "pos1", "pos2", "pos_new", "pri", "prop",
  "sample1", "sample2", "symbol", "total",
  "x", "xend", "y", "yend",
  "sid", "cidx", "is_x", "chrom_id", "aligned",
  "Tot", "k", "between", "rows", "inverted"
))
