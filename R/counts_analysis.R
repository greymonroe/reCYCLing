#' Per-position symbol counts across aligned sequences
#'
#' Computes per-position symbol counts and proportions across a set of aligned
#' sequences of equal length. Positions whose consensus base is a gap character
#' are removed and the remaining positions are renumbered.
#'
#' @param rawseqs Character vector of aligned sequences (all must have the same
#'   \code{nchar()}). Typically \code{monomers_list[[i]]$bponly} from
#'   \code{\link{run_sim_ps}}.
#' @param gap_char Character. The gap symbol. Default \code{"-"}.
#' @param tie_break_levels Character vector giving the priority order for
#'   breaking consensus ties. The gap character is always given lowest priority
#'   regardless of its position in this vector.
#'
#' @return A long-format \code{data.table} with one row per
#'   (non-gap-consensus position, symbol) combination and columns:
#'   \describe{
#'     \item{pos_new}{Renumbered position index (gaps removed).}
#'     \item{pos}{Original aligned position index.}
#'     \item{symbol}{The nucleotide or gap character.}
#'     \item{count}{Number of sequences with this symbol at this position.}
#'     \item{total}{Total sequences at this position.}
#'     \item{prop}{Proportion (\code{count / total}).}
#'     \item{consensus}{Consensus base at this position.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 300, mu_total = 3e-5)
#' mono <- sim$monomers_list[[1]]
#' ct   <- counts_long_nogap(mono$bponly)
#' head(ct)
#' }
counts_long_nogap <- function(rawseqs,
                              gap_char = "-",
                              tie_break_levels = c("A","C","G","T","N","R","Y",
                                                   "S","W","K","M","B","D","H",
                                                   "V","-")) {
  stopifnot(length(rawseqs) >= 2L)
  x  <- as.character(rawseqs)
  Ls <- nchar(x)
  if (length(unique(Ls)) != 1L)
    stop("All sequences must have the same aligned length")
  L <- Ls[1L]
  n <- length(x)

  chars     <- strsplit(x, "", fixed = TRUE)
  mat       <- do.call(rbind, chars)
  syms      <- sort(unique(c(mat)))
  counts_all <- data.table::as.data.table(
    sapply(syms, function(s) colSums(mat == s, na.rm = TRUE))
  )
  counts_all[, pos := seq_len(L)]
  data.table::setcolorder(counts_all, c("pos", setdiff(names(counts_all), "pos")))

  long_counts <- data.table::melt(
    counts_all,
    id.vars      = "pos",
    variable.name = "symbol",
    value.name   = "count"
  )
  long_counts[, total  := sum(count), by = pos]
  long_counts[, prop   := ifelse(total > 0, count / total, NA_real_)]
  long_counts[, symbol := as.character(symbol)]

  pri_map <- seq_along(tie_break_levels)
  names(pri_map) <- tie_break_levels
  if (!gap_char %in% names(pri_map)) {
    tie_break_levels        <- c(tie_break_levels, gap_char)
    pri_map                 <- seq_along(tie_break_levels)
    names(pri_map)          <- tie_break_levels
  }
  pri_map[gap_char] <- -1L
  long_counts[, pri := pri_map[symbol]]

  consensus_dt <- long_counts[
    , {
      m    <- max(count)
      cand <- .SD[count == m]
      cand[which.max(pri)][, .(consensus = symbol)]
    },
    by = pos
  ]

  keep_map <- consensus_dt[consensus != gap_char, .(pos, consensus)]
  keep_map[, pos_new := seq_len(.N)]

  long_counts_nogap <- keep_map[long_counts, on = "pos", nomatch = 0L]
  data.table::setcolorder(long_counts_nogap,
                          c("pos_new", "pos", "symbol", "count",
                            "total", "prop", "consensus", "pri"))
  long_counts_nogap[, pri := NULL][]
}

#' Mutation load of each sequence relative to a consensus
#'
#' Counts the number of positions at which each sequence in \code{rawseqs}
#' differs from the consensus sequence. This gives a simple per-monomer
#' measure of divergence.
#'
#' @param rawseqs Character vector of aligned sequences (all same length).
#' @param consensus A \code{data.table} with columns \code{pos} and
#'   \code{consensus} (one row per alignment position), as returned by
#'   filtering \code{counts_long_nogap()} to rows where
#'   \code{symbol == consensus}.
#'
#' @return An integer vector of the same length as \code{rawseqs}, giving the
#'   number of mismatches to consensus for each sequence.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 300, mu_total = 3e-5)
#' mono <- sim$monomers_list[[1]]
#' ct   <- counts_long_nogap(mono$bponly)
#' cons <- ct[symbol == consensus][order(pos)]
#' load <- mu_load_from_consensus(mono$bponly, cons)
#' hist(load)
#' }
mu_load_from_consensus <- function(rawseqs, consensus) {
  x  <- as.character(rawseqs)
  Ls <- nchar(x)
  if (length(unique(Ls)) != 1L)
    stop("All sequences must have the same aligned length")
  L <- Ls[1L]
  consensus <- consensus[order(pos)]
  if (length(consensus$consensus) != L)
    stop("Consensus length must match sequence length")

  chars   <- strsplit(x, "", fixed = TRUE)
  sapply(chars, function(bases) sum(bases != consensus$consensus))
}
