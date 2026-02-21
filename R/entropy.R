#' Shannon entropy of a character vector
#'
#' Computes the Shannon entropy (in bits, log base 2) of the frequency
#' distribution of elements in \code{x}. Higher values indicate greater
#' diversity among monomer sequences.
#'
#' @param x Character vector of monomer sequences.
#'
#' @return A single non-negative numeric value (entropy in bits).
#' @export
#'
#' @examples
#' shannon_entropy(c("AAA", "AAA", "AAA"))  # 0 — all identical
#' shannon_entropy(c("AAA", "CCC", "GGG"))  # log2(3) — all distinct
shannon_entropy <- function(x) {
  p <- table(x)
  p <- p[p > 0] / sum(p)
  -sum(p * log2(p))
}

#' Extended Shannon entropy summary
#'
#' Returns a \code{data.table} combining Shannon entropy with array-size
#' statistics. Useful for comparing diversity across parameter sweeps.
#'
#' @param x Character vector of monomer sequences.
#'
#' @return A one-row \code{data.table} with columns:
#'   \describe{
#'     \item{H}{Shannon entropy (bits).}
#'     \item{uniqueN}{Number of distinct sequences.}
#'     \item{L}{Total array length (number of monomers).}
#'     \item{H_logL}{H normalised by log(L).}
#'     \item{H_logUniqueN}{H normalised by log(uniqueN).}
#'   }
#' @export
#'
#' @examples
#' shannon_entropyplus(c("AAA", "CCC", "AAA", "GGG"))
shannon_entropyplus <- function(x) {
  p       <- table(x)
  p       <- p[p > 0] / sum(p)
  H       <- -sum(p * log2(p))
  uN      <- uniqueN(x)
  L       <- length(x)
  data.table(H          = H,
             uniqueN    = uN,
             L          = L,
             H_logL     = H / log(L),
             H_logUniqueN = H / log(uN))
}

#' Sliding-window entropy along a monomer array
#'
#' Applies \code{\link{shannon_entropyplus}} over a sliding window of
#' \code{N} consecutive rows of a monomer data.table, returning per-window
#' entropy statistics.
#'
#' @param kb A monomer \code{data.table} with a \code{bponly} column (such as
#'   an element of \code{monomers_list} from \code{\link{run_sim_ps}}).
#' @param N Integer. Window size in rows.
#'
#' @return A \code{data.table} with one row per window position and columns
#'   from \code{\link{shannon_entropyplus}} plus \code{pos} (focal row index).
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 200, mu_total = 3e-5)
#' kb  <- sim$monomers_list[[1]]
#' win <- kb_row_window_entropy(kb, N = 50)
#' }
kb_row_window_entropy <- function(kb, N) {
  kb <- as.data.table(kb)

  out_list <- pbapply::pblapply(seq_len(nrow(kb) - N), function(i) {
    sub_idx <- i:min(nrow(kb), i + N)
    res     <- shannon_entropyplus(kb$bponly[sub_idx])
    res$pos <- i
    res
  })

  rbindlist(out_list, fill = TRUE)
}

#' Summarise entropy across simulation replicates
#'
#' Extracts the final Shannon entropy, total array size, and unique-monomer
#' count for each replicate in a \code{\link{run_sim_ps}} result.
#'
#' @param ps_results Output from \code{\link{run_sim_ps}}.
#'
#' @return A \code{data.table} with one row per replicate and columns:
#'   \describe{
#'     \item{i}{Replicate index.}
#'     \item{H}{Shannon entropy of the final array (bits).}
#'     \item{total}{Total number of monomers in the final array.}
#'     \item{uniqueN}{Number of distinct monomer sequences.}
#'     \item{meanN}{Mean copy number per unique sequence (\code{total / uniqueN}).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 3, max_t = 200, mu_total = 3e-5)
#' get_sim_entropies(sim)
#' }
get_sim_entropies <- function(ps_results) {
  entropies <- rbindlist(lapply(seq_along(ps_results$monomers_list), function(i) {
    mono <- ps_results$monomers_list[[i]]
    data.table(
      i       = i,
      H       = shannon_entropy(mono$bponly),
      total   = nrow(mono),
      uniqueN = uniqueN(mono$bponly)
    )
  }))
  entropies[, meanN := total / uniqueN]
  entropies
}
