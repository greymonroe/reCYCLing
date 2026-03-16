#' Find all pairs of identical repeat units
#'
#' For a monomer data.table (such as an element of \code{monomers_list} from
#' \code{\link{run_sim_ps}}), returns every unordered pair of units that share
#' the same sequence (\code{bponly}). This pairwise table is the foundation
#' for fingerprint plots and distance analyses.
#'
#' @param repeats A \code{data.table} (or coercible object) with at least the
#'   columns \code{bponly}, \code{num}, \code{pos}, \code{chrom}, \code{hap},
#'   \code{sample}, and \code{dir}. Typically an element of
#'   \code{ps_results$monomers_list}.
#' @param max_group_size Integer. Maximum number of identical monomers in a
#'   single group before pairs are sampled instead of fully enumerated. Groups
#'   larger than this threshold contribute at most
#'   \code{max_group_size * (max_group_size - 1) / 2} pairs via random
#'   sampling. Default \code{Inf} (enumerate all pairs). Set to e.g. 500 to
#'   avoid memory issues with very large homogeneous arrays.
#'
#' @return A \code{data.table} with one row per identical pair. Columns are
#'   suffixed \code{1} and \code{2} for the two members of each pair:
#'   \code{bponly}, \code{sample1}/\code{2}, \code{hap1}/\code{2},
#'   \code{chrom1}/\code{2}, \code{num1}/\code{2}, \code{pos1}/\code{2},
#'   \code{dir1}/\code{2}. Returns an empty \code{data.table} when no
#'   duplicate sequences exist.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 300, mu_total = 3e-5)
#' mono <- sim$monomers_list[[1]]
#' ps   <- pairs_identical(mono)
#' head(ps)
#' }
pairs_identical <- function(repeats, max_group_size = Inf) {
  x <- as.data.table(repeats)
  x[, id := .I]

  pair_idx <- x[, {
    if (.N < 2L) {
      NULL
    } else if (.N <= max_group_size) {
      # Generate pairs without materializing full combn matrix
      ids <- id
      n   <- length(ids)
      id1 <- rep(ids, times = (n - 1L):0)
      id2 <- unlist(lapply(seq_len(n - 1L), function(j) ids[(j + 1L):n]),
                     use.names = FALSE)
      data.table(id1 = id1, id2 = id2)
    } else {
      # Sample pairs for very large groups to avoid memory blowup
      ids <- id
      n_sample <- as.integer(max_group_size * (max_group_size - 1L) / 2L)
      s1 <- sample(ids, n_sample, replace = TRUE)
      s2 <- sample(ids, n_sample, replace = TRUE)
      keep <- s1 < s2
      data.table(id1 = s1[keep], id2 = s2[keep])
    }
  }, by = bponly]

  if (nrow(pair_idx) == 0L) return(data.table())

  m1 <- merge(pair_idx,
              x[, .(id1 = id, bponly,
                    sample1 = sample, hap1 = hap, chrom1 = chrom,
                    num1 = num, pos1 = pos, dir1 = dir)],
              by = c("bponly", "id1"), all.x = TRUE, sort = FALSE)
  out <- merge(m1,
               x[, .(id2 = id, bponly,
                     sample2 = sample, hap2 = hap, chrom2 = chrom,
                     num2 = num, pos2 = pos, dir2 = dir)],
               by = c("bponly", "id2"), all.x = TRUE, sort = FALSE)

  out[order(chrom1, chrom2, num1, num2, pos1, pos2)]
}
