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
pairs_identical <- function(repeats) {
  x <- as.data.table(repeats)
  x[, id := .I]

  pair_idx <- x[, {
    if (.N < 2L) NULL else {
      cmb <- utils::combn(id, 2L)
      data.table(id1 = cmb[1, ], id2 = cmb[2, ])
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
