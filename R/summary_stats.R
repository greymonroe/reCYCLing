# ============================================================================
# CANONICAL summary statistics for satellite repeat arrays and genomes.
# THE single source of truth: the SAME functions score simulations (reference
# tables) and real data, so a sim and a real array/genome are guaranteed
# comparable.
#
# Two levels:
#   * ARRAY level (schema v2): knob_summary_stats() on one array
#     (columns bponly, num, dir).
#   * GENOME level (schema g6): knob_genome_stats() on a multi-chromosome
#     single-haplotype genome (columns chrom, num, bponly, dir) -- the output
#     of sim_genome() or assemble_genome_dt().
#
# Design invariants:
#   * MONOMER FILTER: drop monomers whose substitution divergence from the
#     array consensus exceeds div_cutoff (mis-phased / ancient-variant
#     outliers). Applied IDENTICALLY to sim and real; consensus + divergence
#     are recomputed on the retained "main group".
#   * DISTANCE STAT: non-parametric. log10(identical-pair distance) binned
#     into fine equal-width bins over a FIXED range; any coarser binning is
#     re-aggregation, so bin tuning never forces a sim re-run.
#   * FIXED-LENGTH, NA-SAFE outputs: collapsed/empty inputs give NA (or a
#     genuine 0 where absence is informative), never variable length.
# ============================================================================

KNOB_STATS_VERSION <- "v2"
DIV_CUTOFF         <- 22                 # monomer divergence cutoff (subs/monomer)
DIST_NFINE         <- 80L                # fine distance bins stored in the reftable
DIST_LOG_MIN       <- 0                  # log10(1)
DIST_LOG_MAX       <- log10(75000)       # covers the largest observed array (~71k)
DIV_NFINE          <- 22L                # divergence histogram bins (0..DIV_CUTOFF)

# ---------------------------------------------------------------------------
# Array-level (schema v2)
# ---------------------------------------------------------------------------

#' Per-monomer substitution divergence from the array consensus
#'
#' For each monomer, the number of substitutions relative to the column-wise
#' majority consensus of the array. This is the divergence "clock" statistic:
#' with a fixed per-base mutation rate, mean divergence tracks the effective
#' age of the array's internal genealogy.
#'
#' @param bponly Character vector of gapless monomer sequences.
#' @param aligned Optional character vector of gapped (equal-width) MSA
#'   sequences in the same order (real data). When NULL (simulations), all
#'   \code{bponly} must be equal length and a fast matrix path is used.
#' @param min_occ Minimum column occupancy for a gapped column to count.
#' @param gap Gap character in \code{aligned}.
#' @return Numeric vector: substitutions from consensus per monomer.
#' @export
per_monomer_div <- function(bponly, aligned = NULL, min_occ = 0.5, gap = "-") {
  if (is.null(aligned)) {
    L <- nchar(bponly)
    if (length(unique(L)) != 1L)
      stop("per_monomer_div: gapless path needs equal-length sequences")
    mat  <- matrix(unlist(strsplit(bponly, "", fixed = TRUE), use.names = FALSE),
                   nrow = length(bponly), byrow = TRUE)
    cons <- apply(mat, 2, function(col) {
      tb <- tabulate(match(col, c("A","C","G","T")), 4L); c("A","C","G","T")[which.max(tb)] })
    consmat <- matrix(cons, nrow(mat), ncol(mat), byrow = TRUE)
    return(rowSums(mat != consmat))
  }
  mat <- matrix(unlist(strsplit(aligned, "", fixed = TRUE), use.names = FALSE),
                nrow = length(aligned), byrow = TRUE)
  occ  <- colMeans(mat != gap)
  cons <- apply(mat, 2, function(col) { b <- col[col != gap]
                 if (!length(b)) return(gap); names(sort(table(b), decreasing = TRUE))[1] })
  good <- which(occ >= min_occ & cons != gap); cm <- cons[good]
  vapply(seq_len(nrow(mat)), function(i) { s <- mat[i, good]; ng <- s != gap
           sum(s[ng] != cm[ng]) }, numeric(1))
}

#' Sampled distances between identical monomer pairs
#'
#' Uniformly samples pairs of identical monomers within each identical-sequence
#' group (bounding the O(n^2) blow-up for huge groups) and returns their
#' positional distances. The same sampler is used for simulated and real
#' arrays so the resulting distance distributions are comparable.
#'
#' @param mono data.table/data.frame with columns \code{bponly}, \code{num},
#'   \code{dir}.
#' @param max_per_group Maximum sampled pairs per identical-sequence group.
#' @param max_total Maximum total pairs returned.
#' @return data.table with columns \code{dist} (positional distance) and
#'   \code{inverted} (opposite orientations).
#' @export
sample_pair_distances <- function(mono, max_per_group = 200L, max_total = 60000L) {
  x <- as.data.table(mono)[, .(num = as.numeric(num), dir = as.character(dir), bponly)]
  x[, id := .I]
  grps <- x[, .N, by = bponly][N >= 2L]
  if (!nrow(grps)) return(data.table(dist = numeric(0), inverted = logical(0)))
  setkey(x, bponly)
  res <- vector("list", nrow(grps))
  for (gi in seq_len(nrow(grps))) {
    g <- x[.(grps$bponly[gi])]; n <- nrow(g)
    if (n * (n - 1L) / 2L <= max_per_group) {
      cmb <- utils::combn(n, 2L); i1 <- cmb[1, ]; i2 <- cmb[2, ]
    } else {
      i1 <- sample.int(n, max_per_group, replace = TRUE)
      i2 <- sample.int(n, max_per_group, replace = TRUE)
      keep <- i1 != i2; i1 <- i1[keep]; i2 <- i2[keep]
    }
    res[[gi]] <- data.table(dist = abs(g$num[i1] - g$num[i2]), inverted = g$dir[i1] != g$dir[i2])
  }
  out <- rbindlist(res); out <- out[dist > 0]
  if (nrow(out) > max_total) out <- out[sample.int(.N, max_total)]
  out[]
}

# fine log-distance histogram as proportions (internal)
dist_fine_hist <- function(dist, nfine = DIST_NFINE) {
  if (!length(dist)) return(rep(NA_real_, nfine))
  ld <- pmin(pmax(log10(dist), DIST_LOG_MIN), DIST_LOG_MAX - 1e-9)
  br <- seq(DIST_LOG_MIN, DIST_LOG_MAX, length.out = nfine + 1L)
  ct <- tabulate(findInterval(ld, br, rightmost.closed = TRUE), nbins = nfine)
  ct / sum(ct)
}

#' Summary statistics for one repeat array (schema v2)
#'
#' Computes the canonical array-level stat object: retained size, mean
#' divergence from consensus (the clock), identical-pair redundancy, and the
#' fine log-binned identical-pair distance histogram. Monomers more than
#' \code{div_cutoff} substitutions from the consensus are dropped and the
#' consensus recomputed (applied identically to sim and real data).
#'
#' @param mono data.table with columns \code{bponly}, \code{num}, \code{dir}
#'   (sim: from \code{run_sim_ps()$monomers}; real: parsed fasta).
#' @param aligned Optional gapped MSA sequences in the same order (real data).
#' @param div_cutoff Monomer divergence filter (substitutions/monomer).
#' @param nfine Number of fine distance bins.
#' @param max_per_group Pair-sampling cap per identical group.
#' @return A named list ("stat object"); see \code{knob_predictors()} to
#'   flatten it to a numeric vector.
#' @export
knob_summary_stats <- function(mono, aligned = NULL, div_cutoff = DIV_CUTOFF,
                               nfine = DIST_NFINE, max_per_group = 200L) {
  mono <- as.data.table(mono)
  n_in <- nrow(mono)
  div0 <- per_monomer_div(mono$bponly, aligned)
  keep <- is.finite(div0) & div0 <= div_cutoff
  mono2    <- mono[keep]
  aligned2 <- if (is.null(aligned)) NULL else aligned[keep]
  if (nrow(mono2) < 2L)
    return(list(version = KNOB_STATS_VERSION, ok = FALSE, n_input = n_in,
                array_size = nrow(mono2), n_pairs = 0L, frac_kept = mean(keep),
                mean_div_bp = NA_real_, frac_redundant = NA_real_, mean_group_size = NA_real_,
                dist_prop = rep(NA_real_, nfine), div_prop = rep(NA_real_, DIV_NFINE),
                div_vals = numeric(0),
                div_cutoff = div_cutoff))
  div2 <- per_monomer_div(mono2$bponly, aligned2)            # recompute on retained
  ps   <- sample_pair_distances(mono2, max_per_group = max_per_group)
  grp  <- mono2[, .N, by = bponly]$N
  list(
    version         = KNOB_STATS_VERSION, ok = TRUE,
    n_input         = n_in,
    array_size      = nrow(mono2),                            # retained = size context
    frac_kept       = mean(keep),
    n_pairs         = nrow(ps),
    mean_div_bp     = mean(div2),                             # filtered molecular clock
    frac_redundant  = sum(grp[grp >= 2L]) / nrow(mono2),
    mean_group_size = if (any(grp >= 2L)) mean(grp[grp >= 2L]) else NA_real_,
    dist_prop       = dist_fine_hist(ps$dist, nfine),         # length nfine, sums to 1
    div_prop        = { ct <- tabulate(pmin(floor(div2) + 1L, DIV_NFINE), DIV_NFINE); ct / sum(ct) },
    div_vals        = div2,                                   # per-monomer retained divergence
    div_cutoff      = div_cutoff
  )
}

#' Flatten an array stat object to a numeric predictor vector
#'
#' Coarsens the fine distance histogram to \code{K} bins (K must divide the
#' fine bin count) and appends the scalar statistics. This is the ONLY place
#' fit-level binning granularity is chosen; reference tables keep fine bins.
#'
#' @param stat Output of \code{knob_summary_stats()}.
#' @param K Number of coarse distance bins.
#' @return Named numeric vector.
#' @export
knob_predictors <- function(stat, K = 20L) {
  stopifnot(DIST_NFINE %% K == 0L)
  dp <- stat$dist_prop
  cb <- colSums(matrix(dp, nrow = DIST_NFINE %/% K))          # group adjacent fine bins
  v  <- c(setNames(cb, sprintf("d%02d", seq_len(K))),
          mean_div_bp    = stat$mean_div_bp,
          frac_redundant = stat$frac_redundant,
          log_group      = log10(pmax(stat$mean_group_size, 1)),
          log_size       = log10(pmax(stat$array_size, 1)))
  v
}

#' Pool array stat objects into one (weighted)
#'
#' Default weights are \code{n_pairs} (precision). Provided for the
#' combine-inputs alternative and diagnostics; the primary multi-array path
#' pools posteriors, not stats.
#'
#' @param stats List of stat objects from \code{knob_summary_stats()}.
#' @param weights Optional numeric weights.
#' @return A pooled stat object.
#' @export
combine_summary_stats <- function(stats, weights = NULL) {
  stats <- Filter(function(s) isTRUE(s$ok), stats)
  if (!length(stats)) stop("combine_summary_stats: no valid stat objects")
  w <- if (is.null(weights)) vapply(stats, function(s) s$n_pairs, numeric(1)) else weights
  w <- w / sum(w)
  wmean <- function(f) sum(vapply(stats, f, numeric(1)) * w)
  wvec  <- function(f) { M <- vapply(stats, f, numeric(length(stats[[1]]$dist_prop)))
                         as.numeric(M %*% w) }
  list(version = KNOB_STATS_VERSION, ok = TRUE, n_arrays = length(stats),
       array_size = wmean(function(s) s$array_size), n_pairs = sum(vapply(stats, function(s) s$n_pairs, numeric(1))),
       mean_div_bp = wmean(function(s) s$mean_div_bp), frac_redundant = wmean(function(s) s$frac_redundant),
       mean_group_size = wmean(function(s) s$mean_group_size),
       dist_prop = { dp <- wvec(function(s) s$dist_prop); dp / sum(dp) },
       div_prop  = { vp <- wvec(function(s) s$div_prop);  vp / sum(vp) },
       div_cutoff = stats[[1]]$div_cutoff)
}

#' Array stat object straight from a monomer fasta (real data)
#'
#' Parses per-monomer fasta headers by splitting on \code{sep} and reading the
#' positional index and orientation from fixed fields.
#'
#' @param fasta Path to a monomer fasta (possibly gapped/aligned).
#' @param div_cutoff,nfine Passed to \code{knob_summary_stats()}.
#' @param num_field,dir_field 1-based header fields (split on \code{sep})
#'   holding the monomer position and strand.
#' @param sep Header field separator.
#' @return A stat object.
#' @export
knob_stats_from_fasta <- function(fasta, div_cutoff = DIV_CUTOFF, nfine = DIST_NFINE,
                                  num_field = 4L, dir_field = 7L, sep = "_") {
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("need Biostrings")
  dna <- Biostrings::readDNAStringSet(fasta); raw <- as.character(dna)
  nm  <- names(dna); hd <- strsplit(nm, sep, fixed = TRUE)
  mono <- data.table(bponly = gsub("-", "", raw, fixed = TRUE),
                     num = as.numeric(vapply(hd, `[`, "", num_field)),
                     dir = vapply(hd, `[`, "", dir_field))
  knob_summary_stats(mono, aligned = raw, div_cutoff = div_cutoff, nfine = nfine)
}

# ============================================================================
# GENOME-LEVEL, CATEGORY-AWARE summary statistics (multi-chromosome).
#
# Consumes the DATA CONTRACT: a data.table `mono_dt` with one row per monomer
# and columns:
#   chrom  (int 1..K, chromosome / array id)
#   num    (within-chromosome position, 1..n_chrom)
#   bponly (gapless monomer sequence string)
#   dir    ("+"/"-" character orientation)
# The SAME function scores simulated genomes (sim_genome()) AND real
# single-hap genomes (assemble_genome_dt()). Applies the SAME div<cutoff
# filter per chromosome. Output is FIXED LENGTH and NA-safe.
# ============================================================================

GENOME_STATS_VERSION <- "g6"   # g6: + pooled retained-divergence histogram (gdiv_hist) so the
                               #     clock family carries a distribution, not just mean+IQR. ADDITIVE.
                               # g5: + OCCUPANCY stats (per-monomer copies WITHIN own chrom vs BETWEEN
                               #     other chroms, and chromosome-spread), monomer-counts robust to the
                               #     ancestral consensus. g4: + GROUP-SIZE DISTRIBUTION stats (runaway
                               #     identical families). g3: + absolute cross-chrom abundance
                               #     (frac_genome_xchrom, n_xchrom_fam). g2: zero-fill derived-sharing
                               #     block stats. All ADDITIVE.
GX_GBINS <- 7L                          # group-size histogram bins (log10, [1, 5000], monomer-weighted)
GX_OCC   <- 6L                          # occupancy copy-number histogram bins (log10, [1, 2000], monomer-weighted)
GX_SPRB  <- 6L                          # chromosome-spread histogram bins (1,2,3,4,5,>=6), monomer-weighted
GX_NXBINS  <- 20L                       # coarse pooled within-array distance bins
GX_BLKBINS <- 6L                        # shared-block-length histogram bins (log)
GX_JQ      <- c(0.5, 0.9, 0.99)         # inter-chrom Jaccard quantiles reported
GX_DIVBINS <- 11L                       # pooled divergence histogram bins (g6; 0..DIV_CUTOFF)

# ---- ancestral-IBD confound control -------------------------------------------
# In young / low-divergence simulated genomes the K chromosomes still carry many
# copies of the (near-)unmutated ancestral monomer, shared across chromosomes by
# DESCENT, not by exchange. The inter-chromosomal statistics must measure DERIVED
# sharing, so the genome-wide ancestral / major-consensus family (the consensus
# AND its GX_ANC_HAMMING-substitution cloud) is excluded before every cross-chrom
# computation. Contiguous exchanged blocks (runs of cross-chrom matches among
# DERIVED monomers) survive; scattered ancestral singletons do not.
GX_ANC_HAMMING <- 1L                    # ancestral-family substitution band (subs/monomer)
.gx_ancestral_mask <- function(bponly, thr = GX_ANC_HAMMING) {
  n <- length(bponly)
  if (!n) return(logical(0))
  m <- logical(n)
  tb <- sort(table(bponly), decreasing = TRUE)
  if (length(tb)) m <- m | (bponly == names(tb)[1L])    # most common exact monomer
  L <- nchar(bponly)
  Ltab <- sort(table(L), decreasing = TRUE)
  modeL <- as.integer(names(Ltab)[1L])
  sel <- which(L == modeL & modeL > 0L)
  if (length(sel) >= 2L) {                               # modal-length consensus band
    mat <- matrix(unlist(strsplit(bponly[sel], "", fixed = TRUE), use.names = FALSE),
                  nrow = length(sel), byrow = TRUE)
    cons <- apply(mat, 2, function(col) {
      tb2 <- tabulate(match(col, c("A","C","G","T")), 4L)
      if (!sum(tb2)) col[1L] else c("A","C","G","T")[which.max(tb2)] })
    hd <- rowSums(mat != matrix(cons, nrow(mat), ncol(mat), byrow = TRUE))
    m[sel] <- m[sel] | (hd <= thr)
  }
  m
}

#' Assemble an observed single-haplotype genome from per-chromosome fastas
#'
#' Reads per-chromosome monomer fastas for ONE taxon / ONE haplotype and
#' returns the standard genome data contract used by
#' \code{knob_genome_stats()}: one row per monomer with integer chromosome
#' rank, within-chromosome position, gapless sequence, orientation, plus the
#' gapped sequence (column \code{aligned}) for divergence computation.
#'
#' @param files Character vector of per-chromosome fasta paths.
#' @param chrom_field,num_field,dir_field 1-based header fields (split on
#'   \code{sep}) holding the chromosome id, monomer position, and strand.
#' @param sep Header field separator.
#' @return data.table with columns \code{chrom}, \code{num}, \code{bponly},
#'   \code{dir}, \code{aligned}.
#' @export
assemble_genome_dt <- function(files, chrom_field = 3L, num_field = 4L,
                               dir_field = 7L, sep = "_") {
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("need Biostrings")
  parts <- lapply(files, function(f) {
    dna <- Biostrings::readDNAStringSet(f); raw <- as.character(dna)
    nm  <- names(dna); hd <- strsplit(nm, sep, fixed = TRUE)
    chr <- vapply(hd, `[`, "", chrom_field)
    data.table(chrom_id = chr,
               num     = as.numeric(vapply(hd, `[`, "", num_field)),
               dir     = vapply(hd, `[`, "", dir_field),
               bponly  = gsub("-", "", raw, fixed = TRUE),
               aligned = raw)
  })
  dt <- rbindlist(parts)
  # integer chrom rank (gapless 1..K per contract), preserving Chr-id order
  dt[, chrom := as.integer(factor(chrom_id, levels = sort(unique(chrom_id))))]
  dt[, chrom_id := NULL]
  setorder(dt, chrom, num)
  dt[]
}

# fixed-length log-binned histogram (proportions), NA-safe (internal)
.gx_loghist <- function(x, nbins, lo, hi) {
  if (!length(x)) return(rep(NA_real_, nbins))
  lx <- pmin(pmax(log10(pmax(x, 1)), lo), hi - 1e-9)
  br <- seq(lo, hi, length.out = nbins + 1L)
  ct <- tabulate(findInterval(lx, br, rightmost.closed = TRUE), nbins = nbins)
  if (!sum(ct)) return(rep(NA_real_, nbins))
  ct / sum(ct)
}

.gx_qsafe <- function(x, probs) {
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs, names = FALSE, na.rm = TRUE))
}

# within-array blockiness: mean consecutive-identical run length (internal)
.gx_block_runlen <- function(mono2) {
  rl <- mono2[, {
    o <- order(num); s <- bponly[o]
    if (length(s) < 1L) numeric(0) else rle(s)$lengths
  }, by = chrom]$V1
  if (!length(rl)) return(NA_real_)
  mean(rl)
}

# genome-wide identical-family group-size histogram, monomer-weighted (internal)
.gx_grpsize_hist <- function(gsz, nbins = GX_GBINS, lo = 0, hi = log10(5000)) {
  if (!length(gsz)) return(rep(NA_real_, nbins))
  lx  <- pmin(pmax(log10(gsz), lo), hi - 1e-9)
  br  <- seq(lo, hi, length.out = nbins + 1L)
  bin <- findInterval(lx, br, rightmost.closed = TRUE)
  w   <- vapply(seq_len(nbins), function(b) sum(gsz[bin == b]), numeric(1))  # monomer mass per bin
  if (sum(w) == 0) return(rep(NA_real_, nbins))
  w / sum(w)
}

#' Genome-level summary statistics (schema g6)
#'
#' The full statistics catalog for a multi-chromosome single-haplotype repeat
#' genome. Computes, per chromosome and genome-wide: divergence-from-consensus
#' (clock), identical-pair distance structure, redundancy and identical-family
#' group sizes, occupancy (within- vs between-chromosome copy number), and the
#' inter-chromosomal DERIVED-sharing block (Jaccard, cross-chromosome pair
#' fraction, shared-block lengths, absolute cross-chromosome abundance).
#' The same function scores simulated (\code{sim_genome()}) and real
#' (\code{assemble_genome_dt()}) genomes.
#'
#' Use \code{knob_genome_predictors()} for the full flattened catalog and
#' \code{knob_abc_predictors()} for the reduced three-family set designed for
#' ABC fitting.
#'
#' @param mono_dt data.table with columns \code{chrom}, \code{num},
#'   \code{bponly}, \code{dir} (optionally \code{aligned}).
#' @param aligned Optional gapped MSA sequences matching \code{mono_dt} rows;
#'   if \code{mono_dt} carries an \code{aligned} column it is used.
#' @param div_cutoff Monomer divergence filter (substitutions/monomer).
#' @return A named list ("genome stat object").
#' @export
knob_genome_stats <- function(mono_dt, aligned = NULL, div_cutoff = DIV_CUTOFF) {
  mono_dt <- as.data.table(mono_dt)
  stopifnot(all(c("chrom", "num", "bponly", "dir") %in% names(mono_dt)))
  mono_dt[, chrom := as.integer(chrom)]
  mono_dt[, dir   := as.character(dir)]
  # aligned column may already live on mono_dt (assemble_genome_dt); honor arg too
  if (is.null(aligned) && "aligned" %in% names(mono_dt)) aligned <- mono_dt$aligned

  chroms <- sort(unique(mono_dt$chrom))
  n_chrom_in <- length(chroms)

  # ---- per-chromosome within-array stats via knob_summary_stats ---------------
  per <- vector("list", length(chroms))
  keep_all <- logical(nrow(mono_dt))               # genome-wide retained mask
  for (i in seq_along(chroms)) {
    cc  <- chroms[i]
    idx <- which(mono_dt$chrom == cc)
    sub <- mono_dt[idx]
    al  <- if (is.null(aligned)) NULL else aligned[idx]
    st  <- knob_summary_stats(sub, aligned = al, div_cutoff = div_cutoff)
    per[[i]] <- st
    # recompute the same div<cutoff keep mask for the genome-wide identity pass
    div0 <- per_monomer_div(sub$bponly, al)
    keep_all[idx] <- is.finite(div0) & div0 <= div_cutoff
  }
  ok_chr  <- vapply(per, function(s) isTRUE(s$ok), logical(1))
  per_ok  <- per[ok_chr]
  n_chrom <- length(per_ok)

  agg_m  <- function(f) if (n_chrom) mean(vapply(per_ok, f, numeric(1)), na.rm = TRUE) else NA_real_
  agg_iqr<- function(f) if (n_chrom >= 2) as.numeric(diff(quantile(
              vapply(per_ok, f, numeric(1)), c(.25,.75), names=FALSE, na.rm=TRUE))) else NA_real_

  mean_div_bp_m  <- agg_m  (function(s) s$mean_div_bp)
  mean_div_bp_iqr<- agg_iqr(function(s) s$mean_div_bp)
  frac_red_m     <- agg_m  (function(s) s$frac_redundant)
  frac_red_iqr   <- agg_iqr(function(s) s$frac_redundant)
  grpsz_m        <- agg_m  (function(s) s$mean_group_size)
  grpsz_iqr      <- agg_iqr(function(s) s$mean_group_size)
  array_size_m   <- agg_m  (function(s) s$array_size)

  # ---- pooled retained-divergence histogram (g6: the clock DISTRIBUTION) ------
  # Per-monomer retained divergence pooled over all chromosomes, binned over
  # [0, div_cutoff]. Carries the SHAPE of the clock signal (e.g. a saturated
  # vs age-limited genealogy), which mean+IQR alone cannot.
  div_pool <- unlist(lapply(per_ok, function(s) s$div_vals), use.names = FALSE)
  gdiv_hist <- if (length(div_pool)) {
    br <- seq(0, div_cutoff, length.out = GX_DIVBINS + 1L)
    ct <- tabulate(findInterval(pmin(div_pool, div_cutoff - 1e-9), br,
                                rightmost.closed = TRUE), nbins = GX_DIVBINS)
    ct / sum(ct)
  } else rep(NA_real_, GX_DIVBINS)

  # pooled within-array coarse distance histogram (re-aggregate fine -> GX_NXBINS)
  if (n_chrom) {
    fineM <- vapply(per_ok, function(s) {
      dp <- s$dist_prop; np <- s$n_pairs
      if (any(is.na(dp)) || !np) rep(0, DIST_NFINE) else dp * np   # weight by n_pairs
    }, numeric(DIST_NFINE))
    pooled_fine <- if (is.matrix(fineM)) rowSums(fineM) else fineM
    if (sum(pooled_fine) > 0) {
      stopifnot(DIST_NFINE %% GX_NXBINS == 0L)
      wx_dist <- colSums(matrix(pooled_fine, nrow = DIST_NFINE %/% GX_NXBINS))
      wx_dist <- wx_dist / sum(wx_dist)
    } else wx_dist <- rep(NA_real_, GX_NXBINS)
  } else wx_dist <- rep(NA_real_, GX_NXBINS)

  # ---- within-array blockiness (on the div-retained monomers) -----------------
  retained <- mono_dt[keep_all]
  blk_runlen <- if (nrow(retained)) .gx_block_runlen(retained) else NA_real_

  # ---- genome-wide group-size distribution (runaway-family detector) ----------
  grpsize_hist       <- rep(NA_real_, GX_GBINS)
  log_max_grp        <- NA_real_
  log_pairs_per_mono <- NA_real_
  if (nrow(retained)) {
    gsz <- retained[, .N, by = bponly]$N                 # genome-wide identical-family sizes
    grpsize_hist       <- .gx_grpsize_hist(gsz)
    log_max_grp        <- log10(max(gsz))
    log_pairs_per_mono <- log10(pmax(sum(gsz * (gsz - 1) / 2) / nrow(retained), 1e-6))
  }

  # ---- OCCUPANCY: per-monomer copies within own chrom vs on other chroms ------
  # Monomer-counts (not pairs), robust to the ancestral consensus. Separates
  # "smear" sims (few sequences mass-copied across chromosomes; high
  # mean_between_copies) from real genomes (many low-copy inter-chrom seqs).
  mean_within_copies <- NA_real_; mean_between_copies <- NA_real_; frac_interchrom <- NA_real_
  occ_within_hist <- rep(NA_real_, GX_OCC); occ_between_hist <- rep(NA_real_, GX_OCC)
  spread_hist <- rep(NA_real_, GX_SPRB)
  if (nrow(retained)) {
    cellR <- retained[, .(n = .N), by = .(bponly, chrom)]
    stR   <- cellR[, .(Tot = sum(n), k = uniqueN(chrom)), by = bponly]
    cellR <- merge(cellR, stR, by = "bponly"); cellR[, between := Tot - n]
    Mr <- nrow(retained)
    mean_within_copies  <- sum(cellR$n * cellR$n) / Mr                 # monomer-weighted
    mean_between_copies <- sum(cellR$between * cellR$n) / Mr
    frac_interchrom     <- sum(cellR$n[cellR$between > 0L]) / Mr
    occ_within_hist  <- .gx_loghist(rep(cellR$n, cellR$n), GX_OCC, 0, log10(2000))
    bm <- cellR$between > 0L
    occ_between_hist <- if (any(bm)) .gx_loghist(rep(cellR$between[bm], cellR$n[bm]), GX_OCC, 0, log10(2000)) else rep(0, GX_OCC)
    spr <- pmin(rep(stR$k, stR$Tot), GX_SPRB)                          # monomer-weighted spread, capped
    sh  <- tabulate(spr, nbins = GX_SPRB)
    spread_hist <- if (sum(sh)) sh / sum(sh) else rep(NA_real_, GX_SPRB)
  }

  # ============================================================================
  # INTER-CHROMOSOMAL component (DERIVED sharing; ancestral family excluded).
  # ============================================================================
  R <- retained
  anc_m <- if (nrow(R)) .gx_ancestral_mask(R$bponly) else logical(0)
  frac_ancestral_xchrom <- if (nrow(R)) mean(anc_m) else NA_real_
  if (any(anc_m)) R <- R[!anc_m]
  R[, sid := .GRP, by = bponly]                    # global distinct-seq id
  n_seq    <- max(c(R$sid, 0L))
  uchr     <- sort(unique(R$chrom)); K <- length(uchr)
  R[, cidx := match(chrom, uchr)]                  # 1..K chromosome index

  # default NA-safe outputs
  xchrom_jacc_mean <- NA_real_
  xchrom_jacc_q    <- rep(NA_real_, length(GX_JQ))
  frac_pairs_xchrom<- NA_real_
  family_spread    <- NA_real_
  frac_xchrom_inv  <- NA_real_
  blk_xchrom_hist  <- rep(NA_real_, GX_BLKBINS)
  blk_xchrom_mean  <- NA_real_
  # absolute cross-chrom ABUNDANCE: how MUCH of the genome is derived
  # cross-chrom shared, not just the ratio. NA when <2 chroms; a genuine 0 (no
  # derived sharing) is set inside the block below.
  n_xchrom_fam       <- NA_real_     # # distinct derived seqs present on >=2 chroms
  frac_genome_xchrom <- NA_real_     # frac of retained monomers that are derived cross-chrom shared

  if (K >= 2L && n_seq >= 1L && nrow(R)) {
    # (1) pairwise distinct-seq Jaccard via sparse incidence + tcrossprod --------
    if (requireNamespace("Matrix", quietly = TRUE)) {
      inc <- Matrix::sparseMatrix(i = R$cidx, j = R$sid, x = 1, dims = c(K, n_seq))
      inc <- (inc > 0) * 1
      shared <- as.matrix(Matrix::tcrossprod(inc))
      selfd  <- diag(shared)
      unionM <- outer(selfd, selfd, `+`) - shared
      jacc   <- shared / pmax(unionM, 1)
      od     <- jacc[upper.tri(jacc)]
      if (length(od)) {
        xchrom_jacc_mean <- mean(od)
        xchrom_jacc_q    <- .gx_qsafe(od, GX_JQ)
      }
    }

    # (2) identical-monomer pairs: fraction inter-chromosomal + inversion --------
    grp <- R[, .(rows = list(.I)), by = sid]
    grp <- grp[vapply(rows, length, integer(1)) >= 2L]
    if (nrow(grp)) {
      MAXG <- 200L
      tot_pairs <- 0; xchrom_pairs <- 0; xchrom_inv <- 0
      cidx_v <- R$cidx; dir_v <- R$dir
      for (gi in seq_len(nrow(grp))) {
        rs <- grp$rows[[gi]]; n <- length(rs)
        npr <- n * (n - 1L) / 2L
        if (npr <= MAXG) {
          cmb <- utils::combn(n, 2L); a <- rs[cmb[1, ]]; b <- rs[cmb[2, ]]
        } else {
          a <- rs[sample.int(n, MAXG, replace = TRUE)]
          b <- rs[sample.int(n, MAXG, replace = TRUE)]
          keep <- a != b; a <- a[keep]; b <- b[keep]
        }
        if (!length(a)) next
        xc <- cidx_v[a] != cidx_v[b]
        tot_pairs    <- tot_pairs    + length(a)
        xchrom_pairs <- xchrom_pairs + sum(xc)
        xchrom_inv   <- xchrom_inv   + sum(xc & (dir_v[a] != dir_v[b]))
      }
      if (tot_pairs > 0) frac_pairs_xchrom <- xchrom_pairs / tot_pairs
      if (xchrom_pairs > 0) frac_xchrom_inv <- xchrom_inv / xchrom_pairs
    }

    # (3) family spread: mean # chromosomes a shared (>=2-chrom) sequence is on --
    fam <- R[, .(nchr = uniqueN(cidx)), by = sid][nchr >= 2L]
    if (nrow(fam)) family_spread <- mean(fam$nchr)

    # (3b) ABSOLUTE cross-chrom abundance ----------------------------------------
    n_xchrom_fam <- nrow(fam)
    if (n_xchrom_fam > 0L && nrow(retained) > 0L) {
      n_part <- R[sid %in% fam$sid, .N]
      frac_genome_xchrom <- n_part / nrow(retained)
    } else {
      n_xchrom_fam       <- 0L     # genuine zero: no derived cross-chrom sharing
      frac_genome_xchrom <- 0
    }

    # (4) shared-block-length distribution: contiguous runs of inter-chromosomal
    #     matches walked in num order per chromosome = transferred block lengths.
    cross_sid <- fam$sid                            # sequences present on >=2 chroms
    if (length(cross_sid)) {
      R[, is_x := sid %in% cross_sid]
      runs <- R[, {
        o <- order(num); h <- is_x[o]
        if (!any(h)) numeric(0) else { r <- rle(h); as.numeric(r$lengths[r$values]) }
      }, by = chrom]$V1
      if (length(runs)) {
        blk_xchrom_mean <- mean(runs)
        blk_xchrom_hist <- .gx_loghist(runs, GX_BLKBINS, lo = 0, hi = log10(5000))
      }
      R[, is_x := NULL]
    }
  }

  list(
    version          = GENOME_STATS_VERSION,
    ok               = n_chrom >= 1L,
    n_input          = nrow(mono_dt),
    n_chrom_input    = n_chrom_in,
    n_chrom          = n_chrom,
    div_cutoff       = div_cutoff,
    # within-array aggregates
    mean_div_bp_m    = mean_div_bp_m,
    mean_div_bp_iqr  = mean_div_bp_iqr,
    frac_redundant_m = frac_red_m,
    frac_redundant_iqr = frac_red_iqr,
    mean_group_size_m  = grpsz_m,
    mean_group_size_iqr= grpsz_iqr,
    array_size_m     = array_size_m,
    gdiv_hist        = gdiv_hist,                  # length GX_DIVBINS, sums to 1 (g6)
    wx_dist          = wx_dist,                    # length GX_NXBINS, sums to 1
    blk_runlen       = blk_runlen,                 # within-array consecutive run
    # group-size distribution (runaway-family detector)
    grpsize_hist       = grpsize_hist,             # length GX_GBINS, monomer-weighted, sums to 1
    log_max_grp        = log_max_grp,              # log10 largest identical family
    log_pairs_per_mono = log_pairs_per_mono,       # log10 mean identical-pairs / monomer
    # inter-chromosomal (DERIVED sharing: ancestral/major consensus excluded)
    frac_ancestral_xchrom = frac_ancestral_xchrom, # frac of retained monomers dropped as ancestral
    xchrom_jacc_mean = xchrom_jacc_mean,
    xchrom_jacc_q    = xchrom_jacc_q,              # length length(GX_JQ)
    frac_pairs_xchrom= frac_pairs_xchrom,
    frac_xchrom_inv  = frac_xchrom_inv,
    family_spread    = family_spread,
    blk_xchrom_mean  = blk_xchrom_mean,
    blk_xchrom_hist  = blk_xchrom_hist,            # length GX_BLKBINS
    # absolute cross-chrom abundance
    n_xchrom_fam       = n_xchrom_fam,
    frac_genome_xchrom = frac_genome_xchrom,
    # occupancy: per-monomer within-chrom vs between-chrom copy number + spread
    mean_within_copies  = mean_within_copies,
    mean_between_copies = mean_between_copies,
    frac_interchrom     = frac_interchrom,
    occ_within_hist     = occ_within_hist,      # length GX_OCC, monomer-weighted
    occ_between_hist    = occ_between_hist,     # length GX_OCC, monomer-weighted (between>0)
    spread_hist         = spread_hist           # length GX_SPRB, monomer-weighted
  )
}

#' Flatten a genome stat object to the full predictor vector (the catalog)
#'
#' FIXED length, NA-safe. Log-transforms sizes/rates (heavy-tailed / spanning
#' decades). Inter-chromosomal DERIVED-sharing statistics that are undefined
#' because there is NO derived sharing (e.g. a no-exchange simulation) are a
#' genuine ZERO, not missing data, and are zero-filled (log versions floor at
#' the -12 sentinel) so the "no-sharing" region of prior space stays
#' represented in reference tables.
#'
#' @param stat Output of \code{knob_genome_stats()}.
#' @return Named numeric vector (the full statistics catalog).
#' @export
knob_genome_predictors <- function(stat) {
  na0 <- function(x) ifelse(is.finite(x), x, NA_real_)
  z0  <- function(x) ifelse(is.finite(x), x, 0)
  lg  <- function(x) log10(pmax(na0(x), 1e-12))
  lg0 <- function(x) log10(pmax(z0(x), 1e-12))   # -> -12 sentinel floor when absent
  v <- c(
    setNames(stat$wx_dist, sprintf("wxd%02d", seq_len(GX_NXBINS))),
    setNames(na0(stat$gdiv_hist), sprintf("gdv%02d", seq_len(GX_DIVBINS))),
    mean_div_bp_m     = na0(stat$mean_div_bp_m),
    mean_div_bp_iqr   = na0(stat$mean_div_bp_iqr),
    frac_redundant_m  = na0(stat$frac_redundant_m),
    frac_redundant_iqr= na0(stat$frac_redundant_iqr),
    log_group_m       = lg (stat$mean_group_size_m),
    log_group_iqr     = lg (stat$mean_group_size_iqr),
    log_array_size_m  = lg (stat$array_size_m),
    log_n_chrom       = lg (stat$n_chrom),
    blk_runlen        = na0(stat$blk_runlen),
    setNames(na0(stat$grpsize_hist), sprintf("gsz%02d", seq_len(GX_GBINS))),
    log_max_grp        = na0(stat$log_max_grp),
    log_pairs_per_mono = na0(stat$log_pairs_per_mono),
    frac_ancestral_xchrom = na0(stat$frac_ancestral_xchrom),
    log_xchrom_jacc_mean = lg0(stat$xchrom_jacc_mean),
    setNames(lg0(stat$xchrom_jacc_q), sprintf("log_xjacc_q%02d", seq_along(GX_JQ))),
    frac_pairs_xchrom = z0(stat$frac_pairs_xchrom),
    frac_xchrom_inv   = z0(stat$frac_xchrom_inv),
    family_spread     = z0(stat$family_spread),
    log_blk_xchrom_mean = lg0(stat$blk_xchrom_mean),
    setNames(z0(stat$blk_xchrom_hist), sprintf("xblk%02d", seq_len(GX_BLKBINS))),
    frac_genome_xchrom = z0(stat$frac_genome_xchrom),
    log_n_xchrom_fam   = lg0(stat$n_xchrom_fam),
    occ_mean_within     = lg (stat$mean_within_copies),
    occ_mean_between    = lg (stat$mean_between_copies),
    occ_frac_interchrom = z0(stat$frac_interchrom),
    setNames(na0(stat$occ_within_hist),  sprintf("owc%02d", seq_len(GX_OCC))),
    setNames(z0 (stat$occ_between_hist), sprintf("obc%02d", seq_len(GX_OCC))),
    setNames(na0(stat$spread_hist),      sprintf("spr%02d", seq_len(GX_SPRB)))
  )
  v
}

#' The reduced three-family ABC predictor set
#'
#' Selects, from the full catalog, the statistics that inform the ABC fit in
#' the reCYCLing model: one family per scientific question, plus size context
#' for the random forest to condition on.
#'
#' \describe{
#'   \item{Divergence (clock)}{\code{gdv01..gdv11} (pooled retained-divergence
#'     histogram), \code{mean_div_bp_m}, \code{mean_div_bp_iqr}. Pins the age
#'     parameter.}
#'   \item{Identical-pair distance}{\code{wxd01..wxd20} (pooled within-
#'     chromosome log-distance histogram). The near mode tests the local
#'     duplication unit size; the presence of a far mode tests long-range
#'     duplication.}
#'   \item{Cross-chromosome sharing}{\code{frac_pairs_xchrom},
#'     \code{frac_genome_xchrom}, \code{log_n_xchrom_fam},
#'     \code{log_xchrom_jacc_mean}. Tests inter-chromosomal exchange rate and
#'     extent.}
#'   \item{Size context}{\code{log_array_size_m}, \code{log_n_chrom}.}
#' }
#'
#' All other catalog statistics are deliberately EXCLUDED from the fit and
#' remain available as held-out posterior-predictive checks.
#'
#' @param stat A genome stat object (\code{knob_genome_stats()}) or an
#'   already-flattened named vector (\code{knob_genome_predictors()}).
#' @return Named numeric vector of the reduced predictor set.
#' @export
knob_abc_predictors <- function(stat) {
  v <- if (is.list(stat)) knob_genome_predictors(stat) else stat
  keep <- c(sprintf("gdv%02d", seq_len(GX_DIVBINS)),
            "mean_div_bp_m", "mean_div_bp_iqr",
            sprintf("wxd%02d", seq_len(GX_NXBINS)),
            "frac_pairs_xchrom", "frac_genome_xchrom",
            "log_n_xchrom_fam", "log_xchrom_jacc_mean",
            "log_array_size_m", "log_n_chrom")
  missing <- setdiff(keep, names(v))
  if (length(missing))
    stop("knob_abc_predictors: missing stats: ", paste(missing, collapse = ", "))
  v[keep]
}
