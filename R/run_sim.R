# ============================================================
# Internal per-generation simulation helpers
# ============================================================
# These functions are not exported; they are called by run_sim_ps().

# Validate that a value is a probability in [0, 1].
.check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || x < 0 || x > 1)
    stop(name, " must be a single numeric value in [0, 1]")
}

# Validate that a value is a positive integer.
.check_pos_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || x < 1 || x != round(x))
    stop(name, " must be a positive integer")
}

# Sample a chunk size from a distribution specification list.
# dist$type must be one of: "fixed", "poisson", "normal", "geom", "unif", "gamma".
.sample_chunk_size <- function(dist, max_k) {
  type <- dist$type

  s <- switch(type,
    fixed  = dist$value,
    poisson = rpois(1, lambda = dist$lambda),
    normal  = round(rnorm(1, mean = dist$mean, sd = dist$sd)),
    geom    = rgeom(1, prob = dist$prob),
    unif    = {
      if (is.null(dist$min) || is.null(dist$max))
        stop("For type = 'unif', provide dist$min and dist$max")
      round(runif(1, min = dist$min, max = dist$max))
    },
    gamma   = {
      if (is.null(dist$shape) || (is.null(dist$scale) && is.null(dist$rate)))
        stop("For type = 'gamma', provide dist$shape and dist$scale or dist$rate")
      if (!is.null(dist$scale))
        round(rgamma(1, shape = dist$shape, scale = dist$scale))
      else
        round(rgamma(1, shape = dist$shape, rate = dist$rate))
    },
    stop("Unknown chunk size dist$type: ", type)
  )

  s <- max(1L, s)
  min(s, max_k)
}

# Apply local (tandem) duplications to the array.
.local_duplication <- function(units, dirs, p_local_dup, chunk_size_dist) {
  k <- length(units)
  if (k == 0) return(list(units = units, dirs = dirs))

  triggers <- runif(k) < p_local_dup
  idxs     <- which(triggers)
  if (length(idxs) == 0L) return(list(units = units, dirs = dirs))

  for (i in rev(idxs)) {
    current_k <- length(units)
    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- .sample_chunk_size(chunk_size_dist, max_k = max_chunk)
    start <- i
    end   <- start + chunk_size - 1L
    chunk <- units[start:end]
    chunk_dirs <- dirs[start:end]
    tail_units <- if (end < current_k) units[(end + 1L):current_k] else character(0)
    tail_dirs  <- if (end < current_k) dirs[(end + 1L):current_k] else character(0)

    units <- c(units[1:end], chunk, tail_units)
    dirs  <- c(dirs[1:end], chunk_dirs, tail_dirs)
  }
  list(units = units, dirs = dirs)
}

# Apply distal (copy-paste elsewhere) duplications.
# p_distal_dup is a per-ARRAY per-generation probability (not per-unit),
# modeling chromosome-level events like unequal crossover.
.distal_duplication <- function(units, dirs, p_distal_dup, chunk_size_dist,
                                p_invert_distal = 0.5) {
  k <- length(units)
  if (k == 0 || runif(1) >= p_distal_dup) return(list(units = units, dirs = dirs))

  # One distal duplication event per generation (at most)
  i <- sample.int(k, 1L)
  max_chunk <- k - i + 1L
  chunk_size <- .sample_chunk_size(chunk_size_dist, max_k = max_chunk)
  start <- i
  end   <- min(start + chunk_size - 1L, k)
  chunk      <- units[start:end]
  chunk_dirs <- dirs[start:end]

  if (runif(1) < p_invert_distal) {
    chunk      <- rev(chunk)
    chunk_dirs <- ifelse(rev(chunk_dirs) == "+", "-", "+")
  }

  # Insert anywhere in the array
  insert_after <- sample.int(k + 1L, 1L) - 1L
  if (insert_after == 0L) {
    units <- c(chunk, units)
    dirs  <- c(chunk_dirs, dirs)
  } else if (insert_after == k) {
    units <- c(units, chunk)
    dirs  <- c(dirs, chunk_dirs)
  } else {
    units <- c(units[1:insert_after], chunk,
               units[(insert_after + 1L):k])
    dirs  <- c(dirs[1:insert_after], chunk_dirs,
               dirs[(insert_after + 1L):k])
  }
  list(units = units, dirs = dirs)
}

# Apply chunk deletions.
.delete_chunk <- function(units, dirs, p_del_chunk, del_chunk_size_dist) {
  k <- length(units)
  if (k == 0) return(list(units = units, dirs = dirs))

  triggers <- runif(k) < p_del_chunk
  idxs     <- which(triggers)
  if (length(idxs) == 0L) return(list(units = units, dirs = dirs))

  for (i in rev(idxs)) {
    current_k <- length(units)
    if (i > current_k) next
    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- .sample_chunk_size(del_chunk_size_dist, max_k = max_chunk)
    start <- i
    end   <- min(start + chunk_size - 1L, current_k)

    if (start == 1L && end == current_k) {
      units <- character(0)
      dirs  <- character(0)
    } else if (start == 1L) {
      units <- units[(end + 1L):current_k]
      dirs  <- dirs[(end + 1L):current_k]
    } else if (end == current_k) {
      units <- units[1:(start - 1L)]
      dirs  <- dirs[1:(start - 1L)]
    } else {
      units <- c(units[1:(start - 1L)], units[(end + 1L):current_k])
      dirs  <- c(dirs[1:(start - 1L)], dirs[(end + 1L):current_k])
    }
  }
  list(units = units, dirs = dirs)
}

# Apply per-base substitution / insertion / deletion mutations.
.mutate_units <- function(units, mu_total,
                          p_sub = 1, p_ins = 0, p_del_bp = 0,
                          alphabet = c("A", "C", "G", "T")) {
  if (mu_total <= 0) return(units)

  probs <- c(p_sub, p_ins, p_del_bp)
  if (any(probs < 0)) stop("p_sub, p_ins, p_del_bp must be >= 0")
  if (sum(probs) == 0) return(units)
  probs <- probs / sum(probs)

  mutate_one_unit <- function(unit) {
    bases       <- strsplit(unit, "")[[1]]
    L           <- length(bases)
    if (L == 0L) return(unit)
    mutate_mask <- runif(L) < mu_total
    if (!any(mutate_mask)) return(unit)

    # Collect positions to mutate (original positions only)
    mut_positions <- which(mutate_mask)

    # Process in reverse to avoid index shifting issues
    for (j in rev(mut_positions)) {
      if (j > length(bases)) next
      mut_type <- sample(c("sub", "ins", "del"), size = 1, prob = probs)
      if (mut_type == "sub") {
        bases[j] <- sample(setdiff(alphabet, bases[j]), 1)
      } else if (mut_type == "ins") {
        bases <- append(bases, sample(alphabet, 1), after = j)
      } else {
        if (length(bases) > 1L) {
          bases <- bases[-j]
        }
      }
    }
    paste0(bases, collapse = "")
  }

  vapply(units, mutate_one_unit, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

# Advance the array by one generation.
.step_generation <- function(units, dirs,
                             p_local_dup, local_chunk_size_dist,
                             p_distal_dup, distal_chunk_size_dist,
                             p_invert_distal = 0.5,
                             p_del_chunk, del_chunk_size_dist,
                             mu_total,
                             p_sub = 1, p_ins = 0, p_del_bp = 0,
                             alphabet = c("A", "C", "G", "T")) {

  res <- .local_duplication(units, dirs, p_local_dup, local_chunk_size_dist)
  res <- .distal_duplication(res$units, res$dirs, p_distal_dup,
                             distal_chunk_size_dist, p_invert_distal)
  res <- .delete_chunk(res$units, res$dirs, p_del_chunk, del_chunk_size_dist)
  units <- .mutate_units(res$units, mu_total, p_sub, p_ins, p_del_bp, alphabet)
  list(units = units, dirs = res$dirs)
}

# Convert a units vector to a monomer data.table.
.make_unit_dt <- function(units, dirs, hap = 1, chrom = 1, sample = "sim") {
  data.table(
    num    = seq_along(units),
    bponly = units,
    pos    = 1L,
    chrom  = chrom,
    hap    = hap,
    sample = sample,
    dir    = dirs
  )
}

# Run a single replicate (used internally by run_sim_ps).
.run_one_replicate <- function(rep_i, init_l, init_k0, init_sequence_type,
                               ancestor_seq, max_units, max_t, hard_cap,
                               p_local_dup, local_dist,
                               p_distal_dup, distal_dist, p_invert_distal,
                               p_del_chunk, del_dist, p_distal_del = 0,
                               p_conversion = 0, conv_dist = NULL, conv_tract = NULL,
                               mu_total, p_sub, p_ins, p_del_bp,
                               verbose, compute_pairs = TRUE) {
  if (p_distal_del > 0 || p_conversion > 0)
    stop("terminal deletion (p_distal_del) and gene conversion (p_conversion) are ",
         "only implemented in the C++ backend; use backend = 'cpp'.", call. = FALSE)
  units <- init_array(l = init_l, k0 = init_k0,
                      init_sequence_type = init_sequence_type,
                      ancestor_seq = ancestor_seq)
  dirs <- rep("+", length(units))

  # Pre-allocate trajectory storage with a reasonable initial capacity
  capacity <- 1024L
  l_vec <- numeric(capacity)
  H_vec <- numeric(capacity)
  N_vec <- numeric(capacity)
  idx   <- 1L

  l_vec[idx] <- length(units)
  H_vec[idx] <- shannon_entropy(units)
  N_vec[idx] <- uniqueN(units)

  t <- 1L

  while (length(units) < max_units && t < max_t) {
    t <- t + 1L
    if (verbose)
      cat("  Rep", rep_i, "- Gen", t, "- units:", length(units), "\n")

    if (length(units) > hard_cap) {
      warning("Replicate ", rep_i, " exceeded hard_cap (", hard_cap,
              ") at generation ", t, " with ", length(units), " units",
              call. = FALSE)
      break
    }

    if (length(units) == 0L) {
      warning("Replicate ", rep_i, " array went extinct at generation ", t,
              "; re-initializing and resetting clock", call. = FALSE)
      units <- init_array(l = init_l, k0 = init_k0,
                          init_sequence_type = init_sequence_type,
                          ancestor_seq = ancestor_seq)
      dirs <- rep("+", length(units))
      t <- 1L  # reset clock — new knob starts fresh
    }

    res <- .step_generation(
      units, dirs,
      p_local_dup           = p_local_dup,
      local_chunk_size_dist = local_dist,
      p_distal_dup          = p_distal_dup,
      distal_chunk_size_dist = distal_dist,
      p_invert_distal       = p_invert_distal,
      p_del_chunk           = p_del_chunk,
      del_chunk_size_dist   = del_dist,
      mu_total              = mu_total,
      p_sub                 = p_sub,
      p_ins                 = p_ins,
      p_del_bp              = p_del_bp
    )
    units <- res$units
    dirs  <- res$dirs

    idx <- idx + 1L
    if (idx > capacity) {
      capacity <- capacity * 2L
      length(l_vec) <- capacity
      length(H_vec) <- capacity
      length(N_vec) <- capacity
    }
    l_vec[idx] <- length(units)
    H_vec[idx] <- shannon_entropy(units)
    N_vec[idx] <- uniqueN(units)
  }

  # Trim to actual length
  l_vec <- l_vec[seq_len(idx)]
  H_vec <- H_vec[seq_len(idx)]
  N_vec <- N_vec[seq_len(idx)]

  dt_units <- .make_unit_dt(units, dirs)
  ps       <- if (compute_pairs) pairs_identical(dt_units) else data.table()

  list(monomers = dt_units, ps = ps,
       L_vec = l_vec, H_vec = H_vec, N_vec = N_vec)
}

# Run a single replicate using the C++ backend (deferred mutations).
.run_one_replicate_cpp <- function(rep_i, init_l, init_k0, init_sequence_type,
                                   ancestor_seq, max_units, max_t, hard_cap,
                                   p_local_dup, local_dist,
                                   p_distal_dup, distal_dist, p_invert_distal,
                                   p_del_chunk, del_dist, p_distal_del = 0,
                                   p_conversion = 0, conv_dist = NULL, conv_tract = NULL,
                                   p_autocorr_alpha = 0, autocorr_window = 50L,
                                   autocorr_every = 10L,
                                   mu_total, verbose, compute_pairs = TRUE) {
  # Generate ancestor sequence
  units <- init_array(l = init_l, k0 = 1L,
                      init_sequence_type = init_sequence_type,
                      ancestor_seq = ancestor_seq)

  # Convert ancestor to integer vector (1=A, 2=C, 3=G, 4=T)
  base_map <- c(A = 1L, C = 2L, G = 3L, T = 4L)
  ancestor_int <- base_map[strsplit(units[1], "")[[1]]]

  # Cap max_t for C++ (no Inf allowed)
  max_t_int <- if (is.infinite(max_t)) .Machine$integer.max else as.integer(max_t)

  # Call C++ engine
  res <- sim_core_cpp(
    ancestor_seq_r   = ancestor_int,
    init_k0          = as.integer(init_k0),
    max_units        = as.integer(max_units),
    max_t            = max_t_int,
    hard_cap         = as.integer(hard_cap),
    p_local_dup      = p_local_dup,
    local_dist       = local_dist,
    p_distal_dup     = p_distal_dup,
    distal_dist      = distal_dist,
    p_invert_distal  = p_invert_distal,
    p_del_chunk      = p_del_chunk,
    del_dist         = del_dist,
    p_distal_del     = p_distal_del,
    p_conversion     = p_conversion,
    conv_dist        = if (is.null(conv_dist)) list(type="gamma", shape=1, scale=20) else conv_dist,
    conv_tract       = if (is.null(conv_tract)) list(type="gamma", shape=1, scale=1) else conv_tract,
    mu_total         = mu_total,
    p_autocorr_alpha = p_autocorr_alpha,
    autocorr_window  = as.integer(autocorr_window),
    autocorr_every   = as.integer(autocorr_every),
    verbose          = verbose
  )

  if (res$hit_hard_cap) {
    warning("Replicate ", rep_i, " exceeded hard_cap (", hard_cap, ")",
            call. = FALSE)
  }

  seqs <- res$seqs
  dirs <- res$dirs
  n_final <- length(seqs)

  if (n_final == 0L) {
    warning("Replicate ", rep_i, " array went extinct", call. = FALSE)
    seqs <- init_array(l = init_l, k0 = init_k0,
                       init_sequence_type = init_sequence_type,
                       ancestor_seq = ancestor_seq)
    dirs <- rep("+", length(seqs))
    n_final <- length(seqs)
  }

  dt_units <- data.table(
    num    = seq_len(n_final),
    bponly = seqs,
    pos    = 1L,
    chrom  = 1L,
    hap    = 1L,
    sample = "sim",
    dir    = dirs
  )

  ps <- if (compute_pairs) pairs_identical(dt_units) else data.table()

  # L_vec from C++; H_vec and N_vec computed only for the final state
  l_vec <- as.numeric(res$l_vec)
  H_final <- shannon_entropy(seqs)
  N_final <- uniqueN(seqs)
  H_vec <- rep(NA_real_, length(l_vec))
  N_vec <- rep(NA_real_, length(l_vec))
  H_vec[length(H_vec)] <- H_final
  N_vec[length(N_vec)] <- N_final

  list(monomers = dt_units, ps = ps,
       L_vec = l_vec, H_vec = H_vec, N_vec = N_vec)
}

# ============================================================
# Main exported simulation function
# ============================================================

#' Run tandem repeat array evolution simulations
#'
#' Runs a single simulation of tandem repeat array evolution.
#' Starts from an array of \code{init_k0} identical monomers and advances
#' through generations by applying local duplication, distal duplication,
#' chunk deletion, and per-base mutation until an array-size or time cap is
#' reached. For multiple replicates, call this function in a loop or with
#' \code{lapply}.
#'
#' @section Chunk-size distribution specs:
#' \code{local_dist}, \code{distal_dist}, and \code{del_dist} are lists with
#' a \code{type} element and type-specific parameters:
#' \describe{
#'   \item{\code{list(type="fixed", value=10)}}{Always use this chunk size.}
#'   \item{\code{list(type="poisson", lambda=20)}}{Draw from Poisson.}
#'   \item{\code{list(type="normal", mean=20, sd=5)}}{Draw from Normal (rounded).}
#'   \item{\code{list(type="geom", prob=0.1)}}{Draw from Geometric.}
#'   \item{\code{list(type="unif", min=1, max=50)}}{Draw from Uniform (rounded).}
#'   \item{\code{list(type="gamma", shape=2, scale=15)}}{Draw from Gamma (rounded).
#'     Use \code{rate} instead of \code{scale} if preferred.}
#' }
#' Chunk sizes are coerced to an integer \eqn{\geq 1} and clipped to the
#' maximum feasible size from the sampled start position.
#'
#' @param max_units Integer. Stop when the array reaches this many
#'   units. Default 20000.
#' @param max_t Numeric. Maximum number of generations per replicate.
#'   Default \code{Inf} (no time cap).
#' @param hard_cap Integer. Safety stop: break out of the generation loop if
#'   the array exceeds this many units. Default 40000.
#' @param init_l Integer. Length (bp) of one monomer unit. Default 30.
#' @param init_k0 Integer. Number of initial copies in the array. Default 10.
#' @param init_sequence_type Character. \code{"random"} (generate a random
#'   ancestral monomer) or \code{"given"} (use \code{ancestor_seq}).
#' @param ancestor_seq Character scalar. Required when
#'   \code{init_sequence_type = "given"}; must have \code{nchar() == init_l}.
#' @param local_dist List. Chunk-size distribution for local (tandem)
#'   duplication. Default gamma(shape=2, scale=15).
#' @param distal_dist List. Chunk-size distribution for distal duplication.
#'   Default gamma(shape=2, scale=500).
#' @param del_dist List. Chunk-size distribution for chunk deletion.
#'   Default gamma(shape=2, scale=15).
#' @param p_local_dup Numeric. Per-unit per-generation probability of triggering
#'   a local duplication event. Default 0.00015.
#' @param p_distal_dup Numeric. Per-array per-generation probability of
#'   a distal duplication event (e.g., unequal crossover). At most one
#'   event per generation. Default 0.1.
#' @param p_invert_distal Numeric. Probability that a distal-duplicated chunk
#'   is reversed in unit order before insertion. Default 0.5.
#' @param p_del_chunk Numeric. Per-unit per-generation probability of
#'   triggering a (local, internal) chunk deletion of size ~ \code{del_dist}.
#'   Default 0.
#' @param p_distal_del Numeric. Per-array per-generation probability of a
#'   terminal (long-range) deletion: removes a chunk of size ~ \code{distal_dist}
#'   from the END of the array (telomere-proximal truncation; mirror of distal
#'   duplication). At most one event per generation. C++ backend only. Default 0.
#' @param mu_total Numeric. Per-base per-generation mutation probability
#'   (applied independently to every base in every unit). Default 0.0001.
#' @param p_sub Numeric. Relative probability of a substitution mutation.
#'   Default 1.
#' @param p_ins Numeric. Relative probability of a base insertion. Default 0.
#' @param p_del_bp Numeric. Relative probability of a base deletion. Default 0.
#' @param verbose Logical. Print progress messages each generation. Default
#'   \code{FALSE}.
#' @param compute_pairs Logical. If \code{TRUE} (default), compute all pairs
#'   of identical monomers via \code{\link{pairs_identical}} and return them
#'   in \code{ps}. Set to \code{FALSE} to skip this step and return an empty
#'   \code{data.table} for \code{ps}, which avoids the O(n^2) memory and time
#'   cost for large arrays. Plot functions that need pairs will compute them
#'   on the fly when passed a simulation with empty \code{ps}.
#' @param backend Character. Simulation backend to use: \code{"auto"} (default)
#'   selects the C++ backend for substitution-only simulations and falls back to
#'   R for indels; \code{"cpp"} forces the C++ backend (errors if indels are
#'   requested); \code{"r"} forces the pure-R backend. The C++ backend uses
#'   deferred mutations and integer encoding for ~100x speedup, but does not
#'   compute per-generation \code{H_vec} or \code{N_vec} (only the final
#'   values are available).
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{\code{monomers}}{A \code{data.table} describing the final array.
#'       Columns: \code{num} (position index), \code{bponly} (monomer
#'       sequence), \code{pos}, \code{chrom}, \code{hap}, \code{sample},
#'       \code{dir}.}
#'     \item{\code{ps}}{A \code{data.table} of all pairs of identical monomers
#'       (output of \code{\link{pairs_identical}}).}
#'     \item{\code{L_vec}}{Numeric vector of array copy-number (length)
#'       sampled each generation.}
#'     \item{\code{H_vec}}{Numeric vector of Shannon entropy sampled each
#'       generation (final value only with C++ backend).}
#'     \item{\code{N_vec}}{Numeric vector of unique monomer count sampled each
#'       generation (final value only with C++ backend).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(init_l = 30, init_k0 = 10,
#'                   max_t = 500, mu_total = 3e-5)
#'
#' # Inspect the monomer table
#' head(sim$monomers)
#'
#' # Plot a summary
#' plot_ps_summary(sim)
#' }
run_sim_ps <- function(
    max_units = 20000,
    max_t     = Inf,
    hard_cap  = 40000,
    init_l    = 30,
    init_k0   = 10,
    init_sequence_type = "random",
    ancestor_seq       = NULL,
    local_dist  = list(type = "gamma", shape = 2, scale = 15),
    distal_dist = list(type = "gamma", shape = 2, scale = 500),
    del_dist    = list(type = "gamma", shape = 2, scale = 15),
    conv_dist   = list(type = "gamma", shape = 1, scale = 20),  # donor separation (monomers)
    conv_tract  = list(type = "gamma", shape = 1, scale = 1),   # converted tract length (monomers)
    p_local_dup          = 0.00015,
    p_distal_dup         = 0.1,
    p_invert_distal      = 0.5,
    p_del_chunk          = 0.000,
    p_distal_del         = 0.000,
    p_conversion         = 0.000,
    p_autocorr_alpha     = 0.0,
    autocorr_window      = 50L,
    autocorr_every       = 10L,
    mu_total             = 0.0001,
    p_sub                = 1,
    p_ins                = 0,
    p_del_bp             = 0,
    verbose = FALSE,
    backend = c("auto", "cpp", "r"),
    compute_pairs = TRUE
) {
  backend <- match.arg(backend)

  # --- Input validation ---
  .check_pos_int(max_units, "max_units")
  if (!is.infinite(max_t)) .check_pos_int(max_t, "max_t")
  .check_pos_int(hard_cap, "hard_cap")
  .check_pos_int(init_l, "init_l")
  .check_pos_int(init_k0, "init_k0")
  .check_prob(p_local_dup, "p_local_dup")
  .check_prob(p_distal_dup, "p_distal_dup")
  .check_prob(p_invert_distal, "p_invert_distal")
  .check_prob(p_del_chunk, "p_del_chunk")
  .check_prob(p_distal_del, "p_distal_del")
  .check_prob(p_conversion, "p_conversion")
  .check_prob(mu_total, "mu_total")
  if (!is.numeric(p_autocorr_alpha) || length(p_autocorr_alpha) != 1L ||
      p_autocorr_alpha < 0)
    stop("p_autocorr_alpha must be a single numeric >= 0")
  .check_pos_int(autocorr_window, "autocorr_window")
  .check_pos_int(autocorr_every, "autocorr_every")
  if (hard_cap < max_units)
    warning("hard_cap (", hard_cap, ") < max_units (", max_units,
            "); simulations may never reach max_units", call. = FALSE)

  # Select backend
  has_indels <- p_ins > 0 || p_del_bp > 0
  use_cpp <- switch(backend,
    auto = !has_indels,
    cpp  = {
      if (has_indels)
        stop("C++ backend does not support indels (p_ins/p_del_bp > 0)")
      TRUE
    },
    r    = FALSE
  )

  if (use_cpp) {
    .run_one_replicate_cpp(
      rep_i = 1L,
      init_l = init_l, init_k0 = init_k0,
      init_sequence_type = init_sequence_type,
      ancestor_seq = ancestor_seq,
      max_units = max_units, max_t = max_t, hard_cap = hard_cap,
      p_local_dup = p_local_dup, local_dist = local_dist,
      p_distal_dup = p_distal_dup, distal_dist = distal_dist,
      p_invert_distal = p_invert_distal,
      p_del_chunk = p_del_chunk, del_dist = del_dist, p_distal_del = p_distal_del,
      p_conversion = p_conversion, conv_dist = conv_dist, conv_tract = conv_tract,
      p_autocorr_alpha = p_autocorr_alpha, autocorr_window = autocorr_window,
      autocorr_every = autocorr_every,
      mu_total = mu_total, verbose = verbose,
      compute_pairs = compute_pairs
    )
  } else {
    if (p_autocorr_alpha > 0)
      stop("autocatalytic local duplication (p_autocorr_alpha > 0) is only ",
           "implemented in the C++ backend; use backend = 'cpp'.", call. = FALSE)
    .run_one_replicate(
      rep_i = 1L,
      init_l = init_l, init_k0 = init_k0,
      init_sequence_type = init_sequence_type,
      ancestor_seq = ancestor_seq,
      max_units = max_units, max_t = max_t, hard_cap = hard_cap,
      p_local_dup = p_local_dup, local_dist = local_dist,
      p_distal_dup = p_distal_dup, distal_dist = distal_dist,
      p_invert_distal = p_invert_distal,
      p_del_chunk = p_del_chunk, del_dist = del_dist, p_distal_del = p_distal_del,
      p_conversion = p_conversion, conv_dist = conv_dist, conv_tract = conv_tract,
      mu_total = mu_total, p_sub = p_sub, p_ins = p_ins, p_del_bp = p_del_bp,
      verbose = verbose, compute_pairs = compute_pairs
    )
  }
}

# ============================================================
# Multi-chromosome (single-haplotype) genome simulation
# ============================================================

#' Run a multi-chromosome genome simulation with inter-chromosomal translocation
#'
#' Simulates a single haplotype as \code{K} chromosome arrays evolving in one
#' shared per-generation event loop. Each chromosome experiences the usual
#' within-array moves (local/distal duplication, conversion, deletion) and grows
#' toward its own \code{target_sizes[k]}. In addition, with per-generation
#' per-genome probability \code{p_translocation}, a contiguous chunk is copied
#' from a random donor chromosome into a random recipient chromosome (inverted
#' w.p. \code{p_invert_transloc}) -- the inter-chromosomal translocation move.
#'
#' All chromosomes are seeded from the SAME ancestral monomer, so without
#' translocation any cross-chromosome identical-sequence sharing arises only by
#' (vanishing) convergence; translocation is what creates genome-wide sharing.
#'
#' @param K Integer. Number of chromosome arrays.
#' @param target_sizes Integer vector (length K). Target unit count per chrom.
#' @param init_k0 Integer. Initial copies per chromosome. Default 10.
#' @param max_t Numeric. Max generations (Inf allowed -> capped internally).
#' @param hard_caps Integer vector (length K) or NULL. Per-chrom safety cap.
#'   Default 1.5 * target_sizes.
#' @param init_l Integer. Monomer length (bp). Default 178.
#' @param init_sequence_type,ancestor_seq Passed to \code{init_array} to make the
#'   shared ancestral monomer.
#' @param local_dist,distal_dist,del_dist,conv_dist,conv_tract Chunk-size dist
#'   specs (see \code{run_sim_ps}).
#' @param transloc_dist List. Chunk-size dist for translocation (default
#'   gamma(shape=1, scale=transloc_chunk_mean)). If supplied, overrides
#'   \code{transloc_chunk_mean}.
#' @param transloc_chunk_mean Numeric. Mean chunk size for translocation when
#'   \code{transloc_dist} is not supplied. Default 500.
#' @param p_local_dup,p_distal_dup,p_invert_distal,p_del_chunk,p_distal_del,p_conversion
#'   Within-array move probabilities (see \code{run_sim_ps}).
#' @param p_translocation Numeric. Per-genome per-generation probability of an
#'   inter-chromosomal translocation. Default 0.
#' @param p_invert_transloc Numeric. Prob a translocated chunk is inverted.
#'   Default 0.35.
#' @param mu_total Numeric. Per-base per-generation substitution prob.
#' @param verbose Logical.
#'
#' @return A \code{data.table} with one row per monomer and columns
#'   \code{chrom} (int 1..K), \code{num} (int, within-chrom position),
#'   \code{bponly} (gapless sequence string), \code{dir} ("+"/"-").
#'   Carries attributes \code{array_sizes}, \code{total_gens}, \code{n_transloc},
#'   \code{hit_hard_cap}.
#' @export
run_genome_ps <- function(
    K,
    target_sizes,
    init_k0   = 10,
    max_t     = Inf,
    n_generations = NULL,
    size_band     = 0.10,
    hard_caps = NULL,
    init_l    = 178,
    init_sequence_type = "random",
    ancestor_seq       = NULL,
    local_dist  = list(type = "gamma", shape = 2, scale = 15),
    distal_dist = list(type = "gamma", shape = 2, scale = 500),
    del_dist    = list(type = "gamma", shape = 2, scale = 15),
    conv_dist   = list(type = "gamma", shape = 1, scale = 20),
    conv_tract  = list(type = "gamma", shape = 1, scale = 1),
    transloc_dist       = NULL,
    transloc_chunk_mean = 500,
    p_local_dup         = 0.00015,
    p_distal_dup        = 0.1,
    p_invert_distal     = 0.22,
    p_del_chunk         = 0.000,
    p_distal_del        = 0.000,
    p_conversion        = 0.000,
    p_translocation     = 0.000,
    p_invert_transloc   = 0.35,
    p_autocorr_alpha    = 0.0,
    autocorr_window     = 50L,
    autocorr_every      = 10L,
    mu_total            = 5e-5,
    verbose             = FALSE
) {
  .check_pos_int(K, "K")
  if (length(target_sizes) != K)
    stop("target_sizes must have length K (", K, ")")
  target_sizes <- as.integer(round(target_sizes))
  if (any(target_sizes < 1)) stop("all target_sizes must be >= 1")
  .check_pos_int(init_k0, "init_k0")
  .check_pos_int(init_l, "init_l")
  .check_prob(p_local_dup, "p_local_dup")
  .check_prob(p_distal_dup, "p_distal_dup")
  .check_prob(p_invert_distal, "p_invert_distal")
  .check_prob(p_del_chunk, "p_del_chunk")
  .check_prob(p_distal_del, "p_distal_del")
  .check_prob(p_conversion, "p_conversion")
  .check_prob(p_translocation, "p_translocation")
  .check_prob(p_invert_transloc, "p_invert_transloc")
  .check_prob(mu_total, "mu_total")
  if (!is.numeric(p_autocorr_alpha) || length(p_autocorr_alpha) != 1L ||
      p_autocorr_alpha < 0)
    stop("p_autocorr_alpha must be a single numeric >= 0")
  .check_pos_int(autocorr_window, "autocorr_window")
  .check_pos_int(autocorr_every, "autocorr_every")

  if (is.null(hard_caps)) hard_caps <- ceiling(1.5 * target_sizes)
  hard_caps <- as.integer(round(hard_caps))
  if (length(hard_caps) != K) stop("hard_caps must have length K")

  # AGE dial (steady-state mode). If n_generations is given (>0), the genome runs
  # for exactly that many generations: chromosomes grow to their target band then
  # enter net-zero TURNOVER (fresh twins held at steady size via the per-chrom
  # governor), so divergence ~= mu*age and redundancy is maintained. If NULL/<=0,
  # the legacy "grow until all reach target then stop" behavior is used.
  if (is.null(n_generations)) {
    n_generations_int <- -1L
  } else {
    if (!is.numeric(n_generations) || length(n_generations) != 1L ||
        n_generations < 1)
      stop("n_generations must be NULL or a single integer >= 1")
    n_generations_int <- as.integer(round(n_generations))
  }
  if (!is.numeric(size_band) || length(size_band) != 1L ||
      size_band <= 0 || size_band >= 1)
    stop("size_band must be a single numeric in (0, 1)")

  if (is.null(transloc_dist))
    transloc_dist <- list(type = "gamma", shape = 1, scale = transloc_chunk_mean)

  # Shared ancestral monomer (single ancestor for the whole genome)
  anc <- init_array(l = init_l, k0 = 1L,
                    init_sequence_type = init_sequence_type,
                    ancestor_seq = ancestor_seq)
  base_map <- c(A = 1L, C = 2L, G = 3L, T = 4L)
  ancestor_int <- base_map[strsplit(anc[1], "")[[1]]]

  max_t_int <- if (is.infinite(max_t)) .Machine$integer.max else as.integer(max_t)

  res <- sim_genome_cpp(
    ancestor_seq_r    = ancestor_int,
    K                 = as.integer(K),
    target_sizes_r    = target_sizes,
    init_k0           = as.integer(init_k0),
    max_t             = max_t_int,
    hard_caps_r       = hard_caps,
    p_local_dup       = p_local_dup,
    local_dist        = local_dist,
    p_distal_dup      = p_distal_dup,
    distal_dist       = distal_dist,
    p_invert_distal   = p_invert_distal,
    p_del_chunk       = p_del_chunk,
    del_dist          = del_dist,
    p_distal_del      = p_distal_del,
    p_conversion      = p_conversion,
    conv_dist         = conv_dist,
    conv_tract        = conv_tract,
    p_translocation   = p_translocation,
    transloc_dist     = transloc_dist,
    p_invert_transloc = p_invert_transloc,
    mu_total          = mu_total,
    p_autocorr_alpha  = p_autocorr_alpha,
    autocorr_window   = as.integer(autocorr_window),
    autocorr_every    = as.integer(autocorr_every),
    n_generations     = n_generations_int,
    size_band         = size_band,
    verbose           = verbose
  )

  if (isTRUE(res$hit_hard_cap))
    warning("genome sim hit a per-chromosome hard cap", call. = FALSE)

  dt <- data.table(
    chrom  = as.integer(res$chrom),
    num    = as.integer(res$num),
    bponly = res$bponly,
    dir    = res$dir
  )
  data.table::setattr(dt, "array_sizes", as.integer(res$array_sizes))
  data.table::setattr(dt, "total_gens",  res$total_gens)
  data.table::setattr(dt, "n_transloc",  res$n_transloc)
  data.table::setattr(dt, "hit_hard_cap", res$hit_hard_cap)
  dt[]
}
