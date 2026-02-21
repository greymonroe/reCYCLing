# ============================================================
# Internal per-generation simulation helpers
# ============================================================
# These functions are not exported; they are called by run_sim_ps().

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
.local_duplication <- function(units, p_local_dup, chunk_size_dist) {
  k <- length(units)
  if (k == 0) return(units)

  triggers <- runif(k) < p_local_dup
  idxs     <- which(triggers)
  if (length(idxs) == 0L) return(units)

  for (i in rev(idxs)) {
    current_k <- length(units)
    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- .sample_chunk_size(chunk_size_dist, max_k = max_chunk)
    start <- i
    end   <- start + chunk_size - 1L
    chunk <- units[start:end]
    tail  <- if (end < current_k) units[(end + 1L):current_k] else character(0)

    units <- c(units[1:end], chunk, tail)
  }
  units
}

# Apply distal (copy-paste elsewhere) duplications.
.distal_duplication <- function(units, p_distal_dup, chunk_size_dist,
                                p_invert_distal = 0.5) {
  k <- length(units)
  if (k == 0) return(units)

  triggers <- runif(k) < p_distal_dup
  idxs     <- which(triggers)
  if (length(idxs) == 0L) return(units)

  for (i in rev(idxs)) {
    current_k <- length(units)
    if (i > current_k) next
    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size   <- .sample_chunk_size(chunk_size_dist, max_k = max_chunk)
    start        <- i
    end          <- min(start + chunk_size - 1L, current_k)
    chunk        <- units[start:end]

    if (runif(1) < p_invert_distal) chunk <- rev(chunk)

    insert_after <- sample(end:current_k, 1L)

    if (insert_after == current_k) {
      units <- c(units, chunk)
    } else {
      units <- c(units[1:insert_after], chunk,
                 units[(insert_after + 1L):current_k])
    }
  }
  units
}

# Apply chunk deletions.
.delete_chunk <- function(units, p_del_chunk, del_chunk_size_dist) {
  k <- length(units)
  if (k == 0) return(units)

  triggers <- runif(k) < p_del_chunk
  idxs     <- which(triggers)
  if (length(idxs) == 0L) return(units)

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
    } else if (start == 1L) {
      units <- units[(end + 1L):current_k]
    } else if (end == current_k) {
      units <- units[1:(start - 1L)]
    } else {
      units <- c(units[1:(start - 1L)], units[(end + 1L):current_k])
    }
  }
  units
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
    mutate_mask <- runif(L) < mu_total
    if (!any(mutate_mask)) return(unit)

    i <- 1L
    while (i <= length(bases)) {
      if (!mutate_mask[i]) { i <- i + 1L; next }
      mut_type <- sample(c("sub", "ins", "del"), size = 1, prob = probs)
      if (mut_type == "sub") {
        bases[i] <- sample(setdiff(alphabet, bases[i]), 1)
        i <- i + 1L
      } else if (mut_type == "ins") {
        bases        <- append(bases, sample(alphabet, 1), after = i)
        i <- i + 2L
      } else {
        bases        <- bases[-i]
        mutate_mask  <- mutate_mask[-i]
      }
    }
    paste0(bases, collapse = "")
  }

  vapply(units, mutate_one_unit, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

# Advance the array by one generation.
.step_generation <- function(units,
                             p_local_dup, local_chunk_size_dist,
                             p_distal_dup, distal_chunk_size_dist,
                             p_invert_distal = 0.5,
                             p_del_chunk, del_chunk_size_dist,
                             mu_total,
                             p_sub = 1, p_ins = 0, p_del_bp = 0,
                             alphabet = c("A", "C", "G", "T")) {

  units <- .local_duplication(units, p_local_dup, local_chunk_size_dist)
  units <- .distal_duplication(units, p_distal_dup, distal_chunk_size_dist,
                               p_invert_distal)
  units <- .delete_chunk(units, p_del_chunk, del_chunk_size_dist)
  units <- .mutate_units(units, mu_total, p_sub, p_ins, p_del_bp, alphabet)
  units
}

# Convert a units vector to a monomer data.table.
.make_unit_dt <- function(units, hap = 1, chrom = 1, sample = "sim") {
  data.table(
    num    = seq_along(units),
    bponly = units,
    pos    = 1L,
    chrom  = chrom,
    hap    = hap,
    sample = sample,
    dir    = "-"
  )
}

# ============================================================
# Main exported simulation function
# ============================================================

#' Run tandem repeat array evolution simulations
#'
#' Runs \code{n} independent replicates of tandem repeat array evolution.
#' Each replicate starts from an array of \code{init_k0} identical monomers
#' and advances through generations by applying local duplication, distal
#' duplication, chunk deletion, and per-base mutation until an array-size
#' or time cap is reached.
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
#' @param n Integer. Number of independent simulation replicates. Default 1000.
#' @param max_units Integer. Stop a replicate when the array reaches this many
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
#' @param p_distal_dup Numeric. Per-unit per-generation probability of
#'   triggering a distal duplication event. Default 0.0000015.
#' @param p_invert_distal Numeric. Probability that a distal-duplicated chunk
#'   is reversed in unit order before insertion. Default 0.5.
#' @param p_del_chunk Numeric. Per-unit per-generation probability of
#'   triggering a deletion event. Default 0.
#' @param mu_total Numeric. Per-base per-generation mutation probability
#'   (applied independently to every base in every unit). Default 0.0001.
#' @param p_sub Numeric. Relative probability of a substitution mutation.
#'   Default 1.
#' @param p_ins Numeric. Relative probability of a base insertion. Default 0.
#' @param p_del_bp Numeric. Relative probability of a base deletion. Default 0.
#' @param verbose Logical. Print progress messages each generation. Default
#'   \code{FALSE}.
#'
#' @return A named list with five elements, each a list of length \code{n}
#'   (one entry per replicate):
#'   \describe{
#'     \item{\code{monomers_list}}{Each element is a \code{data.table}
#'       describing the final array. Columns: \code{num} (position index),
#'       \code{bponly} (monomer sequence), \code{pos}, \code{chrom},
#'       \code{hap}, \code{sample}, \code{dir}.}
#'     \item{\code{ps_list}}{Each element is a \code{data.table} of all pairs
#'       of identical monomers (output of \code{\link{pairs_identical}}).}
#'     \item{\code{L_vec_list}}{Each element is a numeric vector of array
#'       copy-number (length) sampled each generation.}
#'     \item{\code{H_vec_list}}{Each element is a numeric vector of Shannon
#'       entropy sampled each generation.}
#'     \item{\code{N_vec_list}}{Each element is a numeric vector of unique
#'       monomer count sampled each generation.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Quick single replicate
#' sim <- run_sim_ps(n = 1, init_l = 30, init_k0 = 10,
#'                   max_t = 500, mu_total = 3e-5)
#'
#' # Inspect the monomer table for replicate 1
#' head(sim$monomers_list[[1]])
#'
#' # Plot a summary
#' plot_ps_summary(sim, i = 1)
#' }
run_sim_ps <- function(
    n = 1000,
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
    p_local_dup          = 0.00015,
    p_distal_dup         = 0.0000015,
    p_invert_distal      = 0.5,
    p_del_chunk          = 0.000,
    mu_total             = 0.0001,
    p_sub                = 1,
    p_ins                = 0,
    p_del_bp             = 0,
    verbose = FALSE
) {
  ps_list       <- vector("list", n)
  monomers_list <- vector("list", n)
  L_vec_list    <- vector("list", n)
  H_vec_list    <- vector("list", n)
  N_vec_list    <- vector("list", n)

  for (i in seq_len(n)) {
    cat("Starting replicate", i, "\n")

    units <- init_array(l = init_l, k0 = init_k0,
                        init_sequence_type = init_sequence_type,
                        ancestor_seq = ancestor_seq)

    l_vec <- length(units)
    t     <- 1L
    H_vec <- numeric(0)
    N_vec <- numeric(0)

    while (length(units) < max_units && t < max_t) {
      t <- t + 1L
      if (verbose)
        cat("  Rep", i, "- Gen", t, "- units:", length(units), "\n")
      if (length(units) > hard_cap) break

      if (length(units) == 0L)
        units <- init_array(l = init_l, k0 = init_k0,
                            init_sequence_type = init_sequence_type,
                            ancestor_seq = ancestor_seq)

      units <- .step_generation(
        units,
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

      H_vec <- c(H_vec, shannon_entropy(units))
      l_vec <- c(l_vec, length(units))
      N_vec <- c(N_vec, uniqueN(units))
    }

    L_vec_list[[i]]    <- l_vec
    N_vec_list[[i]]    <- N_vec
    H_vec_list[[i]]    <- H_vec
    dt_units           <- .make_unit_dt(units)
    monomers_list[[i]] <- dt_units
    ps_list[[i]]       <- pairs_identical(dt_units)
  }

  list(monomers_list = monomers_list,
       ps_list       = ps_list,
       L_vec_list    = L_vec_list,
       H_vec_list    = H_vec_list,
       N_vec_list    = N_vec_list)
}
