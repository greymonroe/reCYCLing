#!/usr/bin/env Rscript
# =============================================================================
# SMC-ABC parameter inference for tandem repeat array evolution
# =============================================================================
#
# Uses Sequential Monte Carlo Approximate Bayesian Computation to infer
# simulation parameters that reproduce observed summary statistics.
#
# Usage:
#   Rscript inst/scripts/smc_abc.R \
#     --target_med_near 20 --target_med_far 5000 --target_near_frac 0.5 \
#     --target_load 10 --target_array_size 20000 \
#     --particles 500 --max_generations 10 --outdir results/run1
#
# All --target_* flags control what data you're fitting to.
# All other flags control computational resources and algorithm behavior.
#
# Convergence: stops automatically when the median distance of the top 50%
# improves by less than --conv_tol (default 1%) for --conv_patience
# consecutive generations (default 3).
#
# Runs serially (single core) by default. The C++ simulation backend is fast
# enough that parallelization is not needed for most use cases, and avoids
# resource exhaustion risks.
#
# NOTE: This script mirrors inst/shiny/abc/app.R. Keep them in sync.
# =============================================================================

library(reCYCLing)
library(data.table)

# ---- Configuration with defaults -------------------------------------------

CONFIG <- list(
  # Target summary statistics (from real data)
  target = list(
    med_near          = log10(20),   # log10 of median near-distance
    med_far           = log10(5000), # log10 of median far-distance
    near_frac         = 0.5,         # fraction of pairs below threshold
    mean_load         = 10,          # mean mismatches to consensus
    target_array_size = 20000,       # expected array size (emergent)
    near_far_threshold = 500         # distance threshold separating near/far
  ),

  # Molecular clock
  ancestor_age_my   = 10,          # ancestor age in million years
  gen_time_yr       = 10,          # years per generation
  compression       = 1000,        # time compression factor
  mu_per_base_real  = 5.6e-8,      # real per-base per-gen mutation rate
  init_l            = 178,         # monomer length (bp)
  init_k0           = 10,          # initial copy number

  # SMC-ABC algorithm settings
  n_particles     = 500,
  max_generations = 10,
  retention_frac  = 0.3,
  perturbation_sd = 0.3,

  # Convergence detection
  conv_tol       = 0.01,
  conv_patience  = 3,

  # Simulation
  sim_timeout = 30,

  # Output
  outdir = "results",

  # Which parameters to fit (names from ALL_PARAMS)
  fit_params = c("p_local_dup", "p_distal_dup", "p_del_chunk")
)

# ---- All 9 possible ABC parameters (matching Shiny app) --------------------
# lo/hi define the initial prior sampling range. Perturbation is NOT clamped
# to these bounds — particles can explore beyond them freely.

ALL_PARAMS <- list(
  list(name = "p_local_dup",  lo = -5.5, hi = -3.5, default = -4.5,
       transform = function(x) 10^x, label = "Local dup rate"),
  list(name = "p_distal_dup", lo = -2.5, hi = -1.0, default = -1.3,
       transform = function(x) 10^x, label = "Distal dup rate"),
  list(name = "p_del_chunk",  lo = -6,   hi = -4,   default = -5,
       transform = function(x) 10^x, label = "Deletion rate"),
  list(name = "local_shape",  lo = 0.5,  hi = 5,    default = 2,
       transform = function(x) x,    label = "Local chunk shape"),
  list(name = "local_scale",  lo = 1,    hi = 50,   default = 15,
       transform = function(x) x,    label = "Local chunk scale"),
  list(name = "distal_shape", lo = 0.5,  hi = 5,    default = 2,
       transform = function(x) x,    label = "Distal chunk shape"),
  list(name = "distal_scale", lo = 10,   hi = 5000, default = 500,
       transform = function(x) x,    label = "Distal chunk scale"),
  list(name = "del_shape",    lo = 0.5,  hi = 5,    default = 2,
       transform = function(x) x,    label = "Del chunk shape"),
  list(name = "del_scale",    lo = 1,    hi = 50,   default = 15,
       transform = function(x) x,    label = "Del chunk scale")
)
ALL_PARAM_NAMES <- sapply(ALL_PARAMS, `[[`, "name")

# Fixed values for non-fitted params (populated from defaults, overridable via CLI)
FIXED_VALS <- setNames(
  lapply(ALL_PARAMS, `[[`, "default"),
  ALL_PARAM_NAMES
)

# ---- Parse command-line args ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  v <- if (i < length(args)) args[i + 1] else NULL

  if      (a == "--outdir")          { CONFIG$outdir <- v }
  else if (a == "--particles")       { CONFIG$n_particles <- as.integer(v) }
  else if (a == "--max_generations") { CONFIG$max_generations <- as.integer(v) }
  else if (a == "--sim_timeout")     { CONFIG$sim_timeout <- as.numeric(v) }
  else if (a == "--init_l")          { CONFIG$init_l <- as.integer(v) }
  else if (a == "--init_k0")         { CONFIG$init_k0 <- as.integer(v) }
  else if (a == "--retention_frac")  { CONFIG$retention_frac <- as.numeric(v) }
  else if (a == "--perturbation_sd") { CONFIG$perturbation_sd <- as.numeric(v) }
  else if (a == "--conv_tol")        { CONFIG$conv_tol <- as.numeric(v) }
  else if (a == "--conv_patience")   { CONFIG$conv_patience <- as.integer(v) }
  else if (a == "--ancestor_age_my") { CONFIG$ancestor_age_my <- as.numeric(v) }
  else if (a == "--gen_time_yr")     { CONFIG$gen_time_yr <- as.numeric(v) }
  else if (a == "--compression")     { CONFIG$compression <- as.numeric(v) }
  else if (a == "--mu_per_base_real") { CONFIG$mu_per_base_real <- as.numeric(v) }
  # Target statistics
  else if (a == "--target_med_near")    { CONFIG$target$med_near <- log10(as.numeric(v)) }
  else if (a == "--target_med_far")     { CONFIG$target$med_far <- log10(as.numeric(v)) }
  else if (a == "--target_near_frac")   { CONFIG$target$near_frac <- as.numeric(v) }
  else if (a == "--target_load")        { CONFIG$target$mean_load <- as.numeric(v) }
  else if (a == "--target_array_size")  { CONFIG$target$target_array_size <- as.integer(v) }
  else if (a == "--near_far_threshold") { CONFIG$target$near_far_threshold <- as.numeric(v) }
  # Which params to fit (comma-separated)
  else if (a == "--fit_params") {
    CONFIG$fit_params <- strsplit(v, ",")[[1]]
    for (fp in CONFIG$fit_params) {
      if (!fp %in% ALL_PARAM_NAMES)
        stop("Unknown parameter in --fit_params: ", fp,
             "\nValid: ", paste(ALL_PARAM_NAMES, collapse = ", "))
    }
  }
  # Fixed parameter values (for non-fitted params)
  else if (a == "--fixed_p_local_dup")  { FIXED_VALS$p_local_dup  <- as.numeric(v) }
  else if (a == "--fixed_p_distal_dup") { FIXED_VALS$p_distal_dup <- as.numeric(v) }
  else if (a == "--fixed_p_del_chunk")  { FIXED_VALS$p_del_chunk  <- as.numeric(v) }
  else if (a == "--fixed_local_shape")  { FIXED_VALS$local_shape  <- as.numeric(v) }
  else if (a == "--fixed_local_scale")  { FIXED_VALS$local_scale  <- as.numeric(v) }
  else if (a == "--fixed_distal_shape") { FIXED_VALS$distal_shape <- as.numeric(v) }
  else if (a == "--fixed_distal_scale") { FIXED_VALS$distal_scale <- as.numeric(v) }
  else if (a == "--fixed_del_shape")    { FIXED_VALS$del_shape    <- as.numeric(v) }
  else if (a == "--fixed_del_scale")    { FIXED_VALS$del_scale    <- as.numeric(v) }
  else if (a == "--help") {
    cat("
SMC-ABC parameter inference for reCYCLing tandem repeat simulations.

Target statistics (what your real data looks like):
  --target_med_near NUM      Median near-distance in units (default: 20)
  --target_med_far NUM       Median far-distance in units (default: 5000)
  --target_near_frac NUM     Fraction of pairs below threshold (default: 0.5)
  --target_load NUM          Mean mismatches to consensus (default: 10)
  --target_array_size NUM    Target array size in monomers (default: 20000)
  --near_far_threshold NUM   Distance threshold separating near/far (default: 500)

Molecular clock:
  --ancestor_age_my NUM      Ancestor age in million years (default: 10)
  --gen_time_yr NUM          Years per generation (default: 10)
  --compression NUM          Time compression factor (default: 1000)
  --mu_per_base_real NUM     Per-base per-gen mutation rate (default: 5.6e-8)

Fitted parameters:
  --fit_params LIST          Comma-separated param names to fit
                             (default: p_local_dup,p_distal_dup,p_del_chunk)
                             Valid: p_local_dup, p_distal_dup, p_del_chunk,
                                    local_shape, local_scale, distal_shape,
                                    distal_scale, del_shape, del_scale

Fixed parameter values (for non-fitted params, log10 scale for rates):
  --fixed_p_local_dup NUM    (default: -4.5)
  --fixed_p_distal_dup NUM   (default: -1.3)
  --fixed_p_del_chunk NUM    (default: -5)
  --fixed_local_shape NUM    (default: 2)
  --fixed_local_scale NUM    (default: 15)
  --fixed_distal_shape NUM   (default: 2)
  --fixed_distal_scale NUM   (default: 500)
  --fixed_del_shape NUM      (default: 2)
  --fixed_del_scale NUM      (default: 15)

Algorithm settings:
  --particles NUM            Particles per generation (default: 500)
  --max_generations NUM      Max ABC iterations (default: 10)
  --retention_frac NUM       Fraction kept each gen (default: 0.3)
  --perturbation_sd NUM      Initial exploration width (default: 0.3)
  --conv_tol NUM             Convergence tolerance (default: 0.01)
  --conv_patience NUM        Stall gens before stopping (default: 3)

Simulation settings:
  --sim_timeout NUM          Seconds before killing a slow sim (default: 30)
  --init_l NUM               Monomer length in bp (default: 178)
  --init_k0 NUM              Initial copy number (default: 10)

Output:
  --outdir PATH              Output directory (default: results)
")
    quit(save = "no", status = 0)
  }
  else { stop("Unknown argument: ", a, "\nUse --help for usage.") }

  i <- i + 2
}

dir.create(CONFIG$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Build active PARAMS list from fit_params config -----------------------
# PARAMS = only the parameters being fitted. Others use FIXED_VALS.

PARAMS <- list()
for (pname in CONFIG$fit_params) {
  idx <- which(ALL_PARAM_NAMES == pname)
  PARAMS <- c(PARAMS, list(ALL_PARAMS[[idx]]))
}
param_names <- sapply(PARAMS, `[[`, "name")

cat("Fitting parameters:", paste(param_names, collapse = ", "), "\n")
cat("Fixed values:\n")
for (pname in ALL_PARAM_NAMES) {
  if (!pname %in% param_names) {
    cat(sprintf("  %-15s = %s\n", pname, FIXED_VALS[[pname]]))
  }
}

# ---- Helper functions (matching Shiny app) ---------------------------------

sample_prior_one <- function() {
  row <- as.list(sapply(PARAMS, function(p) runif(1, p$lo, p$hi)))
  names(row) <- param_names
  as.data.table(row)
}

# Build simulation args from a particle row + fixed values.
# Matches Shiny's params_to_args exactly.
params_to_args <- function(row, mu_total) {
  get_val <- function(pname) {
    if (pname %in% names(row)) row[[pname]]
    else if (pname %in% names(FIXED_VALS)) FIXED_VALS[[pname]]
    else {
      idx <- which(ALL_PARAM_NAMES == pname)
      ALL_PARAMS[[idx]]$default
    }
  }
  list(
    p_local_dup  = 10^get_val("p_local_dup"),
    p_distal_dup = 10^get_val("p_distal_dup"),
    p_del_chunk  = 10^get_val("p_del_chunk"),
    mu_total     = mu_total,
    local_dist   = list(type = "gamma", shape = get_val("local_shape"),
                        scale = get_val("local_scale")),
    distal_dist  = list(type = "gamma", shape = get_val("distal_shape"),
                        scale = get_val("distal_scale")),
    del_dist     = list(type = "gamma", shape = get_val("del_shape"),
                        scale = get_val("del_scale"))
  )
}

fit_2gauss <- function(x, max_iter = 100, max_n = 10000) {
  if (length(x) < 10) return(NULL)
  if (length(x) > max_n) x <- x[sample.int(length(x), max_n)]
  k <- tryCatch(kmeans(x, 2, nstart = 5), error = function(e) NULL)
  if (is.null(k)) return(NULL)
  mu1 <- min(k$centers); mu2 <- max(k$centers)
  cl1 <- which.min(k$centers)
  sd1 <- sd(x[k$cluster == cl1]); sd2 <- sd(x[k$cluster != cl1])
  if (is.na(sd1) || sd1 == 0) sd1 <- 0.1
  if (is.na(sd2) || sd2 == 0) sd2 <- 0.1
  w1 <- mean(k$cluster == cl1)
  for (iter in seq_len(max_iter)) {
    d1 <- w1 * dnorm(x, mu1, sd1)
    d2 <- (1 - w1) * dnorm(x, mu2, sd2)
    gamma1 <- d1 / (d1 + d2 + 1e-300)
    w1_new <- mean(gamma1)
    mu1_new <- sum(gamma1 * x) / sum(gamma1)
    mu2_new <- sum((1 - gamma1) * x) / sum(1 - gamma1)
    sd1_new <- sqrt(sum(gamma1 * (x - mu1_new)^2) / sum(gamma1))
    sd2_new <- sqrt(sum((1 - gamma1) * (x - mu2_new)^2) / sum(1 - gamma1))
    converged <- abs(mu1_new - mu1) < 1e-6 && abs(mu2_new - mu2) < 1e-6
    w1 <- w1_new; mu1 <- mu1_new; mu2 <- mu2_new
    sd1 <- max(sd1_new, 0.01); sd2 <- max(sd2_new, 0.01)
    if (converged) break
  }
  if (mu1 > mu2) list(log_mode1 = mu2, log_mode2 = mu1, weight1 = 1 - w1)
  else list(log_mode1 = mu1, log_mode2 = mu2, weight1 = w1)
}

#' Run a single simulation and extract summary statistics.
#' Uses setTimeLimit for wall-clock timeout (no forking).
#' Matches Shiny app's run_and_summarize_one logic.
run_and_summarize <- function(row, config) {
  mu_compressed <- config$mu_per_base_real * config$compression
  sim_gens <- round(config$ancestor_age_my * 1e6 / config$gen_time_yr / config$compression)
  sim_args <- params_to_args(row, mu_compressed)
  target_size <- config$target$target_array_size
  hard_cap <- as.integer(target_size * 3)

  sim <- tryCatch(
    {
      setTimeLimit(elapsed = config$sim_timeout, transient = TRUE)
      suppressWarnings(run_sim_ps(
        max_units = hard_cap + 1L, max_t = sim_gens,
        hard_cap = hard_cap,
        init_l = config$init_l, init_k0 = config$init_k0,
        p_local_dup  = sim_args$p_local_dup,
        p_distal_dup = sim_args$p_distal_dup,
        p_del_chunk  = sim_args$p_del_chunk,
        mu_total     = sim_args$mu_total,
        local_dist   = sim_args$local_dist,
        distal_dist  = sim_args$distal_dist,
        del_dist     = sim_args$del_dist,
        verbose      = FALSE,
        compute_pairs = FALSE
      ))
    },
    error = function(e) NULL,
    finally = setTimeLimit(elapsed = Inf, transient = TRUE)
  )
  if (is.null(sim)) return(NULL)

  mono <- sim$monomers

  # Reject extinct or exploded arrays
  if (nrow(mono) < 50 || nrow(mono) > hard_cap) return(NULL)

  # --- Full summary statistics ---
  # Compute pairwise distances from sampled identical pairs directly,
  # bypassing the full pairs table (which can be millions of rows at 20K units).
  dists <- mono[, {
    if (.N < 2L) NULL
    else {
      positions <- num
      n <- length(positions)
      if (n <= 200L) {
        id1 <- rep(seq_len(n), times = (n - 1L):0)
        id2 <- unlist(lapply(seq_len(n - 1L), function(j) (j + 1L):n),
                       use.names = FALSE)
        list(d = abs(positions[id1] - positions[id2]))
      } else {
        n_sample <- min(5000L, as.numeric(n) * (n - 1L) %/% 2)
        s1 <- sample.int(n, n_sample, replace = TRUE)
        s2 <- sample.int(n, n_sample, replace = TRUE)
        keep <- s1 != s2
        list(d = abs(positions[s1[keep]] - positions[s2[keep]]))
      }
    }
  }, by = bponly]$d
  dists <- dists[dists > 0]
  if (length(dists) < 20) return(NULL)

  # Quantile-based summary stats
  threshold <- config$target$near_far_threshold
  near <- dists[dists <= threshold]
  far  <- dists[dists > threshold]
  near_frac <- length(near) / length(dists)
  med_near  <- if (length(near) >= 5) median(log10(near)) else NA_real_
  med_far   <- if (length(far) >= 5) median(log10(far)) else NA_real_
  if (is.na(med_near) || is.na(med_far)) return(NULL)

  # Mutation load — quick subsample gate first
  sample_n <- min(200, nrow(mono))
  sample_seqs <- mono$bponly[sample.int(nrow(mono), sample_n)]
  ct_quick <- tryCatch(counts_long_nogap(sample_seqs), error = function(e) NULL)
  if (is.null(ct_quick)) return(NULL)
  cons_quick <- ct_quick[symbol == consensus][order(pos)]
  load_quick <- mu_load_from_consensus(sample_seqs, cons_quick)
  mean_load_est <- mean(load_quick)

  # Reject if mutation load wildly off target
  if (!is.na(config$target$mean_load) && config$target$mean_load > 0) {
    if (mean_load_est > config$target$mean_load * 3 ||
        mean_load_est < config$target$mean_load / 3) return(NULL)
  }

  # Full mutation load on all monomers
  ct <- tryCatch(counts_long_nogap(mono$bponly), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  cons <- ct[symbol == consensus][order(pos)]
  load <- mu_load_from_consensus(mono$bponly, cons)

  data.table(
    med_near  = med_near,
    med_far   = med_far,
    near_frac = near_frac,
    mean_load = mean(load),
    n_units   = nrow(mono),
    n_pairs   = length(dists)
  )
}

compute_distance <- function(sim_stats, target, stat_scales = NULL) {
  diffs <- c(
    med_near   = sim_stats$med_near  - target$med_near,
    med_far    = sim_stats$med_far   - target$med_far,
    near_frac  = sim_stats$near_frac - target$near_frac,
    mean_load  = (sim_stats$mean_load - target$mean_load) / max(target$mean_load, 0.1),
    array_size = (log10(sim_stats$n_units) - log10(target$target_array_size))
  )
  if (!is.null(stat_scales)) {
    diffs <- diffs / stat_scales
  }
  # Blended: Euclidean + max deviation (prevents cherry-picking one stat)
  euclidean <- sqrt(mean(diffs^2))
  worst     <- max(abs(diffs))
  0.5 * euclidean + 0.5 * worst
}

compute_stat_scales <- function(particles, target) {
  rel_load <- (particles$mean_load - target$mean_load) / max(target$mean_load, 0.1)
  rel_size <- log10(particles$n_units) - log10(target$target_array_size)
  scales <- c(
    med_near   = sd(particles$med_near, na.rm = TRUE),
    med_far    = sd(particles$med_far, na.rm = TRUE),
    near_frac  = sd(particles$near_frac, na.rm = TRUE),
    mean_load  = sd(rel_load, na.rm = TRUE),
    array_size = sd(rel_size, na.rm = TRUE)
  )
  scales[is.na(scales) | scales < 1e-10] <- 1
  scales
}

perturb_particle <- function(particle, param_sds) {
  new <- copy(particle)
  for (j in seq_along(PARAMS)) {
    pname <- PARAMS[[j]]$name
    new[[pname]] <- new[[pname]] + rnorm(1, 0, param_sds[[pname]])
  }
  # Enforce physical constraints: shape params > 0, scale params > 0
  for (pname in c("local_shape", "distal_shape", "del_shape"))
    if (pname %in% names(new)) new[[pname]] <- max(0.01, new[[pname]])
  for (pname in c("local_scale", "distal_scale", "del_scale"))
    if (pname %in% names(new)) new[[pname]] <- max(0.1, new[[pname]])
  new
}

#' Compute adaptive perturbation SDs from the accepted particle population.
#' Uses 2 * empirical SD per parameter (Beaumont 2009 optimal kernel).
#' Falls back to prior-range-based SD for Gen 0.
compute_param_sds <- function(kept, fallback_scale = 0.3) {
  sds <- list()
  for (j in seq_along(PARAMS)) {
    pname <- PARAMS[[j]]$name
    prior_range <- PARAMS[[j]]$hi - PARAMS[[j]]$lo
    emp_sd <- if (!is.null(kept) && nrow(kept) >= 3)
      sd(kept[[pname]], na.rm = TRUE) else NA_real_
    if (is.na(emp_sd) || emp_sd < 1e-10) {
      sds[[pname]] <- fallback_scale * prior_range
    } else {
      # Adaptive SD, but capped at the prior range to prevent runaway widening
      sds[[pname]] <- min(2 * emp_sd, prior_range * 0.15)
    }
  }
  sds
}

# ---- Fill N valid particle slots (retrying failures) -------------------------

fill_particle_slots <- function(n_slots, propose_fn, config, stat_scales = NULL) {
  particles <- NULL
  n_filled <- 0L
  n_failed <- 0L
  n_timeout <- 0L
  max_attempts <- n_slots * 5L  # safety valve

  for (attempt in seq_len(max_attempts)) {
    if (n_filled >= n_slots) break

    proposed <- propose_fn()
    t0 <- proc.time()[3]
    stats <- run_and_summarize(proposed, config)
    elapsed <- proc.time()[3] - t0

    if (is.null(stats)) {
      if (elapsed >= config$sim_timeout - 0.5) {
        n_timeout <- n_timeout + 1L
      } else {
        n_failed <- n_failed + 1L
      }
      next
    }

    dist <- compute_distance(stats, config$target, stat_scales)
    row <- cbind(proposed, stats, data.table(distance = dist, elapsed = elapsed))

    if (is.null(particles)) {
      particles <- row
    } else {
      particles <- rbind(particles, row, fill = TRUE)
    }
    n_filled <- n_filled + 1L

    if (n_filled %% 10 == 0) {
      cat(sprintf("  filled %d/%d (attempts: %d, failed: %d, timeout: %d)\n",
                  n_filled, n_slots, attempt, n_failed, n_timeout))
    }
  }

  if (n_filled < n_slots) {
    cat(sprintf("  WARNING: only filled %d/%d after %d attempts\n",
                n_filled, n_slots, max_attempts))
  }

  list(particles = particles, n_failed = n_failed, n_timeout = n_timeout)
}

# ---- SMC-ABC main loop ------------------------------------------------------

run_smc_abc <- function(config) {
  cat("=== SMC-ABC Parameter Inference ===\n")
  cat("Particles:", config$n_particles, "\n")
  cat("Max generations:", config$max_generations, "\n")
  cat("Convergence: stop after", config$conv_patience, "gens with <",
      config$conv_tol * 100, "% improvement\n")
  cat("Timeout:", config$sim_timeout, "s per particle\n")
  mu_compressed <- config$mu_per_base_real * config$compression
  sim_gens <- round(config$ancestor_age_my * 1e6 / config$gen_time_yr / config$compression)
  target_size <- config$target$target_array_size
  hard_cap <- as.integer(target_size * 3)
  cat("Molecular clock: ancestor=", config$ancestor_age_my, "My, gen_time=",
      config$gen_time_yr, "yr, compression=", config$compression, "\n")
  cat("Effective mu:", signif(mu_compressed, 3), "| Sim generations:", sim_gens, "\n")
  cat("Hard cap:", hard_cap, "| Target size:", target_size, "\n")
  cat("Targets: med_near=", 10^config$target$med_near,
      " med_far=", 10^config$target$med_far,
      " near_frac=", config$target$near_frac,
      " load=", config$target$mean_load,
      " array_size=", target_size, "\n")
  cat("Output:", config$outdir, "\n\n")

  # Track convergence
  convergence <- data.table(gen = integer(), best = numeric(),
                            median = numeric(), n_failed = integer(),
                            n_timeout = integer())
  stall_count <- 0L

  # -- Generation 0: sample from prior --
  cat("--- Generation 0 (prior) ---\n")
  t0 <- proc.time()
  result <- fill_particle_slots(
    config$n_particles,
    propose_fn = sample_prior_one,
    config = config
  )
  particles <- result$particles
  elapsed <- (proc.time() - t0)[3]

  if (is.null(particles) || nrow(particles) == 0) {
    cat("No valid particles found in generation 0. Stopping.\n")
    return(invisible(NULL))
  }

  best_dist <- min(particles$distance)
  med_dist  <- median(particles$distance)
  convergence <- rbind(convergence,
    data.table(gen = 0L, best = best_dist, median = med_dist,
               n_failed = result$n_failed, n_timeout = result$n_timeout))

  cat(sprintf("  Done: %d particles | Time: %.0fs | Best: %.3f | Median: %.3f\n",
              nrow(particles), elapsed, best_dist, med_dist))
  cat(sprintf("  Failures: %d failed, %d timed out\n",
              result$n_failed, result$n_timeout))

  # Compute stat scales from Gen 0 pilot batch and recompute distances
  stat_scales <- compute_stat_scales(particles, config$target)
  cat("  Stat scales (SD from pilot):",
      paste(names(stat_scales), "=", round(stat_scales, 4), collapse = ", "), "\n")
  for (ri in seq_len(nrow(particles))) {
    particles[ri, distance := compute_distance(
      .SD, config$target, stat_scales),
      .SDcols = c("med_near", "med_far", "near_frac", "mean_load", "n_units")]
  }
  best_dist <- min(particles$distance)
  med_dist  <- median(particles$distance)
  cat(sprintf("  Rescaled: Best: %.3f | Median: %.3f\n", best_dist, med_dist))

  fwrite(particles, file.path(config$outdir, "gen_0_particles.csv"))

  # -- Subsequent generations --
  for (gen in seq_len(config$max_generations)) {
    cat(sprintf("\n--- Generation %d ---\n", gen))

    # Elitist selection: pool this generation with previous survivors,
    # keep the best overall. Good particles survive across generations.
    pool <- particles
    if (exists("kept") && !is.null(kept) && nrow(kept) > 0) {
      pool <- rbind(pool, kept, fill = TRUE)
      pool <- unique(pool, by = intersect(param_names, names(pool)))
    }
    n_keep <- max(3L, as.integer(nrow(particles) * config$retention_frac))
    pool <- pool[order(distance)]
    kept <- pool[seq_len(min(n_keep, nrow(pool)))]

    cat(sprintf("  Kept %d particles (threshold: %.3f)\n",
                nrow(kept), max(kept$distance)))

    # Show current best
    best <- kept[1]
    cat(sprintf("  Best: med_near=%.0f med_far=%.0f frac=%.2f load=%.1f size=%d (dist=%.3f)\n",
                10^best$med_near, 10^best$med_far,
                best$near_frac, best$mean_load, best$n_units, best$distance))

    # Adaptive perturbation: SD per parameter = 2 * empirical SD of survivors
    param_sds <- compute_param_sds(kept, config$perturbation_sd)

    # Fill all slots by perturbing random parents, retrying failures
    t0 <- proc.time()
    result <- fill_particle_slots(
      config$n_particles,
      propose_fn = function() {
        parent <- kept[sample.int(nrow(kept), 1), ..param_names]
        perturb_particle(parent, param_sds)
      },
      config = config,
      stat_scales = stat_scales
    )
    particles <- result$particles
    elapsed <- (proc.time() - t0)[3]

    if (is.null(particles) || nrow(particles) == 0) {
      cat("  No valid particles. Stopping.\n")
      break
    }

    best_dist <- min(particles$distance)
    med_dist  <- median(particles$distance)
    convergence <- rbind(convergence,
      data.table(gen = gen, best = best_dist, median = med_dist,
                 n_failed = result$n_failed, n_timeout = result$n_timeout))

    cat(sprintf("  Done: %d particles | Time: %.0fs | Best: %.3f | Median: %.3f\n",
                nrow(particles), elapsed, best_dist, med_dist))
    cat(sprintf("  Failures: %d failed, %d timed out\n",
                result$n_failed, result$n_timeout))

    fwrite(particles, file.path(config$outdir,
                                sprintf("gen_%d_particles.csv", gen)))

    # -- Convergence check --
    if (nrow(convergence) >= 2) {
      prev_med <- convergence$median[nrow(convergence) - 1]
      improvement <- if (is.finite(prev_med) && prev_med > 0)
        (prev_med - med_dist) / prev_med else 1
      cat(sprintf("  Improvement: %.1f%%\n", improvement * 100))

      if (improvement < config$conv_tol) {
        stall_count <- stall_count + 1L
        cat(sprintf("  Stall %d/%d\n", stall_count, config$conv_patience))
        if (stall_count >= config$conv_patience) {
          cat("  Converged! Stopping.\n")
          break
        }
      } else {
        stall_count <- 0L
      }
    }
  }

  # ---- Final summary --------------------------------------------------------
  cat("\n=== Final Results ===\n")
  cat("Generations completed:", nrow(convergence), "\n")
  cat("Median distance trajectory:",
      paste(round(convergence$median, 3), collapse = " -> "), "\n")
  cat("Failure trajectory:",
      paste(convergence$n_failed + convergence$n_timeout, collapse = " -> "), "\n")

  if (is.null(particles) || nrow(particles) == 0) {
    cat("No valid particles found.\n")
    return(invisible(NULL))
  }

  particles <- particles[order(distance)]
  top <- particles[seq_len(min(20, nrow(particles)))]

  cat("\nTop 20 particles:\n")
  top_display <- copy(top)
  for (j in seq_along(PARAMS)) {
    pn <- PARAMS[[j]]$name
    top_display[[pn]] <- sapply(top_display[[pn]], PARAMS[[j]]$transform)
  }
  top_display[, near := 10^med_near]
  top_display[, far := 10^med_far]
  display_cols <- c(param_names, "near", "far", "near_frac", "mean_load",
                    "n_units", "distance")
  print(top_display[, ..display_cols], digits = 3)

  cat("\nPosterior summary (top 50%):\n")
  n_post <- max(5, as.integer(nrow(particles) * 0.5))
  posterior <- particles[seq_len(n_post)]
  for (pn in param_names) {
    p <- PARAMS[[which(param_names == pn)]]
    raw_vals <- sapply(posterior[[pn]], p$transform)
    cat(sprintf("  %-15s  mean=%-10.4g  sd=%-10.4g  [%-10.4g, %-10.4g]\n",
                pn, mean(raw_vals), sd(raw_vals),
                min(raw_vals), max(raw_vals)))
  }

  fwrite(particles, file.path(config$outdir, "final_accepted.csv"))
  fwrite(convergence, file.path(config$outdir, "convergence.csv"))

  cat("\nFiles saved to:", config$outdir, "\n")
  cat("  gen_*_particles.csv  - all particles per generation\n")
  cat("  final_accepted.csv   - accepted particles, sorted by distance\n")
  cat("  convergence.csv      - convergence + failure history\n")

  invisible(particles)
}

# ---- Run it -----------------------------------------------------------------
run_smc_abc(CONFIG)
