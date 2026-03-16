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
#     --target_mode1 100 --target_mode2 2000 --target_weight1 0.3 \
#     --target_load 10 --particles 200 --max_generations 20 \
#     --max_units 20000 --outdir results/run1
#
# All --target_* and --*_weight flags control what data you're fitting to.
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
# =============================================================================

library(reCYCLing)
library(data.table)

# ---- Configuration with defaults -------------------------------------------

CONFIG <- list(
  # Target summary statistics (what your real data looks like)
  target = list(
    log_mode1  = 2.0,        # log10 of first distance mode (100)
    log_mode2  = 3.3,        # log10 of second distance mode (~2000)
    weight1    = 0.3,        # relative weight of the near-mode peak
    mean_load  = 10          # mean mismatches to consensus per monomer
  ),

  # Distance scaling: computed automatically from Gen 0 pilot batch.
  # Each summary stat is normalized by its SD across the pilot particles,
  # so the distance is unitless and each stat contributes proportionally.

  # SMC-ABC algorithm settings
  n_particles     = 200,     # particles per generation (all slots filled)
  max_generations = 20,      # upper limit on generations
  retention_frac  = 0.5,     # fraction of particles kept each generation
  perturbation_sd = 0.3,     # perturbation scale (used for Gen 0 fallback)

  # Convergence detection
  conv_tol       = 0.01,     # stop if median distance improves < 1%
  conv_patience  = 3,        # must stall for this many consecutive gens

  # Simulation settings
  max_units   = 20000,       # grow array to this many monomers
  max_t       = 50000,       # generation cap per simulation
  sim_timeout = 30,          # seconds before killing a slow sim
  init_l      = 30,          # monomer length (bp)
  init_k0     = 10,          # initial copy number

  # Output
  outdir = "results"
)

# ---- Parse command-line args ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  v <- if (i < length(args)) args[i + 1] else NULL

  matched <- TRUE
  if      (a == "--outdir")          { CONFIG$outdir <- v }
  else if (a == "--particles")       { CONFIG$n_particles <- as.integer(v) }
  else if (a == "--max_generations") { CONFIG$max_generations <- as.integer(v) }
  else if (a == "--max_units")       { CONFIG$max_units <- as.integer(v) }
  else if (a == "--max_t")           { CONFIG$max_t <- as.integer(v) }
  else if (a == "--sim_timeout")     { CONFIG$sim_timeout <- as.numeric(v) }
  else if (a == "--init_l")          { CONFIG$init_l <- as.integer(v) }
  else if (a == "--init_k0")         { CONFIG$init_k0 <- as.integer(v) }
  else if (a == "--retention_frac")  { CONFIG$retention_frac <- as.numeric(v) }
  else if (a == "--perturbation_sd") { CONFIG$perturbation_sd <- as.numeric(v) }
  else if (a == "--conv_tol")        { CONFIG$conv_tol <- as.numeric(v) }
  else if (a == "--conv_patience")   { CONFIG$conv_patience <- as.integer(v) }
  # Target statistics
  else if (a == "--target_mode1")    { CONFIG$target$log_mode1 <- log10(as.numeric(v)) }
  else if (a == "--target_mode2")    { CONFIG$target$log_mode2 <- log10(as.numeric(v)) }
  else if (a == "--target_weight1")  { CONFIG$target$weight1 <- as.numeric(v) }
  else if (a == "--target_load")     { CONFIG$target$mean_load <- as.numeric(v) }
  else if (a == "--help") {
    cat("
SMC-ABC parameter inference for reCYCLing tandem repeat simulations.

Target statistics (what your real data looks like):
  --target_mode1 NUM     First distance mode in units (default: 100)
  --target_mode2 NUM     Second distance mode in units (default: 2000)
  --target_weight1 NUM   Weight of the near-mode peak, 0-1 (default: 0.3)
  --target_load NUM      Mean mismatches to consensus (default: 10)

Distance scaling: Automatically computed from Gen 0 pilot batch. Each
  summary stat is normalized by its SD, so no manual weights needed.

Algorithm settings:
  --particles NUM        Particles per generation (default: 200)
  --max_generations NUM  Max generations before stopping (default: 20)
  --retention_frac NUM   Fraction of particles kept each gen (default: 0.5)
  --perturbation_sd NUM  Perturbation scale for Gen 0 fallback (default: 0.3)
                         After Gen 0, perturbation adapts to empirical spread of survivors.
  --conv_tol NUM         Convergence tolerance, fraction (default: 0.01)
  --conv_patience NUM    Gens without improvement to stop (default: 3)

Simulation settings:
  --max_units NUM        Array size target (default: 20000)
  --max_t NUM            Generation cap per sim (default: 50000)
  --sim_timeout NUM      Seconds before killing a slow sim (default: 30)
  --init_l NUM           Monomer length in bp (default: 30)
  --init_k0 NUM          Initial copy number (default: 10)

Output:
  --outdir PATH          Output directory (default: results)
")
    quit(save = "no", status = 0)
  }
  else { stop("Unknown argument: ", a, "\nUse --help for usage.") }

  i <- i + 2
}

dir.create(CONFIG$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Parameter space definition ---------------------------------------------
# lo/hi define the initial prior sampling range. Perturbation is NOT clamped
# to these bounds — particles can explore beyond them freely.

PARAMS <- list(
  list(name = "p_local_dup",   lo = -4,   hi = -2,   transform = function(x) 10^x),
  list(name = "p_distal_dup",  lo = -6,   hi = -3,   transform = function(x) 10^x),
  list(name = "p_del_chunk",   lo = -6,   hi = -3,   transform = function(x) 10^x),
  list(name = "mu_total",      lo = -5,   hi = -2,   transform = function(x) 10^x),
  list(name = "local_shape",   lo = 0.5,  hi = 5,    transform = function(x) x),
  list(name = "local_scale",   lo = 1,    hi = 50,   transform = function(x) x),
  list(name = "distal_shape",  lo = 0.5,  hi = 5,    transform = function(x) x),
  list(name = "distal_scale",  lo = 10,   hi = 2000, transform = function(x) x),
  list(name = "del_shape",     lo = 0.5,  hi = 5,    transform = function(x) x),
  list(name = "del_scale",     lo = 1,    hi = 50,   transform = function(x) x)
)

param_names <- sapply(PARAMS, `[[`, "name")

# ---- Helper functions -------------------------------------------------------

sample_prior <- function(n) {
  mat <- matrix(NA_real_, nrow = n, ncol = length(PARAMS))
  colnames(mat) <- param_names
  for (j in seq_along(PARAMS)) {
    mat[, j] <- runif(n, PARAMS[[j]]$lo, PARAMS[[j]]$hi)
  }
  as.data.table(mat)
}

sample_prior_one <- function() {
  row <- as.list(sapply(PARAMS, function(p) runif(1, p$lo, p$hi)))
  names(row) <- param_names
  as.data.table(row)
}

params_to_args <- function(row) {
  list(
    p_local_dup  = PARAMS[[1]]$transform(row$p_local_dup),
    p_distal_dup = PARAMS[[2]]$transform(row$p_distal_dup),
    p_del_chunk  = PARAMS[[3]]$transform(row$p_del_chunk),
    mu_total     = PARAMS[[4]]$transform(row$mu_total),
    local_dist   = list(type = "gamma", shape = row$local_shape,
                        scale = row$local_scale),
    distal_dist  = list(type = "gamma", shape = row$distal_shape,
                        scale = row$distal_scale),
    del_dist     = list(type = "gamma", shape = row$del_shape,
                        scale = row$del_scale)
  )
}

fit_2gauss <- function(x, max_iter = 100, max_n = 10000) {
  if (length(x) < 10) return(NULL)
  # Subsample to cap computational cost
  if (length(x) > max_n) x <- x[sample.int(length(x), max_n)]

  k <- tryCatch(kmeans(x, 2, nstart = 5), error = function(e) NULL)
  if (is.null(k)) return(NULL)

  mu1 <- min(k$centers); mu2 <- max(k$centers)
  cl1 <- which.min(k$centers)
  sd1 <- sd(x[k$cluster == cl1])
  sd2 <- sd(x[k$cluster != cl1])
  if (is.na(sd1) || sd1 == 0) sd1 <- 0.1
  if (is.na(sd2) || sd2 == 0) sd2 <- 0.1
  w1 <- mean(k$cluster == cl1)

  for (iter in seq_len(max_iter)) {
    d1 <- w1 * dnorm(x, mu1, sd1)
    d2 <- (1 - w1) * dnorm(x, mu2, sd2)
    denom <- d1 + d2 + 1e-300
    gamma1 <- d1 / denom

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

  if (mu1 > mu2) {
    list(log_mode1 = mu2, log_mode2 = mu1, weight1 = 1 - w1)
  } else {
    list(log_mode1 = mu1, log_mode2 = mu2, weight1 = w1)
  }
}

#' Run a single simulation and extract summary statistics.
#' Uses setTimeLimit for wall-clock timeout (no forking).
#' Includes early rejection checks before expensive summary stats.
run_and_summarize <- function(row, config) {
  sim_args <- params_to_args(row)
  hard_cap <- min(config$max_units * 2L, 40000L)

  sim <- tryCatch(
    {
      setTimeLimit(elapsed = config$sim_timeout, transient = TRUE)
      suppressWarnings(run_sim_ps(
        max_units = config$max_units, max_t = config$max_t,
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
        compute_pairs = FALSE  # skip O(n^2) pairs for ABC; distances sampled below
      ))
    },
    error = function(e) NULL,
    finally = setTimeLimit(elapsed = Inf, transient = TRUE)
  )
  if (is.null(sim)) return(NULL)

  mono <- sim$monomers

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
        n_sample <- min(5000L, n * (n - 1L) %/% 2L)
        s1 <- sample.int(n, n_sample, replace = TRUE)
        s2 <- sample.int(n, n_sample, replace = TRUE)
        keep <- s1 != s2
        list(d = abs(positions[s1[keep]] - positions[s2[keep]]))
      }
    }
  }, by = bponly]$d
  dists <- dists[dists > 0]
  if (length(dists) < 20) return(NULL)

  fit <- fit_2gauss(log10(dists))
  if (is.null(fit)) return(NULL)

  ct <- tryCatch(counts_long_nogap(mono$bponly), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  cons <- ct[symbol == consensus][order(pos)]
  load <- mu_load_from_consensus(mono$bponly, cons)

  data.table(
    log_mode1 = fit$log_mode1,
    log_mode2 = fit$log_mode2,
    weight1   = fit$weight1,
    mean_load = mean(load),
    n_units   = nrow(mono),
    n_pairs   = length(dists)
  )
}

compute_distance <- function(sim_stats, target, stat_scales = NULL) {
  diffs <- c(
    log_mode1 = sim_stats$log_mode1 - target$log_mode1,
    log_mode2 = sim_stats$log_mode2 - target$log_mode2,
    weight1   = sim_stats$weight1   - target$weight1,
    mean_load = (sim_stats$mean_load - target$mean_load) / max(target$mean_load, 0.1)
  )
  # Normalize by SD of each stat (computed from Gen 0 pilot batch)
  if (!is.null(stat_scales)) {
    diffs <- diffs / stat_scales
  }
  sqrt(sum(diffs^2))
}

compute_stat_scales <- function(particles, target) {
  rel_load <- (particles$mean_load - target$mean_load) / max(target$mean_load, 0.1)
  scales <- c(
    log_mode1 = sd(particles$log_mode1, na.rm = TRUE),
    log_mode2 = sd(particles$log_mode2, na.rm = TRUE),
    weight1   = sd(particles$weight1, na.rm = TRUE),
    mean_load = sd(rel_load, na.rm = TRUE)
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
  # Enforce physical constraints only
  for (pname in c("local_shape", "distal_shape", "del_shape"))
    new[[pname]] <- max(0.01, new[[pname]])
  for (pname in c("local_scale", "distal_scale", "del_scale"))
    new[[pname]] <- max(0.1, new[[pname]])
  new
}

compute_param_sds <- function(kept, fallback_scale = 0.3) {
  sds <- list()
  for (j in seq_along(PARAMS)) {
    pname <- PARAMS[[j]]$name
    emp_sd <- if (!is.null(kept) && nrow(kept) >= 3)
      sd(kept[[pname]], na.rm = TRUE) else NA_real_
    if (is.na(emp_sd) || emp_sd < 1e-10) {
      sds[[pname]] <- fallback_scale * (PARAMS[[j]]$hi - PARAMS[[j]]$lo)
    } else {
      sds[[pname]] <- 2 * emp_sd
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
  cat("Max units:", config$max_units, "\n")
  cat("Timeout:", config$sim_timeout, "s per particle\n")
  cat("Targets: mode1=", 10^config$target$log_mode1,
      " mode2=", 10^config$target$log_mode2,
      " w1=", config$target$weight1,
      " load=", config$target$mean_load, "\n")
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
      .SDcols = c("log_mode1", "log_mode2", "weight1", "mean_load")]
  }
  best_dist <- min(particles$distance)
  med_dist  <- median(particles$distance)
  cat(sprintf("  Rescaled: Best: %.3f | Median: %.3f\n", best_dist, med_dist))

  fwrite(particles, file.path(config$outdir, "gen_0_particles.csv"))

  # -- Subsequent generations --
  for (gen in seq_len(config$max_generations)) {
    cat(sprintf("\n--- Generation %d ---\n", gen))

    # Select survivors
    n_keep <- max(3L, as.integer(nrow(particles) * config$retention_frac))
    particles <- particles[order(distance)]
    kept <- particles[seq_len(min(n_keep, nrow(particles)))]

    cat(sprintf("  Kept %d particles (threshold: %.3f)\n",
                nrow(kept), max(kept$distance)))

    # Show current best
    best <- kept[1]
    cat(sprintf("  Best: mode1=%.0f mode2=%.0f w1=%.2f load=%.1f (dist=%.3f)\n",
                10^best$log_mode1, 10^best$log_mode2,
                best$weight1, best$mean_load, best$distance))

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
  top_display[, mode1 := 10^log_mode1]
  top_display[, mode2 := 10^log_mode2]
  print(top_display[, .(p_local_dup, p_distal_dup, p_del_chunk, mu_total,
                         local_shape, local_scale, distal_shape, distal_scale,
                         del_shape, del_scale,
                         mode1, mode2, weight1, mean_load, distance)],
        digits = 3)

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
