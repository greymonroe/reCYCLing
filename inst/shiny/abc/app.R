library(shiny)
library(reCYCLing)
library(ggplot2)
library(data.table)
library(cowplot)
library(scales)

# ============================================================================
# POTENTIAL IMPROVEMENTS (not yet implemented)
# ============================================================================
# 1. [DONE] Adaptive perturbation kernel: Per-parameter SD = 2 * empirical SD
#    of accepted particles (Beaumont 2009). Replaces fixed 0.3 * range * 0.7^gen
#    decay. Future: could extend to full covariance matrix to let correlated
#    parameters move together.
#
# 2. Importance weights: In proper SMC-ABC, each particle carries a weight
#    based on prior(theta) / perturbation_kernel(theta). We currently do
#    unweighted resampling. Adding weights would be more statistically
#    rigorous for posterior inference.
#
# 3. Additional summary statistics: Currently using 4 stats (distance modes,
#    weight, mutation load). Consider adding: Shannon entropy, array size
#    trajectory shape, copy number distribution quantiles.
#
# 4. Multiple replicates per particle: Running each parameter set once adds
#    simulation noise to the distance. Averaging 2-3 replicates per particle
#    would give cleaner distance estimates (at computational cost).
#
# 5. Posterior diagnostics: Effective sample size, parameter correlation
#    plots, prior-vs-posterior overlay, ECDF of distances.
#
# NOTE: When making changes here, also sync to inst/scripts/smc_abc.R
#       (the CLI batch version). Both should stay in sync.
# ============================================================================

# ============================================================================
# ABC helper functions
# ============================================================================

PARAMS <- list(
  # p_local_dup: per-unit per-gen. Vignette used 4e-4.
  # At 10-20K units, want ~1-10 events/gen → 10^-5 to 10^-3
  list(name = "p_local_dup",  lo = -5,  hi = -3,   transform = function(x) 10^x),
  # p_distal_dup: per-array per-gen. Vignette had 1e-5 per-unit at 10K = 0.1/gen
  list(name = "p_distal_dup", lo = -2,  hi = -0.3, transform = function(x) 10^x),
  # p_del_chunk: per-unit per-gen. Vignette had 1e-5.
  list(name = "p_del_chunk",  lo = -6,  hi = -4,   transform = function(x) 10^x),
  list(name = "local_shape",  lo = 0.5, hi = 5,    transform = function(x) x),
  list(name = "local_scale",  lo = 1,   hi = 50,   transform = function(x) x),
  list(name = "distal_shape", lo = 0.5, hi = 5,    transform = function(x) x),
  list(name = "distal_scale", lo = 10,  hi = 2000, transform = function(x) x),
  list(name = "del_shape",    lo = 0.5, hi = 5,    transform = function(x) x),
  list(name = "del_scale",    lo = 1,   hi = 50,   transform = function(x) x)
)
param_names <- sapply(PARAMS, `[[`, "name")

sample_prior_one <- function() {
  row <- as.list(sapply(PARAMS, function(p) runif(1, p$lo, p$hi)))
  names(row) <- param_names
  as.data.table(row)
}

params_to_args <- function(row, mu_total) {
  list(
    p_local_dup  = 10^row$p_local_dup,
    p_distal_dup = 10^row$p_distal_dup,
    p_del_chunk  = 10^row$p_del_chunk,
    mu_total     = mu_total,
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
  # Subsample to cap computational cost — EM on 10K points is fast and gives
  # nearly identical mode estimates as running on millions.
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
    conv <- abs(mu1_new - mu1) < 1e-6 && abs(mu2_new - mu2) < 1e-6
    w1 <- w1_new; mu1 <- mu1_new; mu2 <- mu2_new
    sd1 <- max(sd1_new, 0.01); sd2 <- max(sd2_new, 0.01)
    if (conv) break
  }
  if (mu1 > mu2) list(log_mode1 = mu2, log_mode2 = mu1, weight1 = 1 - w1)
  else list(log_mode1 = mu1, log_mode2 = mu2, weight1 = w1)
}

run_and_summarize_one <- function(row, sim_gens, init_l, init_k0,
                                  sim_timeout, target = NULL, mu_total = 1e-4) {
  sim_args <- params_to_args(row, mu_total)
  t0 <- proc.time()[3]

  # Run for fixed number of generations (molecular clock).
  # max_units set very high — array size is emergent, not a target.
  hard_cap <- 100000L

  sim <- tryCatch(
    suppressWarnings(run_sim_ps(
      max_units = hard_cap, max_t = sim_gens, hard_cap = hard_cap,
      init_l = init_l, init_k0 = init_k0,
      p_local_dup = sim_args$p_local_dup, p_distal_dup = sim_args$p_distal_dup,
      p_del_chunk = sim_args$p_del_chunk, mu_total = sim_args$mu_total,
      local_dist = sim_args$local_dist, distal_dist = sim_args$distal_dist,
      del_dist = sim_args$del_dist, verbose = FALSE,
      compute_pairs = FALSE  # skip O(n^2) pairs for ABC; distances sampled below
    )),
    error = function(e) NULL
  )

  elapsed <- proc.time()[3] - t0

  if (is.null(sim)) return(list(stats = NULL, sim = NULL, elapsed = elapsed))
  mono <- sim$monomers
  n_gens <- length(sim$L_vec)

  # Reject if array went extinct or got absurdly large
  if (nrow(mono) < 50 || nrow(mono) > hard_cap)
    return(list(stats = NULL, sim = NULL, elapsed = elapsed))

  # --- Full summary statistics ---
  # Compute pairwise distances from sampled identical pairs directly,
  # bypassing the full pairs table (which can be millions of rows at 20K units).
  # For each unique sequence group, sample up to max_pairs_per_group pairs.
  dists <- mono[, {
    if (.N < 2L) NULL
    else {
      positions <- num
      n <- length(positions)
      if (n <= 200L) {
        # Small group: enumerate all pairs
        id1 <- rep(seq_len(n), times = (n - 1L):0)
        id2 <- unlist(lapply(seq_len(n - 1L), function(j) (j + 1L):n),
                       use.names = FALSE)
        list(d = abs(positions[id1] - positions[id2]))
      } else {
        # Large group: sample pairs
        n_sample <- min(5000L, as.numeric(n) * (n - 1L) %/% 2)
        s1 <- sample.int(n, n_sample, replace = TRUE)
        s2 <- sample.int(n, n_sample, replace = TRUE)
        keep <- s1 != s2
        list(d = abs(positions[s1[keep]] - positions[s2[keep]]))
      }
    }
  }, by = bponly]$d
  dists <- dists[dists > 0]
  if (length(dists) < 20) return(list(stats = NULL, sim = sim, elapsed = elapsed))

  # Quantile-based summary stats: robust, fast, no distribution assumptions.
  # Split distances into near (<= threshold) and far (> threshold).
  threshold <- if (!is.null(target)) target$near_far_threshold else 100
  near <- dists[dists <= threshold]
  far  <- dists[dists > threshold]
  near_frac <- length(near) / length(dists)
  med_near  <- if (length(near) >= 5) median(log10(near)) else NA_real_
  med_far   <- if (length(far) >= 5) median(log10(far)) else NA_real_

  # Require both components to have enough data
  if (is.na(med_near) || is.na(med_far))
    return(list(stats = NULL, sim = sim, elapsed = elapsed))

  # Mutation load — compute on a quick subsample first to reject wildly off particles
  sample_n <- min(200, nrow(mono))
  sample_seqs <- mono$bponly[sample.int(nrow(mono), sample_n)]
  ct_quick <- tryCatch(counts_long_nogap(sample_seqs), error = function(e) NULL)
  if (is.null(ct_quick)) return(list(stats = NULL, sim = sim, elapsed = elapsed))
  cons_quick <- ct_quick[symbol == consensus][order(pos)]
  load_quick <- mu_load_from_consensus(sample_seqs, cons_quick)
  mean_load_est <- mean(load_quick)

  # Biological gate: if mutation load is wildly off target, reject.
  # This is a hard constraint — we KNOW roughly what load should be.
  if (!is.null(target) && !is.na(target$mean_load) && target$mean_load > 0) {
    if (mean_load_est > target$mean_load * 3 || mean_load_est < target$mean_load / 3)
      return(list(stats = NULL, sim = sim, elapsed = elapsed))
  }

  # Full mutation load on all monomers
  ct <- tryCatch(counts_long_nogap(mono$bponly), error = function(e) NULL)
  if (is.null(ct)) return(list(stats = NULL, sim = sim, elapsed = elapsed))
  cons <- ct[symbol == consensus][order(pos)]
  load <- mu_load_from_consensus(mono$bponly, cons)

  stats <- data.table(med_near = med_near, med_far = med_far,
                      near_frac = near_frac, mean_load = mean(load),
                      n_units = nrow(mono), n_gens = n_gens)
  list(stats = stats, sim = sim, elapsed = elapsed)
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
    new[[pname]] <- max(0.01, new[[pname]])
  for (pname in c("local_scale", "distal_scale", "del_scale"))
    new[[pname]] <- max(0.1, new[[pname]])
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

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well { background: #f8f9fa; }
    h4 { margin-top: 12px; margin-bottom: 4px; color: #555; }
    .btn-success { width: 100%; font-size: 15px; padding: 8px; margin-top: 8px; }
    .btn-danger { width: 100%; font-size: 13px; padding: 6px; margin-top: 4px; }
    .status-text { font-weight: bold; margin: 8px 0; }
    .target-box { background: #e8f4e8; padding: 8px; border-radius: 5px;
                  margin-bottom: 10px; }
  "))),

  titlePanel("reCYCLing ABC Command Center"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      actionButton("start", "Start ABC", class = "btn-success"),
      actionButton("stop", "Stop", class = "btn-danger"),
      htmlOutput("status"),

      h4("Target statistics"),
      div(class = "target-box",
        numericInput("target_med_near", "Median near-distance", 20, min = 1),
        numericInput("target_med_far", "Median far-distance", 5000, min = 10),
        numericInput("target_near_frac", "Fraction near-distance", 0.5,
                     min = 0, max = 1, step = 0.05),
        numericInput("target_load", "Mean mutation load", 10, min = 0,
                     step = 1),
        numericInput("target_array_size", "Target array size", 20000, min = 100,
                     step = 1000),
        numericInput("near_far_threshold", "Near/far threshold (units)", 500,
                     min = 5, step = 50)
      ),

      h4("Molecular clock"),
      div(class = "target-box",
        numericInput("ancestor_age_my", "Ancestor age (million years)", 10,
                     min = 0.1, step = 1),
        numericInput("gen_time_yr", "Generation time (years)", 10,
                     min = 1, step = 1),
        numericInput("compression", "Time compression factor", 1000,
                     min = 1, step = 100),
        numericInput("mu_per_base_real", "Real mutation rate (per base/gen)",
                     5.6e-8, min = 1e-10, max = 1e-6, step = 1e-9),
        numericInput("init_l", "Monomer length (bp)", 178, min = 10, step = 10),
        htmlOutput("clock_info")
      ),

      h4("Algorithm"),
      numericInput("n_particles", "Particles per generation", 500,
                   min = 10, step = 50),
      numericInput("max_generations", "Max generations", 10, min = 1),
      sliderInput("retention_frac", "Retention fraction", 0.3,
                  min = 0.1, max = 0.9, step = 0.05),
      sliderInput("perturbation_sd", "Initial perturbation (Gen 0 only)", 0.3,
                  min = 0.05, max = 0.5, step = 0.05),

      h4("Simulation"),
      numericInput("sim_timeout", "Per-particle timeout (sec)", 30,
                   min = 1, max = 300, step = 5)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("Live Dashboard",
          fluidRow(
            column(9, plotOutput("stats_evolution_plot", height = "300px")),
            column(3, plotOutput("convergence_plot", height = "300px"))
          ),
          fluidRow(
            column(12, plotOutput("param_evolution_plot", height = "350px"))
          )
        ),
        tabPanel("Best Fit Explorer",
          fluidRow(
            column(3,
              h4("Top particles"),
              tableOutput("top_table"),
              numericInput("inspect_rank", "Inspect rank #", 1,
                           min = 1, max = 10),
              actionButton("rerun_best", "Re-run & Visualize",
                           class = "btn-success")
            ),
            column(9,
              plotOutput("best_summary_plot", height = "800px")
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {

  # Reactive state
  state <- reactiveValues(
    running       = FALSE,
    generation    = 0L,
    particle_idx  = 0L,        # counts only VALID particles (filled slots)
    particles     = NULL,       # current gen data.table (params + stats)
    kept          = NULL,       # survivors from selection
    convergence   = data.table(gen = integer(), best = numeric(),
                               median = numeric(), n_valid = integer(),
                               n_failed = integer(), n_timeout = integer(),
                               elapsed_s = numeric()),
    all_particles = list(),     # all gens' particles for posteriors
    best_sim      = NULL,       # stored sim result for best particle
    param_sds     = NULL,       # per-parameter perturbation SDs (adaptive)
    stall_count   = 0L,
    last_particle_time = NA_real_,  # seconds taken by last particle
    n_failed      = 0L,             # failed attempts this generation
    n_timeout     = 0L,             # timed-out attempts this generation
    n_attempts    = 0L,             # total attempts this generation
    gen_start_time = NA_real_,      # wall-clock start of current generation
    stat_scales   = NULL            # SD of each summary stat (from Gen 0 pilot)
  )

  # Derived clock values
  get_clock <- reactive({
    real_gens <- input$ancestor_age_my * 1e6 / input$gen_time_yr
    sim_gens <- round(real_gens / input$compression)
    mu_compressed <- input$mu_per_base_real * input$compression
    expected_load <- mu_compressed * input$init_l * sim_gens
    list(real_gens = real_gens, sim_gens = sim_gens,
         mu_compressed = mu_compressed, expected_load = expected_load)
  })

  output$clock_info <- renderUI({
    ck <- get_clock()
    HTML(paste0(
      '<small style="color: #555;">',
      'Real gens: ', format(ck$real_gens, big.mark = ","),
      '<br>Sim gens: ', format(ck$sim_gens, big.mark = ","),
      '<br>Effective mu: ', signif(ck$mu_compressed, 3),
      '<br>Expected load: ', round(ck$expected_load, 1),
      '</small>'))
  })

  get_target <- reactive({
    list(med_near  = log10(input$target_med_near),
         med_far   = log10(input$target_med_far),
         near_frac = input$target_near_frac,
         mean_load = input$target_load,
         target_array_size = input$target_array_size,
         near_far_threshold = input$near_far_threshold)
  })

  # --- Start / Stop ---
  observeEvent(input$start, {
    state$running <- TRUE
    state$generation <- 0L
    state$particle_idx <- 0L
    state$particles <- NULL
    state$kept <- NULL
    state$convergence <- data.table(gen = integer(), best = numeric(),
                                     median = numeric(), n_valid = integer(),
                                     n_failed = integer(), n_timeout = integer(),
                                     elapsed_s = numeric())
    state$all_particles <- list()
    state$best_sim <- NULL
    # Initial perturbation SDs: use prior range * perturbation_scale for Gen 0
    state$param_sds <- compute_param_sds(NULL, input$perturbation_sd)
    state$stall_count <- 0L
    state$last_particle_time <- NA_real_
    state$n_failed <- 0L
    state$n_timeout <- 0L
    state$n_attempts <- 0L
    state$gen_start_time <- proc.time()[3]
    state$stat_scales <- NULL
  })

  observeEvent(input$stop, { state$running <- FALSE })

  # --- Main cooperative loop: process a batch of particles per cycle ---
  # Runs multiple particles before yielding to Shiny for UI updates.
  # This avoids the overhead of ggplot re-rendering after every single particle.
  observe({
    # Schedule next cycle FIRST so return() inside isolate can't skip it
    invalidateLater(100)
    if (!state$running) return()

    isolate({
      target <- get_target()
      n_particles <- input$n_particles
      timeout <- input$sim_timeout

      # Process a batch of particles. Large budget since plots only update
      # between generations (they depend on all_particles, not particles).
      batch_start <- proc.time()[3]
      batch_budget <- 30  # seconds before yielding to UI

      while (proc.time()[3] - batch_start < batch_budget) {
        gen <- state$generation
        idx <- state$particle_idx

        # --- Need to start a new generation? ---
        if (idx >= n_particles) {
          particles <- state$particles
          n_valid <- nrow(particles)  # all particles are valid (failures retried)

          # After Gen 0: compute stat scales from pilot batch and recompute
          # all distances so they're properly normalized going forward.
          if (gen == 0 && is.null(state$stat_scales)) {
            state$stat_scales <- compute_stat_scales(particles, target)
            cat("Stat scales from Gen 0 pilot:\n")
            print(state$stat_scales)
            # Recompute Gen 0 distances with proper scaling
            for (ri in seq_len(nrow(particles))) {
              particles[ri, distance := compute_distance(
                .SD, target, state$stat_scales),
                .SDcols = c("med_near", "med_far", "near_frac", "mean_load", "n_units")]
            }
            state$particles <- particles
          }

          best_d <- min(particles$distance)
          med_d  <- median(particles$distance)

          gen_elapsed <- if (!is.na(state$gen_start_time))
            round(proc.time()[3] - state$gen_start_time) else NA_real_
          state$convergence <- rbind(state$convergence,
            data.table(gen = gen, best = best_d, median = med_d,
                       n_valid = n_valid,
                       n_failed = state$n_failed,
                       n_timeout = state$n_timeout,
                       elapsed_s = gen_elapsed))
          state$all_particles[[gen + 1]] <- copy(particles)

          # Selection
          n_keep <- max(3L, as.integer(n_valid * input$retention_frac))
          valid <- particles[order(distance)]
          state$kept <- valid[seq_len(min(n_keep, nrow(valid)))]

          # Convergence check
          conv <- state$convergence
          if (nrow(conv) >= 2) {
            prev_med <- conv$median[nrow(conv) - 1]
            curr_med <- med_d
            improvement <- if (is.finite(prev_med) && prev_med > 0)
              (prev_med - curr_med) / prev_med else 1
            if (improvement < 0.01) {
              state$stall_count <- state$stall_count + 1L
            } else {
              state$stall_count <- 0L
            }
            if (state$stall_count >= 3 || gen >= input$max_generations) {
              state$running <- FALSE
              return()
            }
          }

          # Prepare next generation — compute adaptive perturbation SDs
          # from the empirical spread of survivors (Beaumont 2009)
          state$generation <- gen + 1L
          state$particle_idx <- 0L
          state$particles <- NULL
          state$param_sds <- compute_param_sds(state$kept, input$perturbation_sd)
          state$n_failed <- 0L
          state$n_timeout <- 0L
          state$n_attempts <- 0L
          state$gen_start_time <- proc.time()[3]
          next  # Continue processing in this batch
        }

        # --- Propose a particle ---
        if (gen == 0 && is.null(state$kept)) {
          proposed <- sample_prior_one()
        } else {
          kept <- state$kept
          if (is.null(kept) || nrow(kept) == 0) {
            state$running <- FALSE
            return()
          }
          parent <- kept[sample.int(nrow(kept), 1), ..param_names]
          proposed <- perturb_particle(parent, state$param_sds)
        }

        # Run simulation (max_t is the guard rail for slow sims).
        # No setTimeLimit — it can crash R when interrupting C++ code.
        res <- tryCatch(
          run_and_summarize_one(
            proposed,
            sim_gens    = get_clock()$sim_gens,
            init_l      = input$init_l,
            init_k0     = 10L,
            sim_timeout = timeout,
            target      = target,
            mu_total    = get_clock()$mu_compressed
          ),
          error = function(e) {
            list(stats = NULL, sim = NULL, elapsed = NA_real_)
          }
        )

        state$last_particle_time <- res$elapsed
        state$n_attempts <- state$n_attempts + 1L

        timed_out <- is.null(res$stats) && !is.null(res$elapsed) &&
                     res$elapsed >= (timeout - 0.5)

        if (is.null(res$stats)) {
          # Failed particle — track it but DON'T increment particle_idx.
          # We'll retry with a new draw to fill this slot.
          if (timed_out) {
            state$n_timeout <- state$n_timeout + 1L
          } else {
            state$n_failed <- state$n_failed + 1L
          }
          # Safety valve: if too many consecutive failures, give up on this gen
          if (state$n_attempts > n_particles * 5) {
            state$running <- FALSE
            return()
          }
        } else {
          # Valid particle — store it and advance to next slot
          d <- compute_distance(res$stats, target, state$stat_scales)
          row <- cbind(proposed, res$stats, data.table(distance = d,
            elapsed = res$elapsed, status = "ok"))
          if (is.null(state$best_sim) || d < min(state$convergence$best, Inf)) {
            if (!is.null(res$sim)) state$best_sim <- res$sim
          }

          if (is.null(state$particles)) {
            state$particles <- row
          } else {
            state$particles <- rbind(state$particles, row, fill = TRUE)
          }
          state$particle_idx <- idx + 1L
        }
      }  # end batch while loop
    })
  })

  # --- Status display ---
  output$status <- renderUI({
    if (state$running) {
      gen <- state$generation
      idx <- state$particle_idx
      n   <- input$n_particles
      pct <- round(100 * idx / n)
      last_t <- state$last_particle_time
      gen_elapsed <- if (!is.na(state$gen_start_time))
        round(proc.time()[3] - state$gen_start_time) else 0

      timing_info <- if (!is.na(last_t))
        paste0(' | Last: ', round(last_t, 1), 's') else ''

      retry_info <- ''
      if (state$n_failed > 0 || state$n_timeout > 0) {
        parts <- c()
        if (state$n_failed > 0)  parts <- c(parts, paste0(state$n_failed, ' failed'))
        if (state$n_timeout > 0) parts <- c(parts, paste0(state$n_timeout, ' timed out'))
        retry_info <- paste0(' | Retried: ', paste(parts, collapse = ', '))
      }

      HTML(paste0(
        '<div class="status-text" style="color: green;">',
        'Gen ', gen, ' | Filled ', idx, '/', n,
        ' (', pct, '%)', timing_info, retry_info,
        '<br>Attempts: ', state$n_attempts, ' | Elapsed: ', gen_elapsed, 's',
        '</div>'))
    } else if (nrow(state$convergence) > 0) {
      best <- min(state$convergence$best)
      HTML(paste0(
        '<div class="status-text" style="color: #333;">',
        'Done: ', nrow(state$convergence), ' generations | ',
        'Best distance: ', round(best, 3), '</div>'))
    } else {
      HTML('<div class="status-text" style="color: #888;">Ready</div>')
    }
  })

  # --- Convergence plot ---
  output$convergence_plot <- renderPlot({
    conv <- state$convergence
    if (nrow(conv) == 0) {
      # Show progress within current generation
      p <- state$particles
      if (!is.null(p) && nrow(p) > 0) {
        valid <- p[is.finite(distance)]
        if (nrow(valid) > 0) {
          plot_dt <- data.table(idx = seq_len(nrow(valid)),
                                distance = valid$distance)
          return(ggplot(plot_dt, aes(x = idx, y = distance)) +
            geom_point(colour = "steelblue", size = 2, alpha = 0.7) +
            theme_classic(base_size = 13) +
            xlab("Particle") + ylab("Distance") +
            ggtitle("Generation 0 — distances so far"))
        }
      }
      return(ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "Waiting...",
                 size = 6, colour = "grey50"))
    }

    conv <- copy(conv)
    conv[, n_fail_total := n_failed + n_timeout]
    conv[, label := paste0(n_fail_total, "f / ", elapsed_s, "s")]

    ggplot(conv, aes(x = gen)) +
      geom_line(aes(y = median, colour = "Median"), linewidth = 1) +
      geom_point(aes(y = median, colour = "Median"), size = 2) +
      geom_line(aes(y = best, colour = "Best"), linewidth = 1) +
      geom_point(aes(y = best, colour = "Best"), size = 2) +
      geom_text(aes(y = max(median) * 1.08, label = label),
                size = 2, colour = "grey40", angle = 30) +
      scale_colour_manual(values = c("Median" = "steelblue",
                                      "Best" = "firebrick")) +
      theme_classic(base_size = 10) +
      labs(x = "Generation", y = "Distance",
           colour = NULL, title = "Convergence (failures / time)") +
      theme(legend.position = c(0.8, 0.8),
            legend.background = element_rect(fill = "transparent"),
            legend.key = element_rect(fill = "transparent"))
  }, res = 100)

  # --- Summary stat evolution across generations (violins per gen) ---
  output$stats_evolution_plot <- renderPlot({
    all_p <- state$all_particles
    if (length(all_p) == 0) {
      return(ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "Summary stats will appear here",
                 size = 6, colour = "grey50"))
    }

    target <- get_target()

    # Collect stats from all generations
    gen_list <- lapply(seq_along(all_p), function(i) {
      p <- all_p[[i]]
      valid <- p[is.finite(distance)]
      if (nrow(valid) == 0) return(NULL)
      data.table(
        gen        = i - 1L,
        med_near   = 10^valid$med_near,
        med_far    = 10^valid$med_far,
        near_frac  = valid$near_frac,
        mean_load  = valid$mean_load,
        array_size = valid$n_units
      )
    })
    gen_dt <- rbindlist(gen_list[!sapply(gen_list, is.null)])
    if (nrow(gen_dt) == 0) return(ggplot() + theme_void())

    target_vals <- data.table(
      stat = c("med_near", "med_far", "near_frac", "mean_load", "array_size"),
      target = c(10^target$med_near, 10^target$med_far,
                 target$near_frac, target$mean_load, target$target_array_size)
    )

    melt_dt <- melt(gen_dt, id.vars = "gen",
                    measure.vars = c("med_near", "med_far", "near_frac",
                                     "mean_load", "array_size"),
                    variable.name = "stat", value.name = "value")
    melt_dt[, gen_f := factor(gen)]

    ggplot(melt_dt, aes(x = gen_f, y = value)) +
      geom_violin(fill = "steelblue", alpha = 0.4, scale = "width",
                  linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.2) +
      geom_hline(data = target_vals, aes(yintercept = target),
                 colour = "firebrick", linewidth = 0.8, linetype = "dashed") +
      facet_wrap(~ stat, scales = "free_y", nrow = 1) +
      theme_classic(base_size = 11) +
      theme(strip.background = element_blank()) +
      labs(x = "Generation", y = NULL,
           title = "Summary statistics across generations (target = red dashed)")
  }, res = 100)

  # --- Parameter evolution across generations (violins) ---
  output$param_evolution_plot <- renderPlot({
    all_p <- state$all_particles
    if (length(all_p) == 0) {
      return(ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "Parameter posteriors will appear here",
                 size = 6, colour = "grey50"))
    }

    rate_params <- c("p_local_dup", "p_distal_dup", "p_del_chunk")
    shape_params <- c("local_shape", "local_scale", "distal_shape",
                      "distal_scale", "del_shape", "del_scale")

    gen_list <- lapply(seq_along(all_p), function(i) {
      p <- all_p[[i]]
      valid <- p[is.finite(distance)]
      if (nrow(valid) == 0) return(NULL)
      d <- copy(valid)
      # Transform rates to natural scale
      for (j in 1:3) {
        pn <- PARAMS[[j]]$name
        set(d, j = pn, value = 10^d[[pn]])
      }
      d[, gen := i - 1L]
      d
    })
    gen_dt <- rbindlist(gen_list[!sapply(gen_list, is.null)], fill = TRUE)
    if (nrow(gen_dt) == 0) return(ggplot() + theme_void())

    gen_dt[, gen_f := factor(gen)]

    # Rates + diagnostics (log scale for rates)
    diag_params <- c(rate_params, "n_units", "n_gens")
    melt_rates <- melt(gen_dt[, c("gen_f", diag_params), with = FALSE],
                       id.vars = "gen_f", variable.name = "param", value.name = "value")
    p1 <- ggplot(melt_rates, aes(x = gen_f, y = value)) +
      geom_violin(fill = "darkorange", alpha = 0.4, scale = "width", linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.2, alpha = 0.15) +
      facet_wrap(~ param, scales = "free_y", nrow = 1) +
      scale_y_log10() +
      theme_classic(base_size = 9) +
      theme(strip.background = element_blank()) +
      labs(x = "Generation", y = NULL, title = "Rates + diagnostics (n_units, n_gens)")

    # Shape/scale (linear)
    melt_shapes <- melt(gen_dt[, c("gen_f", shape_params), with = FALSE],
                        id.vars = "gen_f", variable.name = "param", value.name = "value")
    p2 <- ggplot(melt_shapes, aes(x = gen_f, y = value)) +
      geom_violin(fill = "darkorange", alpha = 0.4, scale = "width", linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.2, alpha = 0.15) +
      facet_wrap(~ param, scales = "free_y", nrow = 1) +
      theme_classic(base_size = 9) +
      theme(strip.background = element_blank()) +
      labs(x = "Generation", y = NULL, title = "Chunk-size parameters")

    cowplot::plot_grid(p1, p2, ncol = 1)
  }, res = 100)

  # --- Best Fit Explorer ---
  output$top_table <- renderTable({
    all_p <- state$all_particles
    if (length(all_p) == 0) return(NULL)
    last <- all_p[[length(all_p)]]
    valid <- last[is.finite(distance)][order(distance)]
    if (nrow(valid) == 0) return(NULL)
    n_show <- min(10, nrow(valid))
    top <- valid[seq_len(n_show)]
    data.frame(
      rank        = seq_len(n_show),
      distance    = round(top$distance, 3),
      med_near    = round(10^top$med_near),
      med_far     = round(10^top$med_far),
      near_frac   = round(top$near_frac, 2),
      mean_load   = round(top$mean_load, 1),
      p_local     = signif(10^top$p_local_dup, 2),
      p_distal    = signif(10^top$p_distal_dup, 2),
      p_del       = signif(10^top$p_del_chunk, 2),
      local_shp   = round(top$local_shape, 1),
      local_scl   = round(top$local_scale, 1),
      distal_shp  = round(top$distal_shape, 1),
      distal_scl  = round(top$distal_scale, 0),
      del_shp     = round(top$del_shape, 1),
      del_scl     = round(top$del_scale, 1)
    )
  }, digits = 3)

  observeEvent(input$rerun_best, {
    all_p <- state$all_particles
    if (length(all_p) == 0) return()
    last <- all_p[[length(all_p)]]
    valid <- last[is.finite(distance)][order(distance)]
    rank <- input$inspect_rank
    if (rank > nrow(valid)) return()
    row <- valid[rank]

    withProgress(message = "Re-running best particle...", value = 0.3, {
      ck <- get_clock()
      sim_args <- params_to_args(row, ck$mu_compressed)
      sim <- suppressWarnings(run_sim_ps(
        max_units = 100000, max_t = ck$sim_gens, init_l = input$init_l,
        init_k0 = 10,
        p_local_dup = sim_args$p_local_dup,
        p_distal_dup = sim_args$p_distal_dup,
        p_del_chunk = sim_args$p_del_chunk,
        mu_total = ck$mu_compressed,
        local_dist = sim_args$local_dist,
        distal_dist = sim_args$distal_dist,
        del_dist = sim_args$del_dist,
        verbose = FALSE
      ))
      incProgress(0.7)
      state$best_sim <- sim
    })
  })

  output$best_summary_plot <- renderPlot({
    sim <- state$best_sim
    req(sim)
    plot_ps_summary(sim, title = "Best fit simulation")
  }, res = 100)
}

# ============================================================================
shinyApp(ui, server)
