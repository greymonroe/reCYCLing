library(shiny)
library(reCYCLing)
library(ggplot2)
library(data.table)
library(cowplot)
library(scales)

# ============================================================================
# POTENTIAL IMPROVEMENTS (not yet implemented)
# ============================================================================
# 1. Adaptive perturbation kernel: Use the empirical covariance matrix of
#    accepted particles as the perturbation kernel (Beaumont 2009), instead
#    of independent Gaussian perturbation per parameter. This lets correlated
#    parameters move together and improves mixing.
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
  list(name = "p_local_dup",  lo = -4,  hi = -2,   transform = function(x) 10^x),
  list(name = "p_distal_dup", lo = -6,  hi = -3,   transform = function(x) 10^x),
  list(name = "p_del_chunk",  lo = -6,  hi = -3,   transform = function(x) 10^x),
  list(name = "mu_total",     lo = -5,  hi = -2,   transform = function(x) 10^x),
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

params_to_args <- function(row) {
  list(
    p_local_dup  = 10^row$p_local_dup,
    p_distal_dup = 10^row$p_distal_dup,
    p_del_chunk  = 10^row$p_del_chunk,
    mu_total     = 10^row$mu_total,
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

run_and_summarize_one <- function(row, max_units, max_t, init_l, init_k0,
                                  sim_timeout, target = NULL) {
  sim_args <- params_to_args(row)
  t0 <- proc.time()[3]

  # Use a reduced hard_cap to prevent memory blowup in the Shiny context
  hard_cap <- min(max_units * 2L, 20000L)

  sim <- tryCatch(
    suppressWarnings(run_sim_ps(
      max_units = max_units, max_t = max_t, hard_cap = hard_cap,
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
  n_gens <- length(sim$L_vec)  # how many generations the sim actually ran

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
        n_sample <- min(5000L, n * (n - 1L) %/% 2L)
        s1 <- sample.int(n, n_sample, replace = TRUE)
        s2 <- sample.int(n, n_sample, replace = TRUE)
        keep <- s1 != s2
        list(d = abs(positions[s1[keep]] - positions[s2[keep]]))
      }
    }
  }, by = bponly]$d
  dists <- dists[dists > 0]
  if (length(dists) < 20) return(list(stats = NULL, sim = sim, elapsed = elapsed))
  fit <- fit_2gauss(log10(dists))
  if (is.null(fit)) return(list(stats = NULL, sim = sim, elapsed = elapsed))
  ct <- tryCatch(counts_long_nogap(mono$bponly), error = function(e) NULL)
  if (is.null(ct)) return(list(stats = NULL, sim = sim, elapsed = elapsed))
  cons <- ct[symbol == consensus][order(pos)]
  load <- mu_load_from_consensus(mono$bponly, cons)
  stats <- data.table(log_mode1 = fit$log_mode1, log_mode2 = fit$log_mode2,
                      weight1 = fit$weight1, mean_load = mean(load),
                      n_units = nrow(mono), n_gens = n_gens)
  list(stats = stats, sim = sim, elapsed = elapsed)
}

compute_distance <- function(sim_stats, target, stat_scales = NULL) {
  diffs <- c(
    log_mode1 = sim_stats$log_mode1 - target$log_mode1,
    log_mode2 = sim_stats$log_mode2 - target$log_mode2,
    weight1   = sim_stats$weight1   - target$weight1,
    # Use relative error for mean_load so it's on a comparable scale to
    # the log-space mode differences (both are dimensionless ratios)
    mean_load = (sim_stats$mean_load - target$mean_load) / max(target$mean_load, 0.1)
  )
  # Normalize by SD of each stat (computed from Gen 0 pilot batch).
  # Before scaling is available, use raw differences.
  if (!is.null(stat_scales)) {
    diffs <- diffs / stat_scales
  }
  sqrt(sum(diffs^2))
}

#' Compute normalization scales from a set of valid particles.
#' Returns the SD of each summary stat, with a floor to avoid division by zero.
compute_stat_scales <- function(particles, target) {
  # Compute diffs the same way as compute_distance, then take SD
  rel_load <- (particles$mean_load - target$mean_load) / max(target$mean_load, 0.1)
  scales <- c(
    log_mode1 = sd(particles$log_mode1, na.rm = TRUE),
    log_mode2 = sd(particles$log_mode2, na.rm = TRUE),
    weight1   = sd(particles$weight1, na.rm = TRUE),
    mean_load = sd(rel_load, na.rm = TRUE)
  )
  # Floor: if a stat has zero variance, use 1 (no scaling)
  scales[is.na(scales) | scales < 1e-10] <- 1
  scales
}

perturb_particle <- function(particle, sd_scale) {
  new <- copy(particle)
  for (j in seq_along(PARAMS)) {
    pname <- PARAMS[[j]]$name
    range_j <- PARAMS[[j]]$hi - PARAMS[[j]]$lo
    new[[pname]] <- new[[pname]] + rnorm(1, 0, sd_scale * range_j)
    # No hard clamp — allow exploration beyond the initial prior range.
    # Only enforce physical constraints (e.g. shape params must be positive).
  }
  # Enforce physical constraints: shape params > 0, scale params > 0
  for (pname in c("local_shape", "distal_shape", "del_shape"))
    new[[pname]] <- max(0.01, new[[pname]])
  for (pname in c("local_scale", "distal_scale", "del_scale"))
    new[[pname]] <- max(0.1, new[[pname]])
  new
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
        numericInput("target_mode1", "Distance mode 1", 100, min = 1),
        numericInput("target_mode2", "Distance mode 2", 2000, min = 10),
        numericInput("target_weight1", "Weight of mode 1", 0.3,
                     min = 0, max = 1, step = 0.05),
        numericInput("target_load", "Mean mutation load", 10, min = 0,
                     step = 1)
      ),

      h4("Algorithm"),
      numericInput("n_particles", "Particles per generation", 40,
                   min = 10, step = 10),
      numericInput("max_generations", "Max generations", 10, min = 1),
      sliderInput("retention_frac", "Retention fraction", 0.5,
                  min = 0.1, max = 0.9, step = 0.05),
      sliderInput("perturbation_sd", "Perturbation scale", 0.3,
                  min = 0.05, max = 0.5, step = 0.05),

      h4("Simulation"),
      numericInput("max_units", "Max array units", 2000, min = 500,
                   step = 1000),
      numericInput("max_t", "Max generations per sim", 2000,
                   min = 100, step = 500),
      numericInput("sim_timeout", "Per-particle timeout (sec)", 10,
                   min = 1, max = 120, step = 1)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("Live Dashboard",
          fluidRow(
            column(6, plotOutput("convergence_plot", height = "280px")),
            column(6, plotOutput("best_stats_plot", height = "280px"))
          ),
          fluidRow(
            column(12, plotOutput("param_posterior_plot", height = "300px"))
          ),
          fluidRow(
            column(12, h4("Current generation particles"),
                   tableOutput("particles_table"))
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
                               n_failed = integer(), n_timeout = integer()),
    all_particles = list(),     # all gens' particles for posteriors
    best_sim      = NULL,       # stored sim result for best particle
    sd_scale      = 0.3,
    stall_count   = 0L,
    last_particle_time = NA_real_,  # seconds taken by last particle
    n_failed      = 0L,             # failed attempts this generation
    n_timeout     = 0L,             # timed-out attempts this generation
    n_attempts    = 0L,             # total attempts this generation
    gen_start_time = NA_real_,      # wall-clock start of current generation
    stat_scales   = NULL            # SD of each summary stat (from Gen 0 pilot)
  )

  get_target <- reactive({
    list(log_mode1 = log10(input$target_mode1),
         log_mode2 = log10(input$target_mode2),
         weight1   = input$target_weight1,
         mean_load = input$target_load)
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
                                     n_failed = integer(), n_timeout = integer())
    state$all_particles <- list()
    state$best_sim <- NULL
    state$sd_scale <- input$perturbation_sd
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

      # Process a batch of particles (up to 2 seconds of wall time per batch)
      batch_start <- proc.time()[3]
      batch_budget <- 2  # seconds before yielding to UI

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
                .SDcols = c("log_mode1", "log_mode2", "weight1", "mean_load")]
            }
            state$particles <- particles
          }

          best_d <- min(particles$distance)
          med_d  <- median(particles$distance)

          state$convergence <- rbind(state$convergence,
            data.table(gen = gen, best = best_d, median = med_d,
                       n_valid = n_valid,
                       n_failed = state$n_failed,
                       n_timeout = state$n_timeout))
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

          # Prepare next generation
          state$generation <- gen + 1L
          state$particle_idx <- 0L
          state$particles <- NULL
          state$sd_scale <- input$perturbation_sd * (0.7 ^ gen)
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
          proposed <- perturb_particle(parent, state$sd_scale)
        }

        # Run simulation with wall-clock time limit (no forking).
        res <- tryCatch(
          {
            setTimeLimit(elapsed = timeout, transient = TRUE)
            run_and_summarize_one(
              proposed,
              max_units   = input$max_units,
              max_t       = input$max_t,
              init_l      = 30L,
              init_k0     = 10L,
              sim_timeout = timeout,
              target      = target
            )
          },
          error = function(e) {
            list(stats = NULL, sim = NULL, elapsed = timeout)
          },
          finally = {
            setTimeLimit(elapsed = Inf, transient = TRUE)
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

    # Two-panel: convergence + failure rate
    # IMPORTANT: copy() to avoid modifying state$convergence by reference
    conv <- copy(conv)
    conv[, fail_rate := (n_failed + n_timeout) /
           (n_valid + n_failed + n_timeout)]

    p1 <- ggplot(conv, aes(x = gen)) +
      geom_line(aes(y = median, colour = "Median"), linewidth = 1) +
      geom_point(aes(y = median, colour = "Median"), size = 2.5) +
      geom_line(aes(y = best, colour = "Best"), linewidth = 1) +
      geom_point(aes(y = best, colour = "Best"), size = 2.5) +
      scale_colour_manual(values = c("Median" = "steelblue",
                                      "Best" = "firebrick")) +
      theme_classic(base_size = 12) +
      labs(x = "Generation", y = "Distance",
           colour = NULL, title = "Convergence") +
      theme(legend.position = c(0.8, 0.8))

    p2 <- ggplot(conv, aes(x = gen, y = fail_rate)) +
      geom_col(fill = "grey40", alpha = 0.7, width = 0.6) +
      geom_text(aes(label = n_failed + n_timeout), vjust = -0.3, size = 3) +
      scale_y_continuous(labels = scales::percent_format(),
                         limits = c(0, max(0.1, max(conv$fail_rate) * 1.2))) +
      theme_classic(base_size = 12) +
      labs(x = "Generation", y = "Failure rate",
           title = "Sim failures per generation")

    cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1))
  }, res = 100)

  # --- Best particle stats vs targets ---
  output$best_stats_plot <- renderPlot({
    conv <- state$convergence
    all_p <- state$all_particles
    if (length(all_p) == 0) {
      return(ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "No results yet",
                 size = 6, colour = "grey50"))
    }

    target <- get_target()
    last <- all_p[[length(all_p)]]
    valid <- last[is.finite(distance)][order(distance)]
    if (nrow(valid) == 0) return(ggplot() + theme_void())
    best <- valid[1]

    stats_dt <- data.table(
      stat    = c("log(mode1)", "log(mode2)", "weight1", "mean_load"),
      target  = c(target$log_mode1, target$log_mode2,
                  target$weight1, target$mean_load),
      current = c(best$log_mode1, best$log_mode2,
                  best$weight1, best$mean_load)
    )
    stats_long <- melt(stats_dt, id.vars = "stat",
                       variable.name = "type", value.name = "value")

    ggplot(stats_long, aes(x = stat, y = value, fill = type)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c(target = "grey60", current = "firebrick"),
                        labels = c("Target", "Best fit")) +
      theme_classic(base_size = 13) +
      labs(x = NULL, y = "Value", fill = NULL,
           title = paste0("Best fit (dist = ",
                          round(best$distance, 3), ")")) +
      theme(legend.position = c(0.15, 0.85))
  }, res = 100)

  # --- Parameter posteriors ---
  output$param_posterior_plot <- renderPlot({
    all_p <- state$all_particles
    if (length(all_p) == 0) {
      return(ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "Parameter posteriors will appear here",
                 size = 6, colour = "grey50"))
    }

    # Show top 50% of most recent generation
    last <- all_p[[length(all_p)]]
    valid <- last[is.finite(distance)][order(distance)]
    if (nrow(valid) < 3) return(ggplot() + theme_void())
    n_show <- max(3, as.integer(nrow(valid) * 0.5))
    top <- valid[seq_len(n_show)]

    # Transform to natural scale for the 4 rate params
    display <- copy(top)
    for (j in 1:4) {
      pn <- PARAMS[[j]]$name
      set(display, j = pn, value = 10^display[[pn]])
    }

    # Melt for faceted plot — only show the 4 key rates
    key_params <- c("p_local_dup", "p_distal_dup", "p_del_chunk", "mu_total")
    melt_dt <- melt(display[, ..key_params], measure.vars = key_params,
                    variable.name = "param", value.name = "value")

    ggplot(melt_dt, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 20, fill = "steelblue", alpha = 0.6) +
      geom_density(linewidth = 0.5) +
      facet_wrap(~ param, scales = "free", nrow = 1) +
      scale_x_log10() +
      theme_classic(base_size = 11) +
      labs(x = NULL, y = "Density",
           title = paste0("Parameter posteriors (top 50%, gen ",
                          length(all_p) - 1, ")"))
  }, res = 100)

  # --- Particles table (top 8 by distance) ---
  output$particles_table <- renderTable({
    p <- state$particles
    if (is.null(p) || nrow(p) == 0) return(NULL)

    show <- p[order(distance)][seq_len(min(8, nrow(p)))]

    data.frame(
      rank      = seq_len(nrow(show)),
      distance  = round(show$distance, 3),
      mode1     = round(10^show$log_mode1),
      mode2     = round(10^show$log_mode2),
      weight1   = round(show$weight1, 2),
      mean_load = round(show$mean_load, 1),
      units     = show$n_units,
      gens      = show$n_gens,
      time_s    = round(show$elapsed, 1)
    )
  }, digits = 3)

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
      rank      = seq_len(n_show),
      distance  = round(top$distance, 3),
      mode1     = round(10^top$log_mode1),
      mode2     = round(10^top$log_mode2),
      mean_load = round(top$mean_load, 1)
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
      sim_args <- params_to_args(row)
      sim <- suppressWarnings(run_sim_ps(
        max_units = input$max_units, init_l = 30, init_k0 = 10,
        p_local_dup = sim_args$p_local_dup,
        p_distal_dup = sim_args$p_distal_dup,
        p_del_chunk = sim_args$p_del_chunk,
        mu_total = sim_args$mu_total,
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
