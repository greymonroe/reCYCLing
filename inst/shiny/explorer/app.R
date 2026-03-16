library(shiny)
library(reCYCLing)
library(ggplot2)
library(data.table)
library(cowplot)

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .sidebar { background: #f8f9fa; padding: 15px; border-radius: 8px; }
    .well { background: #f8f9fa; }
    h4 { margin-top: 15px; margin-bottom: 5px; color: #555; }
    .btn-primary { width: 100%; margin-top: 10px; margin-bottom: 10px;
                   font-size: 16px; padding: 10px; }
    .sim-time { color: #888; font-style: italic; margin-bottom: 10px; }
    .mode-box { background: #e8f0fe; padding: 8px; border-radius: 5px;
                margin-bottom: 10px; }
  "))),

  titlePanel("reCYCLing Simulation Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      actionButton("run", "Run Simulation", class = "btn-primary"),
      textOutput("sim_time"),

      h4("Stopping rule"),
      div(class = "mode-box",
        radioButtons("stop_mode", NULL,
          choices = c("Grow to size" = "size", "Fixed generations (equilibrium)" = "gens"),
          selected = "gens"),
        conditionalPanel("input.stop_mode == 'size'",
          numericInput("max_units", "Max units", 10000, min = 100, step = 1000)
        ),
        conditionalPanel("input.stop_mode == 'gens'",
          numericInput("max_t_fixed", "Generations", 1000, min = 100, step = 100)
        )
      ),

      h4("Array"),
      numericInput("init_l", tags$span("Monomer length (bp)",
        title = "Length of one repeat unit in base pairs."), 178, min = 10, max = 500),
      numericInput("init_k0", tags$span("Initial copies",
        title = "Starting number of identical monomers in the array."), 10, min = 1, max = 100),

      h4("Duplication rates"),
      sliderInput("p_local_dup", tags$span("Local dup (per unit, log10)",
        title = "Per-unit per-generation probability of tandem duplication. Each unit independently triggers a copy-paste of a nearby chunk."),
                  min = -6, max = -2, value = -4, step = 0.1),
      sliderInput("p_distal_dup", tags$span("Distal dup (per array, log10)",
        title = "Per-array per-generation probability of a distal duplication event (e.g., unequal crossover). At most one event per generation."),
                  min = -3, max = 0, value = log10(0.05), step = 0.1),
      sliderInput("p_invert_distal", tags$span("P(invert distal chunk)",
        title = "Probability that a distal-duplicated chunk is inverted before insertion."),
                  min = 0, max = 1, value = 0.5, step = 0.05),

      h4("Deletion"),
      sliderInput("p_del_chunk", tags$span("Deletion (per unit, log10)",
        title = "Per-unit per-generation probability of deleting a chunk of consecutive monomers."),
                  min = -7, max = -2, value = log10(7e-5), step = 0.1),

      h4("Mutation"),
      sliderInput("mu_total", tags$span("Mutation rate (per base, log10)",
        title = "Per-base per-generation substitution rate. Combined with monomer length and generations, this determines mutation load."),
                  min = -7, max = -3, value = log10(5.6e-5), step = 0.1),

      h4("Chunk-size distributions"),
      helpText("Gamma(shape, scale) for duplication/deletion chunk sizes"),
      fluidRow(
        column(6, numericInput("local_shape", "Local shape", 2,
                               min = 0.1, step = 0.5)),
        column(6, numericInput("local_scale", "Local scale", 15,
                               min = 1, step = 5))
      ),
      fluidRow(
        column(6, numericInput("distal_shape", "Distal shape", 2,
                               min = 0.1, step = 0.5)),
        column(6, numericInput("distal_scale", "Distal scale", 500,
                               min = 10, step = 50))
      ),
      fluidRow(
        column(6, numericInput("del_shape", "Del shape", 2,
                               min = 0.1, step = 0.5)),
        column(6, numericInput("del_scale", "Del scale", 15,
                               min = 1, step = 5))
      ),

      h4("Zoom (fingerprint)"),
      checkboxInput("use_zoom", "Enable zoom", FALSE),
      fluidRow(
        column(6, numericInput("zoom_min", "From", 0, min = 0)),
        column(6, numericInput("zoom_max", "To", 1000, min = 1))
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("Summary",
                 plotOutput("summary_plot", height = "900px")),
        tabPanel("Trajectories & Load",
                 plotOutput("traj_load_plot", height = "700px")),
        tabPanel("Fingerprint",
                 plotOutput("fingerprint_plot", height = "700px")),
        tabPanel("Pair Distances",
                 plotOutput("distances_plot", height = "500px")),
        tabPanel("Copy Numbers",
                 plotOutput("copynums_plot", height = "500px")),
        tabPanel("Mutation Load",
                 plotOutput("mutload_plot", height = "500px")),
        tabPanel("Consensus Support",
                 plotOutput("consensus_plot", height = "500px")),
        tabPanel("Stats",
                 verbatimTextOutput("stats_text"))
      )
    )
  )
)

# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {

  # Helper to run a simulation from the current input values
  run_simulation <- function() {
    t0 <- proc.time()

    # Determine stopping rule
    if (input$stop_mode == "size") {
      max_u <- input$max_units
      max_t <- 1000000L  # effectively unlimited
    } else {
      max_u <- 200000L   # very high ceiling — size is emergent
      max_t <- input$max_t_fixed
    }

    sim <- withProgress(message = "Running simulation...", value = 0, {
      incProgress(0.1, detail = "Evolving array")
      result <- run_sim_ps(
        max_units       = max_u,
        max_t           = max_t,
        hard_cap        = 200000L,
        init_l          = input$init_l,
        init_k0         = input$init_k0,
        p_local_dup     = 10^input$p_local_dup,
        p_distal_dup    = 10^input$p_distal_dup,
        p_invert_distal = input$p_invert_distal,
        p_del_chunk     = 10^input$p_del_chunk,
        mu_total        = 10^input$mu_total,
        local_dist      = list(type = "gamma",
                               shape = input$local_shape,
                               scale = input$local_scale),
        distal_dist     = list(type = "gamma",
                               shape = input$distal_shape,
                               scale = input$distal_scale),
        del_dist        = list(type = "gamma",
                               shape = input$del_shape,
                               scale = input$del_scale),
        verbose         = FALSE
      )
      incProgress(0.7, detail = paste0(
        nrow(result$monomers), " monomers — rendering plots..."))
      result
    })

    elapsed <- round((proc.time() - t0)[3], 2)
    sim_result(sim)
    sim_elapsed(elapsed)
  }

  sim_result <- reactiveVal(NULL)
  sim_elapsed <- reactiveVal(NULL)

  # Run on button click
  observeEvent(input$run, run_simulation())

  # Auto-run once on startup
  auto_ran <- reactiveVal(FALSE)
  observe({
    if (!auto_ran()) {
      auto_ran(TRUE)
      run_simulation()
    }
  })

  output$sim_time <- renderText({
    el <- sim_elapsed()
    if (is.null(el)) return("")
    sim <- sim_result()
    paste0("Completed in ", el, "s | ",
           nrow(sim$monomers), " monomers | ",
           length(sim$L_vec), " generations")
  })

  # --- Zoom helper ---
  get_zoom <- reactive({
    if (input$use_zoom) c(input$zoom_min, input$zoom_max) else NULL
  })

  # --- Plots ---

  output$summary_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    plot_ps_summary(sim, title = "Simulation Summary")
  }, res = 120)

  output$traj_load_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    mono <- sim$monomers

    # Size trajectory
    traj <- data.table(gen = seq_along(sim$L_vec), size = sim$L_vec)
    p1 <- ggplot(traj, aes(x = gen, y = size)) +
      geom_line(colour = "steelblue", linewidth = 0.5) +
      theme_classic(base_size = 12) +
      labs(x = "Generation", y = "Array size", title = "Array size trajectory")

    # Mutation load distribution
    ct <- tryCatch(counts_long_nogap(mono$bponly), error = function(e) NULL)
    if (!is.null(ct)) {
      cons <- ct[symbol == consensus][order(pos)]
      load <- mu_load_from_consensus(mono$bponly, cons)
      p2 <- ggplot(data.frame(load = load), aes(x = load)) +
        geom_histogram(bins = 40, fill = "firebrick", alpha = 0.7, colour = "white") +
        geom_vline(xintercept = mean(load), linetype = "dashed") +
        theme_classic(base_size = 12) +
        labs(x = "Mismatches to consensus", y = "Count",
             title = paste0("Mutation load (mean = ", round(mean(load), 1), ")"))
    } else {
      p2 <- ggplot() + theme_void() + annotate("text", x=.5, y=.5, label="N/A")
    }

    # Copy number distribution
    cn <- mono[, .N, by = bponly]
    p3 <- ggplot(cn, aes(x = N)) +
      geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7, colour = "white") +
      scale_x_log10() +
      theme_classic(base_size = 12) +
      labs(x = "Copies per unique sequence", y = "Count",
           title = paste0("Copy numbers (", nrow(cn), " unique / ",
                          nrow(mono), " total)"))

    # Info panel
    n_gens <- length(sim$L_vec)
    info_text <- paste0(
      "Generations: ", n_gens, "\n",
      "Final size: ", nrow(mono), " monomers\n",
      "Unique sequences: ", uniqueN(mono$bponly), "\n",
      if (!is.null(ct)) paste0("Mean load: ", round(mean(load), 1), "\n") else "",
      "Mean copy number: ", round(mean(cn$N), 1)
    )
    p4 <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = info_text,
               size = 4.5, hjust = 0.5, family = "mono") +
      ggtitle("Summary statistics")

    plot_grid(p1, p2, p3, p4, ncol = 2)
  }, res = 120)

  output$fingerprint_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    zoom <- get_zoom()
    ptsize <- if (!is.null(zoom)) 0.3 else 0.01
    plot_repeat_fingerprint(sim, zoom = zoom, ptsize = ptsize,
                            gridalpha = if (!is.null(zoom)) 0.3 else 0)
  }, res = 120)

  output$distances_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    plot_pair_distances(sim)
  }, res = 120)

  output$copynums_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    plot_copy_numbers(sim)
  }, res = 120)

  output$mutload_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    plot_mutation_load(sim)
  }, res = 120)

  output$consensus_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    plot_consensus_support(sim)
  }, res = 120)

  output$stats_text <- renderPrint({
    sim <- sim_result()
    req(sim)
    mono <- sim$monomers
    cat("=== Simulation Statistics ===\n\n")
    cat("Final array size:", nrow(mono), "monomers\n")
    cat("Generations:", length(sim$L_vec), "\n")
    cat("Unique sequences:", uniqueN(mono$bponly), "\n")
    cat("Shannon entropy:", round(shannon_entropy(mono$bponly), 3), "bits\n")
    cat("Identical pairs:", nrow(sim$ps), "\n\n")

    ct <- counts_long_nogap(mono$bponly)
    cons <- ct[symbol == consensus][order(pos)]
    load <- mu_load_from_consensus(mono$bponly, cons)
    cat("Mean mutation load:", round(mean(load), 2), "mismatches\n")
    cat("Median mutation load:", median(load), "\n\n")

    cn <- as.integer(table(mono$bponly))
    cat("Copy number distribution:\n")
    cat("  Mean:", round(mean(cn), 1), "\n")
    cat("  Max:", max(cn), "\n")
    cat("  Singletons:", sum(cn == 1), "/", length(cn), "unique sequences\n\n")

    dirs <- table(mono$dir)
    cat("Strand orientation:\n")
    for (d in names(dirs)) cat("  ", d, ":", dirs[d], "\n")
  })

}

# ============================================================================
# Launch
# ============================================================================

shinyApp(ui, server)
