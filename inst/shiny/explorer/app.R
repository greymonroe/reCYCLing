library(shiny)
library(reCYCLing)
library(ggplot2)
library(data.table)

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
  "))),

  titlePanel("reCYCLing Simulation Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      actionButton("run", "Run Simulation", class = "btn-primary"),
      textOutput("sim_time"),

      h4("Array"),
      numericInput("max_units", "Max units", 10000, min = 100, step = 1000),
      numericInput("init_l", "Monomer length (bp)", 178, min = 10, max = 500),
      numericInput("init_k0", "Initial copies", 10, min = 1, max = 100),

      h4("Duplication rates"),
      sliderInput("p_local_dup", "Local dup prob (per unit)",
                  min = -5, max = -1, value = log10(0.0004), step = 0.1,
                  post = " (log10)"),
      sliderInput("p_distal_dup", "Distal dup prob (per unit)",
                  min = -7, max = -2, value = log10(1e-5), step = 0.1,
                  post = " (log10)"),
      sliderInput("p_invert_distal", "P(invert distal chunk)",
                  min = 0, max = 1, value = 0.5, step = 0.05),

      h4("Deletion"),
      sliderInput("p_del_chunk", "Deletion prob (per unit)",
                  min = -7, max = -2, value = log10(1e-5), step = 0.1,
                  post = " (log10)"),

      h4("Mutation"),
      sliderInput("mu_total", "Mutation rate (per base)",
                  min = -6, max = -2, value = log10(5e-5), step = 0.1,
                  post = " (log10)"),

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
        tabPanel("Trajectories",
                 plotOutput("traj_plot", height = "500px")),
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

    sim <- withProgress(message = "Running simulation...", value = 0, {
      incProgress(0.1, detail = "Evolving array")
      result <- run_sim_ps(
        max_units       = input$max_units,
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

  output$traj_plot <- renderPlot({
    sim <- sim_result()
    req(sim)
    traj <- data.frame(gen = seq_along(sim$L_vec), size = sim$L_vec)
    ggplot(traj, aes(x = gen, y = size)) +
      geom_line(colour = "steelblue", linewidth = 0.5) +
      theme_classic(base_size = 14) +
      xlab("Generation") +
      ylab("Array size (monomers)")
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
