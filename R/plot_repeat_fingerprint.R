# Internal implementation (not exported) — does the actual ggplot work.
# Takes ps (data.table of pairs) and kb (monomer data.table), both already
# extracted from the simulation object.
.fingerprint_impl <- function(ps, kb = NULL, rotate = FALSE, zoom = NULL,
                              ptsize = 0.001, gridalpha = 1,
                              gbreaks = 20, ptcol = "black") {

  ps[, `:=`(num1 = as.numeric(num1), num2 = as.numeric(num2))]

  if (!is.null(zoom))
    ps <- ps[num1 > zoom[1] & num1 < zoom[2] &
               num2 > zoom[1] & num2 < zoom[2]]

  break_fun  <- int_breaks_pretty(gbreaks)
  lims       <- if (!is.null(zoom)) zoom else range(c(ps$num1, ps$num2), na.rm = TRUE)
  grid_breaks <- break_fun(lims)
  min_lim     <- min(lims)
  max_lim     <- max(lims)

  grid_vert <- data.frame(x = grid_breaks, xend = grid_breaks,
                           y = max_lim,    yend = grid_breaks)
  grid_horiz <- data.frame(x = min_lim, xend = grid_breaks,
                            y = grid_breaks, yend = grid_breaks)

  grid_bp <- seq(ceiling(min_lim), floor(max_lim), by = 1)
  grid_vertbp  <- data.frame(x = grid_bp, xend = grid_bp,
                              y = max_lim, yend = grid_bp)
  grid_horizbp <- data.frame(x = min_lim, xend = grid_bp,
                              y = grid_bp, yend = grid_bp)

  p <- ggplot(ps[num1 != num2], aes(x = num1, y = num2)) +
    geom_segment(data = grid_vertbp,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE,
                 colour = "grey90", linewidth = 0.2, alpha = gridalpha) +
    geom_segment(data = grid_horizbp,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE,
                 colour = "grey90", linewidth = 0.2, alpha = gridalpha) +
    geom_point(size = ptsize, shape = 15, alpha = 1, col = ptcol) +
    coord_fixed() +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(expand = c(0, 0), position = "top",
                       limits = lims, breaks = grid_breaks) +
    scale_y_continuous(expand = c(0, 0), position = "left",
                       limits = lims, breaks = grid_breaks) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 0),
          axis.text.y      = element_text(angle = 45, hjust = 1),
          axis.line        = element_line(),
          legend.position  = "none",
          plot.background  = element_blank(),
          panel.background = element_blank())

  if (is.null(zoom) && !is.null(kb)) {
    p <- p +
      geom_point(data = kb,
                 aes(x = num, y = max(kb$num, na.rm = TRUE) + 100, color = dir),
                 size = 0.1, inherit.aes = FALSE, shape = 15) +
      geom_point(data = kb,
                 aes(x = -100, y = num, color = dir),
                 size = 0.1, inherit.aes = FALSE, shape = 15)
  }

  if (rotate) {
    g <- ggplotGrob(p)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(angle = -45, width = 0.75,
                                     height = 0.75, clip = "off"))
    grid::grid.draw(g)
    grid::popViewport()
  } else {
    p
  }
}

#' Dot-plot fingerprint of identical repeat-unit pairs
#'
#' Draws a dot-matrix ("fingerprint") plot in which each dot represents a pair
#' of identical monomers, placed at their respective positions along the array
#' axes. This visualisation reveals the duplication history of the array:
#' diagonal streaks near the main diagonal indicate local (tandem) duplications,
#' while off-diagonal blocks indicate distal duplications.
#'
#' @param sim Output from \code{\link{run_sim_ps}}.
#' @param rotate Logical. If \code{TRUE}, rotate the plot 45 degrees using grid
#'   graphics (draws directly to the active device and returns \code{NULL}). If
#'   \code{FALSE} (default), return a \code{ggplot2} object.
#' @param zoom Numeric vector of length 2 \code{c(min, max)} defining a
#'   sub-region of the array to zoom into. \code{NULL} (default) shows the
#'   full array.
#' @param ptsize Numeric. Point size. Default 0.01.
#' @param gridalpha Numeric. Alpha for the fine grid lines. Default 0 (hidden).
#' @param gbreaks Integer. Approximate number of major axis breaks. Default 20.
#' @param ptcol Character. Point colour. Default \code{"black"}.
#'
#' @return A \code{ggplot2} object when \code{rotate = FALSE}; draws to the
#'   active graphics device (invisibly \code{NULL}) when \code{rotate = TRUE}.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(max_t = 500, mu_total = 3e-5)
#' plot_repeat_fingerprint(sim)
#' plot_repeat_fingerprint(sim, zoom = c(0, 200))
#' }
plot_repeat_fingerprint <- function(sim,
                                    rotate   = FALSE,
                                    zoom     = NULL,
                                    ptsize   = 0.01,
                                    gridalpha = 0,
                                    gbreaks  = 20,
                                    ptcol    = "black") {
  ps    <- as.data.table(sim$ps)
  monos <- as.data.table(sim$monomers)
  # Compute pairs on the fly if sim was run with compute_pairs = FALSE
  if (nrow(ps) == 0 && nrow(monos) > 1) ps <- pairs_identical(monos)

  .fingerprint_impl(ps, kb = monos,
                    rotate    = rotate,
                    zoom      = zoom,
                    ptsize    = ptsize,
                    gridalpha = gridalpha,
                    gbreaks   = gbreaks,
                    ptcol     = ptcol)
}
