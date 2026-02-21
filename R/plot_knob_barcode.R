#' Barcode visualisation of identical repeat units in a zoom window
#'
#' Draws a linear barcode in which each monomer position in the zoom window
#' is represented by a coloured segment. Monomers sharing the same sequence
#' receive the same colour, making it easy to see the pattern of identical
#' copies within a sub-region of the array.
#'
#' @param ps Pairwise-identical \code{data.table} (such as an element of
#'   \code{ps_results$ps_list} from \code{\link{run_sim_ps}}, or the output
#'   of \code{\link{pairs_identical}}).
#' @param zoom Numeric vector of length 2 \code{c(min, max)} defining the
#'   window of monomer positions to display.
#' @param palette Character. \pkg{RColorBrewer} palette name for colouring
#'   identical-sequence groups. Default \code{"Dark2"}.
#' @param line_width Numeric. Segment line width. Default 2.
#' @param base_size Numeric. Base font size for \code{theme_void}. Default 6.
#' @param label_size Numeric. Size of position-number text labels. Default 2.
#' @param show_labels Logical. Whether to show numeric position labels.
#'   Default \code{TRUE}.
#' @param gapsize Numeric. Fraction of the segment width to use as the gap
#'   between segments. Default 0.9.
#' @param vjusttext Numeric. Vertical justification for position labels.
#'   Default -0.1.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' ps  <- sim$ps_list[[1]]
#' plot_knob_barcode(ps, zoom = c(1, 100))
#' }
plot_knob_barcode <- function(ps,
                              zoom,
                              palette    = "Dark2",
                              line_width = 2,
                              base_size  = 6,
                              label_size = 2,
                              show_labels = TRUE,
                              gapsize = 0.9, vjusttext = -0.1) {
  stopifnot(length(zoom) == 2)

  ps_zoom <- ps[num1 > zoom[1] & num1 < zoom[2] &
                  num2 > zoom[1] & num2 < zoom[2], ]
  ps_zoom <- as.data.frame(ps_zoom)

  if (nrow(ps_zoom) == 0) {
    return(
      ggplot() +
        theme_void(base_size = base_size) +
        ggtitle("No pairs in zoom window")
    )
  }

  ps_zoom$factorlevel <- factor(ps_zoom$bponly)

  range_x <- seq(ceiling(zoom[1]), floor(zoom[2]))
  if (length(range_x) == 0) range_x <- zoom[1]:zoom[2]

  gap <- gapsize * (zoom[2] - zoom[1]) / length(range_x)

  bg <- data.frame(x = range_x, xend = range_x + gap, y = 1, yend = 1)

  lab_nums <- sort(unique(c(ps_zoom$num1, ps_zoom$num2)))
  labs <- data.frame(x = lab_nums + gap / 2, y = 1, lab = lab_nums)

  p <- ggplot() +
    geom_segment(data = bg,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE,
                 linewidth = line_width, lineend = "round",
                 colour = "gray85") +
    geom_segment(data = ps_zoom,
                 aes(x = num1, xend = num1 + gap, y = 1, yend = 1,
                     colour = factorlevel),
                 linewidth = line_width, lineend = "round") +
    geom_segment(data = ps_zoom,
                 aes(x = num2, xend = num2 + gap, y = 1, yend = 1,
                     colour = factorlevel),
                 linewidth = line_width, lineend = "round") +
    scale_color_brewer(palette = palette) +
    coord_cartesian(xlim = zoom, ylim = c(0.95, 1.25), expand = FALSE) +
    theme_void(base_size = base_size) +
    theme(legend.position = "none")

  if (show_labels) {
    p <- p +
      geom_text(data = labs,
                aes(x = x, y = y, label = lab),
                inherit.aes = FALSE,
                size = label_size, angle = 90,
                vjust = vjusttext, hjust = 1.2)
  }
  p
}
