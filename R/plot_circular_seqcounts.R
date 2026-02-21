#' Circular visualisation of per-position base composition
#'
#' Draws a polar (circular) tile plot in which each alignment position is a
#' column of coloured tiles, one tile per allowed nucleotide symbol. Tile
#' transparency (alpha) is proportional to the observed proportion of that
#' symbol at that position, so positions with a strong consensus appear as a
#' single opaque tile while variable positions show a blend of colours.
#'
#' @param seqcounts Output from \code{\link{counts_long_nogap}}: a long-format
#'   \code{data.table} with columns \code{pos_new}, \code{symbol},
#'   \code{prop}, and \code{consensus}.
#' @param allowed Character vector of symbols to display. Default
#'   \code{c("-","A","T","G","C")}.
#' @param title Character. Plot title. Default \code{""}.
#' @param palette Character. \pkg{RColorBrewer} palette name. Default
#'   \code{"Dark2"}.
#' @param show_consensus_labels Logical. If \code{TRUE}, add the consensus
#'   base letter outside each position. Default \code{FALSE}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples
#' \dontrun{
#' sim  <- run_sim_ps(n = 1, max_t = 500, mu_total = 3e-5)
#' mono <- sim$monomers_list[[1]]
#' ct   <- counts_long_nogap(mono$bponly)
#' plot_circular_seqcounts(ct)
#' }
plot_circular_seqcounts <- function(seqcounts,
                                    allowed = c("-", "A", "T", "G", "C"),
                                    title   = "",
                                    palette = "Dark2",
                                    show_consensus_labels = FALSE) {
  dt <- as.data.table(seqcounts)

  needed <- c("pos_new", "symbol", "prop", "consensus")
  stopifnot(all(needed %in% names(dt)))

  dt[, symbol := factor(as.character(symbol), levels = allowed)]

  consensus <- dt[, .(consensus = unique(consensus)), by = pos_new][order(pos_new)]
  max_x     <- max(consensus$pos_new)
  y_levels  <- c(levels(dt$symbol), "CONS")

  cons_lab <- consensus[, {
    ang  <- (pos_new - 0.5) / max_x * 360
    ang2 <- ifelse(ang > 90 & ang < 270, ang + 180, ang)
    hj   <- ifelse(ang > 90 & ang < 270, 1, 0)
    .(pos_new = pos_new,
      symbol  = factor("CONS", levels = y_levels),
      label   = consensus,
      angle   = ang2,
      hjust   = hj)
  }]

  p <- ggplot(dt, aes(x = pos_new, y = symbol)) +
    geom_tile(fill = "white", linewidth = 0.1, width = 1, height = 0.75,
              col = "gray90") +
    geom_tile(fill = "white", linewidth = 0.1, width = 1.1, height = 0.6,
              col = "white") +
    geom_tile(
      aes(fill  = symbol, alpha = prop,
          colour = after_scale(scales::alpha(fill, alpha))),
      linewidth = 0.1, width = 1, height = 0.75
    ) +
    coord_polar() +
    scale_alpha(range = c(0, 0.8), guide = "none") +
    scale_x_continuous(expand = c(0.002, 0.002)) +
    scale_y_discrete(drop = FALSE, limits = y_levels, expand = c(3, 2)) +
    scale_fill_brewer(palette = palette, name = "") +
    scale_color_brewer(palette = palette, name = "") +
    theme_void(base_size = 6) +
    theme(legend.position = c(0.5, 0.5),
          legend.key.size = grid::unit(0.1, "cm"),
          panel.grid      = element_blank(),
          plot.margin     = grid::unit(c(0, 0, 0, 0), "pt")) +
    ggtitle(title)

  if (isTRUE(show_consensus_labels)) {
    p <- p +
      geom_text(data = cons_lab,
                aes(x = pos_new, y = symbol, label = label,
                    angle = angle, hjust = hjust),
                inherit.aes = FALSE, size = 1.5)
  }

  p + geom_segment(x = 0, xend = 0, y = 0, yend = length(y_levels),
                   linewidth = 0.1, col = "gray")
}
