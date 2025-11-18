
library(data.table)
library(patchwork)
library(ggplot2)
# Initialize a tandem array of repeat units
init_array <- function(l      = 178,
                       k0     = 10,
                       init_sequence_type = c("random", "given"),
                       ancestor_seq       = NULL,
                       alphabet = c("A", "C", "G", "T")) {
  init_sequence_type <- match.arg(init_sequence_type)

  if (init_sequence_type == "given") {
    if (is.null(ancestor_seq))
      stop("ancestor_seq must be provided when init_sequence_type = 'given'")
    if (nchar(ancestor_seq) != l)
      stop("ancestor_seq must have length l = ", l)
    unit <- ancestor_seq
  } else {
    unit <- paste0(sample(alphabet, size = l, replace = TRUE), collapse = "")
  }

  # k0 copies of the same starting unit
  rep(unit, k0)
}
sample_chunk_size <- function(dist, max_k) {
  type <- dist$type

  if (type == "fixed") {
    s <- dist$value

  } else if (type == "poisson") {
    s <- rpois(1, lambda = dist$lambda)

  } else if (type == "normal") {
    s <- round(rnorm(1, mean = dist$mean, sd = dist$sd))

  } else if (type == "geom") {
    s <- rgeom(1, prob = dist$prob)

  } else if (type == "unif") {
    if (is.null(dist$min) || is.null(dist$max)) {
      stop("For type = 'unif', provide dist$min and dist$max")
    }
    s <- round(runif(1, min = dist$min, max = dist$max))

  } else if (type == "gamma") {
    # gamma(shape, scale) -> positive real, then round to integer
    if (is.null(dist$shape) || (is.null(dist$scale) && is.null(dist$rate))) {
      stop("For type = 'gamma', provide dist$shape and dist$scale or dist$rate")
    }
    if (!is.null(dist$scale)) {
      s <- round(rgamma(1, shape = dist$shape, scale = dist$scale))
    } else {
      s <- round(rgamma(1, shape = dist$shape, rate = dist$rate))
    }

  } else {
    stop("Unknown chunk size dist$type: ", type)
  }

  s <- max(1L, s)
  s <- min(s, max_k)
  s
}

local_duplication <- function(units,
                              p_local_dup,       # per-unit probability per generation
                              chunk_size_dist) {

  k <- length(units)
  if (k == 0) return(units)

  # which units trigger a local duplication event?
  triggers <- runif(k) < p_local_dup
  idxs <- which(triggers)
  if (length(idxs) == 0L) return(units)

  # apply events from right to left so that earlier indices aren't shifted
  for (i in rev(idxs)) {
    current_k <- length(units)
    # limit chunk size by how much remains from i to end
    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- sample_chunk_size(chunk_size_dist, max_k = max_chunk)

    start <- i
    end   <- start + chunk_size - 1L
    chunk <- units[start:end]

    # safe tail
    if (end < current_k) {
      tail <- units[(end + 1L):current_k]
    } else {
      tail <- character(0)
    }

    units <- c(units[1:end], chunk, tail)
  }

  units
}


## Distal duplication: p_distal_dup = per-unit per-generation
distal_duplication <- function(units,
                               p_distal_dup,     # per-unit probability per generation
                               chunk_size_dist,
                               p_invert_distal = 0.5) {

  k <- length(units)
  if (k == 0) return(units)

  triggers <- runif(k) < p_distal_dup
  idxs <- which(triggers)
  if (length(idxs) == 0L) return(units)

  # process from right to left
  for (i in rev(idxs)) {
    current_k <- length(units)
    if (i > current_k) next

    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- sample_chunk_size(chunk_size_dist, max_k = max_chunk)

    start <- i
    end   <- start + chunk_size - 1L
    end   <- min(end, current_k)
    chunk <- units[start:end]

    # maybe invert: just reverse order of units
    if (runif(1) < p_invert_distal) {
      chunk <- rev(chunk)
    }

    # choose insert position *to the right* of the chunk: between end..current_k
    # positions are between units, so use slots [end .. current_k]
    #  - if insert_after == end-1: just after chunk (tandem-like)
    #  - if insert_after == current_k: after entire array
    insert_after <- sample(end:current_k, 1L)

    if (insert_after == current_k) {
      units <- c(units, chunk)
    } else {
      units <- c(
        units[1:insert_after],
        chunk,
        units[(insert_after + 1L):current_k]
      )
    }
  }

  units
}


## Deletion: p_del_chunk = per-unit per-generation
delete_chunk <- function(units,
                         p_del_chunk,        # per-unit probability per generation
                         del_chunk_size_dist) {

  k <- length(units)
  if (k == 0) return(units)

  triggers <- runif(k) < p_del_chunk
  idxs <- which(triggers)
  if (length(idxs) == 0L) return(units)

  # apply deletions from right to left
  for (i in rev(idxs)) {
    current_k <- length(units)
    if (i > current_k) next  # if array has shrunk past this index

    max_chunk <- current_k - i + 1L
    if (max_chunk <= 0) next

    chunk_size <- sample_chunk_size(del_chunk_size_dist, max_k = max_chunk)

    start <- i
    end   <- start + chunk_size - 1L
    end   <- min(end, current_k)

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

## Mutation: per-base per-generation (unchanged)
mutate_units <- function(units,
                         mu_total,
                         p_sub    = 1,
                         p_ins    = 0,
                         p_del_bp = 0,
                         alphabet = c("A","C","G","T")) {
  if (mu_total <= 0) return(units)

  probs <- c(p_sub, p_ins, p_del_bp)
  if (any(probs < 0)) stop("p_sub, p_ins, p_del_bp must be >= 0")
  if (sum(probs) == 0) return(units)
  probs <- probs / sum(probs)

  mutate_one_unit <- function(unit) {
    bases <- strsplit(unit, "")[[1]]
    L     <- length(bases)

    mutate_mask <- runif(L) < mu_total
    if (!any(mutate_mask)) return(unit)

    i <- 1
    while (i <= length(bases)) {
      if (!mutate_mask[i]) {
        i <- i + 1
        next
      }

      mut_type <- sample(c("sub","ins","del"), size = 1, prob = probs)

      if (mut_type == "sub") {
        current <- bases[i]
        choices <- setdiff(alphabet, current)
        bases[i] <- sample(choices, 1)
        i <- i + 1

      } else if (mut_type == "ins") {
        new_base <- sample(alphabet, 1)
        bases <- append(bases, new_base, after = i)
        i <- i + 2

      } else if (mut_type == "del") {
        bases <- bases[-i]
        mutate_mask <- mutate_mask[-i]
      }
    }

    paste0(bases, collapse = "")
  }

  vapply(units, mutate_one_unit, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

## One-generation step: **args match your original call**
step_generation <- function(units,
                            p_local_dup,
                            local_chunk_size_dist,
                            p_distal_dup,
                            distal_chunk_size_dist,
                            p_invert_distal = 0.5,
                            p_del_chunk,
                            del_chunk_size_dist,
                            mu_total,
                            p_sub    = 1,
                            p_ins    = 0,
                            p_del_bp = 0,
                            alphabet = c("A","C","G","T")) {

  # 1) multiple local dups
  units <- local_duplication(
    units,
    p_local_dup        = p_local_dup,
    chunk_size_dist    = local_chunk_size_dist
  )

  # 2) multiple distal dups
  units <- distal_duplication(
    units,
    p_distal_dup       = p_distal_dup,
    chunk_size_dist    = distal_chunk_size_dist,
    p_invert_distal    = p_invert_distal
  )

  # 3) multiple deletions
  units <- delete_chunk(
    units,
    p_del_chunk        = p_del_chunk,
    del_chunk_size_dist = del_chunk_size_dist
  )

  # 4) per-base mutations (already per bp per gen)
  units <- mutate_units(
    units,
    mu_total = mu_total,
    p_sub    = p_sub,
    p_ins    = p_ins,
    p_del_bp = p_del_bp,
    alphabet = alphabet
  )

  units
}





make_unit_dt <- function(units, hap=1, chrom=1, sample="test") {
  data.table(
    num    = seq_along(units),
    bponly    =units,
    pos    = 1,
    chrom  = chrom,
    hap    =hap,
    sample = sample,
    dir="-"
  )
}

run_sim_ps <- function(
    n = 1000,
    max_units = 20000,
    max_t=Inf,
    hard_cap  = 40000,
    init_l = 30,
    init_k0 = 10,
    init_sequence_type = "random",
    ancestor_seq=NULL,
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
  # store one ps per replicate (list in case ps is not a scalar)
  ps_list <- vector("list", n)
  monomers_list <- vector("list", n)
  L_vec_list<-vector("list", n)
  H_vec_list<-vector("list", n)
  N_vec_list<-vector("list", n)
  for (i in seq_len(n)) {
    cat("Starting replicate", i, "\n")

    # initial array
    units <- init_array(
      l = init_l,
      k0 = init_k0,
      init_sequence_type = init_sequence_type,
      ancestor_seq=ancestor_seq
    )

    l_vec <- 20L
    t     <- 1L
    H_vec<-c()
    N_vec<-c()
    while (length(units) < max_units & t<max_t) {
      t <- t + 1L

      if (verbose) {
        cat("  Rep", i, "- Generation", t,
            "- number of units:", length(units), "\n")
      }



      if (length(units) > hard_cap) break()

      if (length(units) == 0L) {
        units <- init_array(
          l = init_l,
          k0 = init_k0,
          init_sequence_type = init_sequence_type,
          ancestor_seq=ancestor_seq
        )
      }

      units <- step_generation(
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
      H<-shannon_entropy(units)
      H_vec<-c(H_vec, H)
      l_vec <- c(l_vec, length(units))
      N_vec <- c(N_vec, uniqueN(units))
    }

    L_vec_list[[i]]<-l_vec
    N_vec_list[[i]]<-N_vec
    H_vec_list[[i]]<-H_vec
    dt_units <- make_unit_dt(units)
    monomers_list[[i]]<-dt_units
    ps_list[[i]] <- pairs_identical(dt_units)
  }

  list(monomers_list=monomers_list, ps_list=ps_list, L_vec_list=L_vec_list, H_vec_list=H_vec_list, N_vec_list=N_vec_list)
}

plot.repeat.fingerprint <- function(ps, kb=NULL, rotate=T, zoom=NULL, ptsize=0.001, gridalpha=1,gbreaks=20, ptcol="black") {


  # assume you already have ps, zoom, etc.

  ps[, `:=`(num1 = as.numeric(num1), num2 = as.numeric(num2))]

  if (!is.null(zoom)) {
    ps <- ps[num1 > zoom[1] & num1 < zoom[2] &
               num2 > zoom[1] & num2 < zoom[2]]
  }

  break_fun <- int_breaks_pretty(gbreaks)
  if (!is.null(zoom)) {
    lims <- zoom
  } else {
    lims <- range(c(ps$num1, ps$num2), na.rm = TRUE)
  }
  grid_breaks <- break_fun(lims)

  min_lim <- min(lims)
  max_lim <- max(lims)

  # vertical lines: top edge -> diagonal y = x
  grid_vert <- data.frame(
    x    = grid_breaks,
    xend = grid_breaks,
    y    = max_lim,
    yend = grid_breaks            # intersection with diagonal at y = x = break
  )

  # horizontal lines: left edge -> diagonal y = x
  grid_horiz <- data.frame(
    x    = min_lim,
    xend = grid_breaks,
    y    = grid_breaks,
    yend = grid_breaks
  )

  # choose integer bp positions between limits
  grid_bp <- seq(ceiling(min_lim), floor(max_lim), by = 1)

  # vertical lines: top edge -> diagonal y = x
  grid_vertbp <- data.frame(
    x    = grid_bp,
    xend = grid_bp,
    y    = max_lim,
    yend = grid_bp          # intersection with diagonal at y = x = break
  )

  # horizontal lines: left edge -> diagonal y = x
  grid_horizbp <- data.frame(
    x    = min_lim,
    xend = grid_bp,
    y    = grid_bp,
    yend = grid_bp
  )

  p1_base <- ggplot(ps[num1 != num2], aes(x = num1, y = num2)) +
    geom_segment(
      data = grid_vertbp,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey90", linewidth = 0.2, alpha = gridalpha
    ) +
    geom_segment(
      data = grid_horizbp,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey90", linewidth = 0.2, alpha = gridalpha
    ) +
    geom_point(size = ptsize, shape = 15, alpha = 1, col=ptcol) +
    coord_fixed() +
    theme_classic(base_size = 6) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(
      expand   = c(0, 0),
      position = "top",
      limits   = lims,
      breaks   = grid_breaks
    ) +
    scale_y_continuous(
      expand   = c(0, 0),
      position = "left",
      limits   = lims,
      breaks   = grid_breaks
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.text.y = element_text(angle = 45, hjust = 1),
      axis.line   = element_line(),
      legend.position = "none",
      plot.background = element_blank(),
      panel.background = element_blank()
    )

  if(is.null(zoom)){

    p1_base<-p1_base+geom_point(data = kb, aes(x = num, y = max(kb$num, na.rm = TRUE) + 100, color = dir),
                                size = 0.1, inherit.aes = FALSE, shape=15) +
      geom_point(data = kb, aes(x = -100, y = num, color = dir),
                 size = 0.1, inherit.aes = FALSE, shape=15)
  }


  if(rotate){
    g <- ggplotGrob(p1_base)

    grid.newpage()
    pushViewport(viewport(angle = -45, width = 0.75, height = 0.75, clip = "off"))
    grid.draw(g)
    popViewport()
  } else {
    return(p1_base)
  }
}


pairs_identical <- function(repeats) {
  x <- as.data.table(repeats)
  x[, id := .I]

  # for each bponly, output all id pairs; emit NULL for singleton groups
  pair_idx <- x[, {
    if (.N < 2L) NULL else {
      cmb <- utils::combn(id, 2L)
      data.table(id1 = cmb[1,], id2 = cmb[2,])
    }
  }, by = bponly]

  if (nrow(pair_idx) == 0L) return(data.table())

  # join metadata for id1 and id2 with base-merge to avoid i.* NSE
  m1 <- merge(pair_idx,
              x[, .(id1 = id, bponly, sample1 = sample, hap1 = hap, chrom1 = chrom, num1 = (num), pos1 = (pos), dir1=dir)],
              by = c("bponly", "id1"), all.x = TRUE, sort = FALSE)
  out <- merge(m1,
               x[, .(id2 = id, bponly, sample2 = sample, hap2 = hap, chrom2 = chrom, num2 = (num), pos2 = (pos),  dir2=dir)],
               by = c("bponly", "id2"), all.x = TRUE, sort = FALSE)

  # order nicely and drop helper ids if you want
  out[order(chrom1, chrom2, num1=(num1), num2=(num2), pos1=(pos1), pos2=(pos2))]
}


shannon_entropy <- function(x) {
  p <- table(x)
  p <- p[p > 0] / sum(p)   # convert to probabilities, drop zeros if any
  -sum(p * log2(p))
}


shannon_entropyplus <- function(x) {
  p <- table(x)
  p <- p[p > 0] / sum(p)   # convert to probabilities, drop zeros if any
  H=-sum(p * log2(p))
  uniqueN=uniqueN(x)
  L=length(x)
  H_logL=H/log(L)
  H_logUniqueN=H/log(uniqueN)
  data.table(H, uniqueN, L, H_logL, H_logUniqueN)
}


kb_row_window_entropy <- function(kb, N) {
  kb <- as.data.table(kb)

  out_list <- pblapply(seq_len(nrow(kb)-N), function(i) {
    # subset: row i plus next N rows (clipped at end)
    sub_idx <- i:min(nrow(kb), i + N)
    sub     <- kb[sub_idx]

    # your existing function
    res <- shannon_entropyplus(sub$bponly)

    # position = focal row index
    res$pos <- i

    res
  })

  rbindlist(out_list, fill = TRUE)
}

plot_knob_barcode <- function(ps,
                              zoom,
                              palette    = "Dark2",
                              line_width = 2,
                              base_size  = 6,
                              label_size = 2,
                              show_labels = TRUE,
                              gapsize=0.9, vjusttext=-0.1) {
  stopifnot(length(zoom) == 2)

  # subset ps in the zoom window and coerce to data.frame
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

  # integer range across zoom
  range_x <- seq(ceiling(zoom[1]), floor(zoom[2]))
  if (length(range_x) == 0) {
    range_x <- zoom[1]:zoom[2]
  }

  # segment length scaled to window width
  gap <- gapsize * (zoom[2] - zoom[1]) / length(range_x)

  # background gray segments at each bp
  bg <- data.frame(
    x    = range_x,
    xend = range_x + gap,
    y    = 1,
    yend = 1
  )

  # label positions only where we actually have segments in ps
  lab_nums <- sort(unique(c(ps_zoom$num1, ps_zoom$num2)))
  labs <- data.frame(
    x   = lab_nums + gap / 2,  # center the label over the segment
    y   = 1,
    lab = lab_nums
  )

  p <- ggplot() +
    # gray background "grid" segments
    geom_segment(
      data = bg,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = line_width,
      lineend   = "round",
      colour    = "gray85"
    ) +
    # coloured segments for num1 / num2 from ps_zoom
    geom_segment(
      data = ps_zoom,
      aes(x = num1, xend = num1 + gap, y = 1, yend = 1, colour = factorlevel),
      linewidth = line_width,
      lineend   = "round"
    ) +
    geom_segment(
      data = ps_zoom,
      aes(x = num2, xend = num2 + gap, y = 1, yend = 1, colour = factorlevel),
      linewidth = line_width,
      lineend   = "round"
    ) +
    scale_color_brewer(palette = palette) +
    coord_cartesian(xlim = zoom, ylim = c(0.95, 1.25), expand = FALSE) +
    theme_void(base_size = base_size) +
    theme(
      legend.position = "none"
    )

  if (show_labels) {
    p <- p +
      geom_text(
        data = labs,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size   = label_size,
        angle  = 90,
        vjust  = vjusttext,
        hjust  = 1.2
      )
  }

  p
}



counts_long_nogap <- function(rawseqs,
                              gap_char = "-",
                              tie_break_levels = c("A","C","G","T","N","R","Y","S","W","K","M","B","D","H","V","-")) {
  stopifnot(length(rawseqs) >= 2L)
  x <- as.character(rawseqs)
  Ls <- nchar(x)
  if (length(unique(Ls)) != 1L) stop("All sequences must have the same aligned length")
  L <- Ls[1L]
  n <- length(x)

  # build matrix (n x L)
  chars <- strsplit(x, "", fixed = TRUE)
  mat <- do.call(rbind, chars)

  # count every symbol present (including gaps/ambiguities)
  syms <- sort(unique(c(mat)))
  counts_all <- data.table::as.data.table(
    sapply(syms, function(s) colSums(mat == s, na.rm = TRUE))
  )
  counts_all[, pos := seq_len(L)]
  data.table::setcolorder(counts_all, c("pos", setdiff(names(counts_all), "pos")))

  # long format with per-position totals & proportions
  long_counts <- data.table::melt(
    counts_all,
    id.vars = "pos",
    variable.name = "symbol",
    value.name = "count"
  )
  long_counts[, total := sum(count), by = pos]
  long_counts[, prop  := ifelse(total > 0, count / total, NA_real_)]
  long_counts[, symbol := as.character(symbol)]

  # consensus with tie-break (prefer bases over gap using provided priority)
  pri_map <- seq_along(tie_break_levels)
  names(pri_map) <- tie_break_levels
  # ensure gap has the lowest priority
  if (!gap_char %in% names(pri_map)) {
    tie_break_levels <- c(tie_break_levels, gap_char)
    pri_map <- seq_along(tie_break_levels); names(pri_map) <- tie_break_levels
  }
  pri_map[gap_char] <- -1L

  long_counts[, pri := pri_map[symbol]]

  consensus_dt <- long_counts[
    , {
      m <- max(count)
      cand <- .SD[count == m]
      cand[which.max(pri)][, .(consensus = symbol)]
    },
    by = pos
  ]

  # keep positions where consensus != gap, and renumber
  keep_map <- consensus_dt[consensus != gap_char, .(pos, consensus)]
  keep_map[, pos_new := seq_len(.N)]

  # filter long table and attach new numbering
  long_counts_nogap <- keep_map[long_counts, on = "pos", nomatch = 0L]
  data.table::setcolorder(long_counts_nogap, c("pos_new", "pos", "symbol", "count", "total", "prop", "consensus", "pri"))

  # drop helper
  long_counts_nogap[, pri := NULL][]  # return data.table
}


get_sim_entropies<-function(ps_results){
  entropies<-rbindlist(lapply(1:length(ps_results$monomers_list), function(i){
    H<-shannon_entropy(ps_results$monomers_list[[i]]$bponly)
    return(data.table(i, H, total=nrow(ps_results$monomers_list[[i]]), uniqueN=uniqueN(ps_results$monomers_list[[i]]$bponly)))
  }))
  entropies$meanN<-entropies$total/entropies$uniqueN
  return(entropies)
}

int_breaks <- function(step = 1) {
  function(x) seq(floor(min(x, na.rm=TRUE)),
                  ceiling(max(x, na.rm=TRUE)), by = step)
}
int_breaks_pretty <- function(n = 5) {
  function(x) {
    b <- pretty(range(x, na.rm = TRUE), n)
    b[b == round(b)]
  }
}


plot_circular_seqcounts <- function(seqcounts,
                                    allowed = c("-", "A","T","G","C"),
                                    title   = "",
                                    palette = "Dark2",
                                    show_consensus_labels = FALSE) {
  require(data.table)
  require(ggplot2)

  dt <- as.data.table(seqcounts)

  # checks
  needed <- c("pos_new","symbol","prop","consensus")
  stopifnot(all(needed %in% names(dt)))

  # symbol factor (inner rings)
  dt[, symbol := factor(as.character(symbol), levels = allowed)]

  # consensus by position (outer ring)
  consensus <- dt[, .(consensus = unique(consensus)), by = pos_new][order(pos_new)]
  max_x <- max(consensus$pos_new)

  # y-levels include an extra "CONS" outer band
  y_levels <- c(levels(dt$symbol), "CONS")

  # optional consensus labels arranged around the ring
  cons_lab <- consensus[, {
    ang  <- (pos_new - 0.5) / max_x * 360
    ang2 <- ifelse(ang > 90 & ang < 270, ang + 180, ang)  # keep upright
    hj   <- ifelse(ang > 90 & ang < 270, 1, 0)
    .(pos_new = pos_new,
      symbol  = factor("CONS", levels = y_levels),
      label   = consensus,
      angle   = ang2,
      hjust   = hj)
  }]

  # base plot
  p <- ggplot(dt, aes(x = pos_new, y = symbol)) +
    # subtle background tiles for spacing/stroke
    geom_tile(fill = "white", linewidth = 0.1, width = 1, height = 0.75, col = "gray90") +
    geom_tile(fill = "white", linewidth = 0.1, width = 1.1, height = 0.6, col = "white") +

    # data tiles with alpha = proportion; stroke follows fill alpha
    geom_tile(
      aes(
        fill   = symbol,
        alpha  = prop,
        colour = after_scale(scales::alpha(fill, alpha))
      ),
      linewidth = 0.1, width = 1, height = 0.75
    ) +

    coord_polar() +
    scale_alpha(range = c(0, 0.8), guide = "none") +
    scale_x_continuous(expand = c(0.002, 0.002)) +
    scale_y_discrete(drop = FALSE, limits = y_levels, expand = c(3, 2)) +
    scale_fill_brewer(palette = palette, name = "") +
    scale_color_brewer(palette = palette, name = "") +
    theme_void(base_size = 6) +
    theme(
      legend.position = c(0.5, 0.5),
      legend.key.size = grid::unit(0.1, "cm"),
      panel.grid      = element_blank(),
      plot.margin     = grid::unit(c(0, 0, 0, 0), "pt")
    ) +
    ggtitle(title)

  # optional outer consensus text labels
  if (isTRUE(show_consensus_labels)) {
    p <- p +
      geom_text(
        data = cons_lab,
        aes(x = pos_new, y = symbol, label = label, angle = angle, hjust = hjust),
        inherit.aes = FALSE,
        size = 1.5
      )
  }

  # a faint radial guide line (center)
  p <- p + geom_segment(x = 0, xend = 0, y = 0, yend = length(y_levels), linewidth = 0.1, col = "gray")

  p
}

mu_load_from_consensus <- function(rawseqs, consensus) {
  x  <- as.character(rawseqs)
  Ls <- nchar(x)

  if (length(unique(Ls)) != 1L) {
    stop("All sequences must have the same aligned length")
  }

  L <- Ls[1L]
  consensus<-consensus[order(pos)]
  if (length(consensus$consensus) != L) {
    stop("Consensus length must match sequence length")
  }

  # build matrix of characters and count mismatches to consensus
  chars   <- strsplit(x, "", fixed = TRUE)
  mu_load <- sapply(chars, function(i) sum(i != consensus$consensus))

  return(mu_load)
}

plot_ps_summary <- function(ps_results, i, title = NULL,
                            rel_widths_bottom = c(1, 1, 1, 1),
                            rel_heights = c(2, 1)) {
  library(data.table)
  library(ggplot2)
  library(cowplot)

  ## pull out data
  ps    <- ps_results$ps_list[[i]]
  monos <- as.data.table(ps_results$monomers_list[[i]])

  ## counts + consensus
  counts <- as.data.table(counts_long_nogap(monos$bponly))
  consensus <- counts[symbol == consensus][order(pos)]

  ## add mutation load to monomers (in place)
  monos[, load := mu_load_from_consensus(bponly, consensus)]
  ps_results$monomers_list[[i]] <- monos

  base_theme <- theme_classic(base_size = 6)

  ## p1: distance between pairs
  p1 <- ggplot(ps, aes(x = abs(num1 - num2))) +
    geom_histogram(bins=50) +
    scale_x_log10() +
    base_theme +
    xlab("|distance|") +
    ylab("Count")

  ## p2: fingerprint (big, on its own row)
  p2 <- plot.repeat.fingerprint(ps, kb = monos,
                                rotate = FALSE,
                                ptsize = 0.01,
                                gridalpha = 0) +
    theme(panel.grid       = element_blank(),
          panel.background = element_blank())

  if (!is.null(title)) {
    p2 <- p2 +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  ## p3: copy number of each unique monomer sequence
  p3 <- ggplot(data.table(table(monos$bponly)), aes(x = N)) +
    geom_histogram(bins=50) +
    scale_x_log10() +
    base_theme +
    xlab("Monomer copy number") +
    ylab("Count")

  ## p4: mutation load histogram
  p4 <- ggplot(monos, aes(x = load)) +
    geom_histogram(bins=50) +
    base_theme +
    xlab("Mismatch load to consensus") +
    ylab("Count")

  ## p5: consensus prop histogram
  p5 <- ggplot(consensus, aes(x = prop)) +
    geom_histogram(bins=50) +
    base_theme +
    xlab("Consensus prop") +
    ylab("Count")

  ## bottom row: all small panels in one row, adjustable widths
  bottom_row <- plot_grid(p1, p3, p4, p5,
                          nrow = 1,
                          rel_widths = rel_widths_bottom)

  ## top: p2 big; bottom: others
  combined <- plot_grid(p2, bottom_row,
                        ncol = 1,
                        rel_heights = rel_heights)

  combined
}

