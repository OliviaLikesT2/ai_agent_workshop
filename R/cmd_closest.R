# mytools closest -a <file|-> -b <file> [-d]
# For each -a feature, the nearest -b feature(s) on the same chromosome.
# Both inputs must be sorted. Ties are all reported (bedtools -t all).

closest_usage <- "usage: mytools closest -a <file|-> -b <file> [-d]"

closest_parse_args <- function(args) {
  a <- NULL; b <- NULL; d <- FALSE; i <- 1L
  while (i <= length(args)) {
    if (args[i] == "-a" && i < length(args)) { a <- args[i + 1L]; i <- i + 2L }
    else if (args[i] == "-b" && i < length(args)) { b <- args[i + 1L]; i <- i + 2L }
    else if (args[i] == "-d") { d <- TRUE; i <- i + 1L }
    else bed_fail(paste0("mytools closest: unexpected argument '", args[i], "'\n", closest_usage), 2L)
  }
  if (is.null(a) || is.null(b)) bed_fail(closest_usage, 2L)
  if (identical(b, "-")) bed_fail("mytools closest: -b must be a file, not stdin", 2L)
  list(a = a, b = b, d = d)
}

bed_lines <- function(df) {
  out <- paste(df$chrom, df$start, df$end, sep = "\t")
  has_extra <- nzchar(df$extra)
  out[has_extra] <- paste(out[has_extra], df$extra[has_extra], sep = "\t")
  out
}

bed_ncol <- function(df) {
  if (nrow(df) == 0 || !nzchar(df$extra[1])) return(3L)
  3L + nchar(df$extra[1]) - nchar(gsub("\t", "", df$extra[1], fixed = TRUE)) + 1L
}

# What bedtools prints for -b when nothing is on the chromosome. The shape
# follows -b's column count: BED4-6 get "." / -1 / "." for name/score/strand,
# anything wider gets "." throughout. (bedtools 2.31.1.)
closest_null_row <- function(ncol) {
  cols <- if (ncol <= 6L) c(".", "-1", "-1", ".", "-1", ".")[seq_len(ncol)]
          else c(".", "-1", "-1", rep(".", ncol - 3L))
  paste(cols, collapse = "\t")
}

closest_index <- function(b, source = "-b") {
  check_sorted(b, source)
  span <- bedtools_span(b$start, b$end)
  lines <- bed_lines(b)
  by_chrom <- lapply(split(seq_len(nrow(b)), b$chrom), function(i) {
    e <- span$end[i]
    eo <- order(e)
    list(idx = i, start = span$start[i], end = e, cummax_end = cummax(e),
         end_sorted = e[eo], end_order = eo)
  })
  list(chroms = by_chrom, lines = lines, null_row = closest_null_row(bed_ncol(b)))
}

closest_format <- function(a_lines, b_lines, dist, report_dist) {
  if (report_dist) paste(a_lines, b_lines, dist, sep = "\t") else paste(a_lines, b_lines, sep = "\t")
}

# Pure: output lines for the rows of `a`, in `a` order.
# Vectorised for the common shapes (one overlapping -b row, or a single nearest
# neighbour); only rows with ties fall through to the per-row loop.
closest_rows <- function(a, index, report_dist = FALSE) {
  n <- nrow(a)
  if (n == 0) return(character())
  span <- bedtools_span(a$start, a$end)
  a_lines <- bed_lines(a)
  out <- vector("list", n)

  for (chrom in unique(a$chrom)) {
    rows <- which(a$chrom == chrom)
    ix <- index$chroms[[chrom]]
    if (is.null(ix)) {
      out[rows] <- as.list(closest_format(a_lines[rows], index$null_row, "-1", report_dist))
      next
    }
    as <- span$start[rows]; ae <- span$end[rows]
    nb <- length(ix$start)
    at <- function(v, i, ok) ifelse(ok, v[pmax(pmin(i, length(v)), 1L)], NA_integer_)

    k <- findInterval(ae - 1L, ix$start)                    # -b rows with start < a.end
    has_hit <- k >= 1L & at(ix$cummax_end, k, k >= 1L) > as
    has_hit[is.na(has_hit)] <- FALSE
    prev_cm <- at(ix$cummax_end, k - 1L, k >= 2L)
    single_hit <- has_hit & at(ix$end, k, k >= 1L) > as & (is.na(prev_cm) | prev_cm <= as)

    nl <- findInterval(as, ix$end_sorted)                   # -b rows with end <= a.start
    maxend <- at(ix$end_sorted, nl, nl >= 1L)
    lo <- findInterval(maxend - 1L, ix$end_sorted) + 1L
    dl <- as - maxend + 1L
    nr <- k + 1L                                            # first -b row with start >= a.end
    minstart <- at(ix$start, nr, nr <= nb)
    hi <- findInterval(minstart, ix$start)
    dr <- minstart - ae + 1L

    use_left  <- !is.na(dl) & (is.na(dr) | dl < dr)
    use_right <- !is.na(dr) & (is.na(dl) | dr < dl)
    single_left  <- !has_hit & use_left  & (nl - lo) == 0L
    single_right <- !has_hit & use_right & (hi - nr) == 0L

    sel <- rep(NA_integer_, length(rows)); dist <- rep(NA_integer_, length(rows))
    sel[single_hit]   <- k[single_hit];                  dist[single_hit]   <- 0L
    sel[single_left]  <- ix$end_order[nl[single_left]];  dist[single_left]  <- dl[single_left]
    sel[single_right] <- nr[single_right];               dist[single_right] <- dr[single_right]
    simple <- !is.na(sel)
    if (any(simple)) {
      out[rows[simple]] <- as.list(closest_format(
        a_lines[rows[simple]], index$lines[ix$idx[sel[simple]]], dist[simple], report_dist))
    }

    for (r in which(!simple)) {
      i <- rows[r]
      if (has_hit[r]) {
        hits <- integer(); j <- k[r]
        while (j >= 1L && ix$cummax_end[j] > as[r]) {
          if (ix$end[j] > as[r]) hits <- c(hits, j)
          j <- j - 1L
        }
        chosen <- sort(hits); d <- 0L
      } else {
        left  <- if (!is.na(dl[r])) ix$end_order[lo[r]:nl[r]] else integer()
        right <- if (!is.na(dr[r])) nr[r]:hi[r] else integer()
        if (use_left[r])       { chosen <- left;  d <- dl[r] }
        else if (use_right[r]) { chosen <- right; d <- dr[r] }
        else                   { chosen <- c(left, right); d <- dl[r] }
        chosen <- sort(chosen)
      }
      out[[i]] <- closest_format(a_lines[i], index$lines[ix$idx[chosen]], d, report_dist)
    }
  }
  unlist(out, use.names = FALSE)
}

cmd_closest <- function(args) {
  opt <- closest_parse_args(args)
  index <- closest_index(read_bed(opt$b), opt$b)
  r <- bed_reader(opt$a)
  on.exit(r$close())
  state <- NULL
  source <- if (identical(opt$a, "-")) "<stdin>" else opt$a
  while (!is.null(chunk <- r$read(50000L))) {
    state <- check_sorted(chunk, source, state)
    lines <- closest_rows(chunk, index, opt$d)
    if (length(lines)) writeLines(lines)
  }
  invisible(NULL)
}
