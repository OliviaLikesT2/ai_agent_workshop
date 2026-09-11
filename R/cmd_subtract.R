# mytools subtract -a <file|-> -b <file>
# Removes from each -a feature every region covered by any -b feature.

subtract_usage <- "usage: mytools subtract -a <file|-> -b <file>"

subtract_parse_args <- function(args) {
  a <- NULL; b <- NULL; i <- 1L
  while (i <= length(args)) {
    if (args[i] == "-a" && i < length(args)) { a <- args[i + 1L]; i <- i + 2L }
    else if (args[i] == "-b" && i < length(args)) { b <- args[i + 1L]; i <- i + 2L }
    else bed_fail(paste0("mytools subtract: unexpected argument '", args[i], "'\n", subtract_usage), 2L)
  }
  if (is.null(a) || is.null(b)) bed_fail(subtract_usage, 2L)
  if (identical(b, "-")) bed_fail("mytools subtract: -b must be a file, not stdin", 2L)
  list(a = a, b = b)
}

# Merge sorted-by-start intervals into disjoint, ascending, non-touching runs.
union_intervals <- function(start, end) {
  n <- length(start)
  if (n == 0) return(list(start = integer(), end = integer()))
  o <- order(start, end); start <- start[o]; end <- end[o]
  cm <- cummax(end)
  first <- c(TRUE, start[-1] > cm[-n])
  last  <- c(first[-1], TRUE)
  list(start = start[first], end = cm[last])
}

# Per-chromosome union of -b using bedtools' zero-length expansion.
subtract_index <- function(b) {
  span <- bedtools_span(b$start, b$end)
  idx <- split(seq_len(nrow(b)), b$chrom)
  lapply(idx, function(i) union_intervals(span$start[i], span$end[i]))
}

# Pure: returns the rows of `a` with `index` regions removed, in `a` order.
subtract_rows <- function(a, index) {
  n <- nrow(a)
  if (n == 0) return(a)
  span <- bedtools_span(a$start, a$end)
  zero <- a$start == a$end

  keep_idx <- integer(); keep_start <- integer(); keep_end <- integer()
  out_idx <- list(); out_start <- list(); out_end <- list()

  for (chrom in unique(a$chrom)) {
    rows <- which(a$chrom == chrom)
    u <- index[[chrom]]
    if (is.null(u) || length(u$start) == 0) {
      keep_idx <- c(keep_idx, rows); keep_start <- c(keep_start, a$start[rows]); keep_end <- c(keep_end, a$end[rows])
      next
    }
    first <- findInterval(span$start[rows], u$end) + 1L        # first union with end > a.start
    last  <- findInterval(span$end[rows] - 1L, u$start)        # last union with start < a.end
    untouched <- first > last
    r <- rows[untouched]
    keep_idx <- c(keep_idx, r); keep_start <- c(keep_start, a$start[r]); keep_end <- c(keep_end, a$end[r])

    for (k in which(!untouched)) {
      i <- rows[k]
      cur <- span$start[i]; a_end <- span$end[i]
      ps <- integer(); pe <- integer()
      for (j in first[k]:last[k]) {
        if (u$start[j] > cur) { ps <- c(ps, cur); pe <- c(pe, u$start[j]) }
        cur <- max(cur, u$end[j])
      }
      if (cur < a_end) { ps <- c(ps, cur); pe <- c(pe, a_end) }
      if (zero[i]) {
        # bedtools prints a zero-length feature unchanged unless fully covered
        if (length(ps)) { ps <- a$start[i]; pe <- a$end[i] }
      }
      if (length(ps)) {
        out_idx[[length(out_idx) + 1L]] <- rep.int(i, length(ps))
        out_start[[length(out_start) + 1L]] <- ps
        out_end[[length(out_end) + 1L]] <- pe
      }
    }
  }

  idx   <- c(keep_idx, unlist(out_idx))
  start <- c(keep_start, unlist(out_start))
  end   <- c(keep_end, unlist(out_end))
  o <- order(idx, start)
  res <- a[idx[o], , drop = FALSE]
  res$start <- as.integer(start[o]); res$end <- as.integer(end[o])
  rownames(res) <- NULL
  res
}

cmd_subtract <- function(args) {
  opt <- subtract_parse_args(args)
  index <- subtract_index(read_bed(opt$b))
  r <- bed_reader(opt$a)
  on.exit(r$close())
  while (!is.null(chunk <- r$read(50000L))) write_bed(subtract_rows(chunk, index))
  invisible(NULL)
}
