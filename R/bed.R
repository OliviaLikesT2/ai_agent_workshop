# Shared BED I/O and interval helpers. Functions only; sourcing this has no side effects.
# BED is 0-based, half-open: chr1 100 200 covers bases 100..199.

bed_fail <- function(msg, status) {
  cat(msg, "\n", file = stderr(), sep = "")
  quit(save = "no", status = status)
}

bed_open <- function(path) {
  if (identical(path, "-")) {
    con <- file("stdin")
  } else {
    if (!file.exists(path)) bed_fail(sprintf("%s: no such file", path), 2L)
    con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  }
  open(con, "r")
  con
}

bed_parse_error <- function(source, line, msg) {
  stop(structure(
    class = c("bed_parse_error", "error", "condition"),
    list(message = sprintf("%s:%d: %s", source, line, msg), call = NULL)
  ))
}

# Parse raw lines into a data.frame(chrom, start, end, extra).
# `extra` is everything after column 3, verbatim, or "" if there is nothing.
# Skips blank, "#", "track" and "browser" lines. Errors name source:line.
parse_bed <- function(lines, source = "<input>", first_line = 1L) {
  lineno <- seq_along(lines) + first_line - 1L
  keep <- !(grepl("^\\s*$", lines) |
            startsWith(lines, "#") |
            startsWith(lines, "track") |
            startsWith(lines, "browser"))
  lines <- lines[keep]
  lineno <- lineno[keep]

  empty <- data.frame(chrom = character(), start = integer(), end = integer(),
                      extra = character(), stringsAsFactors = FALSE)
  if (length(lines) == 0) return(empty)

  ntab <- nchar(lines) - nchar(gsub("\t", "", lines, fixed = TRUE))
  bad <- which(ntab < 2)
  if (length(bad)) {
    bed_parse_error(source, lineno[bad[1]],
      sprintf("expected at least 3 tab-separated fields, got %d", ntab[bad[1]] + 1L))
  }

  chrom   <- sub("\t.*$", "", lines)
  rest    <- sub("^[^\t]*\t", "", lines)
  start_s <- sub("\t.*$", "", rest)
  rest    <- sub("^[^\t]*\t", "", rest)
  end_s   <- sub("\t.*$", "", rest)
  extra   <- ifelse(ntab >= 3, sub("^[^\t]*\t", "", rest), "")

  is_int <- function(s) grepl("^[0-9]+$", s)
  bad <- which(!is_int(start_s))
  if (length(bad)) {
    bed_parse_error(source, lineno[bad[1]],
      sprintf("start is not a non-negative integer ('%s')", start_s[bad[1]]))
  }
  bad <- which(!is_int(end_s))
  if (length(bad)) {
    bed_parse_error(source, lineno[bad[1]],
      sprintf("end is not a non-negative integer ('%s')", end_s[bad[1]]))
  }

  start <- as.integer(start_s)
  end   <- as.integer(end_s)
  bad <- which(start > end)
  if (length(bad)) {
    bed_parse_error(source, lineno[bad[1]],
      sprintf("start > end (%d > %d)", start[bad[1]], end[bad[1]]))
  }

  data.frame(chrom = chrom, start = start, end = end, extra = extra,
             stringsAsFactors = FALSE)
}

# Chunked reader for subcommands that stream. Usage:
#   r <- bed_reader(path); while (!is.null(df <- r$read(n))) { ... }; r$close()
bed_reader <- function(path) {
  con <- bed_open(path)
  source <- if (identical(path, "-")) "<stdin>" else path
  next_line <- 1L
  list(
    read = function(n = 10000L) {
      lines <- readLines(con, n = n, warn = FALSE)
      if (length(lines) == 0) return(NULL)
      df <- tryCatch(parse_bed(lines, source, next_line),
                     bed_parse_error = function(e) bed_fail(conditionMessage(e), 1L))
      next_line <<- next_line + length(lines)
      df
    },
    close = function() close(con)
  )
}

read_bed <- function(path) {
  r <- bed_reader(path)
  on.exit(r$close())
  chunks <- list()
  while (!is.null(df <- r$read(100000L))) chunks[[length(chunks) + 1L]] <- df
  if (length(chunks) == 0) return(parse_bed(character()))
  do.call(rbind, chunks)
}

# Two half-open intervals overlap iff a.start < b.end AND b.start < a.end. Strict.
# Bookended intervals (a.end == b.start) do not overlap. Vectorised.
overlaps <- function(a_start, a_end, b_start, b_end) {
  a_start < b_end & b_start < a_end
}

# bedtools treats a zero-length feature at p as spanning [p-1, p+1) for overlap
# and distance purposes, while printing it unchanged. Use these coordinates for
# any comparison; print the originals. (Verified against bedtools 2.31.1.)
bedtools_span <- function(start, end) {
  z <- start == end
  list(start = ifelse(z, start - 1L, start), end = ifelse(z, end + 1L, end))
}

# Error unless rows are sorted the way bedtools requires for merge/closest:
# no chromosome reappears after another, and starts never decrease within one.
# `state` carries chrom/start/seen across chunks; returns the updated state.
check_sorted <- function(df, source, state = NULL) {
  if (is.null(state)) state <- list(chrom = NA_character_, start = -1L, seen = character())
  n <- nrow(df)
  if (n == 0) return(state)
  chrom <- df$chrom; start <- df$start
  prev_chrom <- c(state$chrom, chrom[-n])
  prev_start <- c(state$start, start[-n])
  same <- !is.na(prev_chrom) & chrom == prev_chrom
  bad <- which(same & start < prev_start)
  if (length(bad)) {
    bed_fail(sprintf("%s: unsorted input at %s:%d (previous start %d)",
                     source, chrom[bad[1]], start[bad[1]], prev_start[bad[1]]), 1L)
  }
  runs <- rle(chrom)$values
  if (identical(runs[1], state$chrom)) runs <- runs[-1]
  dup <- runs[duplicated(runs) | runs %in% state$seen]
  if (length(dup)) {
    bed_fail(sprintf("%s: unsorted input, %s reappears after other chromosomes", source, dup[1]), 1L)
  }
  list(chrom = chrom[n], start = start[n], seen = c(state$seen, runs))
}

write_bed <- function(df, con = stdout()) {
  if (nrow(df) == 0) return(invisible(NULL))
  out <- paste(df$chrom, as.integer(df$start), as.integer(df$end), sep = "\t")
  has_extra <- nzchar(df$extra)
  out[has_extra] <- paste(out[has_extra], df$extra[has_extra], sep = "\t")
  writeLines(out, con)
  invisible(NULL)
}
