#!/usr/bin/env Rscript
# Unit tests for R/bed.R. Runs without bedtools: Rscript tests/test_bed.R
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "harness.R"))

# --- overlaps(): strict on both sides ---------------------------------------

check("overlapping pair overlaps",          overlaps(100L, 200L, 150L, 250L))
check("nested interval overlaps",           overlaps(100L, 400L, 200L, 300L))
check("identical intervals overlap",        overlaps(100L, 200L, 100L, 200L))
check("bookended pair does NOT overlap",   !overlaps(100L, 200L, 200L, 300L))
check("bookended, other order",            !overlaps(200L, 300L, 100L, 200L))
check("disjoint pair does not overlap",    !overlaps(100L, 200L, 300L, 400L))
check("interval at position 0 overlaps",    overlaps(0L, 100L, 50L, 150L))
check("one-base overlap counts",            overlaps(100L, 200L, 199L, 300L))
check("overlaps is vectorised",
      identical(overlaps(c(100L, 100L), c(200L, 200L), c(150L, 200L), c(250L, 300L)),
                c(TRUE, FALSE)))

# --- parse_bed(): shape and skipping ----------------------------------------

df <- parse_bed(c("chr1\t100\t200\ta\t0\t+", "chr1\t300\t400"))
check("parse returns chrom/start/end/extra",
      identical(names(df), c("chrom", "start", "end", "extra")))
check("start and end are integer", is.integer(df$start) && is.integer(df$end))
check("extra keeps trailing columns verbatim", df$extra[1] == "a\t0\t+")
check("extra is empty string when only BED3", df$extra[2] == "")
check("columns past 6 survive",
      parse_bed("chr1\t1\t2\tn\t0\t+\tseven\teight")$extra == "n\t0\t+\tseven\teight")

df <- parse_bed(c("# comment", "track name=x", "browser position chr1", "",
                  "   ", "chr1\t0\t100"))
check("comment/track/browser/blank lines skipped", nrow(df) == 1 && df$start == 0L)
check("empty input gives zero rows", nrow(parse_bed(character())) == 0)
check("all-comment input gives zero rows", nrow(parse_bed("# only")) == 0)

# --- edge cases the fixtures are built around --------------------------------

df <- parse_bed("chr1\t500\t500")
check("zero-length interval is legal", nrow(df) == 1 && df$start == df$end)
df <- parse_bed("chr2\t0\t0")
check("zero-length interval at position 0 is legal", nrow(df) == 1 && df$end == 0L)
df <- parse_bed("chr1\t0\t100")
check("interval at position 0 parses", df$start == 0L)

# --- parse_bed(): rejections, with file:line in the message -----------------

e <- check_error("start > end rejected", parse_bed("chr1\t500\t400", "x.bed"), "bed_parse_error")
check("start > end message names file:line", grepl("^x.bed:1: start > end \\(500 > 400\\)", conditionMessage(e)))

e <- check_error("non-integer start rejected", parse_bed("chr1\tabc\t400", "x.bed"), "bed_parse_error")
check("non-integer message names the value", grepl("'abc'", conditionMessage(e)))

check_error("negative coordinate rejected", parse_bed("chr1\t-5\t400"), "bed_parse_error")
check_error("float coordinate rejected",    parse_bed("chr1\t1.5\t400"), "bed_parse_error")
check_error("too few fields rejected",      parse_bed("chr1\t100"), "bed_parse_error")

e <- check_error("line number counts skipped lines",
                 parse_bed(c("# hdr", "chr1\t1\t2", "chr1\t9\t3"), "y.bed"), "bed_parse_error")
check("error is on line 3", grepl("^y.bed:3:", conditionMessage(e)))

e <- check_error("first_line offset applied",
                 parse_bed("chr1\t9\t3", "z.bed", first_line = 41L), "bed_parse_error")
check("offset line number is 41", grepl("^z.bed:41:", conditionMessage(e)))

# --- write_bed(): exact bytes -----------------------------------------------

out <- capture.output(write_bed(parse_bed(c("chr1\t100\t200\ta\t0\t+", "chr1\t300\t400"))))
check("write_bed emits tab-separated rows", identical(out, c("chr1\t100\t200\ta\t0\t+", "chr1\t300\t400")))
check("write_bed on zero rows prints nothing",
      length(capture.output(write_bed(parse_bed(character())))) == 0)
check("write_bed never uses scientific notation",
      capture.output(write_bed(parse_bed("chr1\t100000\t2000000"))) == "chr1\t100000\t2000000")

tmp <- tempfile(); on.exit(unlink(tmp))
con <- file(tmp, "wb"); write_bed(parse_bed("chr1\t1\t2"), con); close(con)
check("write_bed uses LF with trailing newline", identical(readBin(tmp, "raw", 100), charToRaw("chr1\t1\t2\n")))

done()
