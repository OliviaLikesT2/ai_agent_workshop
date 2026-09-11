#!/usr/bin/env Rscript
# Unit tests for R/cmd_subtract.R. Runs without bedtools: Rscript tests/test_subtract.R
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "harness.R"))
source(file.path(.tests_dir, "..", "R", "cmd_subtract.R"))

bed <- function(...) parse_bed(c(...))
subtr <- function(a, b) subtract_rows(a, subtract_index(b))
rows <- function(df) paste(df$chrom, df$start, df$end, df$extra, sep = "\t")

# --- union_intervals ---------------------------------------------------------

u <- union_intervals(c(100L, 150L, 400L), c(200L, 300L, 500L))
check("union merges overlapping runs", identical(u$start, c(100L, 400L)) && identical(u$end, c(300L, 500L)))
u <- union_intervals(c(100L, 200L), c(200L, 300L))
check("union joins touching runs", identical(u$start, 100L) && identical(u$end, 300L))
u <- union_intervals(c(300L, 100L), c(400L, 200L))
check("union sorts its input", identical(u$start, c(100L, 300L)))
u <- union_intervals(c(100L, 120L), c(500L, 130L))
check("union handles nested input", identical(u$start, 100L) && identical(u$end, 500L))
check("union of nothing is empty", length(union_intervals(integer(), integer())$start) == 0)

# --- the splitting logic, half-open bounds -----------------------------------

r <- subtr(bed("chr1\t100\t400\tA"), bed("chr1\t200\t300\tB"))
check("middle overlap splits into two", identical(rows(r), c("chr1\t100\t200\tA", "chr1\t300\t400\tA")))

r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t50\t150\tB"))
check("left-edge trim", identical(rows(r), "chr1\t150\t200\tA"))
r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t150\t250\tB"))
check("right-edge trim", identical(rows(r), "chr1\t100\t150\tA"))

r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t100\t200\tB"))
check("identical -b removes the row", nrow(r) == 0)
r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t0\t500\tB"))
check("fully covered row vanishes", nrow(r) == 0)
r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t100\t150\tB1", "chr1\t150\t200\tB2"))
check("covered by union of two vanishes", nrow(r) == 0)

# --- bookended: strict overlap, nothing removed ------------------------------

r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t0\t100\tL", "chr1\t200\t300\tR"))
check("bookended -b removes nothing", identical(rows(r), "chr1\t100\t200\tA"))

# --- position 0 ---------------------------------------------------------------

r <- subtr(bed("chr1\t0\t100\tA"), bed("chr1\t0\t50\tB"))
check("position 0 left trim", identical(rows(r), "chr1\t50\t100\tA"))

# --- zero-length: oracle behaviour (bedtools 2.31.1), see tests/golden/subtract.sh
# bedtools treats a zero-length feature at p as [p-1, p+1) for the overlap test.

r <- subtr(bed("chr1\t0\t100\tA"), bed("chr1\t100\t100\tz"))
check("zero-length -b at 100 removes base 99 (oracle: chr1 0 99)", identical(rows(r), "chr1\t0\t99\tA"))
r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t100\t100\tz"))
check("zero-length -b at 100 removes base 100 (oracle: chr1 101 200)", identical(rows(r), "chr1\t101\t200\tA"))
r <- subtr(bed("chr1\t500\t500\tz"), bed("chr1\t500\t500\tzb"))
check("zero-length -a fully covered by zero-length -b vanishes", nrow(r) == 0)
r <- subtr(bed("chr1\t300\t300\tz"), bed("chr1\t200\t300\tL"))
check("zero-length -a partially hit is printed unchanged", identical(rows(r), "chr1\t300\t300\tz"))
r <- subtr(bed("chr1\t300\t300\tz"), bed("chr1\t200\t300\tL", "chr1\t300\t400\tR"))
check("zero-length -a covered by union vanishes", nrow(r) == 0)
r <- subtr(bed("chr2\t0\t0\tz"), bed("chr2\t0\t10\tB"))
check("zero-length -a at 0 survives (oracle: chr2 0 0)", identical(rows(r), "chr2\t0\t0\tz"))

# --- order, columns, absent chromosomes --------------------------------------

r <- subtr(bed("chr2\t100\t200\tX", "chr1\t100\t400\tA", "chr1\t50\t60\tB"),
         bed("chr1\t200\t300\tB"))
check("-a order preserved, pieces adjacent", identical(rows(r), c("chr2\t100\t200\tX", "chr1\t100\t200\tA", "chr1\t300\t400\tA", "chr1\t50\t60\tB")))
r <- subtr(bed("chr9\t1\t2\tA\t0\t+\textra"), bed("chr1\t1\t2"))
check("chrom absent from -b passes through with all columns", identical(rows(r), "chr9\t1\t2\tA\t0\t+\textra"))
r <- subtr(bed("chr1\t100\t200\tA"), bed("chr1\t120\t130\tB1", "chr1\t150\t160\tB2", "chr1\t125\t155\tB3"))
check("overlapping -b regions subtract as their union", identical(rows(r), c("chr1\t100\t120\tA", "chr1\t160\t200\tA")))
check("empty -a gives empty result", nrow(subtr(parse_bed(character()), bed("chr1\t1\t2"))) == 0)

# --- argument parsing --------------------------------------------------------

check("parse -a -b", identical(subtract_parse_args(c("-a", "x", "-b", "y")), list(a = "x", b = "y")))
check("parse order-independent", identical(subtract_parse_args(c("-b", "y", "-a", "-")), list(a = "-", b = "y")))

done()
