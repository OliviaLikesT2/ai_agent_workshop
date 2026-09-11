#!/usr/bin/env Rscript
# Unit tests for R/cmd_closest.R. Runs without bedtools: Rscript tests/test_closest.R
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "harness.R"))
source(file.path(.tests_dir, "..", "R", "cmd_closest.R"))

bed <- function(...) parse_bed(c(...))
clo <- function(a, b, d = TRUE) closest_rows(a, closest_index(b), d)

# --- distance rule: 0 on overlap, else gap + 1 (bookended = 1) ----------------

check("overlap is distance 0",
      identical(clo(bed("chr1\t100\t200\tA"), bed("chr1\t150\t250\tB")), "chr1\t100\t200\tA\tchr1\t150\t250\tB\t0"))
check("nested is distance 0",
      identical(clo(bed("chr1\t100\t400\tA"), bed("chr1\t200\t300\tB")), "chr1\t100\t400\tA\tchr1\t200\t300\tB\t0"))
check("identical is distance 0",
      identical(clo(bed("chr1\t100\t200\tA"), bed("chr1\t100\t200\tB")), "chr1\t100\t200\tA\tchr1\t100\t200\tB\t0"))
check("bookended on the left is distance 1 (oracle)",
      identical(clo(bed("chr1\t300\t400\tA"), bed("chr1\t200\t300\tL")), "chr1\t300\t400\tA\tchr1\t200\t300\tL\t1"))
check("bookended on the right is distance 1 (oracle)",
      identical(clo(bed("chr1\t300\t400\tA"), bed("chr1\t400\t500\tR")), "chr1\t300\t400\tA\tchr1\t400\t500\tR\t1"))
check("gap of 20 is distance 21 (oracle: a14/b11)",
      identical(clo(bed("chr2\t100\t120\ta14"), bed("chr2\t140\t160\tb11")), "chr2\t100\t120\ta14\tchr2\t140\t160\tb11\t21"))

# --- nearer side wins; ties report both, in -b order -------------------------

check("nearer right neighbour wins",
      identical(clo(bed("chr1\t100\t200\tA"), bed("chr1\t0\t50\tL", "chr1\t210\t220\tR")), "chr1\t100\t200\tA\tchr1\t210\t220\tR\t11"))
check("nearer left neighbour wins",
      identical(clo(bed("chr1\t100\t200\tA"), bed("chr1\t80\t95\tL", "chr1\t210\t220\tR")), "chr1\t100\t200\tA\tchr1\t80\t95\tL\t6"))
check("equidistant left and right both reported, left first",
      identical(clo(bed("chr1\t10\t20\tA"), bed("chr1\t0\t5\tL", "chr1\t25\t35\tR")),
                c("chr1\t10\t20\tA\tchr1\t0\t5\tL\t6", "chr1\t10\t20\tA\tchr1\t25\t35\tR\t6")))
check("two left ties with the same end both reported",
      identical(clo(bed("chr1\t500\t600\tA"), bed("chr1\t100\t400\tL1", "chr1\t300\t400\tL2")),
                c("chr1\t500\t600\tA\tchr1\t100\t400\tL1\t101", "chr1\t500\t600\tA\tchr1\t300\t400\tL2\t101")))
check("all overlapping -b rows reported, distance 0 each",
      identical(clo(bed("chr1\t190\t195\tA"), bed("chr1\t100\t200\tB1", "chr1\t150\t250\tB2", "chr1\t180\t300\tB3")),
                c("chr1\t190\t195\tA\tchr1\t100\t200\tB1\t0", "chr1\t190\t195\tA\tchr1\t150\t250\tB2\t0", "chr1\t190\t195\tA\tchr1\t180\t300\tB3\t0")))
check("a farther overlapping -b is not hidden by a nearer non-overlapping one",
      identical(clo(bed("chr1\t190\t195\tA"), bed("chr1\t100\t300\tB1", "chr1\t196\t197\tnear")),
                "chr1\t190\t195\tA\tchr1\t100\t300\tB1\t0"))

# --- zero-length: bedtools spans [p-1, p+1) for overlap and distance ---------

check("zero-length -a at 300 vs 400-500 is 100, not 101 (oracle)",
      identical(clo(bed("chr1\t300\t300\tz"), bed("chr1\t400\t500\tR")), "chr1\t300\t300\tz\tchr1\t400\t500\tR\t100"))
check("zero-length -b at 400 right of 300-350 is 50 (oracle)",
      identical(clo(bed("chr1\t300\t350\tA"), bed("chr1\t400\t400\tz")), "chr1\t300\t350\tA\tchr1\t400\t400\tz\t50"))
check("zero-length -b at 400 left of 500-600 is 100 (oracle)",
      identical(clo(bed("chr1\t500\t600\tA"), bed("chr1\t400\t400\tz")), "chr1\t500\t600\tA\tchr1\t400\t400\tz\t100"))
check("zero-length at 100 overlaps 0-100 (oracle: a01/b02 distance 0)",
      identical(clo(bed("chr1\t0\t100\ta01"), bed("chr1\t100\t100\tb02")), "chr1\t0\t100\ta01\tchr1\t100\t100\tb02\t0"))
check("zero-length at 300 overlaps 200-300 (oracle: a16/b12 distance 0)",
      identical(clo(bed("chr2\t300\t300\ta16"), bed("chr2\t200\t300\tb12", "chr2\t400\t500\tb13")), "chr2\t300\t300\ta16\tchr2\t200\t300\tb12\t0"))
check("zero-length at position 0 overlaps 0-10 (oracle: a12/b10)",
      identical(clo(bed("chr2\t0\t0\ta12"), bed("chr2\t0\t10\tb10")), "chr2\t0\t0\ta12\tchr2\t0\t10\tb10\t0"))

# --- position 0 ---------------------------------------------------------------

check("interval at position 0 overlapping",
      identical(clo(bed("chr1\t0\t100\tA"), bed("chr1\t0\t50\tB")), "chr1\t0\t100\tA\tchr1\t0\t50\tB\t0"))

# --- null rows ---------------------------------------------------------------

check("null row BED3", closest_null_row(3L) == ".\t-1\t-1")
check("null row BED4", closest_null_row(4L) == ".\t-1\t-1\t.")
check("null row BED5", closest_null_row(5L) == ".\t-1\t-1\t.\t-1")
check("null row BED6", closest_null_row(6L) == ".\t-1\t-1\t.\t-1\t.")
check("null row 7 columns is all dots (oracle)", closest_null_row(7L) == ".\t-1\t-1\t.\t.\t.\t.")
check("chrom absent from -b gives null row and -1 distance",
      identical(clo(bed("chr9\t1\t2\tq"), bed("chr1\t5\t6\tn\t0\t+")), "chr9\t1\t2\tq\t.\t-1\t-1\t.\t-1\t.\t-1"))
check("without -d no distance column",
      identical(clo(bed("chr9\t1\t2\tq"), bed("chr1\t5\t6"), d = FALSE), "chr9\t1\t2\tq\t.\t-1\t-1"))
check("bed_ncol counts columns from the first row",
      bed_ncol(bed("chr1\t1\t2")) == 3L && bed_ncol(bed("chr1\t1\t2\tn\t0\t+")) == 6L)

# --- order and shape ---------------------------------------------------------

check("-a order preserved across chromosomes",
      identical(clo(bed("chr2\t1\t2\tA", "chr1\t1\t2\tB"), bed("chr2\t5\t6\tX", "chr1\t5\t6\tY")),
                c("chr2\t1\t2\tA\tchr2\t5\t6\tX\t4", "chr1\t1\t2\tB\tchr1\t5\t6\tY\t4")))
check("empty -a gives no lines", length(clo(parse_bed(character()), bed("chr1\t1\t2"))) == 0)

# --- sortedness --------------------------------------------------------------

check("parse -a -b -d", identical(closest_parse_args(c("-d", "-a", "x", "-b", "y")), list(a = "x", b = "y", d = TRUE)))
check("parse without -d", identical(closest_parse_args(c("-a", "-", "-b", "y")), list(a = "-", b = "y", d = FALSE)))
st <- check_sorted(bed("chr1\t1\t2", "chr1\t5\t6", "chr2\t0\t1"), "t")
check("check_sorted accepts sorted input and carries state", identical(st$seen, c("chr1", "chr2")) && st$start == 0L)
st2 <- check_sorted(bed("chr2\t3\t4"), "t", st)
check("check_sorted continues the same chromosome across chunks", identical(st2$seen, c("chr1", "chr2")))

done()
