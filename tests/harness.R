# Minimal unit-test harness, base R only. Source this, then:
#   check("name", condition)      -- records pass/fail
#   check_error("name", expr, class = "bed_parse_error")  -- expects expr to signal
#   done()                        -- prints summary, exits non-zero on any failure
# Sources R/bed.R relative to the tests/ directory so tests run from anywhere.

.tests_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
})
source(file.path(.tests_dir, "..", "R", "bed.R"))

.results <- new.env()
.results$pass <- 0L
.results$fail <- 0L

check <- function(name, cond) {
  ok <- isTRUE(cond)
  cat(if (ok) "ok   " else "FAIL ", name, "\n", sep = "")
  if (ok) .results$pass <- .results$pass + 1L else .results$fail <- .results$fail + 1L
  invisible(ok)
}

check_error <- function(name, expr, class = "error") {
  caught <- tryCatch({ force(expr); NULL }, condition = function(e) e)
  check(name, !is.null(caught) && inherits(caught, class))
  invisible(caught)
}

done <- function() {
  cat("---\n", .results$pass, " passed, ", .results$fail, " failed\n", sep = "")
  quit(save = "no", status = if (.results$fail == 0L) 0L else 1L)
}
