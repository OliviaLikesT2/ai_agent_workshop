# SPEC.md — mytools

A small reimplementation of a subset of bedtools, in R. The governing rule for every
decision below: **do what bedtools does.** Real `bedtools` on the files in `data/` is
the oracle; where this document and bedtools disagree, bedtools is right and this
document is wrong.

## 1. Scope

v1 ships five subcommands:

| Subcommand  | Flags in v1              | Rationale |
|-------------|--------------------------|-----------|
| `sort`      | (none)                   | prerequisite for `merge` and `closest` |
| `merge`     | `-d <int>`               | the one flag that changes bookended behaviour |
| `intersect` | `-u`, `-v`, `-wa`        | the three most-used output modes |
| `subtract`  | (none)                   | shares the overlap predicate with `intersect` |
| `closest`   | `-d`                     | distance reporting is what makes it useful |

Explicitly **not** in v1: `-s`/`-S` strand awareness, `-f`/`-r` minimum-overlap
fractions, `-t` tie handling for `closest` (bedtools' default `all` applies), BED12,
GFF3, VCF, `-header`, `-g` genome files, `-sorted`/`-g` for `intersect`, any
subcommand not listed above. Say no to all of it.

## 2. Invocation

    mytools sort      -i <file|->
    mytools merge     -i <file|-> [-d N]
    mytools intersect -a <file|-> -b <file> [-u | -v | -wa]
    mytools subtract  -a <file|-> -b <file>
    mytools closest   -a <file|-> -b <file> [-d]
    mytools --version

Flag names and meanings match bedtools exactly. `-` means stdin; at most one input
may be `-`. `-b` must be a real file (it is read fully; see §6).

## 3. Input format

BED3 through BED6, tab-separated:

    chrom  start  end  [name  score  strand]

- Column count may vary between lines; treat missing trailing columns as absent.
  Columns beyond the sixth are carried through untouched, never interpreted.
- Files ending in `.gz` are read through gzip transparently.
- Lines beginning with `#`, `track`, or `browser` are skipped silently.
- Blank lines are skipped.
- `start` and `end` are non-negative integers. `start > end` is an error (§7).
- `start == end` (zero-length) is **legal**. See §4.

## 4. Interval semantics

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199.

- Overlap: `a.start < b.end AND b.start < a.end`. Strict `<` on both sides.
- Bookended intervals (`a.end == b.start`) do **not** overlap.
- Bookended intervals **do** merge at `-d 0` (bedtools' default), because `merge`
  joins features whose gap is `<= d`, and the gap here is 0.
- Minimum overlap to count: 1 base. No fractional thresholds in v1.
- Zero-length intervals are legal and real `bedtools` treats them in ways you will
  not predict. Do not reason about them from first principles — run bedtools on
  `data/a.bed` and match whatever it prints. Two of them are in the fixtures
  specifically so you find this out early.
- `closest` distance: 0 for overlapping features; otherwise the gap between the
  nearer ends, as bedtools computes it. Ties are all reported (bedtools `-t all`).
  No feature on the same chromosome: bedtools prints `.` fields and `-1`; match it.

## 5. Output

- Tab-separated, LF line endings, trailing newline on the final line.
- Empty result: print nothing, exit 0.
- `sort`: all input columns preserved. Order is by chrom (lexicographic — `chr17`
  before `chr7`), then `start`, then `end`.
- `merge`: BED3 only (`chrom start end`). Input columns are dropped.
- `intersect`: default prints the intersected region carrying `-a`'s trailing
  columns. `-wa` prints `-a`'s original interval, once per overlapping `-b` feature.
  `-u` prints each `-a` feature at most once. `-v` prints `-a` features with no
  overlap. `-u`, `-v` and `-wa` are mutually exclusive.
- `subtract`: `-a` features with `-b` regions removed. One feature may become two,
  or vanish entirely.
- `closest`: each `-a` row followed by its nearest `-b` row on the same line; with
  `-d`, a final integer distance column.
- Input order is preserved for `intersect` and `subtract` (bedtools does not sort
  `-a` for you, and neither do we).

## 6. Memory model

Like bedtools: stream where sorted input makes it possible, hold in memory where it
doesn't.

- `sort` holds the whole input in memory. It is the only subcommand allowed to.
- `merge` streams, assuming sorted input. It does **not** sort for you; unsorted
  input is an error (§7), same as bedtools.
- `intersect` and `subtract` load `-b` into memory (grouped by chrom, sorted by
  start) and stream `-a` line by line.
- `closest` requires both inputs sorted and streams them together; unsorted input is
  an error (§7), same as bedtools.
- Target: inputs up to ~10^6 intervals — the ~500,000-interval `reads.bed` from
  `bedtools bamtobed` on `/data/HG002.neighbourhoods.bam` must run to completion in
  reasonable time and memory. Slower than bedtools is fine; quadratic is not.
- No mmap, no index files, no threads, no third-party packages.

## 7. Errors and exit codes

Errors go to **stderr**. stdout carries data only.

| Situation                              | exit |
|----------------------------------------|------|
| Success (including empty output)       | 0    |
| Malformed line / non-integer coords    | 1    |
| `start > end`                          | 1    |
| Unsorted input to `merge` or `closest` | 1    |
| Unknown flag / missing required arg    | 2    |
| Mutually exclusive flags together      | 2    |
| Input file does not exist              | 2    |
| No arguments at all                    | 2 (prints usage) |
| `--version`, `--help`                  | 0    |

Data problems are 1. Caller problems are 2. Messages name the file and line number:
`a.bed:14: start > end (500 > 400)`.

## 8. Correctness

Real `bedtools` is the oracle. Every one of these must produce byte-identical output
to its bedtools equivalent on the files in `data/`:

    mytools sort -i data/a.bed                       == bedtools sort -i data/a.bed
    mytools sort -i data/b.bed                       == bedtools sort -i data/b.bed
    mytools merge -i <sorted a.bed>                  == bedtools merge -i <sorted a.bed>
    mytools merge -d 10 -i <sorted a.bed>            == bedtools merge -d 10 -i <sorted a.bed>
    mytools intersect -a data/a.bed -b data/b.bed    == bedtools intersect -a ... -b ...
    mytools intersect -u  ...                        == bedtools intersect -u ...
    mytools intersect -v  ...                        == bedtools intersect -v ...
    mytools intersect -wa ...                        == bedtools intersect -wa ...
    mytools subtract -a data/a.bed -b data/b.bed     == bedtools subtract -a ... -b ...
    mytools closest -a <sorted a> -b <sorted b>      == bedtools closest -a ... -b ...
    mytools closest -d -a <sorted a> -b <sorted b>   == bedtools closest -d -a ... -b ...

Each case also runs once reading `-a`/`-i` from stdin via `-`, and compares exit
codes, not just stdout. Also required: `mytools --version` prints a version and
exits 0.

Accepted deviations from bedtools: none. If you find one you cannot fix, write it
down here with the reason.

## 9. Language and layout

- **Implemented in R**, invoked via `Rscript`. Every subcommand and every test. No
  second language without asking (see `CLAUDE.md`).
- Entry point: the executable `mytools` script in the repo root. The golden tests
  invoke it as `./mytools`; put the repo root on `PATH` if you want the bare name.
- Base R only — no CRAN dependencies at runtime or in tests.
- Tests live in `tests/`: golden tests driven by `tests/run_golden.sh`, unit tests
  runnable without bedtools installed.
