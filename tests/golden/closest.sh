# Golden cases for `mytools closest`. Sourced by tests/run_golden.sh.
# closest needs sorted input and a.bed/b.bed are deliberately unsorted, so the
# fixtures are sorted into $tmp first (with bedtools sort: fixture prep, not the
# thing under test).

bedtools sort -i "$DATA/a.bed" > "$tmp/a.sorted.bed"
bedtools sort -i "$DATA/b.bed" > "$tmp/b.sorted.bed"

check       "closest a b"                 -- closest    -a "$tmp/a.sorted.bed" -b "$tmp/b.sorted.bed"
check       "closest -d a b"              -- closest -d -a "$tmp/a.sorted.bed" -b "$tmp/b.sorted.bed"
check       "closest -d b a (chr3 has no -b)" -- closest -d -a "$tmp/b.sorted.bed" -b "$tmp/a.sorted.bed"
check_stdin "closest -d stdin b" "$tmp/a.sorted.bed" -- closest -d -a - -b "$tmp/b.sorted.bed"
check       "closest unsorted -a: exit code" -- closest -a "$DATA/a.bed" -b "$tmp/b.sorted.bed"
check       "closest unsorted -b: exit code" -- closest -a "$tmp/a.sorted.bed" -b "$DATA/b.bed"

# Null-row shape depends on how many columns -b has.
cut -f1-3 "$tmp/b.sorted.bed" > "$tmp/b3.bed"
cut -f1-4 "$tmp/b.sorted.bed" > "$tmp/b4.bed"
cut -f1-5 "$tmp/b.sorted.bed" > "$tmp/b5.bed"
printf 'chr9\t1\t2\tq\n' > "$tmp/clo_nochrom.bed"
check "closest null row, BED3 -b"      -- closest    -a "$tmp/clo_nochrom.bed" -b "$tmp/b3.bed"
check "closest -d null row, BED3 -b"   -- closest -d -a "$tmp/clo_nochrom.bed" -b "$tmp/b3.bed"
check "closest -d null row, BED4 -b"   -- closest -d -a "$tmp/clo_nochrom.bed" -b "$tmp/b4.bed"
check "closest -d null row, BED5 -b"   -- closest -d -a "$tmp/clo_nochrom.bed" -b "$tmp/b5.bed"
check "closest -d null row, BED6 -b"   -- closest -d -a "$tmp/clo_nochrom.bed" -b "$tmp/b.sorted.bed"
printf 'chr1\t1\t2\tn\t0\t+\tseven\n' > "$tmp/b7.bed"
check "closest -d null row, 7-col -b"  -- closest -d -a "$tmp/clo_nochrom.bed" -b "$tmp/b7.bed"

# One small case per edge the fixtures were built around.
printf 'chr1\t300\t400\tA\n' > "$tmp/clo_a.bed"
printf 'chr1\t200\t300\tL\nchr1\t400\t500\tR\n' > "$tmp/clo_bookended_b.bed"
check "closest -d bookended both sides: distance 1, both reported" -- closest -d -a "$tmp/clo_a.bed" -b "$tmp/clo_bookended_b.bed"

printf 'chr1\t320\t330\tinner\n' > "$tmp/clo_nested_b.bed"
check "closest -d nested: distance 0"  -- closest -d -a "$tmp/clo_a.bed" -b "$tmp/clo_nested_b.bed"
printf 'chr1\t300\t400\tsame\n' > "$tmp/clo_identical_b.bed"
check "closest -d identical: distance 0" -- closest -d -a "$tmp/clo_a.bed" -b "$tmp/clo_identical_b.bed"

printf 'chr1\t300\t300\tz\n' > "$tmp/clo_zero_a.bed"
printf 'chr1\t400\t500\tR\n' > "$tmp/clo_zero_b.bed"
check "closest -d zero-length -a"      -- closest -d -a "$tmp/clo_zero_a.bed" -b "$tmp/clo_zero_b.bed"
printf 'chr1\t400\t400\tz\n' > "$tmp/clo_zero_b2.bed"
check "closest -d zero-length -b right" -- closest -d -a "$tmp/clo_a.bed" -b "$tmp/clo_zero_b2.bed"
printf 'chr1\t500\t600\tA\n' > "$tmp/clo_a2.bed"
check "closest -d zero-length -b left" -- closest -d -a "$tmp/clo_a2.bed" -b "$tmp/clo_zero_b2.bed"
printf 'chr1\t0\t100\tA\n' > "$tmp/clo_pos0_a.bed"
printf 'chr1\t100\t100\tz\n' > "$tmp/clo_pos0_b.bed"
check "closest -d position 0 vs zero-length at 100: overlap" -- closest -d -a "$tmp/clo_pos0_a.bed" -b "$tmp/clo_pos0_b.bed"
printf 'chr1\t0\t0\tz\n' > "$tmp/clo_pos0z_a.bed"
printf 'chr1\t0\t10\tB\n' > "$tmp/clo_pos0z_b.bed"
check "closest -d zero-length at position 0" -- closest -d -a "$tmp/clo_pos0z_a.bed" -b "$tmp/clo_pos0z_b.bed"

printf 'chr1\t100\t400\tL1\nchr1\t300\t400\tL2\n' > "$tmp/clo_tie_b.bed"
check "closest -d two left ties, same end" -- closest -d -a "$tmp/clo_a2.bed" -b "$tmp/clo_tie_b.bed"
printf 'chr1\t0\t5\tL\nchr1\t25\t35\tR\n' > "$tmp/clo_tie2_b.bed"
printf 'chr1\t10\t20\tA\n' > "$tmp/clo_tie2_a.bed"
check "closest -d equidistant left and right" -- closest -d -a "$tmp/clo_tie2_a.bed" -b "$tmp/clo_tie2_b.bed"
printf 'chr1\t100\t200\tB1\nchr1\t150\t250\tB2\nchr1\t180\t300\tB3\n' > "$tmp/clo_multi_b.bed"
printf 'chr1\t190\t195\tA\n' > "$tmp/clo_multi_a.bed"
check "closest -d three overlapping: all reported" -- closest -d -a "$tmp/clo_multi_a.bed" -b "$tmp/clo_multi_b.bed"

printf 'chr2\t1\t2\tA\nchr1\t1\t2\tB\n' > "$tmp/clo_order_a.bed"
printf 'chr2\t5\t6\tA\nchr1\t5\t6\tB\n' > "$tmp/clo_order_b.bed"
check "closest -d non-lexicographic but consistent chrom order" -- closest -d -a "$tmp/clo_order_a.bed" -b "$tmp/clo_order_b.bed"

# Usage errors are not golden-tested: bedtools exits 1 for them, SPEC.md says 2.
