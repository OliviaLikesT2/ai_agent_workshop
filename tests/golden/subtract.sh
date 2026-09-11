# Golden cases for `mytools subtract`. Sourced by tests/run_golden.sh.

check       "subtract a.bed - b.bed"        -- subtract -a "$DATA/a.bed" -b "$DATA/b.bed"
check_stdin "subtract stdin - b.bed" "$DATA/a.bed" -- subtract -a - -b "$DATA/b.bed"

# bedtools 2.31.1 crashes ("illegal bin number -1") when -b holds a zero-length
# feature at position 0 (a12: chr2 0 0), so the swapped and self cases use a.bed
# with that one line removed. a07 (chr1 500 500) is still in there.
grep -v $'^chr2\t0\t0\t' "$DATA/a.bed" > "$tmp/a_no_a12.bed"
check "subtract b.bed - a.bed (minus a12)"  -- subtract -a "$DATA/b.bed" -b "$tmp/a_no_a12.bed"
check "subtract a.bed - a.bed (minus a12)"  -- subtract -a "$DATA/a.bed" -b "$tmp/a_no_a12.bed"

# One small case per edge the fixtures were built around, so a failure names it.
printf 'chr1\t100\t200\tA\n'                 > "$tmp/sub_bookended_a.bed"
printf 'chr1\t0\t100\tL\nchr1\t200\t300\tR\n' > "$tmp/sub_bookended_b.bed"
check "subtract bookended: nothing removed" -- subtract -a "$tmp/sub_bookended_a.bed" -b "$tmp/sub_bookended_b.bed"

printf 'chr1\t100\t400\tA\n'                 > "$tmp/sub_nested_a.bed"
printf 'chr1\t200\t300\tB\n'                 > "$tmp/sub_nested_b.bed"
check "subtract nested: split in two"       -- subtract -a "$tmp/sub_nested_a.bed" -b "$tmp/sub_nested_b.bed"
check "subtract nested, swapped: vanishes"  -- subtract -a "$tmp/sub_nested_b.bed" -b "$tmp/sub_nested_a.bed"

printf 'chr1\t100\t200\tA\nchr1\t100\t200\tA2\n' > "$tmp/sub_identical_a.bed"
printf 'chr1\t100\t200\tB\n'                 > "$tmp/sub_identical_b.bed"
check "subtract identical: both vanish"     -- subtract -a "$tmp/sub_identical_a.bed" -b "$tmp/sub_identical_b.bed"

printf 'chr1\t500\t500\tz\nchr1\t500\t600\tA\nchr1\t400\t500\tL\n' > "$tmp/sub_zero_a.bed"
printf 'chr1\t500\t500\tzb\n'                > "$tmp/sub_zero_b.bed"
check "subtract zero-length in -b"          -- subtract -a "$tmp/sub_zero_a.bed" -b "$tmp/sub_zero_b.bed"
printf 'chr1\t300\t300\tz\n'                 > "$tmp/sub_zero2_a.bed"
printf 'chr1\t200\t300\tL\nchr1\t300\t400\tR\n' > "$tmp/sub_zero2_b.bed"
check "subtract zero-length -a covered by union" -- subtract -a "$tmp/sub_zero2_a.bed" -b "$tmp/sub_zero2_b.bed"
printf 'chr1\t250\t260\tB\n'                 > "$tmp/sub_zero3_b.bed"
check "subtract zero-length -a untouched"   -- subtract -a "$tmp/sub_zero2_a.bed" -b "$tmp/sub_zero3_b.bed"

printf 'chr1\t0\t100\tA\n'                   > "$tmp/sub_pos0_a.bed"
printf 'chr1\t0\t50\tB\n'                    > "$tmp/sub_pos0_b.bed"
check "subtract position 0: left trim"      -- subtract -a "$tmp/sub_pos0_a.bed" -b "$tmp/sub_pos0_b.bed"

printf 'chr1\t100\t200\tA\n'                 > "$tmp/sub_union_a.bed"
printf 'chr1\t100\t150\tB1\nchr1\t150\t200\tB2\n' > "$tmp/sub_union_b.bed"
check "subtract covered by union of two"    -- subtract -a "$tmp/sub_union_a.bed" -b "$tmp/sub_union_b.bed"

printf 'chr9\t1\t2\tA\n'                     > "$tmp/sub_nochrom_a.bed"
check "subtract chrom absent from -b"       -- subtract -a "$tmp/sub_nochrom_a.bed" -b "$DATA/b.bed"

# Usage errors are not golden-tested: bedtools exits 1 for them, SPEC.md says 2.
