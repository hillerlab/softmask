#!/usr/bin/env bash
# merge_repeatmasker_repeats.sh
# ---------------------------------------------------------------------------
# Merge per-chunk RepeatMasker .out files into a single tab-separated table.
#
#   * keeps ONE header (a clean, flat one)
#   * renumbers the ID column GLOBALLY, preserving within-chunk fragment
#     grouping (rows that shared an ID inside a chunk keep a shared ID)
#   * adds chrom (chunk name with the :START-END stripped)
#   * adds abs_begin / abs_end -> genome coordinates
#
# Input: ONE file-of-filenames (fofn), one .out path per line, in the order
# you want IDs assigned. This avoids "argument list too long" with many chunks.
#
# Chunk sequences are assumed faidx-style "seqname:START-END" (1-based,
# inclusive START). Names without that suffix -> chrom == query, offset 0.
#
# Usage:
#   merge_repeatmasker_repeats.sh paths.fofn > merged.tsv
# ---------------------------------------------------------------------------
set -euo pipefail

fofn="${1:?usage: merge_repeatmasker_repeats.sh <file-of-paths.fofn>}"

# pick whatever awk is present (gawk in a conda env, mawk in a slim image)
awk_bin="$(command -v gawk || command -v mawk || command -v awk)"

# Stream every listed file into ONE awk process. We do NOT put the paths on a
# command line (ARG_MAX-safe) and we keep a single awk process so the global
# ID offset survives across files. File boundaries are recovered from the
# RepeatMasker header line rather than from separate file arguments.
while IFS= read -r f || [ -n "$f" ]; do
    [ -n "$f" ] || continue            # tolerate blank lines in the fofn
    cat -- "$f"
done < "$fofn" | "$awk_bin" '
BEGIN {
    OFS = "\t"
    id_offset = 0          # added to every ID of the current file
    cur_max   = 0          # max original ID seen in the current file
    print "SW_score","perc_div","perc_del","perc_ins","query",      \
          "begin","end","q_left","strand","repeat","class_family",   \
          "r_begin","r_end","r_left","ID","overlap","chrom","abs_begin","abs_end"
}

# A RepeatMasker header line marks the start of a new file: roll the offset
# forward by the previous file max, then skip the header line itself.
$1 == "SW"    { id_offset += cur_max; cur_max = 0; next }
$1 == "score" { next }
NF == 0       { next }

# Ambiguity-only and repeat-free inputs contain a diagnostic sentence instead
# of the normal table. A real RepeatMasker result row always starts with its
# numeric Smith-Waterman score.
$1 !~ /^[0-9]+$/ { next }

{
    # --- chrom + absolute coordinates ---------------------------------------
    q = $5
    if (match(q, /:[0-9]+-[0-9]+$/)) {            # chunked sequence?
        chrom  = substr(q, 1, RSTART - 1)         # name before the ":"
        split(substr(q, RSTART + 1), c, "-")      # "START-END"
        offset = c[1] - 1                          # 1-based inclusive -> shift
    } else {                                       # un-chunked
        chrom  = q
        offset = 0
    }
    abs_begin = $6 + offset
    abs_end   = $7 + offset

    # --- global ID, preserving fragment grouping -----------------------------
    old_id = $15 + 0
    if (old_id > cur_max) cur_max = old_id
    new_id = old_id + id_offset

    # --- optional trailing "*" (higher-scoring overlapping match) -------------
    overlap = ($16 == "*") ? "*" : ""

    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14, \
          new_id, overlap, chrom, abs_begin, abs_end
}
'
