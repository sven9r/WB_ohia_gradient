#!/usr/bin/env bash

# Primer and inline-barcode trimming engine.
#
# This script processes sequences that already have the forward and reverse
# reads merged and oriented as:
#
#   5' [inline barcode] [forward primer] [insert] [reverse-primer match] 3'
#
# It retains a record only when the complete barcode/primer pattern matches
# the start and end of the sequence. It removes those matching flanks and
# filters on the LENGTH OF THE REMAINING INSERT. This is exact, case-sensitive
# regular-expression matching, without substitutions, indels, alignment scores,
# or an allowed-mismatch setting. A '.' in a primer matches any one character.
# No record is deduplicated, clustered, denoised, or quality-filtered here.
#
# Runtime: Bash >= 4 (the barcode lookup uses an associative array), awk, grep,
# sed, tr, rev, basename, date, and gzip for APSCALE mode. Run with this script's
# shebang or a Bash >= 4 interpreter, not with `sh`. Linux installations usually
# provide a suitable Bash. macOS's system /bin/bash is often version 3.2; use a
# separately installed Bash >= 4 on PATH or invoke its absolute path explicitly.
# `#!/usr/bin/env bash` uses the first Bash on PATH; verify with `bash --version`.
#
# UNIX mode (default): reads *.fasta in the CURRENT WORKING DIRECTORY, makes or
# reuses a nonl_fasta/ cache, and writes trimmed FASTA with replacement headers.
# Example, with primers supplied in their usual oligonucleotide orientations:
#   bash trim_fasta.sh -i ATCACG -f ACGTACGT -r TGCATGCA \
#     -o 1_marker.fas -l 50 -m 150
# The example primers are illustrative, not study-specific primer sequences.
#
# APSCALE mode: reads exactly one -s file named *_PE.fastq.gz and writes gzipped
# FASTQ while retaining the record identifiers and slicing qualities with bases.
# It expects ordinary four-line FASTQ records; wrapped FASTQ is not parsed.
#   bash trim_fasta.sh -p APSCALE -s sample_PE.fastq.gz \
#     -i ATCACG -f ACGTACGT -r TGCATGCA -o trimmed.fastq.gz -l 50 -m 150
#
# Use -c ONLY if the supplied reverse primer is already written as it appears
# at the right end of the merged sequence. Otherwise the script calculates its
# reverse complement. Use -w to translate uppercase ambiguity codes to '.',
# understanding that this accepts ANY character rather than the code's IUPAC
# subset. For example, R becomes '.', not '[AG]'.
#
# For reads without an inline barcode, omit -i. The special value X only selects
# the '_nobarcode' HEADER LABEL below; X remains literal in the matching pattern
# and counts as one character when trimming. An unrecognized barcode is matched
# as supplied and receives the label UNKNOWN; it is not automatically rejected.
#
# Output files already present trigger an interactive overwrite prompt. Existing
# FASTA cache files are never refreshed automatically: use a separate working
# directory or refresh the cache after changing source FASTA files. No output
# directory is created for -o; its parent directory must already exist.
#
# The executable commands are preserved from the source engine; the explanatory
# comments below describe that engine's actual data handling.


# CHANGELOG:
# 12/16/2024:
# - Added functionality to map inline barcodes to predefined abbreviations (A, B, C, etc.).
# - Updated suffix generation for headers to reflect the mapped abbreviation.

# 12/28/2024:
# - Reverted to original logic for matching sequences and pattern generation.
# - Retained fixes for header naming conventions and .gz removal.
# - No log file is generated unless explicitly required.
# - Added special handling for 'X' as no inline barcode, assigned '_nobarcode' in headers.
# - Fixed bad array subscript error when inline_barcode is invalid or empty.
# 02/19/2026:
# - Added explicit pipeline switching with `-p UNIX|APSCALE` so the same
#   trimming engine can support both legacy UNIX FASTA output and APSCALE
#   single-file FASTQ output workflows.
# - Added APSCALE input binding via `-s <input_file>` and strict filename
#   validation requiring `_PE.fastq.gz`, preventing accidental processing of
#   incompatible files.
# - Implemented FASTQ-native APSCALE trimming path:
#   reads gzipped 4-line FASTQ records, matches
#   `inline + forward_primer + amplicon + reverse_primer`,
#   trims both sequence and quality strings consistently, and writes
#   gzipped FASTQ output.
# - Kept UNIX mode behavior unchanged for backward compatibility:
#   scanning `*.fasta`, converting to `nonl_fasta/*_nonl.fasta`, and writing
#   trimmed FASTA records.
# - Fixed output creation bug in `check_file_exists()` where non-existing
#   targets previously returned non-zero and skipped writing entirely.
#   New behavior now allows creation of new output files by default.
# - Extended usage help text to document new APSCALE flags and constraints.
# - Clarified wrapper compatibility: callers must use this updated script
#   version when passing `-p`, otherwise older copies without APSCALE options
#   will fail with "invalid option -p".

# An absent destination is allowed immediately. For an existing file, Enter
# accepts the default overwrite, Y/y accepts, and N/n skips the write block.
# A different answer repeats the prompt. This makes existing destinations
# interactive; the trailing completion message is printed even if writing was
# declined, so check the actual output when running the script manually.
# Function to check if a file exists and prompt for overwrite if necessary
check_file_exists() {
    local file="$1"
    if [ -f "$file" ]; then
        read -p "$file already exists. Do you want to overwrite it? (Y/n): " choice
        case "$choice" in
            y|Y|"") return 0;;
            n|N) return 1;;
            *) echo "Invalid choice. Please enter Y or N."; check_file_exists "$file";;
        esac
    else
        return 0
    fi
}

# Only A/T/C/G (in either case) are complemented by tr; all other characters
# pass through unchanged before rev reverses the string. With -w, ambiguous
# uppercase primer bases have already become '.', which remains a wildcard
# after reversal. Input primers/read sequences must use matching letter case.
# Function to calculate reverse complement of a DNA sequence
reverse_complement() {
    local sequence="$1"
    echo "$sequence" | tr 'ATCGatcg' 'TAGCtagc' | rev
}

# This is deliberately a permissive wildcard translation, NOT a full IUPAC
# expansion. Each uppercase R, Y, S, W, K, M, B, D, H, V, or N becomes one '.',
# which the later regular expression treats as any single character. Lowercase
# ambiguity codes are not translated. The barcode is not passed through here.
# Function to translate ambiguous bases to dot (.)
translate_wobble_bases() {
    local sequence="$1"
    echo "$sequence" | tr 'RYSWKMBDHVN' '.'
}

# Flags accepted by getopts below:
#   -i optional prefix barcode; only one barcode is matched in this invocation.
#   -f/-r required primer strings; -r is reverse-complemented unless -c is set.
#   -o required output path; FASTA in UNIX mode, gzip FASTQ in APSCALE mode.
#   -l required minimum insert length, INCLUSIVE after both flanks are removed.
#   -m optional maximum insert length, also INCLUSIVE; omitted means no maximum.
#   -c a switch with no argument: -r already has the observed read orientation.
#   -w a switch with no argument: turn uppercase primer ambiguity codes into '.'.
#   -p UNIX or APSCALE, case-insensitive; UNIX is used when -p is omitted.
#   -s required only for APSCALE: one existing file ending in _PE.fastq.gz.
# Primers and barcodes are interpolated directly into regular expressions;
# provide literal DNA characters and optional single-position '.' wildcards.
# Other regex syntax would change matching while trim lengths still use the
# literal string lengths, so this is not an arbitrary-regex primer interface.
# Usage instructions
usage() {
    echo "Usage: $0 -i <inline_barcode> -f <forward_primer> -r <reverse_primer> -o <output_file> -l <min_amplicon_size> [-m <max_amplicon_size>] [-c] [-w] [-p <pipeline>] [-s <input_file>]"
    echo "Options:"
    echo "  -i <inline_barcode>       Inline barcode sequence (default: none)"
    echo "  -f <forward_primer>       Forward primer sequence (A, T, C, G, .)"
    echo "  -r <reverse_primer>       Reverse primer sequence (A, T, C, G, .)"
    echo "  -o <output_file>          Output file"
    echo "  -l <min_amplicon_size>    Minimum amplicon size (e.g., 50)"
    echo "  -m <max_amplicon_size>    Maximum amplicon size (optional)"
    echo "  -c                        Flag to indicate reverse primer is already in reverse complement"
    echo "  -w                        Allow wobble bases (ambiguous bases) and translate them to '.'"
    echo "  -p <pipeline>             Pipeline mode: UNIX (default) or APSCALE"
    echo "  -s <input_file>           Input file for APSCALE mode (must end with _PE.fastq.gz)"
    echo ""
    exit 1
}

# The Bash getopts specification distinguishes arguments (a following colon)
# from switches. Missing required option values and unknown options show help.
# Parse command line options
reverse_complement_flag=false
wobble_flag=false
pipeline="UNIX"
while getopts ":i:f:r:o:l:m:cwp:s:" opt; do
    case $opt in
        i) inline_barcode="$OPTARG";;
        f) forward_primer="$OPTARG";;
        r) reverse_primer="$OPTARG";;
        o) output_file="$OPTARG";;
        l) min_amplicon_size="$OPTARG";;
        m) max_amplicon_size="$OPTARG";;
        c) reverse_complement_flag=true;;
        w) wobble_flag=true;;
        p) pipeline="$OPTARG";;
        s) input_file="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; usage;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    esac
done

# Presence checks below do not validate sequence alphabets, numeric bounds,
# gzip integrity, or complete FASTQ structure. Supply integer insert lengths
# and input in the formats described above.
# Check if required arguments are provided
if [ -z "$forward_primer" ] || [ -z "$reverse_primer" ] || [ -z "$output_file" ] || [ -z "$min_amplicon_size" ]; then
    echo "Missing required arguments."
    usage
fi

# Normalize the mode only; barcode, primer, and read case are left unchanged.
pipeline=$(echo "$pipeline" | tr '[:lower:]' '[:upper:]')
if [ "$pipeline" != "UNIX" ] && [ "$pipeline" != "APSCALE" ]; then
    echo "Invalid pipeline mode: $pipeline"
    usage
fi

if [ "$pipeline" = "APSCALE" ]; then
    if [ -z "$input_file" ]; then
        echo "APSCALE mode requires -s <input_file>."
        usage
    fi
    if [[ "$input_file" != *_PE.fastq.gz ]]; then
        echo "Invalid input file '$input_file'. APSCALE requires files ending with _PE.fastq.gz."
        exit 1
    fi
    if [ ! -f "$input_file" ]; then
        echo "Input file not found: $input_file"
        exit 1
    fi
fi

# UNIX replacement headers use the output basename as their marker label:
# strip a terminal .fas and then a leading 1_. For example, 1_COI.fas gives
# COI. Different output suffixes remain part of that label.
# Extract primer suffix from the output file
primer_suffix=$(basename "$output_file" .fas | sed 's/^1_//')

# Translate BEFORE reversing the reverse primer, so all ambiguity positions
# become orientation-independent dots. Without -w, ambiguity letters are literal
# letters in the matching pattern (not degenerate nucleotide alternatives).
# Translate wobble bases if wobble flag is set
if [ "$wobble_flag" = true ]; then
    forward_primer=$(translate_wobble_bases "$forward_primer")
    reverse_primer=$(translate_wobble_bases "$reverse_primer")
fi

# Reverse complement the reverse primer if the -c flag is not set
if [ "$reverse_complement_flag" = false ]; then
    reverse_primer=$(reverse_complement "$reverse_primer")
fi

# These six predefined barcode strings affect naming through this lookup. Matching
# still uses the complete -i sequence, including for values outside this table.
# An associative array requires Bash >= 4; this is why POSIX sh is unsuitable.
# Map inline barcodes to abbreviations
declare -A barcode_map
barcode_map=(
    ["ATCACG"]="A"
    ["CGATGT"]="B"
    ["TTAGGC"]="C"
    ["TGACCA"]="D"
    ["GCCAAT"]="E"
    ["CAGATC"]="F"
)

# This chooses a FASTA header abbreviation, not a pass/fail barcode check.
# No -i (empty) yields an empty prefix match and the _nobarcode label. X selects
# the same label but, as noted above, still participates in matching literally.
# The UNIX header adds '_' before this abbreviation, so a barcode-free label
# has two underscores before nobarcode.
# Validate inline barcode and assign abbreviation
if [ -z "$inline_barcode" ] || [ "$inline_barcode" = "X" ]; then
    barcode_abbreviation="_nobarcode"
elif [ -n "${barcode_map[$inline_barcode]}" ]; then
    barcode_abbreviation="${barcode_map[$inline_barcode]}"
else
    barcode_abbreviation="UNKNOWN"
fi

if [ "$pipeline" = "APSCALE" ]; then
    if check_file_exists "$output_file"; then
        # Read gzip-compressed FASTQ as groups of four lines: identifier, sequence,
        # separator, qualities. The pattern is anchored with ^ and $; .* allows any
        # insert between the flanks. Only forward-oriented records matching BOTH ends
        # are retained. This block does not search the opposite strand or tolerate
        # primer mismatches (apart from any '.' wildcards supplied/created above).
        #
        # total_forward_len is barcode length + forward-primer length; trim_len is the
        # full sequence length minus both flanks. Keep min <= trim_len <= max, with
        # the upper comparison skipped when -m is omitted. Identical substr offsets
        # remove barcode/primers from the sequence and its quality string together.
        # The original identifier and '+' line are printed unchanged; no new FASTQ
        # identifiers or quality values are made. gzip -c compresses the output even
        # if the chosen destination name does not itself end in .gz.
        gzip -dc "$input_file" | awk \
            -v inline="$inline_barcode" \
            -v fwd="$forward_primer" \
            -v rev="$reverse_primer" \
            -v min_size="$min_amplicon_size" \
            -v max_size="${max_amplicon_size:-}" \
            -v total_forward_len="$(( ${#inline_barcode} + ${#forward_primer} ))" \
            -v len_reverse_primer="${#reverse_primer}" \
            '
            NR % 4 == 1 { h = $0 }
            NR % 4 == 2 { seq = $0 }
            NR % 4 == 3 { p = $0 }
            NR % 4 == 0 {
                q = $0
                if (seq ~ ("^" inline fwd ".*" rev "$")) {
                    trim_len = length(seq) - total_forward_len - len_reverse_primer
                    if (trim_len >= min_size && (max_size == "" || trim_len <= max_size)) {
                        trimmed_seq = substr(seq, total_forward_len + 1, trim_len)
                        trimmed_qual = substr(q, total_forward_len + 1, trim_len)
                        print h
                        print trimmed_seq
                        print p
                        print trimmed_qual
                    }
                }
            }' | gzip -c > "$output_file"
    fi
else
    # UNIX mode assumes a working directory containing plain-text *.fasta.
    # Its cache is shared by every invocation from that directory, allowing
    # repeated barcode/primer calls to reuse the same flattened sequences.
    # Create folder for nonl_fasta
    mkdir -p nonl_fasta

    # Flatten multiline FASTA: keep each >header on its own line and join all
    # following sequence lines until the next header. The output remains FASTA,
    # with one full sequence per line for the later line-based grep.
    # Conversion runs ONLY if that cache filename does not already exist.
    # Thus editing a source FASTA alone does not update a previous cache file.
    # Convert .fasta files
    for file in *.fasta; do
        nonl_file="nonl_fasta/$(basename "$file" .fasta)_nonl.fasta"
        if [ ! -f "$nonl_file" ]; then
            echo "Processing file: $file"
            awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$file" > "$nonl_file"
        fi
    done


    # Check if output file already exists
    if check_file_exists "$output_file"; then
        > "$output_file"

        # This UNIX block always creates/overwrites name_transformation.log
        # once output writing is accepted. It records filename-to-sample-label
        # transformations, not per-read matches, counts, or processing quality.
        # This current behavior supersedes the older changelog's no-log remark.
        # Create or clear the logfile at the start
        echo "Name Transformation Log - $(date)" > name_transformation.log

        # Process ALL cached *_nonl.fasta files, including cached files whose
        # original source FASTA may no longer be present in the directory.
        # Sample labels derive from cache basenames, not original FASTA IDs:
        # remove _nonl.fasta and any .gz text, then turn '-' into '_'.
        for file in nonl_fasta/*_nonl.fasta; do
            base_name=$(basename "$file" | sed 's/_nonl\.fasta//; s/\.gz//; s/-/_/g')

            # Log the transformation
            echo "Original: $file -> Transformed: $base_name" >> name_transformation.log

            echo "Searching file: $file"
            # grep -E matches the complete sequence line, anchored at both
            # ends. There is no mismatch allowance, strand search, or minimum
            # overlap setting. Unmatched records are silently discarded.
            # Since only sequence lines are carried forward, the original
            # FASTA record IDs are not available in the output below.
            matching_sequences=$(grep -E "^${inline_barcode}${forward_primer}.*${reverse_primer}$" "$file")

            # For each match, remove the prefix barcode + forward primer and
            # suffix reverse primer using their string lengths. The bounds
            # apply to the remaining insert, inclusively; omitting -m removes
            # only the upper bound. Each accepted input occurrence is written:
            # this is trimming/filtering, with no dereplication step.
            #
            # Replacement header: >sample__marker_barcode, for example
            # >sample_01__COI_A from sample-01.fasta, 1_COI.fas, and ATCACG.
            # Reads from the same source/barcode/marker receive the SAME header;
            # neither their original IDs nor unique per-read numbers are kept.
            if [ -n "$matching_sequences" ]; then
                echo "$matching_sequences" | awk -v min_size="$min_amplicon_size" -v max_size="${max_amplicon_size:-}" \
                    -v sample="$base_name" -v primer_suffix="$primer_suffix" -v barcode="$barcode_abbreviation" \
                    -v total_forward_len="$(( ${#inline_barcode} + ${#forward_primer} ))" -v len_reverse_primer="${#reverse_primer}" \
                    '{
                        trimmed_seq = substr($0, total_forward_len + 1, length($0) - total_forward_len - len_reverse_primer);
                        seq_length = length(trimmed_seq);
                        if (seq_length >= min_size && (max_size == "" || seq_length <= max_size)) {
                            header = ">" sample "__" primer_suffix "_" barcode;
                            print header "\n" trimmed_seq;
                        }
                    }' >> "$output_file"
            fi
        done
    fi
fi

echo "Processing completed. Results written to ${output_file}."
