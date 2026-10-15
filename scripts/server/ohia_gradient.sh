#!/usr/bin/env bash
# =============================================================================
# ANML sequence processing: annotated record of the February 2025 study workflow
# =============================================================================
#
# Provenance
# ----------
# This is a documented copy of the study-specific `ohia_gradient.sh` command log
# supplied by the author. The original file is unchanged. The publication archive
# retains the ANML branch used for the paper; commands for the other marker in the
# original working log are omitted.
#
# This file records a MANUAL workflow. Read and run the relevant sections in an
# interactive shell with the appropriate tools and input files. It is not a
# tested, end-to-end pipeline, and the interactive trimming step requires user
# selections. In particular, it cannot reconstruct the complete trimming setup
# from the command log alone.
#
# Changes made for publication
# ----------------------------
# * Added this explanation and comments around the original ANML commands.
# * Converted the original free-text title and interactive-selection notes into
#   comments so that the file has valid Bash syntax.
# * Replaced machine-specific executable/database paths with configurable
#   environment variables, and sample-file paths with relative paths.
# * Retained the original tools, ANML filenames, options, and numerical settings.
#   No filtering, denoising, or trimming settings have been inferred or added.
#
# External requirements
# ---------------------
# PEAR, VSEARCH, USEARCH, BLAST+, GNU-compatible `sed -r`, standard Unix utilities,
# and the custom programs `trimflow.py` and `blast2taxonomy_v1.4.2.py` are used.
# trimflow.py, its matching trim_fasta.sh, configuration, trimming environment
# and blast2taxonomy_v1.4.2.py are bundled in this directory.
# Taxonomy-tool dependencies and documentation: https://github.com/yschoeneberg/blast2taxonomy
# The original USEARCH binary was named usearch11.0.667.
#
# The original log activated these existing mamba environments at the indicated
# stages: pear, trimming_env, and blast2tax. These names do not specify the exact
# package versions. Change environment names only to match your own installation.
#
# Before running commands, set the external database path, e.g.:
#   export NT_DATABASE="/path/to/the/nt/database/prefix"
# TRIMFLOW_PY and BLAST2TAXONOMY_PY default to the bundled scripts; set either
# variable explicitly to use another copy.
# Optionally set USEARCH_BIN if the USEARCH executable is not called `usearch`.
# The database prefix must refer to the February 2025 nt snapshot used for the
# study if the original reference searches are to be reproduced. The original
# database directory was labelled `25_02_04_ncbi_database`.
#
# Inputs and working directories
# ------------------------------
# Begin sections 1 and 2 in the directory containing the demultiplexed paired-end
# FASTQ.gz files. The filename-generation command expects the original Illumina
# naming pattern, including _R1_/_R2_ and _S<number>_L001_R2_001.fastq, and filenames
# without whitespace. Inspect the generated pear.sh and fastq.sh before running
# them. Section 3 uses the bundled interactive trimming program. Sections 4–6
# run from the sequence directory containing its `2_merged_ANML.fas` output.
#
# Locate the bundled wrapper independently of the directory containing reads.
# This leaves the working directory unchanged for the commands below.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TRIMFLOW_PY="${TRIMFLOW_PY:-$SCRIPT_DIR/trimflow.py}"
BLAST2TAXONOMY_PY="${BLAST2TAXONOMY_PY:-$SCRIPT_DIR/blast2taxonomy_v1.4.2.py}"
if [ ! -f "$TRIMFLOW_PY" ] || [ ! -x "$TRIMFLOW_PY" ]; then
    printf '%s\n' "Missing or non-executable trimming program: $TRIMFLOW_PY" >&2
    exit 1
fi

USEARCH_BIN="${USEARCH_BIN:-usearch}"

# -----------------------------------------------------------------------------
# 1. Merge paired-end reads with PEAR
# -----------------------------------------------------------------------------
mamba activate pear

# Generate one PEAR command for each reverse-read filename. The substitutions
# construct the matching forward-read filename and the output prefix.
# Original settings: minimum overlap 50 bp (-v), quality threshold 20 (-q),
# and 20 threads (-j). The generated script is kept as an intermediate record.
ls *fastq.gz | awk '{print "pear -v 50 -q 20 -j 20 -f "$0" -r "$0" -o "$0}' | grep _R2_ | sed 's/_R2_/_R1_/' | sed -r 's/_S[0-9]{1,3}_L001_R2_001.fastq//2' > pear.sh
chmod u+x pear.sh
./pear.sh

# -----------------------------------------------------------------------------
# 2. Filter merged reads and convert FASTQ to FASTA
# -----------------------------------------------------------------------------
# PEAR's assembled FASTQ files are the inputs. VSEARCH applies the recorded
# maximum expected-error threshold of 1 and relabels reads with the OTU prefix.
# The second occurrence of `assembled.fastq` in each generated command is changed
# to `fasta`, leaving the input filename intact and naming the FASTA output.
ls *assembled.fastq | awk '{print "vsearch -fastq_filter "$0" -fastaout "$0" -relabel OTU -fastq_maxee 1"}' | sed 's/assembled.fastq/fasta/2' > fastq.sh
chmod u+x fastq.sh
./fastq.sh

# -----------------------------------------------------------------------------
# 3. Select the ANML primer pair and inline barcodes interactively
# -----------------------------------------------------------------------------
mamba activate trimming_env

# Launch the bundled wrapper from the directory containing merged FASTA files.
# At its prompts, choose ANML, inline barcodes A,B,C, and pipeline mode UNIX.
# Accept the current directory for output so the following commands can find it.
# These are interactive selections, not command-line arguments. Primer sequences
# and expected insert lengths come from the configured primer table; the code
# comments explain the matching rules and the 95–105% insert-length window.
# The ANML output used below is named `2_merged_ANML.fas`.
"$TRIMFLOW_PY"

# -----------------------------------------------------------------------------
# 4. Dereplicate ANML reads and infer denoised sequence variants
# -----------------------------------------------------------------------------
# Work in the sequence directory containing `2_merged_ANML.fas`.
# Dereplication writes unique sequences and their counts (`-sizeout`). The
# original log then uses USEARCH UNOISE3; no extra UNOISE options are added here.
vsearch -fastx_uniques 2_merged_ANML.fas -sizeout -fastaout 3_ANML_derep.fa
"$USEARCH_BIN" -unoise3 3_ANML_derep.fa -zotus 4_zotu_ANML.fas

# Map the merged reads to the inferred variants. Preserve both the abundance
# table (`-otutabout`) and read-mapping output (`-mapout`) from the original run.
# The command records 12 threads and otherwise uses the original tool defaults.
"$USEARCH_BIN" -otutab 2_merged_ANML.fas \
  -zotus 4_zotu_ANML.fas \
  -otutabout 5_ANML_zotu_table.txt \
  -mapout 5_map_ANML.txt \
  -threads 12

# -----------------------------------------------------------------------------
# 5. Search the inferred ANML sequences against the recorded nt reference
# -----------------------------------------------------------------------------
mamba activate blast2tax

: "${NT_DATABASE:?Set NT_DATABASE to the February 2025 nt database prefix.}"

# Retain at most 10 target sequences per query, using 12 threads. The explicitly
# recorded tabular columns end with `staxids`, which supplies taxonomic IDs for
# the next step. Relative filenames replace the original server-specific paths.
blastn \
    -db "$NT_DATABASE" \
    -query 4_zotu_ANML.fas \
    -num_threads 12 \
    -max_target_seqs 10 \
    -out 6_ANML_blast.out \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"

# -----------------------------------------------------------------------------
# 6. Attach taxonomy to BLAST hits using the bundled version 1.4.2 utility
# -----------------------------------------------------------------------------
if [ ! -f "$BLAST2TAXONOMY_PY" ] || [ ! -x "$BLAST2TAXONOMY_PY" ]; then
    printf '%s\n' "Missing or non-executable taxonomy program: $BLAST2TAXONOMY_PY" >&2
    exit 1
fi

# Column 13 is the `staxids` field in the BLAST output above. The original
# invocation's `-c 13 -t 20 -f` settings are retained exactly. Consult the original
# utility's own documentation for its other options; they are not inferred here.
"$BLAST2TAXONOMY_PY" -i 6_ANML_blast.out -o 7_ANML.tax.out -c 13 -t 20 -f
