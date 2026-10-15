#!/usr/bin/env python3

"""Build and run primer/barcode trimming commands with the paired shell engine.

This file coordinates the workflow; ``trim_fasta.sh`` performs sequence matching,
primer/barcode removal, and insert-length filtering. Run the wrapper from the
folder containing the input sequences. UNIX mode uses FASTA files and combines
barcode outputs by marker; APSCALE mode uses individual ``*_PE.fastq.gz`` files
and retains FASTQ sequence/quality records.

Configuration comes from ``config.yaml`` beside this file unless the
``TRIMFLOW_CONFIG`` environment variable names another file. At startup the
archive configuration reads the bundled ANML primer CSV beside this script.
Custom configurations without ``use_bundled_primer_database: true`` retain the
original Google Sheet download to ``primer_database_path`` in the launch working
directory. Primer pairs and expected insert lengths come from the selected CSV;
amplicons must match its ``amplicon_name`` values exactly, including case.

The output-directory prompt controls the trimmed sequence files and UNIX merged
files. ``run_trimfasta.sh`` and ``trimflow_log.txt`` are written in the launch
working directory; the shell engine also places its UNIX cache and name log
there. See the engine's comments for matching details and FASTA header format.

Source history:
- December 2024: interactive primer/barcode selection, command generation,
  naming adjustments, and merging of nonempty barcode files.
- 19 February 2026: explicit UNIX/APSCALE paths, output-directory selection,
  script-local configuration with an environment override, and local engine
  fallback. The source's processing logic is retained here with expanded comments.
"""

import glob
import os
import shlex
import subprocess
import csv
import pandas as pd
import yaml
import re
from datetime import datetime  # Timestamp entries in trimflow_log.txt.

# Locate the default configuration beside this wrapper, independently of the
# directory containing the reads. Input discovery below still uses the launch CWD.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# A caller may supply a different configuration without editing the script.
CONFIG_PATH = os.environ.get("TRIMFLOW_CONFIG", os.path.join(SCRIPT_DIR, "config.yaml"))

# =============================================================================
# Utility Functions
# =============================================================================

def load_config(config_path):
    """Parse the selected YAML file and return its contents.

    The configuration supplies primer-database and engine paths, inline-barcode
    mappings, and an optional Google Sheet ID for online mode. With
    ``use_bundled_primer_database: true``,
    a relative primer-database path is resolved beside this script at startup.
    Other relative paths retain their launch-working-directory interpretation.
    Loading errors terminate the wrapper before command generation.
    """
    try:
        with open(config_path, "r") as file:
            return yaml.safe_load(file)
    except Exception as e:
        print(f"Failed to load configuration: {e}")
        exit(1)

def fetch_database(config):
    """Download the configured Google Sheet and replace the local primer CSV.

    Called before ``load_database`` when the configuration does not enable the
    bundled primer database. pandas reads the sheet's CSV export directly;
    ``primer_database_path`` names the file written by ``to_csv`` and then read
    by the wrapper. A failed download or write terminates execution. The archive
    configuration skips this function and reads its bundled ANML table instead.
    """
    print("Fetching primer database from Google Sheets...")
    try:
        sheet_id = config["google_sheet_id"]
        csv_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv"
        df = pd.read_csv(csv_url)
        df.to_csv(config["primer_database_path"], index=False)
        print(f"Primer database saved to '{config['primer_database_path']}'.")
    except Exception as e:
        print(f"Failed to fetch the primer database: {e}")
        exit(1)

def load_database(config):
    """Group CSV rows into forward/reverse primer pairs by amplicon name.

    Required columns are ``amplicon_name``, ``type``, and ``sequence``. Amplicon
    names and primer sequences have surrounding whitespace removed; primer type
    is compared in lowercase. The CSV's ``sequence`` column is used for both
    primers, so reverse primers should be supplied in their synthesis orientation:
    the shell engine obtains the reverse complement used for matching.

    ``amplicon_length_wo_primer`` is the expected insert length in bases. Missing,
    empty, or nonnumeric values use 200. The first row encountered for an amplicon
    establishes its length; later rows of each primer type replace that primer's
    sequence. A missing forward or reverse primer is skipped later when commands
    are built. The function itself does not impose a particular marker list.
    """
    try:
        with open(config["primer_database_path"], mode="r") as file:
            reader = csv.DictReader(file)
            amplicons = {}
            for row in reader:
                amplicon = row["amplicon_name"].strip()
                primer_type = row["type"].strip().lower()
                sequence = row["sequence"].strip()

                # Spreadsheet exports may express integer lengths as "64.0".
                # Convert through float, then int; use the existing 200-bp default
                # for blank values or strings that cannot be parsed.
                try:
                    length = int(float(row.get("amplicon_length_wo_primer", 200) or 200))
                except ValueError:
                    length = 200

                # Length is stored when this marker first appears in the CSV.
                if amplicon not in amplicons:
                    amplicons[amplicon] = {"Forward": None, "Reverse": None, "Length": length}

                if primer_type == "forward":
                    amplicons[amplicon]["Forward"] = sequence
                elif primer_type == "reverse":
                    amplicons[amplicon]["Reverse"] = sequence
            return amplicons
    except Exception as e:
        print(f"Failed to load the primer database: {e}")
        exit(1)

def log_changes(log_file, description):
    """Append one timestamped message to the supplied log-file path.

    Calls below pass ``trimflow_log.txt`` as a relative path, so this text log
    lives in the launch working directory. Entries describe skipped markers and
    successful UNIX file concatenations; they are not per-read processing logs.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"[{timestamp}] {description}\n")


def resolve_trim_fasta_path(config):
    """Find the configured shell engine, falling back to the adjacent copy.

    The configured path is used if it exists. Otherwise, look beside this Python
    file for ``trim_fasta.sh`` and print the fallback selected. The wrapper exits
    if neither location exists. The engine is invoked directly by the generated
    shell commands, so its executable permission and shell dependencies apply.
    """
    configured = config.get("trim_fasta_path", "")
    if configured and os.path.exists(configured):
        return configured

    fallback = os.path.join(SCRIPT_DIR, "trim_fasta.sh")
    if os.path.exists(fallback):
        print(
            f"Configured trim script path '{configured}' not found. "
            f"Using local fallback '{fallback}'."
        )
        return fallback

    print(
        f"Failed to locate trim_fasta.sh. Checked configured path '{configured}' "
        f"and fallback '{fallback}'."
    )
    exit(1)

# =============================================================================
# Core Workflow
# =============================================================================

def run_unix_pipeline(trim_fasta_path, config, amplicons, selected_amplicons, selected_barcodes, output_dir):
    """Run the FASTA route and concatenate nonempty barcode files per marker.

    Write one shell command for each selected amplicon/barcode combination. Each
    engine invocation scans ``*.fasta`` in the current working directory; the
    wrapper does not enumerate individual FASTA input files. Trimmed files go to
    ``output_dir`` as ``1_<lowercase-marker>[_<barcode>].fas``.

    The command script is always ``run_trimfasta.sh`` in the launch directory.
    After it runs, consider only expected output files that exist and have a
    nonzero byte size. More than one such file is concatenated into
    ``2_merged_<marker>.fas`` in barcode-selection order. If exactly one is
    nonempty, it remains a ``1_...`` file without creating a ``2_merged_...`` copy.
    This step concatenates records; it does not deduplicate sequences.
    """
    print("\nGenerating trimming commands...")
    # A configured label (for example A) resolves to its DNA barcode. An
    # unrecognized token is passed through as a custom sequence; X is handled
    # separately below by omitting the engine's -i option entirely.
    full_barcode_map = config["inline_barcodes"]
    barcode_map = {b.strip().upper(): full_barcode_map.get(b.strip().upper(), b.strip()) for b in selected_barcodes}
    os.makedirs(output_dir, exist_ok=True)

    # Recreate the command file in the launch CWD on each run. The sequence
    # outputs can be elsewhere, as selected by the caller.
    script_path = os.path.join(os.getcwd(), "run_trimfasta.sh")
    primer_output_files = {}

    with open(script_path, "w") as script:
        script.write("#!/bin/bash\n\n")
        for amplicon in selected_amplicons:
            amplicon = amplicon.strip()
            if amplicon not in amplicons:
                log_changes("trimflow_log.txt", f"Amplicon {amplicon} not found. Skipping...")
                continue

            details = amplicons[amplicon]
            forward, reverse, length = details["Forward"], details["Reverse"], details["Length"]

            if not forward or not reverse:
                log_changes("trimflow_log.txt", f"Missing primers for {amplicon}. Skipping...")
                continue

            # Filter the insert after removing the barcode and both primers.
            # Bounds are inclusive in the engine and truncated with int here:
            # an expected length of 64 gives 60..67; 188 gives 178..197.
            min_length, max_length = int(length * 0.95), int(length * 1.05)
            # The engine turns each uppercase ambiguity symbol into regex '.'.
            # This enables a wildcard at that position, rather than an IUPAC base
            # set or general mismatch tolerance. Literal primer bases stay exact.
            wobble_flag = "-w" if re.search(r"[RYSWKMBDHVN]", forward + reverse) else ""

            if amplicon not in primer_output_files:
                primer_output_files[amplicon] = []

            for barcode, sequence in barcode_map.items():
                # Lowercasing applies to these intermediate filenames only.
                # CSV lookup and the merged filename retain the selected key.
                primer_name = amplicon.lower()
                if barcode == "X":
                    output_file = os.path.join(output_dir, f"1_{primer_name}.fas")
                    inline_flag = ""
                else:
                    output_file = os.path.join(output_dir, f"1_{primer_name}_{barcode}.fas")
                    inline_flag = f"-i {shlex.quote(sequence)}"

                # Remember the expected files for the later concatenation step.
                primer_output_files[amplicon].append(output_file)

                # Quote each shell argument so paths and primer strings are
                # passed as single arguments. UNIX is explicitly selected.
                command = (
                    f"{shlex.quote(trim_fasta_path)} -p UNIX -f {shlex.quote(forward)} "
                    f"-r {shlex.quote(reverse)} -o {shlex.quote(output_file)} "
                    f"{inline_flag} -l {min_length} -m {max_length} {wobble_flag}\n"
                )
                script.write(command)

    print(f"Trimming commands written to '{script_path}'. Running the script...")
    os.chmod(script_path, 0o755)
    # Execution is synchronous. The existing wrapper does not request
    # check=True or inspect the shell return code; engine output appears on
    # the terminal. The UNIX route then assesses the expected output files.
    subprocess.run(script_path, shell=True)

    for primer, files in primer_output_files.items():
        # "Nonempty" here means a nonzero byte count, not FASTA validation.
        non_empty_files = [f for f in files if os.path.exists(f) and os.path.getsize(f) > 0]

        if len(non_empty_files) > 1:
            merged_file = os.path.join(output_dir, f"2_merged_{primer}.fas")
            print(f"Merging files for primer {primer} into {merged_file}...")
            with open(merged_file, "w") as outfile:
                for file in non_empty_files:
                    with open(file, "r") as infile:
                        outfile.writelines(infile.readlines())
            log_changes("trimflow_log.txt", f"Merged files for primer {primer} into {merged_file}.")
        else:
            if len(non_empty_files) == 1:
                print(f"Only one file contains data for primer {primer}. No merging required.")
            else:
                print(f"No valid data found for primer {primer}. Skipping merging.")


def run_apscale_pipeline(trim_fasta_path, config, amplicons, selected_amplicons, selected_barcodes, output_dir):
    """Run the per-file FASTQ route without a final marker-level merge.

    Discover ``*_PE.fastq.gz`` files in the launch working directory and generate
    one engine command for every input/amplicon/barcode combination. The paired
    engine trims sequence and quality strings together. Each output is named
    ``<sample_prefix>_<amplicon>_<barcode>_PE_trimmed.fastq.gz`` inside
    ``output_dir``; the sample prefix is the input basename without the exact
    ``_PE.fastq.gz`` suffix. Selecting X retains ``_X_`` in the output filename
    but supplies no inline barcode to the matching engine.

    As in UNIX mode, ``run_trimfasta.sh`` is written in the launch directory and
    executed synchronously. This route does not concatenate the resulting FASTQ
    files or perform paired-end assembly.
    """
    print("\nGenerating APSCALE trimming commands...")
    # Resolve configured labels to sequences; retain custom tokens as DNA.
    # X is handled separately below and suppresses the inline-barcode flag.
    full_barcode_map = config["inline_barcodes"]
    barcode_map = {b.strip().upper(): full_barcode_map.get(b.strip().upper(), b.strip()) for b in selected_barcodes}
    os.makedirs(output_dir, exist_ok=True)

    # Discovery is deliberately limited to this suffix in the current folder.
    # Sorting fixes the order of the generated per-input commands.
    input_files = sorted(glob.glob("*_PE.fastq.gz"))
    if not input_files:
        print("No input files ending with '_PE.fastq.gz' found in the current directory.")
        return

    script_path = os.path.join(os.getcwd(), "run_trimfasta.sh")
    with open(script_path, "w") as script:
        script.write("#!/bin/bash\n\n")

        for amplicon in selected_amplicons:
            amplicon = amplicon.strip()
            if amplicon not in amplicons:
                log_changes("trimflow_log.txt", f"Amplicon {amplicon} not found. Skipping...")
                continue

            details = amplicons[amplicon]
            forward, reverse, length = details["Forward"], details["Reverse"], details["Length"]

            if not forward or not reverse:
                log_changes("trimflow_log.txt", f"Missing primers for {amplicon}. Skipping...")
                continue

            # Filter the insert after removing the barcode and both primers.
            # Bounds are inclusive in the engine and truncated with int here:
            # an expected length of 64 gives 60..67; 188 gives 178..197.
            min_length, max_length = int(length * 0.95), int(length * 1.05)
            # The engine turns each uppercase ambiguity symbol into regex '.'.
            # This enables a wildcard at that position, rather than an IUPAC base
            # set or general mismatch tolerance. Literal primer bases stay exact.
            wobble_flag = "-w" if re.search(r"[RYSWKMBDHVN]", forward + reverse) else ""

            for input_file in input_files:
                if not input_file.endswith("_PE.fastq.gz"):
                    continue

                # Keep the complete sample name preceding the required suffix.
                sample_prefix = os.path.basename(input_file)[: -len("_PE.fastq.gz")]
                for barcode, sequence in barcode_map.items():
                    if barcode == "X":
                        barcode_suffix = "X"
                        inline_flag = ""
                    else:
                        barcode_suffix = barcode
                        inline_flag = f"-i {shlex.quote(sequence)}"

                    output_file = os.path.join(
                        output_dir, f"{sample_prefix}_{amplicon}_{barcode_suffix}_PE_trimmed.fastq.gz"
                    )
                    command = (
                        f"{shlex.quote(trim_fasta_path)} -p APSCALE -s {shlex.quote(input_file)} "
                        f"-f {shlex.quote(forward)} -r {shlex.quote(reverse)} "
                        f"-o {shlex.quote(output_file)} {inline_flag} -l {min_length} -m {max_length} {wobble_flag}\n"
                    )
                    script.write(command)

    print(f"APSCALE commands written to '{script_path}'. Running the script...")
    os.chmod(script_path, 0o755)
    # Execution is synchronous. The existing wrapper does not request
    # check=True or inspect the shell return code; engine output appears on
    # the terminal. The UNIX route then assesses the expected output files.
    subprocess.run(script_path, shell=True)

# =============================================================================
# Main Execution
# =============================================================================

# Importing this module defines functions only. Running it directly selects
# the primer table, requests selections, and executes the generated commands.
if __name__ == "__main__":
    config = load_config(CONFIG_PATH)
    trim_fasta_path = resolve_trim_fasta_path(config)
    if config.get("use_bundled_primer_database", False):
        # Keep the publication table independent of the directory holding reads.
        # os.path.join preserves an explicitly configured absolute path.
        config["primer_database_path"] = os.path.join(SCRIPT_DIR, config["primer_database_path"])
    else:
        fetch_database(config)
    amplicons = load_database(config)

    print("\nAvailable Primers:")
    for amplicon, details in amplicons.items():
        print(f"Amplicon: {amplicon}")
        print(f"  Forward Primer: {details['Forward']}")
        print(f"  Reverse Primer: {details['Reverse']}")
        print(f"  Expected Length: {details['Length']} bp\n")

    # Marker names are case-sensitive CSV keys. For this study select ANML;
    # whitespace is stripped when commands are generated.
    selected_amplicons = input("Enter the amplicons to process (comma-separated): ").split(",")
    print("\nAvailable inline barcodes: A, B, C, D, E, F, X (no barcode)")
    selected_barcodes = input("Enter the inline barcodes or custom sequences (comma-separated): ").split(",")
    # Only the literal APSCALE selection takes the FASTQ route. Every other
    # response, including an empty response, follows the UNIX FASTA route.
    pipeline_mode = input("Select pipeline mode (UNIX or APSCALE): ").strip().upper()
    # This setting chooses sequence-output placement; input discovery and
    # wrapper logs/command script remain relative to the launch directory.
    default_output_dir = os.getcwd()
    output_dir = input(
        f"Enter output directory (press Enter for default: {default_output_dir}): "
    ).strip() or default_output_dir

    if pipeline_mode == "APSCALE":
        run_apscale_pipeline(trim_fasta_path, config, amplicons, selected_amplicons, selected_barcodes, output_dir)
    else:
        run_unix_pipeline(trim_fasta_path, config, amplicons, selected_amplicons, selected_barcodes, output_dir)
