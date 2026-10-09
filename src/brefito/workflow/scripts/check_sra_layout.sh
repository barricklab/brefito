#!/usr/bin/env bash
#
# Check that an SRA run actually has the read layout that data.csv claims it has,
# BEFORE spending time dumping it.
#
# Usage: check_sra_layout.sh <run_accession> <paired|single> [error|warn] [num_spots]
#
# This reads only run metadata with vdb-dump (a few seconds, no read data is
# downloaded). For each of the first <num_spots> spots it counts how many reads are
# biological (not technical barcodes/adapters) AND have a nonzero length. A run whose
# spots have two such reads is paired; one is single-end.
#
# Spots are sampled rather than just checking the first one because runs exist whose
# first spots are unpaired even though the run as a whole is paired. This is the same
# property that makes 'fastq-dump --split-files' produce R1/R2 files with unequal
# numbers of reads, so the maximum over the sampled spots is what decides the layout.
#
# If the metadata cannot be read (no network, an unusual VDB configuration, ...) this
# prints a warning and exits 0: a probe failure must never block a download.

set -uo pipefail

run_accession="${1:?usage: check_sra_layout.sh <run_accession> <paired|single> [error|warn] [num_spots]}"
expected_layout="${2:?usage: check_sra_layout.sh <run_accession> <paired|single> [error|warn] [num_spots]}"
mismatch_mode="${3:-error}"
num_spots="${4:-100}"

layout_counts=$(vdb-dump -R "1-${num_spots}" -C READ_LEN,READ_TYPE "$run_accession" 2>/dev/null | awk '
    BEGIN { max_reads = 0; min_reads = -1; num_spots = 0 }
    /READ_LEN:/ {
        sub(/^[^:]*: */, "")
        num_lengths = split($0, read_length, /, */)
        have_lengths = 1
        next
    }
    /READ_TYPE:/ {
        if (!have_lengths) next
        sub(/^[^:]*: */, "")
        num_types = split($0, read_type, /, */)
        biological_reads = 0
        for (i = 1; i <= num_types; i++)
            if (read_type[i] ~ /BIOLOGICAL/ && read_length[i] + 0 > 0)
                biological_reads++
        if (biological_reads > max_reads) max_reads = biological_reads
        if (min_reads < 0 || biological_reads < min_reads) min_reads = biological_reads
        num_spots++
        have_lengths = 0
        next
    }
    END { print max_reads, min_reads, num_spots }
')

read -r max_reads min_reads spots_examined <<< "$layout_counts"

if [ -z "${spots_examined:-}" ] || [ "${spots_examined:-0}" -eq 0 ]; then
    echo "WARNING: could not read layout metadata for SRA run ${run_accession}." >&2
    echo "         Skipping the paired/single-end check and downloading anyway." >&2
    exit 0
fi

if [ "$max_reads" -ge 2 ]; then
    observed_layout="paired"
else
    observed_layout="single"
fi

echo "SRA run ${run_accession}: ${observed_layout}-end (${max_reads} biological read(s) per spot in the first ${spots_examined} spots)."

if [ "$observed_layout" = "paired" ] && [ "$min_reads" -lt 2 ]; then
    echo "         Note: some spots have only one read. Those reads cannot be paired and"
    echo "         will be written to sra-downloads/${run_accession}.unpaired.fastq.gz."
fi

if [ "$observed_layout" = "$expected_layout" ]; then
    exit 0
fi

if [ "$expected_layout" = "paired" ]; then
    problem="data.csv asks for paired-end reads, but SRA run ${run_accession} is single-end."
    suggestion="Use type 'illumina-SE' (or 'nanopore') for this accession instead of 'illumina-PE'."
else
    problem="data.csv asks for single-end reads, but SRA run ${run_accession} is paired-end."
    suggestion="Use type 'illumina-PE' for this accession. Reading paired data as single-end joins both mates into one read."
fi

if [ "$mismatch_mode" = "warn" ]; then
    echo "WARNING: ${problem}" >&2
    echo "         ${suggestion}" >&2
    echo "         Continuing anyway because SRA_IGNORE_LAYOUT_MISMATCH is set." >&2
    exit 0
fi

echo "ERROR: ${problem}" >&2
echo "       ${suggestion}" >&2
echo "       Pass --config SRA_IGNORE_LAYOUT_MISMATCH=1 to downgrade this to a warning." >&2
exit 1
