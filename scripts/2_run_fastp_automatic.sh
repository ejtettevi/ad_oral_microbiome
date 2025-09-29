#!/bin/bash
# run_fastp_automatic.sh
# Automated fastp QC for all paired-end FASTQ files in RAW_DIR

set -euo pipefail

RAW_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/raw"
OUT_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/fastp"
THREADS=8

mkdir -p "$OUT_DIR"

echo "fastp version:"
fastp --version

for R1 in "$RAW_DIR"/*_R1.fastq*; do
    SAMPLE=$(basename "$R1" | sed 's/_R1\.fastq.*//')
    R2="$RAW_DIR/${SAMPLE}_R2.fastq"
    # If gzipped, adjust R2 filename
    if [[ ! -f "$R2" ]]; then
        R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"
    fi

    if [[ -f "$R2" ]]; then
        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$OUT_DIR/${SAMPLE}_R1.fastq" \
            -O "$OUT_DIR/${SAMPLE}_R2.fastq" \
            --detect_adapter_for_pe \
            --html "$OUT_DIR/${SAMPLE}_fastp.html" \
            --json "$OUT_DIR/${SAMPLE}_fastp.json" \
            --thread $THREADS

        echo "Processed $SAMPLE"
    else
        echo "Warning: Paired file for $R1 not found. Skipping."
    fi
done

# Run MultiQC for summary
multiqc "$OUT_DIR" -o "$OUT_DIR"
