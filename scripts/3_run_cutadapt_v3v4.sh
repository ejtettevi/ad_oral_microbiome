#!/bin/bash
# run_cutadapt_illumina_16s.sh
# Process Illumina 16S V3-V4 sequencing data
# Focus on adapter removal and quality trimming rather than primer removal

set -euo pipefail

INPUT_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/fastp"
OUTPUT_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/cutadapt"
THREADS=8

# Illumina adapter sequences (more relevant than primers for Illumina data)
ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Optional: Partial primer sequences (in case some primer remnants remain)
PRIMER_FWD_PARTIAL="CCTACGGG"  # First 8bp of 341F
PRIMER_REV_PARTIAL="GACTACHV"  # First 8bp of 806R

mkdir -p "$OUTPUT_DIR"

echo "cutadapt version:"
cutadapt --version
echo ""

# Check a few input files first
echo "=== INPUT FILE CHECK ==="
echo "Checking first few input files..."

# Fix: Use a counter to limit to first 3 files
counter=0
for R1 in "$INPUT_DIR"/*_R1.fastq*; do
    if [[ $counter -ge 3 ]]; then
        break
    fi
    
    if [[ -s "$R1" ]]; then
        SAMPLE=$(basename "$R1" | sed 's/_R1\.fastq.*//')
        echo "Sample $SAMPLE: $(stat -f%z "$R1") bytes"
        echo "  First few sequences:"
        if [[ "$R1" == *.gz ]]; then
            zcat "$R1" | head -8
        else
            head -8 "$R1"
        fi
        echo "  ---"
        ((counter++))
    fi
done
echo ""

# Process all _R1.fastq and _R1.fastq.gz files
echo "=== PROCESSING SAMPLES ==="
processed_count=0
failed_count=0

for R1 in "$INPUT_DIR"/*_R1.fastq*; do
    # Only process if file exists and is not empty
    [[ -s "$R1" ]] || continue

    SAMPLE=$(basename "$R1" | sed 's/_R1\.fastq.*//')
    echo "Processing sample: $SAMPLE"
    
    # Find matching R2 (either .fastq or .fastq.gz)
    R2="$INPUT_DIR/${SAMPLE}_R2.fastq"
    if [[ ! -f "$R2" ]]; then
        R2="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"
    fi

    if [[ -f "$R2" && -s "$R2" ]]; then
        # Set output file extensions to match input (fastq or fastq.gz)
        if [[ "$R1" == *.gz ]]; then
            OUT_R1="$OUTPUT_DIR/${SAMPLE}_R1.trimmed.fastq.gz"
            OUT_R2="$OUTPUT_DIR/${SAMPLE}_R2.trimmed.fastq.gz"
        else
            OUT_R1="$OUTPUT_DIR/${SAMPLE}_R1.trimmed.fastq"
            OUT_R2="$OUTPUT_DIR/${SAMPLE}_R2.trimmed.fastq"
        fi

        # More appropriate parameters for Illumina 16S data
        cutadapt \
            -a "$ADAPTER_FWD" \
            -A "$ADAPTER_REV" \
            -g "$PRIMER_FWD_PARTIAL" \
            -G "$PRIMER_REV_PARTIAL" \
            --minimum-length 100 \
            --maximum-length 300 \
            --quality-cutoff=20,20 \
            --error-rate=0.1 \
            --overlap=3 \
            --pair-filter=any \
            -j $THREADS \
            -o "$OUT_R1" \
            -p "$OUT_R2" \
            "$R1" "$R2" > "$OUTPUT_DIR/${SAMPLE}_cutadapt.log" 2>&1

        # Check if output files were created and have content
        if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
            echo "  ✅ Success: R1=$(stat -f%z "$OUT_R1") bytes, R2=$(stat -f%z "$OUT_R2") bytes"
            ((processed_count++))
        else
            echo "  ❌ Failed: Empty output files"
            echo "    Check log: $OUTPUT_DIR/${SAMPLE}_cutadapt.log"
            ((failed_count++))
        fi
    else
        echo "  ❌ Warning: Paired file for $R1 not found or empty. Skipping."
        ((failed_count++))
    fi
done

echo ""
echo "=== PROCESSING SUMMARY ==="
echo "Successfully processed: $processed_count samples"
echo "Failed: $failed_count samples"
echo ""

# Count output files
non_empty=$(find "$OUTPUT_DIR" -name "*trimmed.fastq*" -size +0 2>/dev/null | wc -l)
empty=$(find "$OUTPUT_DIR" -name "*trimmed.fastq*" -size 0 2>/dev/null | wc -l)

echo "Non-empty trimmed files: $non_empty"
echo "Empty trimmed files: $empty"
echo ""

# Show some example log content for diagnosis
echo "=== EXAMPLE LOG CONTENT ==="
log_file=$(find "$OUTPUT_DIR" -name "*cutadapt.log" | head -1)
if [[ -f "$log_file" ]]; then
    echo "Content from $(basename "$log_file"):"
    head -20 "$log_file"
    echo ""
    echo "Summary from same log:"
    tail -10 "$log_file"
else
    echo "No log files found"
fi
echo ""

# Run MultiQC for summary
if [[ $processed_count -gt 0 ]]; then
    echo "Running MultiQC..."
    multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" --force --quiet
    echo "MultiQC report generated: $OUTPUT_DIR/multiqc_report.html"
else
    echo "Skipping MultiQC - no samples were successfully processed"
fi

echo ""
echo "=== NEXT STEPS ==="
if [[ $failed_count -gt 0 ]]; then
    echo "To diagnose failures, check individual log files:"
    echo "ls -la $OUTPUT_DIR/*.log | head -5"
    echo ""
    echo "Or examine a specific log:"
    echo "cat $OUTPUT_DIR/[SAMPLE_NAME]_cutadapt.log"
fi