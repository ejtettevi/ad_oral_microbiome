#!/bin/bash
# run_bowtie2_decontam.sh
set -euo pipefail
set -x

INPUT_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/cutadapt"
OUT_DIR="/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/qc/hostfree"
BOWTIE2_INDEX="/Users/edkwami/bowtie2_db/hg_39"
THREADS=4
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "=== INITIAL CHECKS ==="
# Check if Bowtie2 index exists
if [[ ! -f "${BOWTIE2_INDEX}.1.bt2" ]]; then
    echo "ERROR: Bowtie2 index not found at ${BOWTIE2_INDEX}"
    echo "Please check if the index files exist:"
    ls -la "${BOWTIE2_INDEX}"* 2>/dev/null || echo "No index files found"
    exit 1
fi

echo "Bowtie2 index found: ${BOWTIE2_INDEX}"
echo "Available disk space:"
df -h "$OUT_DIR"

echo "=== SOFTWARE VERSIONS ==="
echo "bowtie2 version:"
bowtie2 --version | head -3
echo "samtools version:"
samtools --version | head -3

echo "=== INPUT FILES ==="
echo "Looking for files in: $INPUT_DIR"
shopt -s nullglob
R1_FILES=("$INPUT_DIR"/*_R1.trimmed.fastq*)
echo "Found ${#R1_FILES[@]} R1 files to process"

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 files found matching pattern *_R1.trimmed.fastq*"
    exit 1
fi

# Process files one by one
PROCESSED=0
FAILED=0

for R1 in "${R1_FILES[@]}"; do
    echo "=========================================="
    echo "Processing file: $R1"
    
    if [[ ! -s "$R1" ]]; then
        echo "WARNING: File $R1 is empty or does not exist, skipping"
        ((FAILED++))
        continue
    fi

    SAMPLE=$(basename "$R1" | sed 's/_R1\.trimmed\.fastq.*//')
    echo "Sample name: $SAMPLE"

    # Find corresponding R2 file
    R2=""
    for ext in ".fastq.gz" ".fastq"; do
        potential_R2="$INPUT_DIR/${SAMPLE}_R2.trimmed${ext}"
        if [[ -f "$potential_R2" && -s "$potential_R2" ]]; then
            R2="$potential_R2"
            break
        fi
    done

    if [[ -z "$R2" ]]; then
        echo "WARNING: Paired file for $R1 not found or empty. Skipping."
        ((FAILED++))
        continue
    fi

    echo "Found pair: $R1 and $R2"
    
    # Check file sizes
    echo "File sizes:"
    ls -lh "$R1" "$R2"

    # Output file names
    if [[ "$R1" == *.gz ]]; then
        OUT_R1="$OUT_DIR/${SAMPLE}_R1.hostfree.fastq.gz"
        OUT_R2="$OUT_DIR/${SAMPLE}_R2.hostfree.fastq.gz"
    else
        OUT_R1="$OUT_DIR/${SAMPLE}_R1.hostfree.fastq"
        OUT_R2="$OUT_DIR/${SAMPLE}_R2.hostfree.fastq"
    fi

    echo "Output files will be: $OUT_R1 and $OUT_R2"

    # Skip if output already exists
    if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
        echo "Output files already exist, skipping $SAMPLE"
        ((PROCESSED++))
        continue
    fi

    # Bowtie2 alignment with better error handling
    echo "Running bowtie2 alignment for $SAMPLE..."
    BOWTIE_LOG="$LOG_DIR/${SAMPLE}_bowtie2.log"
    SAM_FILE="$OUT_DIR/${SAMPLE}_contam.sam"
    
    if ! bowtie2 -x "$BOWTIE2_INDEX" \
                  -1 "$R1" \
                  -2 "$R2" \
                  --very-sensitive \
                  -p $THREADS \
                  -S "$SAM_FILE" \
                  2> "$BOWTIE_LOG"; then
        echo "ERROR: Bowtie2 alignment failed for $SAMPLE"
        echo "Check log file: $BOWTIE_LOG"
        cat "$BOWTIE_LOG"
        ((FAILED++))
        continue
    fi

    # Check if SAM file was created and is not empty
    if [[ ! -s "$SAM_FILE" ]]; then
        echo "ERROR: SAM file was not created or is empty for $SAMPLE"
        ((FAILED++))
        continue
    fi

    echo "Bowtie2 alignment completed. SAM file size:"
    ls -lh "$SAM_FILE"

    # Extract unmapped reads (host-free)
    echo "Extracting unmapped reads for $SAMPLE..."
    UNMAPPED_BAM="$OUT_DIR/${SAMPLE}_unmapped.bam"
    UNMAPPED_SORTED_BAM="$OUT_DIR/${SAMPLE}_unmapped_sorted.bam"
    
    if ! samtools view -b -f 12 -F 256 "$SAM_FILE" > "$UNMAPPED_BAM"; then
        echo "ERROR: Failed to extract unmapped reads for $SAMPLE"
        ((FAILED++))
        rm -f "$SAM_FILE"
        continue
    fi

    if ! samtools sort -n "$UNMAPPED_BAM" -o "$UNMAPPED_SORTED_BAM"; then
        echo "ERROR: Failed to sort unmapped reads for $SAMPLE"
        ((FAILED++))
        rm -f "$SAM_FILE" "$UNMAPPED_BAM"
        continue
    fi

    echo "Unmapped reads extracted and sorted. BAM file size:"
    ls -lh "$UNMAPPED_SORTED_BAM"

    # Output host-free FASTQ
    echo "Converting to FASTQ for $SAMPLE..."
    if [[ "$OUT_R1" == *.gz ]]; then
        if ! samtools fastq -1 >(gzip > "$OUT_R1") -2 >(gzip > "$OUT_R2") "$UNMAPPED_SORTED_BAM"; then
            echo "ERROR: Failed to convert to gzipped FASTQ for $SAMPLE"
            ((FAILED++))
            rm -f "$SAM_FILE" "$UNMAPPED_BAM" "$UNMAPPED_SORTED_BAM"
            continue
        fi
    else
        if ! samtools fastq -1 "$OUT_R1" -2 "$OUT_R2" "$UNMAPPED_SORTED_BAM"; then
            echo "ERROR: Failed to convert to FASTQ for $SAMPLE"
            ((FAILED++))
            rm -f "$SAM_FILE" "$UNMAPPED_BAM" "$UNMAPPED_SORTED_BAM"
            continue
        fi
    fi

    # Clean up intermediate files
    rm -f "$SAM_FILE" "$UNMAPPED_BAM" "$UNMAPPED_SORTED_BAM"

    # Verify output files
    if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
        echo "✅ Successfully decontaminated $SAMPLE"
        echo "Output files:"
        ls -lh "$OUT_R1" "$OUT_R2"
        ((PROCESSED++))
    else
        echo "ERROR: Output files are empty for $SAMPLE"
        ((FAILED++))
    fi
    
    echo "---"
done

echo "=========================================="
echo "=== FINAL SUMMARY ==="
echo "Total files processed successfully: $PROCESSED"
echo "Total files failed: $FAILED"
echo "Output directory contents:"
ls -lh "$OUT_DIR"/*.fastq* 2>/dev/null || echo "No output files found"

if [[ $FAILED -gt 0 ]]; then
    echo "WARNING: Some files failed to process. Check log files in $LOG_DIR"
    ls -la "$LOG_DIR"
fi
