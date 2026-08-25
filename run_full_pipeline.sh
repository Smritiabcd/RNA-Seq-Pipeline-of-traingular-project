#!/bin/bash
set -euo pipefail

# ==========================================================================
# Bee Virome Pipeline — Full Combined Workflow
# Raw FastQC -> fastp -> Trimmed FastQC -> MultiQC
#   -> Host removal (Apis cerana -> Apis mellifera -> Varroa)
#   -> Alignment to new virus -> stats/coverage/depth
# Loops over ALL paired samples in raw_fastq/
# ==========================================================================

# ---- CONFIG ----
PROJECT_DIR="/mnt/d/bee_immune_analysis"
RAW_DIR="$PROJECT_DIR/raw_fastq"
QC_DIR="$PROJECT_DIR/qc"
TRIM_DIR="$PROJECT_DIR/trimmed"
HOST_REMOVED_DIR="$PROJECT_DIR/host_removed"
ALIGN_DIR="$PROJECT_DIR/aligned/newvirus"
LOG_DIR="$PROJECT_DIR/logs"

APIS_CERANA="$PROJECT_DIR/references/apis_cerana_genome.fa"
APIS_MELLIFERA="$PROJECT_DIR/references/apis_mellifera.fa"
VARROA="$PROJECT_DIR/references/varroa_genome.fa"
NEW_VIRUS="$PROJECT_DIR/references/newvirus.fasta"

THREADS=8
FASTQC_THREADS=4

# ---- SETUP DIRECTORIES ----
mkdir -p "$QC_DIR/fastqc_raw" "$QC_DIR/fastqc_trimmed" "$QC_DIR/fastp_reports" \
         "$QC_DIR/multiqc_final" "$TRIM_DIR" "$HOST_REMOVED_DIR" "$ALIGN_DIR" "$LOG_DIR"

# ==========================================================================
# STEP 1: FastQC on RAW reads (all samples together)
# ==========================================================================
echo "[$(date)] STEP 1: FastQC (raw) — all samples"
cd "$RAW_DIR"
fastqc *.fq.gz -o "$QC_DIR/fastqc_raw" -t "$FASTQC_THREADS" \
    2> "$LOG_DIR/fastqc_raw.log"

# ==========================================================================
# STEP 2: Indexing host + virus reference genomes (once, before the loop)
# ==========================================================================
echo "[$(date)] STEP 2: Indexing reference genomes (if needed)"
[ -f "${APIS_CERANA}.bwt" ]    || bwa index "$APIS_CERANA"    2>> "$LOG_DIR/bwa_index.log"
[ -f "${APIS_MELLIFERA}.bwt" ] || bwa index "$APIS_MELLIFERA" 2>> "$LOG_DIR/bwa_index.log"
[ -f "${VARROA}.bwt" ]         || bwa index "$VARROA"         2>> "$LOG_DIR/bwa_index.log"
[ -f "${NEW_VIRUS}.bwt" ]      || bwa index "$NEW_VIRUS"      2>> "$LOG_DIR/bwa_index.log"

# ==========================================================================
# MAIN LOOP — per-sample: fastp -> trimmed FastQC -> host removal -> virus alignment
# ==========================================================================
cd "$RAW_DIR"
for r1 in *_1.fq.gz; do
    SAMPLE=$(basename "$r1" _1.fq.gz)
    R1="${SAMPLE}_1.fq.gz"
    R2="${SAMPLE}_2.fq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: mate file $R2 not found for $R1 — skipping $SAMPLE"
        continue
    fi

    echo ""
    echo "=========================================================="
    echo "  SAMPLE: $SAMPLE"
    echo "=========================================================="

    # ---- STEP 3: fastp ----
    echo "[$(date)] STEP 3: fastp — $SAMPLE"
    fastp \
        -i "$RAW_DIR/$R1" -I "$RAW_DIR/$R2" \
        -o "$TRIM_DIR/${SAMPLE}_1.fq.gz" \
        -O "$TRIM_DIR/${SAMPLE}_2.fq.gz" \
        --html "$QC_DIR/fastp_reports/${SAMPLE}.fastp.html" \
        --json "$QC_DIR/fastp_reports/${SAMPLE}.fastp.json" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --low_complexity_filter \
        --complexity_threshold 30 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --cut_tail \
        2> "$LOG_DIR/fastp_${SAMPLE}.log"

    # ---- STEP 4: FastQC on trimmed reads ----
    echo "[$(date)] STEP 4: FastQC (trimmed) — $SAMPLE"
    fastqc -t "$FASTQC_THREADS" -o "$QC_DIR/fastqc_trimmed" \
        "$TRIM_DIR/${SAMPLE}_1.fq.gz" \
        "$TRIM_DIR/${SAMPLE}_2.fq.gz" \
        2> "$LOG_DIR/fastqc_trimmed_${SAMPLE}.log"

    # ---- STEP 5a: Remove Apis cerana reads ----
    echo "[$(date)] STEP 5a: Host removal — Apis cerana — $SAMPLE"
    bwa mem -t "$THREADS" "$APIS_CERANA" \
        "$TRIM_DIR/${SAMPLE}_1.fq.gz" \
        "$TRIM_DIR/${SAMPLE}_2.fq.gz" \
        2> "$LOG_DIR/bwa_cerana_${SAMPLE}.log" \
      | samtools view -b -@ "$THREADS" \
      | samtools sort -n -@ "$THREADS" \
      | samtools fastq -@ "$THREADS" \
        -f 12 -F 256 \
        -1 "$HOST_REMOVED_DIR/${SAMPLE}_no_cerana_1.fq.gz" \
        -2 "$HOST_REMOVED_DIR/${SAMPLE}_no_cerana_2.fq.gz" \
        -0 /dev/null -s /dev/null \
        2> "$LOG_DIR/fastq_cerana_${SAMPLE}.log"

    echo "  Reads remaining after Apis cerana removal ($SAMPLE):"
    zcat "$HOST_REMOVED_DIR/${SAMPLE}_no_cerana_1.fq.gz" | awk 'NR%4==1' | wc -l

    # ---- STEP 5b: Remove Apis mellifera reads ----
    echo "[$(date)] STEP 5b: Host removal — Apis mellifera — $SAMPLE"
    bwa mem -t "$THREADS" "$APIS_MELLIFERA" \
        "$HOST_REMOVED_DIR/${SAMPLE}_no_cerana_1.fq.gz" \
        "$HOST_REMOVED_DIR/${SAMPLE}_no_cerana_2.fq.gz" \
        2> "$LOG_DIR/bwa_mellifera_${SAMPLE}.log" \
      | samtools view -b -@ "$THREADS" \
      | samtools sort -n -@ "$THREADS" \
      | samtools fastq -@ "$THREADS" \
        -f 12 -F 256 \
        -1 "$HOST_REMOVED_DIR/${SAMPLE}_no_mellifera_1.fq.gz" \
        -2 "$HOST_REMOVED_DIR/${SAMPLE}_no_mellifera_2.fq.gz" \
        -0 /dev/null -s /dev/null \
        2> "$LOG_DIR/fastq_mellifera_${SAMPLE}.log"

    echo "  Reads remaining after Apis mellifera removal ($SAMPLE):"
    zcat "$HOST_REMOVED_DIR/${SAMPLE}_no_mellifera_1.fq.gz" | awk 'NR%4==1' | wc -l

    # ---- STEP 5c: Remove Varroa destructor reads ----
    echo "[$(date)] STEP 5c: Host removal — Varroa destructor — $SAMPLE"
    bwa mem -t "$THREADS" "$VARROA" \
        "$HOST_REMOVED_DIR/${SAMPLE}_no_mellifera_1.fq.gz" \
        "$HOST_REMOVED_DIR/${SAMPLE}_no_mellifera_2.fq.gz" \
        2> "$LOG_DIR/bwa_varroa_${SAMPLE}.log" \
      | samtools view -b -@ "$THREADS" \
      | samtools sort -n -@ "$THREADS" \
      | samtools fastq -@ "$THREADS" \
        -f 12 -F 256 \
        -1 "$HOST_REMOVED_DIR/${SAMPLE}_clean_1.fq.gz" \
        -2 "$HOST_REMOVED_DIR/${SAMPLE}_clean_2.fq.gz" \
        -0 /dev/null -s /dev/null \
        2> "$LOG_DIR/fastq_varroa_${SAMPLE}.log"

    echo "  Reads remaining after Varroa removal — final clean reads ($SAMPLE):"
    zcat "$HOST_REMOVED_DIR/${SAMPLE}_clean_1.fq.gz" | awk 'NR%4==1' | wc -l

    # ---- STEP 6: Align clean reads to new virus ----
    echo "[$(date)] STEP 6: BWA-MEM -> new virus — $SAMPLE"
    bwa mem -t "$THREADS" "$NEW_VIRUS" \
        "$HOST_REMOVED_DIR/${SAMPLE}_clean_1.fq.gz" \
        "$HOST_REMOVED_DIR/${SAMPLE}_clean_2.fq.gz" \
        2> "$LOG_DIR/bwa_newvirus_${SAMPLE}.log" \
      | samtools view -b -F 4 -@ "$THREADS" \
      | samtools sort -@ "$THREADS" \
      -o "$ALIGN_DIR/${SAMPLE}_vs_newvirus.bam"

    samtools index "$ALIGN_DIR/${SAMPLE}_vs_newvirus.bam"

    echo "[$(date)] Alignment stats — $SAMPLE:"
    samtools flagstat "$ALIGN_DIR/${SAMPLE}_vs_newvirus.bam" \
        | tee "$ALIGN_DIR/${SAMPLE}_flagstat.txt"

    samtools coverage "$ALIGN_DIR/${SAMPLE}_vs_newvirus.bam" \
        | tee "$ALIGN_DIR/${SAMPLE}_coverage.txt"

    samtools depth -a "$ALIGN_DIR/${SAMPLE}_vs_newvirus.bam" \
        > "$ALIGN_DIR/${SAMPLE}_depth.txt"

    echo "  Zero-depth positions ($SAMPLE):"
    awk '$3==0' "$ALIGN_DIR/${SAMPLE}_depth.txt" | wc -l

    echo "  Significant virus hits (coverage>10%, depth>1) — $SAMPLE:"
    awk 'NR>1 && $6>10 && $7>1' "$ALIGN_DIR/${SAMPLE}_coverage.txt" \
        | sort -k7 -rn \
        | column -t

done

# ==========================================================================
# STEP 7: MultiQC — combined report across raw FastQC, trimmed FastQC, fastp
# ==========================================================================
echo ""
echo "[$(date)] STEP 7: MultiQC — combined report"
cd "$QC_DIR"
multiqc fastqc_raw fastqc_trimmed fastp_reports \
    -o multiqc_final \
    -n multiqc_final_report

echo ""
echo "=========================================================="
echo "PIPELINE COMPLETE"
echo "  Raw FastQC:       $QC_DIR/fastqc_raw"
echo "  fastp reports:    $QC_DIR/fastp_reports"
echo "  Trimmed FastQC:   $QC_DIR/fastqc_trimmed"
echo "  Trimmed reads:    $TRIM_DIR"
echo "  Host-removed:     $HOST_REMOVED_DIR"
echo "  Virus alignments: $ALIGN_DIR"
echo "  MultiQC report:   $QC_DIR/multiqc_final/multiqc_final_report.html"
echo "=========================================================="
