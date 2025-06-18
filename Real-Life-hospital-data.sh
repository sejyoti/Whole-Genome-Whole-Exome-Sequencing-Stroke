#!/bin/bash
set -euo pipefail

# Define paths
REF="/home/sejyoti/Downloads/Biocompin/Stroke_WES/hg38/hg38.fa"
READS_DIR="/home/sejyoti/Downloads/Biocompin/Stroke_WES/reads"
OUT_DIR="/home/sejyoti/Downloads/Biocompin/Stroke_WES/alignment_output"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Step 1: Check and index reference genome if not already indexed
if [ ! -f "${REF}.bwt" ]; then
    echo "[INFO] Indexing reference genome using BWA..."
    bwa index "$REF"
else
    echo "[INFO] Reference genome already indexed."
fi

# Suvendu (Recurrent)
echo "[INFO] Aligning Recurrent sample: Sumitra..."
bwa mem -t 8 -R "@RG\tID:RE2\tSM:RecurrentSumitra\tPL:ILLUMINA\tLB:Lib3" "$REF" \
    "$READS_DIR/Recurrent_21204200682_Suvendu_R1.fastp.fastq.gz" \
    "$READS_DIR/Recurrent_21204200682_Suvendu_R2.fastp.fastq.gz" \
    > "$OUT_DIR/Recurrent_Suvendu.sam"

echo "[INFO] Processing Suvendu Recurrent alignment..."
samtools view -@ 8 -bS "$OUT_DIR/Recurrent_Suvendu.sam" | \
    samtools sort -@ 8 -o "$OUT_DIR/Recurrent_Suvendu.sorted.bam"

samtools index "$OUT_DIR/Recurrent_Suvendu.sorted.bam"
rm "$OUT_DIR/Recurrent_Suvendu.sam"

# Rabi (Non-Recurrent)
echo "[INFO] Aligning Non-Recurrent sample: Rabi..."
bwa mem -t 8 -R "@RG\tID:NR2\tSM:NonRecurrentRabi\tPL:ILLUMINA\tLB:Lib4" "$REF" \
    "$READS_DIR/NR_30404200661_Rabi_R1.fastq.gz" \
    "$READS_DIR/NR_30404200661_RABI_R2.fastq.gz" \
    > "$OUT_DIR/Non_Recurrent_Rabi.sam"

echo "[INFO] Processing Rabi Non-Recurrent alignment..."
samtools view -@ 8 -bS "$OUT_DIR/Non_Recurrent_Rabi.sam" | \
    samtools sort -@ 8 -o "$OUT_DIR/Non_Recurrent_Rabi.sorted.bam"

samtools index "$OUT_DIR/Non_Recurrent_Rabi.sorted.bam"
rm "$OUT_DIR/Non_Recurrent_Rabi.sam"
