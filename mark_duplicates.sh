#!/bin/bash

# Define your environment name
ENV_NAME="gatk-env"  # ← change this if needed

# Define input/output directory
ALIGN_DIR="/your directory/alignment_output"

# Activate conda environment
source ~/anaconda3/etc/profile.d/conda.sh  # Adjust path if needed
conda activate "$ENV_NAME"

# Mark duplicates for Non Recurrent sample
gatk MarkDuplicates \
  -I "$ALIGN_DIR/Non_Recurrent.sorted.bam" \
  -O "$ALIGN_DIR/de_dup.Non_Recurrent.sorted.bam" \
  -M "$ALIGN_DIR/de_dup.Non_Recurrent.metrics.txt"

# Index the deduplicated BAM
samtools index "$ALIGN_DIR/de_dup.Non_Recurrent.sorted.bam"

# Mark duplicates for Recurrent Event sample
gatk MarkDuplicates \
  -I "$ALIGN_DIR/Recurrent_event.sorted.bam" \
  -O "$ALIGN_DIR/de_dup.Recurrent_event.sorted.bam" \
  -M "$ALIGN_DIR/de_dup.Recurrent_event.metrics.txt"

# Index the deduplicated BAM
samtools index "$ALIGN_DIR/de_dup.Recurrent_event.sorted.bam"

echo "✅ Duplicate marking and indexing complete!"

