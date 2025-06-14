#!/bin/bash

echo "STEP 5: Collect Alignment & Insert Size Metrics"

# Define paths
ref="/your directory/hg38/hg38.fa"
aligned_reads="/your directory/alignment_output"

# ▶️ Non_Recurrent Sample

gatk CollectAlignmentSummaryMetrics \
  R=${ref} \
  I=${aligned_reads}/de_dup_bqsr.Non_Recurrent.sorted.bam \
  O=${aligned_reads}/alignment_metrics.Non_Recurrent.txt

gatk CollectInsertSizeMetrics \
  INPUT=${aligned_reads}/de_dup_bqsr.Non_Recurrent.sorted.bam \
  OUTPUT=${aligned_reads}/insert_size_metrics.Non_Recurrent.txt \
  HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.Non_Recurrent.pdf


# ▶️ Recurrent_event Sample

gatk CollectAlignmentSummaryMetrics \
  R=${ref} \
  I=${aligned_reads}/de_dup_bqsr.Recurrent_event.sorted.bam \
  O=${aligned_reads}/alignment_metrics.Recurrent_event.txt

gatk CollectInsertSizeMetrics \
  INPUT=${aligned_reads}/de_dup_bqsr.Recurrent_event.sorted.bam \
  OUTPUT=${aligned_reads}/insert_size_metrics.Recurrent_event.txt \
  HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.Recurrent_event.pdf

echo "✅ Metrics collected for both samples!"
