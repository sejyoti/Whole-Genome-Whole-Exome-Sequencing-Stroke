#!/bin/bash

echo "STEP 4: Base quality recalibration"

# Define paths
ref="/your directory/hg38/hg38.fa"
known_sites="/your directory/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/your directory/alignment_output"
data="/your directory/alignment_output"

# ▶️ Non-Recurrent Sample

# 1. Build recalibration model
gatk BaseRecalibrator \
  -I ${aligned_reads}/de_dup.Non_Recurrent.sorted.bam \
  -R ${ref} \
  --known-sites ${known_sites} \
  -O ${data}/Non_Recurrent.recal_data.table

# 2. Apply the model to adjust base quality scores
gatk ApplyBQSR \
  -I ${aligned_reads}/de_dup.Non_Recurrent.sorted.bam \
  -R ${ref} \
  --bqsr-recal-file ${data}/Non_Recurrent.recal_data.table \
  -O ${aligned_reads}/de_dup_bqsr.Non_Recurrent.sorted.bam


# ▶️ Recurrent Event Sample

# 1. Build recalibration model
gatk BaseRecalibrator \
  -I ${aligned_reads}/de_dup.Recurrent_event.sorted.bam \
  -R ${ref} \
  --known-sites ${known_sites} \
  -O ${data}/Recurrent_event.recal_data.table

# 2. Apply the model to adjust base quality scores
gatk ApplyBQSR \
  -I ${aligned_reads}/de_dup.Recurrent_event.sorted.bam \
  -R ${ref} \
  --bqsr-recal-file ${data}/Recurrent_event.recal_data.table \
  -O ${aligned_reads}/de_dup_bqsr.Recurrent_event.sorted.bam

echo "✅ BQSR completed for both Non-Recurrent and Recurrent Event samples!"
