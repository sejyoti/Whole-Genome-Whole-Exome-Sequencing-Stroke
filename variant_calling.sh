#!/bin/bash

echo "STEP 6: Call Variants - GATK HaplotypeCaller"

# Define directories and files
ref="/your directory/hg38/hg38.fa"
aligned_reads="/your directory/alignment_output"
results="/your directory/alignment_output"  # output in same dir

# 🧪 Call variants for Non_Recurrent sample
gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/de_dup_bqsr.Non_Recurrent.sorted.bam \
  -O ${results}/raw_variants.Non_Recurrent.vcf

# 🧪 Call variants for Recurrent_event sample
gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/de_dup_bqsr.Recurrent_event.sorted.bam \
  -O ${results}/raw_variants.Recurrent_event.vcf

echo "✅ Variant calling completed. Raw VCFs generated!"
