#!/bin/bash

echo "STEP 6: Call Variants - GATK HaplotypeCaller"

# Define directories and files
ref="/home/sejyoti/Downloads/Biocompin/Stroke_WES/hg38/hg38.fa"
aligned_reads="/home/sejyoti/Downloads/Biocompin/Stroke_WES/alignment_output"
results="/home/sejyoti/Downloads/Biocompin/Stroke_WES/alignment_output"  # output in same dir

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
