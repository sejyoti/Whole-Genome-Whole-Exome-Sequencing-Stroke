# WGS/WES Variant Calling and Annotation Pipeline

This repository contains a customizable and modular pipeline for processing Whole Genome Sequencing (WGS) or Whole Exome Sequencing (WES) data using Unix shell scripting. The pipeline includes quality control, alignment, variant calling, annotation, and comprehensive report generation using MultiQC and ANNOVAR.
Analysis done on the basis of patient data.
---

## 📁 Directory Structure

project/
│
├── raw_data/ # Input FASTQ files
├── scripts/ # Shell scripts for each pipeline step
├── reference/ # Reference genome, index files, and annotation databases
├── results/ # Output files (BAM, VCF, QC reports, etc.)
├── multiqc_reports/ # Aggregated QC reports
└── annovar_annotation/ # Annotated variants


---

## 🔧 Requirements

- Unix/Linux system
- Bash (Shell scripting)
- BWA
- Samtools
- Picard
- GATK (v4+)
- FastQC
- MultiQC
- ANNOVAR
- bcftools
- Python (for MultiQC and auxiliary scripts)

---

## 🔄 Pipeline Steps

1. **Quality Control**
   - Tool: `FastQC`, `MultiQC`
   - Command:
     ```bash
     fastqc raw_data/*.fastq.gz -o results/fastqc/
     multiqc results/fastqc/ -o multiqc_reports/
     ```

2. **Reference Indexing**
   - Tool: `BWA`, `Samtools`
   - Command:
     ```bash
     bwa index reference/genome.fa
     samtools faidx reference/genome.fa
     ```

3. **Alignment**
   - Tool: `BWA-MEM`
   - Command:
     ```bash
     bwa mem -t 8 reference/genome.fa sample_R1.fastq.gz sample_R2.fastq.gz > results/sample.sam
     ```

4. **Conversion, Sorting, and Marking Duplicates**
   - Tool: `Samtools`, `Picard`
   - Commands:
     ```bash
     samtools view -Sb results/sample.sam > results/sample.bam
     samtools sort results/sample.bam -o results/sample.sorted.bam
     picard MarkDuplicates I=results/sample.sorted.bam O=results/sample.de_dup.bam M=results/sample.metrics.txt
     ```

5. **Base Recalibration & Variant Calling**
   - Tool: `GATK`
   - Command:
     ```bash
     gatk HaplotypeCaller -R reference/genome.fa -I results/sample.de_dup.bam -O results/sample.g.vcf.gz -ERC GVCF
     ```

6. **Joint Genotyping (Optional for Cohorts)**
   - Tool: `GATK GenotypeGVCFs`
   - Command:
     ```bash
     gatk GenotypeGVCFs -R reference/genome.fa -V results/sample.g.vcf.gz -O results/sample.vcf
     ```

7. **Variant Filtering**
   - Tool: `GATK VariantFiltration` / `bcftools`
   - Command:
     ```bash
     bcftools filter -s LOWQUAL -e '%QUAL<30 || DP<10' results/sample.vcf -o results/sample.filtered.vcf
     ```

8. **Variant Annotation**
   - Tool: `ANNOVAR`
   - Command:
     ```bash
     ./annovar/table_annovar.pl results/sample.filtered.vcf annovar/humandb/ -buildver hg38 \
     -out annovar_annotation/sample_annotated -remove -protocol refGene,clinvar_20220320,dbnsfp42a \
     -operation g,f,f -nastring . -vcfinput
     ```

9. **MultiQC Report Aggregation**
   - Tool: `MultiQC`
   - Command:
     ```bash
     multiqc results/ -o multiqc_reports/
     ```

---

## 📌 Notes

- Update file paths and filenames in the scripts according to your directory structure.
- Ensure reference genome and ANNOVAR databases are pre-downloaded and indexed.
- Scripts are modular and can be adapted for both WGS and WES depending on coverage and target regions.

---

## 📜 License

MIT License

---

## 👤 Author

Developed by Sejyoti Chakraborty, Lead Analyst at BiocompIn.  
Contact:sejyotichakraborty@gmail.com

---

