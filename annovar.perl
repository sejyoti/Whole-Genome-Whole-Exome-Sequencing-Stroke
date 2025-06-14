./annovar/table_annovar.pl results/sample.filtered.vcf annovar/humandb/ -buildver hg38 \
-out annovar_annotation/sample_annotated -remove -protocol refGene,clinvar_20220320,dbnsfp42a \
-operation g,f,f -nastring . -vcfinput
