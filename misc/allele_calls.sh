for i in `seq 1 22`; do awk 'BEGIN{OFS="\t"} NR>1{print $1,($2-1),$2, $22}' chr${i}.calls.txt > greatapes_mrca_allele_chr${i}.bed; done
parallel 'bcftools annotate -a greatapes_mrca_allele_sorted.bed.gz -h header_line.txt -c CHROM,FROM,TO,ANC_ALLELE {} > {= s/.vcf.gz/_anc.vcf/ =}' ::: ../processed/*vcf.gz
