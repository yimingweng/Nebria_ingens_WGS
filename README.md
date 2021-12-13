# Nebria_ingens_WGS
This work is to search selective sweep using three different methods. The first analysis is OmegaPlus v2.2.2 which detect the highest ω statist in a giving region. The second approach is PCAdapt which detects outlier SNPs by calculating the z-score of regression using genotype of focal SNP as response variable and principal components (PCs) representing the population structure as explanatory variable (Luu, Bazin, & Blum, 2017). Finally, the putative selective sweep will be scanned by calculating the composite likelihood, the μ statistic, from the signature of reduction of genetic diversity, shift of the site allele frequency, and the change of linkage disequilibrium using RAiSD (Alachiotis & Pavlidis, 2018). 

## OmegaPlus
1. The program OmegaPlus v2.2.2 can take vcf as input file, but it doesn't take locus with more than 2 alleles, so those multiallelic loci needs to be trimmed. This is actually good thing, because it has been shown that the multiallelic loci tend to show more genotype errors. The other operation is to remove the loci with minor allele less than 0.05 (maf=0.05). Additionally, the SNPs on the contigs smaller than 10kbp were excluded from the analyses as the SNPs on those small contigs might not be accepted by the programs or could potentially bias the results.
2. OmegaPlus is calculating the maximum ω statistic in a giving window of genome, could be chromosome or scaffold or contig. The main input parameter is the grid size, which is the number of *ω* statistic being calculate in this genome region. Since we have many contigs and they are different in size, it is better to fix the region size than grid size, so that the interval between two location to assess omega can be fixed throughout the genome.
```
mkdir omegaplus_out
IFS=$'\n'
for contig in $(cat contig_10k_list)
do
  name=$(echo ${contig} | cut -d $'\t' -f 1)
  size=$(echo ${contig} | cut -d $'\t' -f 2)
  grid=$(( ${size}/100 ))
  echo -e "working on contig ${name}"
  vcftools --vcf lcs_final_maf005.vcf --chr ${name} --out ${name} --recode
  bcftools view --max-alleles 2 --exclude-types indels ${name}.recode.vcf >> ${name}.recode_biallic.vcf
  ../OmegaPlus-F  -name ${name}  -input ${name}.recode_biallic.vcf -ld RSQUARE -grid ${grid} -minwin 1000 -maxwin 10000 -seed 12345 -threads 24
  rm ${name}.recode.vcf ${name}.log ${name}.recode_biallic.vcf
  mv *.${name} ./omegaplus_out/
done
```

