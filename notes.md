## OmegaPlus
1. Run omegaplus, for each contig individually, as those contigs are different in size and I will have to fix the grid size for them. The maximum omega statistics will be calculated every 100 bps within the range of 1,000 to 10,000 bps. 
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
2. The omega statistics are calculated contig by contig, so now put them together to find the highest omega statistics.

3. here is the example of the individual output:
```
head OmegaPlus_Report.328
//1
125.00  0.000000
222.77  0.000000
320.54  0.000000
418.31  0.000000
516.07  0.000000
613.84  1.612176
711.61  1.612176
809.38  1.612176
```
4. And there are 1,448 contigs (I removed the contigs smaller than 10Kbp)
```
ls OmegaPlus_Report* | wc -l
1448
```
5. Merge those files
```
IFS=$'\n'
for file in $(ls OmegaPlus_Report*)
do
  echo "working on ${file}..."
  contig=$(echo ${file} | cut -d "." -f 2)
  for line in $(cat ${file})
  do
    echo -e "${contig}\t${line}" | sed 's/\.[0-9]*//1' >> omegaPlus_out
  done
done
cat omegaPlus_out | grep -v "/" | sort -n -k1,1 -k2,2 >> omegaPlus_out_final
```

2. After merging the omegaplus outputs of those contigs, use R to find the largest 0.999 quantile of omega statistics. The R script is written in "**omegaplus_outliers.r**" in scripts folder.
```
# R environment
setwd("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/omegaplus")
omegaplus.out <- read.table("omegaPlus_out_final")
# define outlier and write them into a bed file
upper_bound <- quantile(omegaplus.out$V3, 0.999)
outliers <- omegaplus.out[which(omegaplus.out$V3 > upper_bound),]
# 1447 SNPs being defined as outlier
write.csv(outliers, "omegaplus_999quantile.csv", quote=F, row.names = F)
bed <- cbind(outliers$V1, outliers$V2, outliers$V2)
write.table(bed, "omegaplus_999quantile.bed", quote=F, col.names = F, row.names = F, sep="\t")
```

3. This outputs 1,477 sites, so use these 1,477 outliers to find the candidate genes
```
working dir: sean@Fuji:/media/Jade/YMW/gff/WGS/
# convert the bed files from DOS line endings to Unix line endings
# I ran this because I encounted this error from bedtools: ERROR: Received illegal bin number 4294967295 from getBin call
dos2unix omegaplus_999quantile.bed
# find the intersection of omegaplus outlier and the genes (gff.bed)
bedtools intersect -a gff.bed -b omegaplus_999quantile.bed -wa | sour -n -k1,1 -k2,2 >> omegaplus_gff.bed
```

4. Get the gene list from the bed file
    - generate bed file with omegaplus coordinates:
```
working dir: sean@Fuji:/media/Jade/YMW/gff/WGS/omegaplus_out/
dos2unix omegaplus_999quantile.bed
bedtools intersect -a omegaplus_999quantile.bed -b ../gff.bed -wa >> omegaplus_gff_tmp
```

 - get the omega statistic for the regions
```
IFS=$'\n'
for site in $(cat omegaplus_gff_tmp)
do
  echo "working on ${site}..."
  contig=$(echo ${site} | cut -d $'\t' -f 1)
  pos=$(echo ${site} | cut -d $'\t' -f 2)
  cat /media/Jade/YMW/OmegaPlus_v2.2.2_Linux/ingens_scan/omegaPlus_out_final | grep -Pw "${contig}"$'\t'"${pos}" >> omegaplus_gff_coordinate
done
bedtools intersect -a ../gff.bed -b omegaplus_999quantile.bed -wa >> omegaplus_gff.bed.tmp
paste omegaplus_gff.bed.tmp omegaplus_gff_coordinate | cut -d $'\t' -f 1,2,3,6 >> omegaplus_gene_position
rm omegaplus_gff.bed.tmp omegaplus_gff_coordinate
```
  - get the mean of omega statistic for the same sites
```
IFS=$'\n'
for gene in $(cat omegaplus_gene_position)
do
  contig=$(echo ${gene} | cut -d $'\t' -f 1)
  start=$(echo ${gene} | cut -d $'\t' -f 2)
  stop=$(echo ${gene} | cut -d $'\t' -f 3)
  omega=$(echo ${gene} | cut -d $'\t' -f 4)
  mean_omega=$(cat omegaplus_gene_position | grep -Pw "${contig}"$'\t'"${start}"$'\t'"${stop}" | cut -d $'\t' -f 4 | tr '\n' ' ' |awk '{s=0; for (i=1;i<=NF;i++)s+=$i; print s/NF;}')
  echo -e "${contig}\t${start}\t${stop}\t${mean_omega}" >> omegaplus_gene_position_mean
done
```


  - grap the gene from the gff file
```
# working dir: sean@Fuji:/media/Jade/YMW/gff/WGS/omegaplus_out/
IFS=$'\n'
for site in $(cat omegaplus_gene_position_mean | sort -u -n -k1,1 -k2,2)
do
  echo "working on ${site}..."
  contig=$(echo ${site} | cut -d $'\t' -f 1)
  start=$(echo ${site} | cut -d $'\t' -f 2)
  stop=$(echo ${site} | cut -d $'\t' -f 3)
  omega=$(echo ${site} | cut -d $'\t' -f 4)
  gene=$(cat ../../Nebria_riversi_contig_annotation.gff3 | grep -Pw "gene" | grep -Pw "^${contig}" | grep -Pw "${start}" | grep -Pw "${stop}" | cut -d $'\t' -f 9 | cut -d ";" -f 6 | cut -d " " -f 1 | cut -d "=" -f 2)
  ID=$(cat ../../Nebria_riversi_contig_annotation.gff3 | grep -Pw "gene" | grep -Pw "^${contig}" | grep -Pw "${start}" | grep -Pw "${stop}" | cut -d $'\t' -f 9 | cut -d ";" -f 1 | cut -d "=" -f 2)
  echo -e "${contig}\t${start}\t${omega}\t${ID}\t${gene}" >> omegaplus_gene_list
done
cat omegaplus_gene_list | grep "XP" | wc -l 
# in 121 genes, 98 genes have gene annotation
# sort the final table
cat omegaplus_gene_list | sort -n -k3,3 >> omegaplus_gene_final_list
rm *tmp* *position*
```



