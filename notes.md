## OmegaPlus
[OmegaPlus](https://cme.h-its.org/exelixis/web/software/omegaplus/index.html) is a very user-friendly tool to detect the selective sweep. It finds the elevated omega statistic (the ratio of LD (linkage disequilibrium) across site to LD from the same site, where the site doesn't need to be the variation site). For more detail, please read [Kim and Nielsen 2004](https://academic.oup.com/genetics/article/167/3/1513/6050636) and [Alachiotis et al, 2012](https://academic.oup.com/bioinformatics/article/28/17/2274/245799), or this [manual](https://cme.h-its.org/exelixis/resource/download/software/OmegaPlus_Manual.pdf).

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
2. The omega statistics are calculated contig by contig, so now put them together to find the sites with highest omega statistic.

  - here is the example of the individual output:
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
  - And there are 1,448 contigs with this kind of output (Originally there are 2,137 contig, but I removed the contigs smaller than 10Kbp)
```
ls OmegaPlus_Report* | wc -l
1448
```
3. Merge the output of omegaplus from each contig, as add the contig name to each site
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

4. After merging the omegaplus outputs of those contigs, use R to find the largest 0.999 quantile of omega statistics. The R script is written in "**omegaplus_outliers.r**" in the [script folder](https://github.com/yimingweng/Nebria_ingens_WGS/tree/main/scripts). There are another twos methods I considered to use to define outliers but I decided not use them because they gave me similar results to the quantile method, and everything here is about how we want to set up the threshold, if no simulation result can be used to calculate the p-value with reasonable cutoff.
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

5. The quantile 0.999 gave an output with top 1,477 sites. Let's use these 1,477 outliers to find the candidate genes
```
working dir: sean@Fuji:/media/Jade/YMW/gff/WGS/
# convert the bed files from DOS line endings to Unix line endings
# I ran this because I encountered this error from bedtools: ERROR: Received illegal bin number 4294967295 from getBin call
dos2unix omegaplus_999quantile.bed
# find the intersection of omegaplus outlier and the genes (gff.bed)
bedtools intersect -a gff.bed -b omegaplus_999quantile.bed -wa | sour -n -k1,1 -k2,2 >> omegaplus_gff.bed
```

4. Because I wanted to have omega statistic of those genes, so I found the mean omega of sites being identified from the same gene to represent the omega statistic of the gene. With this information, I will be able to search gene function from the gene with highest mean omega statistic.
  - First, get the gene list from the bed file by generating bed file with omegaplus coordinates:
```
working dir: sean@Fuji:/media/Jade/YMW/gff/WGS/omegaplus_out/
dos2unix omegaplus_999quantile.bed
bedtools intersect -a omegaplus_999quantile.bed -b ../gff.bed -wa >> omegaplus_gff_tmp
```

  - get the omega statistic for the sites
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
  - check the final product: **omegaplus_gene_final_list**. Note that some (98 in 121) of genes do not have protein/function annotation. 
```
cat omegaplus_gene_final_list  | sort -nr -k3,3 | head
FORMAT: CONTIG  START_POS OMEGA_STAT GENE_NAME GENE_ID
773     72683   627.243 Nriv.00g078820  XP_975582.3
79      7684    407.265 Nriv.00g009670
1944    84176   346.908 Nriv.00g174690  XP_018578822.1
1335    15429   315.833 Nriv.00g135090  XP_018332723.1
1137    176     232.358 Nriv.00g114980  XP_003702090.1
989     25749   231.529 Nriv.00g100680  XP_019873439.1
1514    35749   219.906 Nriv.00g147950  XP_018577451.1
755     207314  199.232 Nriv.00g076550  XP_018579891.1
755     183745  199.232 Nriv.00g076490  XP_021942915.1
862     592     183.663 Nriv.00g084070  XP_022910755.1
```

5. Use this table to search gene function through internet. The results are manually edited in the excel file called "**omegaplus_gene_function.xlsx**".

6. This work will need to cite these papers:
    - N. Alachiotis, A. Stamatakis, and P. Pavlidis. OmegaPlus: A Scalable Tool for Rapid Detection of Selective Sweeps in Whole-Genome Datasets. Bioinformatics, 2012.

    - N. Alachiotis, P. Pavlidis, and A. Stamatakis. Exploiting Multi-grain Parallelism for  Efficient Selective Sweep Detection. Algorithms and Architectures for Parallel Processing, 56-68, 2012.

## RAiSD
[RAiSD](https://github.com/alachins/raisd) is another tools for Genome Scan developed more recently by [Alachiotis & Pavlidis, 2018](https://www.nature.com/articles/s42003-018-0085-8). This tools considers a composite likelihood called μ statistic three different signals of selective sweep, including the flattening site allele frequency spectrum (afs), the decrease of nucleotide diversity, and also the omega statistic. Because RAiSD takes afs information, it makes no sense to perform this analysis on the whole dataset fulfilling the genetic structure, so I use this tool to detect selective sweep population by population. And because RAiSD was performed on individual population, it will generate gene list for each population. The genes being detected in most populations will be considered a possible genes that is under positive selection. Additionally I will use the union of genes from all the populations to find the intersection with the other two lists from Omegaplus and PCAdapt.

1. Run RAiSD for each population, here I need lists of sample from each population, and use this list to extract the sample from the all-sample vcf file. The population vcf files will then be used to run RAiSD.
  - An example of the lists:
```
ls /media/Jade/YMW/fst/samples_of_populations | sed 's/\t/\n/g'
Army.txt
Conness.txt
Crabtree.txt
Donohue.txt
HungryPacker.txt
Italy.txt
Lamarck.txt
Lyell.txt
Millys.txt
Monarch.txt
Pear.txt
Piute.txt
Recess.txt
Ritter.txt
Ruby.txt
SamMack.txt
Selden.txt
SForester.txt
Taboose.txt
Treasure.txt
Wright.txt

head /media/Jade/YMW/fst/samples_of_populations/Army.txt
SDS06-306A
SDS06-306B
SDS06-322B
SDS06-322C
DS06-322D
SDS06-332A
SDS06-332B
SDS06-332C
SDS06-332D
SDS06-340A
```
```
# working dir: sean@Fuji:/media/Data1/Yiming/software/RAiSD/raisd-master
# select populations with sufficient sample size to run RAiSD by vcftools
mkdir ingens
cd ingens
# working dir: sean@Fuji:/media/Data1/Yiming/software/RAiSD/raisd-master/ingens
for pop in $(ls /media/Jade/YMW/fst/samples_of_populations/*)
do
  name=$(echo "$pop" | cut -d "/" -f 7 |cut -d "." -f 1)
  vcftools --gzvcf /media/Jade/YMW/imputation/lcs_final_maf005.vcf.gz --keep $pop --recode --out ./${name}
  rm *log
  mkdir ${name}_out
  /media/Data1/Yiming/software/RAiSD/raisd-master/bin/release/RAiSD -n ${name} -I ${name}.recode.vcf -f -c 3 -R -O -m 0.05 -P -D -A 0.995
  mv *${name}* ./${name}_out
done

# -n output prefix
# -I input vcf file (can be gzipped)
# -f overwight the existing output
# -m threshold for cutting the minor allele loci
# -P Generates plots in pdf
# -D Generates a site report
# -A Provides a probability value to be used for the quantile function in R and generates a Manhattan plot for the final mu-statistic score using Rscript (activates -s, -t, and -R).

# for unknown reason the Manhattan plot wasn't generated with this command

## for each population, collect the mu statistics from the output files
# working dir: sean@Fuji:/media/Data1/Yiming/software/RAiSD/raisd-master/ingens
for pop in $(ls ./)
do
  name=$(echo "${pop}" | cut -d "_" -f 1)
  echo "working on ${name}..."
  for contig in $(ls ${pop}/RAiSD_Report.${name}* | cut -d "/" -f 2 | cut -d "." -f 3 | sort -n)
  do
    awk -F '\t' -v OFS='\t' '{ $(NF+1) = '$contig'; print }' ./${pop}/RAiSD_Report.${name}.${contig} >> ./${name}_all_mu_stat
  done
done
```

2. Now all the populations have their own μ statistic outputs (I put this data in my own dropbox [here](https://www.dropbox.com/home/Data/Nebria_ingens/RAiSD_out)). Let's use R script called [**raisd_outliers.r**](https://github.com/yimingweng/Nebria_ingens_WGS/blob/main/scripts/raisd_outliers.r) to find the SNPs with highest μ statistic (again, I used 0.999 quantile) in all variation sites of the population. 
```
## R environment
library(openintro)

setwd("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/RAiSD/RAiSD_all_mu_statistics")
pop <- list.files(path = "C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/RAiSD/RAiSD_all_mu_statistics")

q=0.999
for (i in 1:length(pop)){
  population <- pop[i]
  print(paste("working on",population, " ..."))
  tiff(file=paste(population, ".tiff", sep=""))
  raisd.out <- read.table(population)
  raisd.out$V9 <- c(1:length(raisd.out$V7))
  upper_bound <- quantile(raisd.out$V7, q)
  raisd.out$Colour[raisd.out$V7>=upper_bound]="red"
  raisd.out$Colour[raisd.out$V7<=upper_bound]="black"
  plot(raisd.out$V9, raisd.out$V7, pch=20, cex=0.2, col=raisd.out$Colour)
  dev.off()
}

pdf("all_distribution.pdf")
for (i in 1:length(pop)){
  population <- pop[i]
  print(paste("working on",population, "..."))
  raisd.out <- read.table(population)
  upper_bound <- quantile(raisd.out$V7, q)
  densityPlot(raisd.out$V7)
  abline(v=upper_bound, col="red")
  outliers <- raisd.out[which(raisd.out$V7 > upper_bound),]
  bed <- cbind(outliers$V8, outliers$V1, outliers$V1)
  write.table(bed, file=paste(population, "999quantile.bed", sep=""), quote=F, col.names = F, row.names = F, sep="\t")
}
dev.off()
```


3. This R script generates outputs of top 0.001 SNPs with highest μ statistic. Now let's use bedtools to find the gene location for those SNPS.

```
# use bedtools to find the selective genes for each population
# working dir: sean@Fuji:/media/Data1/Yiming/software/RAiSD/raisd-master/ingens/intersection

mkdir gene_list
for file in $(ls ./*stat999quantile.bed)
do
  name=$(echo ${file} | cut -d "/" -f 2 | cut -d "_" -f 1)
  echo "working on ${name}..."
  bedtools intersect -a gff.bed -b ${name}_raisd.bed -wa | sort | uniq  >> ${name}_raisd_gff.bed
  IFS=$'\n'
  for site in $(cat ${name}_raisd_gff.bed)
  do
    contig=$(echo ${site} | cut -d $'\t' -f 1)
    start=$(echo ${site} | cut -d $'\t' -f 2)
    stop=$(echo ${site} | cut -d $'\t' -f 3)
    gene=$(cat /media/Jade/YMW/gff/Nebria_riversi_contig_annotation.gff3 | grep "gene"  | grep -Pw "^${contig}" | grep "${start}" | grep "${stop}" | grep -o "ID=.*" | cut -d ";" -f 1 | cut -d "=" -f 2) 
    echo "${gene}" >> ./gene_list/${name}_raisd_gene_list
  done
done

# working dir: sean@Fuji:/media/Data1/Yiming/software/RAiSD/raisd-master/ingens/intersection/gene_list

IFS=$'\n'
for gene in $(cat * | sort | uniq)
do
  count=$(cat *_gene_list | grep "${gene}" | wc -l)
  echo -e "${gene}\t${count}" >> all_genes
done
cat all_genes | sort -n -k3,3 >> all_raisd_genes_sorted
cat all_raisd_genes_sorted | wc -l 
# 949 genes are found in at least one population 
```

4. In the union 949 genes, there are 60 genes being detected at least in 4 populations. These genes will be pulled out and check their function.
```
for gene in $(cat all_raisd_genes_sorted | grep "XP")
  do 
  count=$(echo $gene | cut -d $'\t' -f 3)
  if [[ $count -ge 4 ]]
  then 
    echo $gene >> raisd_final_genes
  fi
done
cat  raisd_final_genes | wc -l 
# 60 genes
```
5. Use this table to search gene function through internet. The results are manually edited in the excel file called "raisd_gene_function.xlsx".

## Pcadapt
The third methods is Pcadapt which is quite different from OmegaPlus and RAiSD. It detects the SNPs that constructs abnormal genetic structure of fiven individual then the rest. This work was done following this awesome [vignettes](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html) by the authors of the original Pcadapt [paper](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12592). Note that I have simplified some of the steps. 

1. This program is based on R, so I will need to convert the input vcf file to the format that can be imported to R. Bed file of plink is one of it. So,convert the vcf file to bed file. Note that the contig name could make me in trouble, it has limitation on its format. See the script for more details.
```
working dir: sean@Fuji:/media/Jade/YMW/imputation/
zcat lcs_final.vcf.gz | grep "#" >> tmp.vcf
zcat lcs_final.vcf.gz | grep -v "#" | sed "s/^/Contig_/g" >> tmp.vcf
bgzip tmp.vcf
tabix tmp.vcf.gz
bcftools view -H tmp.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > lcs_final.txt
vcftools --gzvcf tmp.vcf.gz --plink --chrom-map lcs_final.txt --out tmp
../../../Data1/Yiming/software/plink --file tmp --allow-extra-chr --make-bed --out lcs_final 
```

2. When the input file in in bed file (map file is not required), it is ready to be imported to R and run Pcadapt. The converted bed file is now in my [dropbox](https://www.dropbox.com/preview/Data/Nebria_ingens/Nebria_ingens_GEA_data/lcs_final.bed?context=browse&role=personal) 
```
## R environment

#install.packages("pcadapt")
#install.packages("vcfR")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("qvalue")
library(qvalue)
library(pcadapt)

setwd("E:/YiMing/pcadapt/")
# convert the vcf file to bed file
variants <- read.pcadapt("E:/YiMing/pcadapt/lcs_final.bed", type = "bed", type.out = "bed")

# Choosing the number K of Principal Components (here we I use 10, see result of sNMF)
# Conduct pcadapt analysis, with minor allele threashold to be 0.05 (374*2*0.05= only loci with > 38 alt alleles will be considered)
# the LD thinning parameters are following the default setting (500 window size, r suqare=0.1)
x <- pcadapt(input = variants, K = 10, min.maf = 0.05 , LD.clumping = list(size = 500, thr = 0.1))
# see how many sites are removed from doing pca 
summary(x$loadings) # 1,642,161 sites are removed by maf and LD thining

# Scree plot to see the contribution of each PC:
x11()
plot(x, option = "screeplot")
```

2. Check how many PCs to be included to calculate the Mahalanobis distance. Here I took 10, as the published data has shown that this species complex roughly has good clusters from K=5 to K=10.
### **Scree plot here**

3. I can also view the population by providing the population informatin for those individuals.
```
# Score plot
Pop <- read.table("E:/YiMing/pcadapt/ingens_lcs_final_list", header = T, sep = "\t")
poplist.sp <- Pop$sp
poplist.names <- Pop$pop
x11()
plot(x, option = "scores", i = 1, j = 2, pop = poplist.sp)
x11()
plot(x, option = "scores", i = 9, j = 10, pop = poplist.names)
# Use up to K10
```
![](@attachment/Clipboard_2021-05-27-15-06-24.png)
![](@attachment/Clipboard_2021-05-27-15-07-57.png)



4. Check the result with Manhattan Plot
```
# Manhattan Plot
x11()
plot(x , option = "manhattan")
```
![](@attachment/Clipboard_2021-05-27-15-14-14.png)


5. Detect outliers with different threshold. I expected to see more outliers than the previous two methods so the main goal of this analysis to to provide another filter to get the final intersection of candidate genes (OmegaPLus+RAiSD+Pcadapt), so unlike thre previous two methods, I will NOT find the most extreme genes in small number and study their fuctions.
```
### detect outliers with different methods
# import and sort the loci position
pos <- read.table("E:/YiMing/pcadapt/lcs_final.bim", header=F, sep="\t")
p <- data.frame(pos=pos$V2,pvalue=x$pvalues)

# write p-value result into disk first, note that no cutoff is applied yet
p$pos <- gsub('Contig_', '', p$pos)
p.final.table <- data.frame(do.call('rbind', strsplit(as.character(p$pos),':',fixed=TRUE)), p$pvalue)
colnames(p.final.table) <- c("contig","pos", "pvalue")
p.final.table$contig <- as.numeric(p.final.table$contig)
p.final.table$pos <- as.numeric(p.final.table$pos)
p.final.table <- p.final.table[order(p.final.table$contig, p.final.table$pos),]

# remove NA (loci didn't pass the 0.05 maf filter)
p.final.table <- na.omit(p.final.table) 
write.table(p.final.table, "E:/YiMing/pcadapt/pvalue_outliers.txt", row.names=FALSE, sep="\t", quote = FALSE)

# Choosing 0.1 a cutoff for outlier detection
##### Q value
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
q.p <- data.frame(pos=pos$V2,pvalue=x$pvalues,qval=qval)
q.outliers <- q.p$pos[which(q.p$qval < alpha)]
length(q.outliers)

# create final table by adding the contig and position information
q.out <- data.frame(pos=q.outliers)
q.out.table <- merge(q.out, q.p, by.x = "pos")
q.out.table$pos <- gsub('Contig_', '', q.out.table$pos)
q.final.table <- data.frame(do.call('rbind', strsplit(as.character(q.out.table$pos),':',fixed=TRUE)), q.out.table$pvalue, q.out.table$qval)
colnames(q.final.table) <- c("contig","pos", "pvalue", "qvalue")
q.final.table$contig <- as.numeric(q.final.table$contig)
q.final.table$pos <- as.numeric(q.final.table$pos)
q.final.table <- q.final.table[order(q.final.table$contig, q.final.table$pos),]
write.table(q.final.table, "E:/YiMing/pcadapt/qvalue_outliers.txt", row.names=FALSE, sep="\t", quote = FALSE)

##### Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
padj.p <- data.frame(pos=pos$V2,pvalue=x$pvalues,padj=padj)
padj.outliers <- padj.p$pos[which(padj.p$padj < alpha)]
length(padj.outliers)
padj.out <- data.frame(pos=padj.outliers)
padj.out.table <- merge(padj.out, padj.p, by.x = "pos")
padj.out.table$pos <- gsub('Contig_', '', padj.out.table$pos)
padj.final.table <- data.frame(do.call('rbind', strsplit(as.character(padj.out.table$pos),':',fixed=TRUE)), padj.out.table$pvalue, padj.out.table$padj)
colnames(padj.final.table) <- c("contig","pos", "pvalue", "BH-p")
padj.final.table$contig <- as.numeric(padj.final.table$contig)
padj.final.table$pos <- as.numeric(padj.final.table$pos)
padj.final.table <- padj.final.table[order(padj.final.table$contig, padj.final.table$pos),]
write.table(padj.final.table, "E:/YiMing/pcadapt/BH_outliers.txt", row.names=FALSE, sep="\t", quote = FALSE)


##### Bonferroni correction
bonadj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
bonadj.p <- data.frame(pos=pos$V2,pvalue=x$pvalues,bonadj=bonadj)
bonadj.outliers <- bonadj.p$pos[which(bonadj.p$bonadj < alpha)]
length(bonadj.outliers)
bonadj.out <- data.frame(pos=bonadj.outliers)
bonadj.out.table <- merge(bonadj.out, bonadj.p, by.x = "pos")
bonadj.out.table$pos <- gsub('Contig_', '', bonadj.out.table$pos)
bonadj.final.table <- data.frame(do.call('rbind', strsplit(as.character(bonadj.out.table$pos),':',fixed=TRUE)), bonadj.out.table$pvalue, bonadj.out.table$bonadj)
colnames(bonadj.final.table) <- c("contig","pos", "pvalue", "adjust")
bonadj.final.table$contig <- as.numeric(bonadj.final.table$contig)
bonadj.final.table$contig <- as.numeric(bonadj.final.table$contig)
bonadj.final.table <- bonadj.final.table[order(bonadj.final.table$contig, bonadj.final.table$pos),]
write.table(bonadj.final.table, "E:/YiMing/pcadapt/bonferroni_outliers.txt", row.names=FALSE, sep="\t", quote = FALSE)
```

6. Finally, I used the most strict method (Bonferroni) at alpha=0.05 to define the outliers. Even with this strict threshold, I still got 14,532 SNPs passing this threshold, which hits 2,564 genes (see next section: Intersection of Three).