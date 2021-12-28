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

# Choosing the number K of Principal Components
# Scree plot:
x <- pcadapt(input = variants, K = 10, min.maf = 0.05 , LD.clumping = list(size = 500, thr = 0.1))
summary(x$loadings)

x11()
plot(x, option = "screeplot")

# Score plot
Pop <- read.table("E:/YiMing/pcadapt/ingens_lcs_final_list", header = T, sep = "\t")
poplist.sp <- Pop$sp
poplist.names <- Pop$pop
x11()
plot(x, option = "scores", i = 1, j = 2, pop = poplist.sp)
x11()
plot(x, option = "scores", i = 9, j = 10, pop = poplist.names)
# Use up to K10


# Manhattan Plot
x11()
plot(x , option = "manhattan")

# QQ plot
x11()
plot(x, option = "qqplot")

# Histograms of the test statistic and of the p-values
x11()
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# histogram of the test statistic D
x11()
plot(x, option = "stat.distribution")


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

