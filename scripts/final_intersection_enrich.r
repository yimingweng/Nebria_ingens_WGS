if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(clusterProfiler)

setwd("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/")

term2gene <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/enrichment/term2gene", header=FALSE, sep="\t")
colnames(term2gene) <- c("GO", "gene")


WGS_intersect.genes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/data/WGS_final_enrich_test", header=FALSE, sep="\t")
WGS_intersect.genes <- as.vector(WGS_intersect.genes$V1)

# set BH adjusted p-value and q-value to be 0.05
WGS_intersect.005 <-enricher(WGS_intersect.genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1, qvalueCutoff = 0.05, TERM2GENE=term2gene, TERM2NAME=NA)
WGS_intersect.005
result.WGS_intersect.005 <- as.data.frame(WGS_intersect.005)
result.WGS_intersect.005
write.table(result.WGS_intersect.005, "WGS_intersect_enriched005_go.txt", sep="\t", quote=F, row.name=F)


# omegaplus+pcadapt
omepca_intersect.genes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/data/omegaplus_pcadapt_enrich_test", header=FALSE, sep="\t")
omepca_intersect.genes <- as.vector(omepca_intersect.genes$V1)

# set BH adjusted p-value and q-value to be 0.2
# this results in 36 Go terms
omepca_intersect.005 <-enricher(omepca_intersect.genes, pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 1, qvalueCutoff = 0.1, TERM2GENE=term2gene, TERM2NAME=NA)
omepca_intersect.005
result.omepca_intersect.005 <- as.data.frame(omepca_intersect.005)
result.omepca_intersect.005
write.table(result.omepca_intersect.005, "omepca_intersect_enriched005_go.txt", sep="\t", quote=F, row.name=F)


# omegaplus
omegaplus_intersect.genes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/data/omegaplus_enrich_test", header=FALSE, sep="\t")
omegaplus_intersect.genes <- as.vector(omegaplus_intersect.genes$V1)

# set BH adjusted p-value and q-value to be 0.2
# this results in 36 Go terms
omegaplus_intersect.005 <-enricher(omegaplus_intersect.genes, pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 1, qvalueCutoff = 0.1, TERM2GENE=term2gene, TERM2NAME=NA)
omegaplus_intersect.005
result.omegaplus_intersect.005 <- as.data.frame(omegaplus_intersect.005)
result.omegaplus_intersect.005
write.table(result.omegaplus_intersect.005, "omegaplus_intersect_enriched005_go.txt", sep="\t", quote=F, row.name=F)




# omegaplus
omegaplus_intersect.genes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/data/omegaplus_k7_log_resist_enrich_test", header=FALSE, sep="\t")
omegaplus_intersect.genes <- as.vector(omegaplus_intersect.genes$V1)

# set BH adjusted p-value and q-value to be 0.2
# this results in 36 Go terms
omegaplus_intersect.005 <-enricher(omegaplus_intersect.genes, pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 1, qvalueCutoff = 0.1, TERM2GENE=term2gene, TERM2NAME=NA)
omegaplus_intersect.005
result.omegaplus_intersect.005 <- as.data.frame(omegaplus_intersect.005)
result.omegaplus_intersect.005
write.table(result.omegaplus_intersect.005, "omegaplus_enriched005_go.txt", sep="\t", quote=F, row.name=F)

