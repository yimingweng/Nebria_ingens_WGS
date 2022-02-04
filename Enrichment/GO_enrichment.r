if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO", force = TRUE)
library(topGO)

# read the go annotation database from Trinotate
gene2go <- readMappings("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/Nriv_trinotate_BP_GOs", sep = "\t", IDsep = ",")
gene2go[lapply(gene2go,length)==0] <-  NA

# import the final gene list from the three way intersection
geneSel <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/all_way_intersection.txt", header = F)
geneSel <- geneSel$V1
allGenes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/all_gene", header = F)
allGenes <- allGenes$V1
geneList <- factor(as.integer(allGenes %in% geneSel))
names(geneList) <- allGenes
str(geneList)


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go)
GOdata 
# among 17834 genes, 9348 genes have Go annotation
# in the 9348 genes, 115 genes are significant

# enrichment tests
# take GO hierarchy into account
resultFisher.weit <- runTest(GOdata, algorithm="weight01", statistic="fisher")
allRes.weit <- GenTable(GOdata, classicFisher = resultFisher.weit, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = length(resultFisher.weit@score))
allRes.weit

# Do not take the GO hierarchy into account 
resultFisher.cls <- runTest(GOdata, algorithm="classic", statistic="fisher")
allRes.cls <- GenTable(GOdata, classicFisher = resultFisher.cls, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = length(resultFisher.weit@score))
allRes.cls

# Put the results of the two together
pValue.cls <- score(resultFisher.cls)
pValue.weit <- score(resultFisher.weit)[names(pValue.cls)]
sel.go <- names(pValue.weit)[pValue.weit < 0.05]
sel.fun <- allRes.weit[allRes.weit$GO.ID %in% as.vector(sel.go),]
sel.genes <- genesInTerm(GOdata, sel.go)
gene_list <- NULL
for (i in 1:length(sel.go))
{
  go <- sel.go[i]
  genesforterm <- sel.genes[go][[1]]
  genesforterm <- paste(genesforterm, collapse=';')
  df <- data.frame(gene_list=genesforterm)
  gene_list <- rbind(gene_list,df)
}
sel.go.stat <- cbind(sel.fun, weight = pValue.weit[sel.go], classic = pValue.cls[sel.go], gene=gene_list$gene_list)
write.csv(sel.go.stat, "C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/WGS_weight_GO_enrich_p005.csv", row.names = F, quote = F)


######################  RAiSD+OmegaPlus #######################
# read the go annotation database from Trinotate
gene2go <- readMappings("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/Nriv_trinotate_BP_GOs", sep = "\t", IDsep = ",")
gene2go[lapply(gene2go,length)==0] <-  NA

# import the final gene list from the three way intersection
geneSel <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/RAiSD_OmegaPlus_intersection.txt", header = F)
geneSel <- geneSel$V1
allGenes <- read.table("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/all_gene", header = F)
allGenes <- allGenes$V1
geneList <- factor(as.integer(allGenes %in% geneSel))
names(geneList) <- allGenes
str(geneList)


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go)
GOdata 
# among 17834 genes, 9348 genes have Go annotation
# in the 9348 genes, 115 genes are significant

# enrichment tests
# take GO hierarchy into account
resultFisher.weit <- runTest(GOdata, algorithm="weight01", statistic="fisher")
allRes.weit <- GenTable(GOdata, classicFisher = resultFisher.weit, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = length(resultFisher.weit@score))
allRes.weit

# Do not take the GO hierarchy into account 
resultFisher.cls <- runTest(GOdata, algorithm="classic", statistic="fisher")
allRes.cls <- GenTable(GOdata, classicFisher = resultFisher.cls, orderBy = "resultFisher", ranksOf = "classicFisher", numChar=1000, topNodes = length(resultFisher.weit@score))
allRes.cls

# Put the results of the two together
pValue.cls <- score(resultFisher.cls)
pValue.weit <- score(resultFisher.weit)[names(pValue.cls)]
sel.go <- names(pValue.weit)[pValue.weit < 0.05]
sel.fun <- allRes.weit[allRes.weit$GO.ID %in% as.vector(sel.go),]
sel.genes <- genesInTerm(GOdata, sel.go)
gene_list <- NULL
for (i in 1:length(sel.go))
{
  go <- sel.go[i]
  genesforterm <- sel.genes[go][[1]]
  genesforterm <- paste(genesforterm, collapse=';')
  df <- data.frame(gene_list=genesforterm)
  gene_list <- rbind(gene_list,df)
}
sel.go.stat <- cbind(sel.fun, weight = pValue.weit[sel.go], classic = pValue.cls[sel.go], gene=gene_list$gene_list)
write.csv(sel.go.stat, "C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/Enrichment/RDOMG_weight_GO_enrich_p005.csv", row.names = F, quote = F)



