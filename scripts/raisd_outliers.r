library(openintro)

setwd("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/RAiSD/RAiSD_outliers")
pop <- list.files(path = "C:/Users/wengz/Dropbox/Data/Nebria_ingens/RAiSD_out", full.names=TRUE)

q=0.995
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
  write.table(bed, file=paste(population, "995quantile.bed", sep=""), quote=F, col.names = F, row.names = F, sep="\t")
}
dev.off()

