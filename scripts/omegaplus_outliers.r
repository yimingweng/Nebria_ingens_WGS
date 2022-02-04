library(openintro)
library(ggplot2)
library(moments)
library(robustbase)
library(univOutl)
setwd("C:/Users/wengz/Dropbox/Chapter 3/Nebria_ingens_WGS/omegaplus")
omegaplus.out <- read.table("omegaPlus_out_final")



# using quantile to define outlier and write them into a bed file
quantile=0.995
upper_bound <- quantile(omegaplus.out$V3, quantile)
outliers <- omegaplus.out[which(omegaplus.out$V3 > upper_bound),]
length(outliers$V3)
# 1447 SNPs being defined as outlier
write.csv(outliers, "omegaplus_995quantile.csv", quote=F, row.names = F)
bed <- cbind(outliers$V1, outliers$V2, outliers$V2)
write.table(bed, "omegaplus_995quantile.bed", quote=F, col.names = F, row.names = F, sep="\t")

# visualize the distribution of omega statistic
x11()
ggplot(omegaplus.out, aes(x=log(V3))) +
  geom_density(alpha=.2, fill="lightblue", color="darkblue") +
  geom_vline(xintercept = log(upper_bound), linetype="dotted", color = "red", size=1) +
  #geom_text(aes(x=log(upper_bound), label="cutoff at 99.5% quantile", y=0.6), col="red", angle=90, vjust = 1.2, text=element_text(size=14)) +
  xlab(expression(log(omega))) + ylab("density") +
  xlim(0, max(log(omegaplus.out$V3)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# manhattan plot
omegaplus.out$V4 <- (1:length(omegaplus.out$V3))
x11()
ggplot(omegaplus.out, aes(x=V4, y=log(V3))) + 
  geom_point(size=1, col="black") +
  xlab("position") + 
  ylab(expression(log(omega))) +
  geom_hline(yintercept=log(upper_bound), linetype="dashed", color = "red", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


### other possible ways to get the extreme omega statistics (for comparison only)
#1) considering the distribution of omega statistic and get the most extreme omegas
densityPlot(log10(omegaplus.out$V3), xlim=c(0,3), col="red")
# define the cutoff to be log10(omegaplus.out$V3) > 1.5 (omega > 31.6)
abline(v=1.5, lty=2)
out <- which(log10(omegaplus.out$V3) > 1.5) # omega > 31.6
length(out) # 3876
outliers <- omegaplus.out[which(log10(omegaplus.out$V3) > 1.5), ]
bed <- cbind(outliers$V1, outliers$V2, outliers$V2)
write.table(bed, "omegaplus_316.bed", quote=F, col.names = F, row.names = F, sep="\t")


# 2) find the outlier considering the skewness of the distribution
# visualize the distribution of omega statistic
ggplot(omegaplus.out, aes(x=V3)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")

# calculate the skewness of the distribution (symmetry skewness should be 0)
skewness(omegaplus.out$V3) # 54.81566 (positive, right skewed)

# calculate the medcouple 
mc(omegaplus.out$V3)[1] # 0.3344475
# ???0.6 < mc < 0.6 so use adjusted boxplot by Hubert and Vandervieren (2008)

# use resistant with log transformation distribution
# use k=7 as the extension of the whiskers (Q3+7IQR)
adj_boxplot <- boxB(omegaplus.out$V3,  k = 8, method="resistant", logt=TRUE)

omegaplus.out$no <- c(1:length(omegaplus.out$V3))

out <- as.data.frame(adj_boxplot[["outliers"]])
colnames(out) <- "no"
adj_boxplot.out <- match_df(omegaplus.out, out, on = "no")
bed <- cbind(adj_boxplot.out$V1, adj_boxplot.out$V2, adj_boxplot.out$V2)
write.table(bed, "omegaplus_k7_log_resist.bed", quote=F, col.names = F, row.names = F, sep="\t")
