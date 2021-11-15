#Step 1. Read in file of subcompartment designations.
setwd()
head(file)
colnames(file)<-c("chr","start","end","id")
library(GenomicRanges)
all <- with(file, GRanges(chr, IRanges(start, end)))
clustera1 <- with(subset(file,file[,4]=="a1") , GRanges(chr, IRanges(start, end)))
clustera2 <- with(subset(file,file[,4]=="a2") , GRanges(chr, IRanges(start, end)))
clusterb1 <- with(subset(file,file[,4]=="b1") , GRanges(chr, IRanges(start, end)))
clusterb2 <- with(subset(file,file[,4]=="b2") , GRanges(chr, IRanges(start, end)))

#Step 2. Bin the mm10 genome into 250kb blocks.
library(GenomicRanges)
chrSizes<- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")
bins   <- tileGenome(chrSizes, tilewidth=2.5e5, cut.last.tile.in.chrom=T)
bins

#Step 3. Read in and format HOMER Hi-C TRANS files into [list].
##Note: Calculate interactions per 250kb by subcompartment or whole genome. Example below:
setwd()
library(stringr)

cluster <- NULL
cluster <- data.frame(matrix("", nrow = length(all)))
for (i in 1:length(list)) {
list <- read.table(list[i],header=FALSE,sep=" ")
head(list)
list[,5:6] <- str_split_fixed(file[,5], "=", 2)
list[,6] <- as.numeric(file[,6])
trans.range.rearrange <- list[,2:4]
trans.range.rearrange <- unique(trans.range.rearrange)
colnames(trans.range.rearrange) <- c("CHR", "START", "END")
circos.interactions <-trans.range.rearrange
head(circos.interactions)
circos.interactions <- with(circos.interactions, GRanges(CHR, IRanges(START, END)))
circos <- as.data.frame(countOverlaps(all,circos.interactions))
cluster <- cbind.data.frame(cluster, circos)
}

#Step 4. Calculate Pearson's correlation matrix on samples of interest.
cormat <- round(cor(cluster[,c()]),2)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

library(ggplot2)
col.plot<-c("#FF0000","#FFFFFF","#0070C0")
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(aes(fill = value))+scale_fill_gradient2(breaks=c(0.2,0.6,1),midpoint=0.6,low="#0070C0", mid="#FFFFFF",high="#FF0000")
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(aes(fill = value))+scale_fill_gradient(low="white",high="#0070C0")
