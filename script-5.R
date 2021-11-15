#Step 1. Read in [file] containing subcompartment designations.
##Create separate variables for each subcompartment.
setwd()
head(file)
library(GenomicRanges)
all <- with(file, GRanges(chr, IRanges(start, end)))
clustera1 <- with(subset(file,file[,4]=="a1") , GRanges(chr, IRanges(start, end)))
clustera2 <- with(subset(file,file[,4]=="a2") , GRanges(chr, IRanges(start, end)))
clusterb1 <- with(subset(file,file[,4]=="b1") , GRanges(chr, IRanges(start, end)))
clusterb2 <- with(subset(file,file[,4]=="b2") , GRanges(chr, IRanges(start, end)))

#Step 2. Read in mm10 Repeatmasker .bed [repeats] file.
setwd()
head(repeats)
colnames(repeats) <- c("chr","start","end","name","type","family")

repeats <- with(repeats, GRanges(chr, IRanges(start, end)))
ERV1 <- with(repeats[grep("ERV1",repeats$family),] , GRanges(chr, IRanges(start, end)))
ERVK <- with(repeats[grep("ERVK",repeats$family),] , GRanges(chr, IRanges(start, end)))
ERVL <- with(repeats[grep("ERVL",repeats$family),] , GRanges(chr, IRanges(start, end)))
LINE <- with(repeats[grep("LINE",repeats$type),] , GRanges(chr, IRanges(start, end)))
SINE <- with(repeats[grep("SINE",repeats$type),] , GRanges(chr, IRanges(start, end)))
Simple_repeat <- with(repeats[grep("Simple_repeat",repeats$type),] , GRanges(chr, IRanges(start, end)))
Satellite <- with(repeats[grep("Satellite",repeats$type),] , GRanges(chr, IRanges(start, end)))
DNA <- with(repeats[grep("DNA",repeats$type),] , GRanges(chr, IRanges(start, end)))
tRNA <- with(repeats[grep("tRNA",repeats$type),] , GRanges(chr, IRanges(start, end)))

test <- repeats[!repeats %in% (ERV1),]
test <- test[!test %in% (ERVK),]
test <- test[!test %in% (ERVL),]
test <- test[!test %in% (LINE),]
test <- test[!test %in% (SINE),]
test <- test[!test %in% (Simple_repeat),]
test <- test[!test %in% (Satellite),]
test <- test[!test %in% (DNA),]
Other <- test[!test %in% (tRNA),]

#Step 3. Bin the mm10 genome into 250kb blocks.
library(GenomicRanges)
chrSizes<- c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566)
names(chrSizes) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")
bins   <- tileGenome(chrSizes, tilewidth=2.5e5, cut.last.tile.in.chrom=T)
bins

#Step 4. Perform Fisher's 2x2 testing for repeat densities per subcompartment (>99th%ile)
library(dplyr)
library(exact2x2)

list <- list(ERV1,ERVK,ERVL,LINE,SINE,Simple_repeat,Satellite,DNA,tRNA,Other)

for (m in 1:length(list)) {
  gr <- list[[m]]
  overlaps <- countOverlaps(bins,gr)
  df <- NULL
  clustera1.yes <- countOverlaps(clustera1,gr)
  clustera1.not <- countOverlaps(bins[!bins %in% (clustera1),],gr)
  clustera2.yes <- countOverlaps(clustera2,gr)
  clustera2.not <- countOverlaps(bins[!bins %in% (clustera2),],gr)
  clusterb1.yes <- countOverlaps(clusterb1,gr)
  clusterb1.not <- countOverlaps(bins[!bins %in% (clusterb1),],gr)
  clusterb2.yes <- countOverlaps(clusterb2,gr)
  clusterb2.not <- countOverlaps(bins[!bins %in% (clusterb2),],gr)
  samples.y <- mget(ls(pattern="yes"))
  samples.n <- mget(ls(pattern="not"))
  for (i in 1:length(samples.y)) {
    for (j in 1:length(samples.n)) {
      if (i != j) {
        next
      }
      else {
        yy <- sum(as.numeric(unlist(samples.y[i])) > quantile(overlaps, .99))
        yn <- sum(as.numeric(unlist(samples.y[i])) <= quantile(overlaps, .99))
        ny <- sum(as.numeric(unlist(samples.n[j])) > quantile(overlaps, .99))
        nn <- sum(as.numeric(unlist(samples.n[j])) <= quantile(overlaps, .99))
        data <- matrix(c(nn,ny,yn,yy), 2, 2, byrow=TRUE)
        df <- cbind(df,fisher.exact(data, y = NULL, or = 1, alternative = "greater",
                                                      tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95,
                                                      tol = 0.00001, midp=FALSE)$p.value)
      }
    }
  }
  assign(paste("fish",p,sep=""),df)
}

all <- rbind(fish1,fish2,fish3,fish4,fish5,fish6,fish7,fish8,fish9,fish10)
all.final <- -log10(final)
all.final