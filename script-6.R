#Step 1. Set working directory and read/format subcompartment designations.
head(file)
colnames(file) <- c("chr","start","end","id")
library(GenomicRanges)
all <- with(file, GRanges(chr, IRanges(start, end)))
all$id <- file$id
all$rowid <- c(1:9853)
head(all)

#Step 2. Read/format UCSC Repeatmasker [reps] file of mm10 ERVKs.
setwd()
head(reps)
colnames(reps) <- c("chr","start","end","name","length","strand")
reps.gr <- with(reps, GRanges(chr, IRanges(start, end)))
all$ervk <- countOverlaps(all,reps.gr)

#Step 3. Select coordinates with high ERVK densities (>99%ile).
perc <- all[(elementMetadata(all)[, "ervk"] >= quantile(all$ervk,0.99))]

#Step 4. Read in SPRET/EiJ and matched C57BL/6J Hi-C matrix files.
library(scales)
setwd()
hicmat <- dir()
names <- c("female.glia.c57",
           "female.glia.spretus",
           "female.neuron.c57",
           "female.neuron.spretus",
           "male.glia.c57",
           "male.glia.spretus",
           "male.neuron.c57",
           "male.neuron.spretus")

#Step 5. Create abbreviated matrices of Hi-C values at high ERVK densities.
for (n in 1:8) {
  setwd("~/Desktop/Akbarian\ Lab/spretus-hicmat/")
  infile <- paste(hicmat[n])
  final <- read.table(infile, header=FALSE ,sep=" ")
  arr <- NULL
  ids <- perc$rowid
  for (p in 1:length(ids)) {
    row <- final[ids[p],]
    col <- row[,ids]
    vec <- as.vector(as.numeric(as.matrix(col)))
    arr <- rbind(arr,vec)
  }
  arr <- as.matrix(arr)
  assign(paste(names[n]), arr)
}  

#Step 6. Sum related matrices.
neuron.c57 <- female.neuron.c57 + male.neuron.c57
neuron.spretus <- female.neuron.spretus + male.neuron.spretus
glia.c57 <- female.glia.c57 + male.glia.c57
glia.spretus <- female.glia.spretus + male.glia.spretus

#Step 7. Generate heatmaps.
my_colors <- colorRampPalette(c("cyan", "deeppink3"))
heatmap(neuron.c57, Colv = NA, Rowv = NA,col = my_colors(100))
heatmap(neuron.spretus, Colv = NA, Rowv = NA,col = my_colors(100))
heatmap(glia.c57, Colv = NA, Rowv = NA,col = my_colors(100))
heatmap(glia.spretus, Colv = NA, Rowv = NA,col = my_colors(100))

#Step 8. Perform statistical comparisons.
wilcox.test(neuron.spretus,glia.spretus,paired=TRUE,alternative="two.sided")$p.value
wilcox.test(neuron.c57,glia.c57,paired=TRUE,alternative="two.sided")$p.value
wilcox.test(neuron.c57,neuron.spretus,paired=TRUE,alternative="two.sided")$p.value
wilcox.test(glia.c57,glia.spretus,paired=TRUE,alternative="two.sided")$p.value
