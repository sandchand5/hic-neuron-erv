#Step 0. Use .hicdump to create input files (UNIX).

#Step 1. Set working directory to folder with input files.
setwd()

#Step 2. Create matrix with odd rows and even columns.
##Do this step for each Hi-C sample.
for (j in c(2,4,6,8,10,12,14,16,18)) {
  var <- NULL
  for (i in c(1,3,5,7,9,11,13,15,17,19)) {
    infile <- paste(i, "_" , j, ".txt", sep="")
    d.frame <- read.table(infile, header=FALSE ,sep="\t")
    if (i<j){
      d.frame[,ncol(d.frame)] <- NULL
    } else {
      d.frame <- t(d.frame)
    }
    var <- rbind.data.frame(var, d.frame)
    var <- var[rowSums(is.na(var)) != ncol(var),]
    assign(paste("chr", j, sep = ""), var)
  }
}

hicgenome <- cbind.data.frame(chr2, chr4, chr6, chr8, chr10, chr12, chr14, chr16, chr18)

##Note: Lengths correspond to chromosomal lengths (mm10) in 250kb bins.
rownames(hicgenome) <- c(paste("chr1:", seq(from = 1, to = 195500001-250000, by = 250000), "-", seq(from = 250000, to = 195500000, by = 250000) , sep=""),
                                    paste("chr3:", seq(from = 1, to = 160000001-250000, by = 250000), "-", seq(from = 250000, to = 160000000, by = 250000) , sep=""),
                                    paste("chr5:", seq(from = 1, to = 151750001-250000, by = 250000), "-", seq(from = 250000, to = 151750000, by = 250000) , sep=""),
                                    paste("chr7:", seq(from = 1, to = 145500001-250000, by = 250000), "-", seq(from = 250000, to = 145500000, by = 250000) , sep=""),
                                    paste("chr9:", seq(from = 1, to = 124500001-250000, by = 250000), "-", seq(from = 250000, to = 124500000, by = 250000) , sep=""),
                                    paste("chr11:", seq(from = 1, to = 122000001-250000, by = 250000), "-", seq(from = 250000, to = 122000000, by = 250000) , sep=""),
                                    paste("chr13:", seq(from = 1, to = 120500001-250000, by = 250000), "-", seq(from = 250000, to = 120500000, by = 250000) , sep=""),
                                    paste("chr15:", seq(from = 1, to = 104000001-250000, by = 250000), "-", seq(from = 250000, to = 104000000, by = 250000) , sep=""),
                                    paste("chr17:", seq(from = 1, to = 95000001-250000, by = 250000), "-", seq(from = 250000, to = 95000000, by = 250000) , sep=""),
                                    paste("chr19:", seq(from = 1, to = 61500001-250000, by = 250000), "-", seq(from = 250000, to = 61500000, by = 250000) , sep=""))
colnames(hicgenome.fem) <- c(paste("chr2:", seq(from = 1, to = 182250001-250000, by = 250000), "-", seq(from = 250000, to = 182250000, by = 250000) , sep=""),
                                    paste("chr4:", seq(from = 1, to = 156500001-250000, by = 250000), "-", seq(from = 250000, to = 156500000, by = 250000) , sep=""),
                                    paste("chr6:", seq(from = 1, to = 149750001-250000, by = 250000), "-", seq(from = 250000, to = 149750000, by = 250000) , sep=""),
                                    paste("chr8:", seq(from = 1, to = 129500001-250000, by = 250000), "-", seq(from = 250000, to = 129500000, by = 250000) , sep=""),
                                    paste("chr10:", seq(from = 1, to = 130750001-250000, by = 250000), "-", seq(from = 250000, to = 130750000, by = 250000) , sep=""),
                                    paste("chr12:", seq(from = 1, to = 120250001-250000, by = 250000), "-", seq(from = 250000, to = 120250000, by = 250000) , sep=""),
                                    paste("chr14:", seq(from = 1, to = 125000001-250000, by = 250000), "-", seq(from = 250000, to = 125000000, by = 250000) , sep=""),
                                    paste("chr16:", seq(from = 1, to = 98250001-250000, by = 250000), "-", seq(from = 250000, to = 98250000, by = 250000) , sep=""),
                                    paste("chr18:", seq(from = 1, to = 90750001-250000, by = 250000), "-", seq(from = 250000, to = 90750000, by = 250000) , sep=""))

hicgenome.fem[is.na(hicgenome.fem)] <- 0

#Step 3. Add all Hi-C matrices under name "sum".
head(sum)

#Step 4. Perform k-means clustering with cluster numbers ranging from 2 to 10.
library(tidyverse)
library(cluster)
library(factoextra)

for (k in c(2:10)) {
  var <- NULL
  var <- kmeans(sum, centers = k, nstart = 10, iter.max = 25)
  assign(paste("o", k, sep = ""), var)
}

for (k in c(2:10)) {
  var <- NULL
  var <- kmeans(t(sum), centers = k, nstart = 10, iter.max = 25)
  assign(paste("e", k, sep = ""), var)
}

#Step 5. Create a table with all cluster designations for k=2 to k=10.
##Note: "o" is for odds and "e" is for evens.
odds <- cbind.data.frame(o2$cluster,o3$cluster,o4$cluster,o5$cluster,o6$cluster,o7$cluster,o8$cluster,o9$cluster,o10$cluster)
evens <- cbind.data.frame(e2$cluster,e3$clustereo4$cluster,e5$cluster,e6$cluster,e7$cluster,e8$cluster,e9$cluster,e10$cluster)