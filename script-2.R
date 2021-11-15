#Step 1. Read in [file] with subcompartment designations.
setwd()
colnames(file) <- c("chr","start","end","id")

chr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")
subcomps <- c("a1","a2","b1","b2")

#####Hi-C TRANS
#Step 2. Apply the following commands on the [sum] Hi-C matrix used to generate the subcompartments.
##Note: Manually substitute subcompartment combinations as desired. Example below:
arr <- NULL
for (m in c(1:18)) {
  sel <- sum[which(ids$chr == chr[m] & ids$id == "a1"),]
  sel <- sel[,(max(which(ids$chr == chr[m]))+1):ncol(sel)]
  ext <- sel[,(which((ids$id == "a1"))-(max(which(ids$chr == chr[m]))))[(which((ids$id == "a1"))-(max(which(ids$chr == chr[m])))) > 0]]
  vec <- as.vector(as.numeric(as.matrix(ext)))
  arr <- append(arr,vec)
}
mean(arr)

#####Hi-C CIS
#Step 3. Apply the following commands on the [sum] Hi-C matrix used to generate the subcompartments.
##Note: Manually substitute subcompartment combinations as desired. Example below:
##Note: This block is for comparisons between differing subcompartments.
arr <- NULL
for (m in 1:19) {
  sel <- sum[which(ids$chr == chr[m] & ids$id == "a1"),]
  sel <- sel[,(min(which(ids$chr == chr[m]))):(max(which(ids$chr == chr[m])))]
  ext <- sel[,which(ids$chr == chr[m] & ids$id == "a2")-(min(which(ids$chr == chr[m])))+1]
  vec <- as.vector(as.numeric(as.matrix(ext)))
  arr <- append(arr,vec)
}
mean(arr)

#Step 3. Apply the following commands on the [sum] Hi-C matrix used to generate the subcompartments.
##Note: Manually substitute subcompartment combinations as desired. Example below:
##Note: This block is for comparisons between identical subcompartments.
arr <- NULL
for (m in c(1:19)) {
  sel <- sum[which(ids$chr == chr[m] & ids$id == "a1"),]
  sel <- sel[,(min(which(ids$chr == chr[m]))):(max(which(ids$chr == chr[m])))]
  ext <- sel[,which(ids$chr == chr[m] & ids$id == "a1")-(min(which(ids$chr == chr[m])))+1]
  mat <- NULL
  for (r in 1:nrow(ext)) {
    tri <- ext[r,r:ncol(ext)]
    mat <- append(mat,tri)
  }
  vec <- as.vector(as.numeric(as.matrix(mat)))
  arr <- append(arr,vec)
}
mean(arr)