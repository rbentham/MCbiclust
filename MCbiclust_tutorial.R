# New package testing.
library(MCbiclust)


# Load example CCLE data and mitochondrial genes
data(CCLE_data)
data(Mitochondrial_genes)

# Select mitochondrial genes in CCLE
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]

# Generate possible seeds in such a way that seeds are unlikely to contain many identical elements.
# This function is useful when running FindSeed multiple times on a HPC
random.seed.list <- SeedGenerator(seed.size = 10,numbers = 10,sample.length = dim(CCLE.mito)[2],break.num = 1,attempts = 100)

# Test intersections between seeds
# Make list of every possible pair of the 1000 generated seeds
combinations <- lapply(apply(combn(1000,2),2,list),unlist)

# Function to compare the length of intersection
len.inter <- function(x,y) length(intersect(x,y))

intersections <- unlist(lapply(combinations,function(x) do.call(len.inter,random.seed.list[x])))
max(intersections)
sum(intersections)

# In comparison, simply generate seeds purely randomly
random.seed.matrix <- t(replicate(1000, sample(c(1:dim(CCLE.mito)[2]),10)))
random.seed.list2 <- lapply(apply(random.seed.matrix, 1, list),unlist)
intersections2 <- unlist(lapply(combinations,function(x) do.call(len.inter,random.seed.list2[x])))
max(intersections2)
sum(intersections)

# Find bicluster seed for first of these initial seeds
CCLE.seed <- FindSeed(gem = CCLE.mito,seed.size = 10,iterations = 1000, initial.seed = random.seed.list[[1]])

# Can also calculate the correlation score from a specific seed as follows:
CorScoreCalc(CCLE.mito,CCLE.seed)

CCLE.mito.cor <- cor(t(CCLE.mito[,CCLE.seed]))
heatmap.2(CCLE.mito.cor,trace = "none")

CCLE.hicor.genes <- HclustGenesHiCor(CCLE.mito,CCLE.seed,cuts = 8)
CCLE.mito.cor2 <- cor(t(CCLE.mito[as.numeric(CCLE.hicor.genes),CCLE.seed]))
CCLE.heat <- heatmap.2(CCLE.mito.cor2,trace = "none")

# Silhoette hclust groupings

# Order Samples by strength of correlation
CCLE.samp.sort <- SampleSort(CCLE.mito[as.numeric(CCLE.hicor.genes),],seed = CCLE.seed,num.cores = 3)

# Plot correlation score vector across ordered samples
CCLE.cor.score.vec <- CalcCorScoreVec(CCLE.samp.sort,CCLE.mito[as.numeric(CCLE.hicor.genes),])
plot(c(10:967),CCLE.cor.score.vec,xlab="Sample size",ylab="Correaltion Score")

# There are three main ways of plotting the "fork"
CCLE.cor.hclust <- hclust(dist(CCLE.mito.cor2))
mito.gene.group1 <- labels(as.dendrogram(CCLE.cor.hclust)[[1]])
mito.gene.group2 <- labels(as.dendrogram(CCLE.cor.hclust)[[2]])

# 1 Average gene expression
Average.mito.g1 <- colMeans(CCLE.mito[as.numeric(CCLE.hicor.genes), CCLE.samp.sort][mito.gene.group1,])
Average.mito.g2 <- colMeans(CCLE.mito[as.numeric(CCLE.hicor.genes), CCLE.samp.sort][mito.gene.group2,])

plot(Average.mito.g1, ylab = "Average gene group 1", xlab = "Index")
plot(Average.mito.g2, ylab = "Average gene group 2", xlab = "Index")

# 2 Point score calculation

## Input
## 1: Gene expression matrix
## 2: location or names of gene group 1
## 3: location or names of gene group 2
CCLE.point.score <- PointScoreCalc(gem = CCLE.mito[as.numeric(CCLE.hicor.genes), CCLE.samp.sort],
                                   gloc1 = mito.gene.group1,gloc2 = mito.gene.group2)

plot(CCLE.point.score, ylab = "Point Score", xlab = "Index")

# Principal component analysis

# Run pca on first 10 samples
pca.matrix <- CCLE.mito[as.numeric(CCLE.hicor.genes), CCLE.samp.sort[c(1:10)]]

pca.results <- prcomp(t(pca.matrix),scores=TRUE,cor=TRUE,center=TRUE)
pc1.percent.var <- (pca.results$sdev/sum(pca.results$sdev))[1]

pca.loadings <- pca.results$rotation

pc1.vec <- seq(length = dim(CCLE.mito)[2])
sum.res.vec <- seq(length = dim(CCLE.mito)[2])

hi.cor.matrix <- CCLE.mito[as.numeric(CCLE.hicor.genes), ]

for(i1 in seq(length = dim(CCLE.mito)[2])){
  pca.l <- as.matrix(pca.loadings)
  hi.cor.mat.samp <- as.matrix(hi.cor.matrix[,CCLE.samp.sort[i1]])
  pc1.vec[i1]<-lsfit(pca.l,hi.cor.mat.samp)$coef[2]
  sum.res.vec[i1]<-sum(lsfit(pca.l,hi.cor.mat.samp)$residuals^2)
}

plot(pc1.vec, ylab="PC1",xlab="Index")
# Calculate correlation vector



# Gene set enrichment analysis

# Mann-Whitney

# Z test