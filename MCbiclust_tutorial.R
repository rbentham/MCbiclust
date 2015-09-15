# New package testing.
library(MCbiclust)
library(gplots)
library(ggplot2)

# 1. Load example CCLE data and mitochondrial genes ----
data(CCLE_data)
data(Mitochondrial_genes)

# Select mitochondrial genes in CCLE
mito.loc <- which(as.character(CCLE_data[,2]) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_data[mito.loc,-c(1,2)]
row.names(CCLE.mito) <- CCLE_data[mito.loc,2]


# 2. Find bicluster seed for first of these initial seeds ----
set.seed(102)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 10000,
                      initial.seed = random.seed.list[[1]])

# Can also calculate the correlation score from a specific seed as follows:
CorScoreCalc(CCLE.mito,CCLE.seed)

# Calculate the correlation matrix
CCLE.mito.cor <- cor(t(CCLE.mito[,CCLE.seed]))
heatmap.2(CCLE.mito.cor,trace = "none")

# Only select the most highly correlating genes. 
# See help(HclustGenesHiCor) for details
CCLE.hicor.genes <- HclustGenesHiCor(CCLE.mito,CCLE.seed,cuts = 8)
CCLE.mito.cor2 <- cor(t(CCLE.mito[as.numeric(CCLE.hicor.genes),CCLE.seed]))
CCLE.heat <- heatmap.2(CCLE.mito.cor2,trace = "none")

# Can examine the two anti-correlating gene lists
CCLE.groups <- list(labels(CCLE.heat$rowDendrogram[[1]]),
                    labels(CCLE.heat$rowDendrogram[[2]]))

# 3. Calculate the Correlation Vector ----

# More generally can calculate a correlation vector to measure how 
# each gene measured in the gene expression matrix matches this
# pattern.

CCLE.gene.vec <- GeneVecFun(CCLE.mito, CCLE.seed, 10)
CCLE.cor.vec <- CalcCorVector(gene.vec = CCLE.gene.vec,
                              gem = CCLE_data[,-c(1,2)][,CCLE.seed])

# 4. Gene Set Enrichment ----

# Mann-Whitney is used as distributions are not normal.
# Note there are many online methods of finding significant gene sets (e.g. gprofiler)

GSE.MW <- GOEnrichmentAnalysis(gene.names = as.character(CCLE_data[,2]),
                               gene.values = CCLE.cor.vec,
                               sig.rate = 0.05)

                      

# 5. Sample ordering ----

# Can calculate ordering with only the highly correlating genes
# or the original gene list, recommended to do so with the highly
# correlating genes.

# Can run on multiple cores to speed up, sort.length = NULL (default) 
# will sort all the samples. Still takes awhile to run. Consider
# running on HPC or leaving for long time.

CCLE.samp.sort <- SampleSort(CCLE.mito[as.numeric(CCLE.hicor.genes),],
                             seed = CCLE.seed,num.cores = 3,
                             sort.length = 100)

# 6. PCA and ggplot2 ----

# function takes input as the submatrix of highly correlated genes,
# sample ordering and number of samples to calculate the principal components with

top.mat <- CCLE.mito[as.numeric(CCLE.hicor.genes),]
pc1.vec <- PC1VecFun(top.mat, CCLE.samp.sort, 10)

# Also can calculate average gene set values
# can use correlation vector values, or actual gene set values
# from heatmap

# rowmeans(CCLE.mito[CCLE.groups[[1]],CCLE.samp.sort])
# rowmeans(CCLE.mito[CCLE.groups[[2]],CCLE.samp.sort])

# rowmeans(CCLE_data[order(CCLE.cor.vec,decreasing = T)[seq(1000)],CCLE.samp.sort])
# rowmeans(CCLE_data[order(CCLE.cor.vec,decreasing = F)[seq(1000)],CCLE.samp.sort])

# Can use kmeans clustering to separate samples into upper/lower fork
fork.class <- ForkClassifier(pc1.vec,100)
fork.status <- rep("Normal",100)
fork.status[fork.class$Upper] <- "Upper"
fork.status[fork.class$Lower] <- "Lower"

CCLE.df <- data.frame(CCLE.name = colnames(CCLE_data)[-c(1,2)][CCLE.samp.sort],
                      PC1 = pc1.vec,
                      Order = seq(length = length(pc1.vec)),
                      Fork = fork.status)

ggplot(CCLE.df, aes(Order,PC1)) +
  geom_point(aes(colour=factor(Fork))) + ylab("PC1")

# 7. Dealing with multiple runs ----

# First calculate up 100 runs with 500 iterations each, 
# multiple runs will typically be calculated with many more
# iterations, but requires high-performance computing to do
# efficiently.

CCLE.multi.seed <- list()
initial.seed1 <- list()

for(i in seq(100)){
  set.seed(i)
  initial.seed1[[i]] <- sample(seq(length = dim(CCLE.mito)[2]),10)
  CCLE.multi.seed[[i]] <- FindSeed(gem = CCLE.mito,
                                   seed.size = 10,
                                   iterations = 500,
                                   initial.seed = initial.seed1[[i]])
}

# Second calculate the corresponding correlation vector for each run
CCLE.cor.vec.multi <- list()

for(i in seq(100)){
  CCLE.gene.vec <- GeneVecFun(CCLE.mito, CCLE.multi.seed[[i]], 10)
  CCLE.cor.vec.multi[[i]] <- CalcCorVector(gene.vec = CCLE.gene.vec,
                                           gem = CCLE_data[,-c(1,2)][,CCLE.multi.seed[[i]]])
}

# Convert correlation vector list to a data matrix
multi.run.cor.vec.mat <- matrix(0,length(CCLE.cor.vec.multi[[1]]),length(CCLE.cor.vec.multi))
for(i in 1:100){
  multi.run.cor.vec.mat[,i] <- CCLE.cor.vec.multi[[i]]
}
rm(CCLE.cor.vec.multi)

# Visualise as a heatmap
routput.corvec.matrix.cor.heat <- heatmap.2(abs(cor((multi.run.cor.vec.mat))),trace="none",
                                            distfun = function(c){as.dist(1 - abs(c))})

# Use Silhouette clustering function to find 
# add random cor vector as comparison for use in clustering

multi.clust.groups <- SilhouetteClustGroups(cor.vec.mat = multi.run.cor.vec.mat,
                                            max.clusters = 20,
                                            plots = T)

# Can visualsise the different patterns found
microarray.genes <- as.character(CCLE_data[,2])
average.corvec <- lapply(X = multi.clust.groups,
                         FUN = function(x) rowMeans(multi.run.cor.vec.mat[,x]))

CVPlot(cv.df = as.data.frame(average.corvec),
        geneset.loc = mito.loc,
        geneset.name = "Mitochondrial",
        alpha1 = 0.1)

# Can also calculate the gene set enrichment from the average correlation vectors.
corvec.gsea <- lapply(X = average.corvec,
                      FUN = function(x) GOEnrichmentAnalysis(gene.names = microarray.genes,
                                                             gene.values = x,
                                                             sig.rate = 0.05))
                        
# Can sort the pattern found, with function that selects the top genes in the average
# ordering and then selects the initial seed with the highest correlation score
# corresponding to those top genes.

CCLE.samp.multi.sort <-  MultiSampleSort(CCLE_data[,-c(1,2)], average.corvec, 750,
                            multi.clust.groups, initial.seed1,2,50)

            
# 8. Comparing results with clinical/other data ----

# In analysis of the results it is important to compare results with clinical data
# this data must be loaded separately into R and then tests of significance must take
# place

# If data is downloaded, load into R in standard way
# CCLE_clinical <- read.delim("~/Documents/PhD/Datasets/CCLE/CCLE_sample_info_file_2012-04-06.txt", dec=",")
# CCLE_copy <- read.delim("~/Documents/PhD/Datasets/CCLE/CCLE_copynumber/CCLE_copynumber_byGene_2012-04-06.txt")

# 1) Tissue of origin
data(CCLE_samples)

# often other data does not have exact same colnames or is missing/has addition values
# this bit of analysis is often messy, part of data wrangling 

CCLE.samples.names <- as.character(CCLE_samples[,1])
CCLE.data.names <- colnames(CCLE_data)[-c(1,2)]

# In this case some samples have an additional X not present in some CCLE_samples data
# so it is necessary to add it for consistency.
CCLE.samples.names[c(1:15)] <- paste("X",CCLE.samples.names[c(1:15)], sep="")
CCLE_samples$CCLE.name <- CCLE.samples.names

library(dplyr)
CCLE.df.samples <- inner_join(CCLE.df,CCLE_samples,by="CCLE.name")

ggplot(CCLE.df.samples, aes(Order,PC1)) +
  geom_point(aes(colour=factor(Site.Primary))) + ylab("PC1")

# 2) Copynumber alterations
data(CCLE_copy)

# This dataset is bigger and more complicated, containing the copynumber alterations (CNA) across
# the entire genome, it can be analysed by looking for differences in two groups.

# remove NA genes
CCLE_copy <- CCLE_copy[-which(is.na(CCLE_copy[,5]) == T),]

CCLE.copy.names <- colnames(CCLE_copy)[-c(1:4)]
which(CCLE.copy.names %in% CCLE.data.names)

a2 <- which(CCLE.copy.names %in% CCLE.data.names)
a3 <- seq(length = length(CCLE.copy.names))[-a2]

a4 <- which(CCLE.df$Fork == "Upper")
a5 <- which(CCLE.df$Fork == "Lower")

av.copy.upper <- rowMeans(CCLE_copy[,-c(1:4)][,-(a3)][,CCLE.samp.sort[a4]])
av.copy.lower <- rowMeans(CCLE_copy[,-c(1:4)][,-(a3)][,CCLE.samp.sort[a5]])

copy.change <- av.copy.upper - av.copy.lower
plot(copy.change,type = "l")

# Find p-values by a permutation test
max.change2 <- list()
for(i in seq(100)){
  set.seed(i)
  rand.g1 <- sample(seq(100),length(a4))
  rand.g2 <- seq(100)[-rand.g1]
  rand.copy1 <- rowMeans(CCLE_copy[,-c(1:4)][,-(a3)][,CCLE.samp.sort[rand.g1]])
  rand.copy2 <- rowMeans(CCLE_copy[,-c(1:4)][,-(a3)][,CCLE.samp.sort[rand.g2]])
  
  max.change2[[i]] <- (abs(rand.copy1 - rand.copy2))
}

max.change2all <- unlist(max.change2)
all.pvalue <- seq(length = length(copy.change))
for(i in 1:length(copy.change)){
  all.pvalue[i] <- length(which(max.change2all > abs(copy.change[i])))/2312400
}

all.pvalue.adj <- p.adjust(all.pvalue)
sig.copy.change <- which(all.pvalue.adj < 0.05)

sig.copy.df <- cbind(CCLE_copy[sig.copy.change,c(1:4)],
                     Average_copy_g1=av.copy.upper[sig.copy.change],
                     Average_copy_g2=av.copy.lower[sig.copy.change],
                     copy_change=copy.change[sig.copy.change],
                     max_pvalue = max.pvalue)

# 9. Initial seed generation for HPC ----
# If running on a HPC need multiple initial seeds, possible to generate
# them such that they are unlikely to contain many identical elements.
# See help(SeedGenerator) for information on the algorithm.

seed.numbers <- 1000
random.seed.list <- SeedGenerator(seed.size = 10,
                                  numbers = seed.numbers,
                                  sample.length = dim(CCLE.mito)[2],
                                  break.num = 1,
                                  attempts = 100)

# Note: for a long list of seeds takes awhile to generate.

# Can demonstrate the effectiveness of this seed generation vs randomly
# generated seeds

# a) Make list of every possible pair of the 1000 generated seeds

combinations <- lapply(apply(combn(seed.numbers,2),2,list),unlist)

# b) Function to compare the length of intersection
len.inter <- function(x,y) length(intersect(x,y))


# c) Generate seeds purely randomly
set.seed(101)
random.seed.matrix <- t(replicate(1000, sample(c(1:dim(CCLE.mito)[2]),10)))
random.seed.list2 <- lapply(apply(random.seed.matrix, 1, list),unlist)

# d) Calculate intersections for both methods
intersections <- unlist(lapply(combinations,function(x) do.call(len.inter,random.seed.list[x])))
intersections2 <- unlist(lapply(combinations,function(x) do.call(len.inter,random.seed.list2[x])))

max(intersections)
sum(intersections)
max(intersections2)
sum(intersections2)
sum(intersections2) - sum(intersections)

rm(combinations)


# 10. Alternative methods based on known gene regulation groups ----

# Point Scoring algorithm - alternative way of scoring samples based on know regulation groups.

CCLE.point.score <- PointScoreCalc(gem = CCLE.mito[as.numeric(CCLE.hicor.genes), CCLE.samp.sort],
                                   gloc1 = CCLE.groups[[1]],gloc2 = CCLE.groups[[2]])

plot(CCLE.point.score, ylab = "Point Score", xlab = "Index")

# FindSeedGroups - adaption of FindSeed algorithm that seeks for a known pattern
# of two gene groups that are strongly correlated with eachother, and anti-correlated
# with the other.
CCLE.seed.groups <- FindSeedGroups(gem = CCLE.mito,seed.size = 10,iterations = 10000,
                                   group1.loc = CCLE.groups[[1]], group2.loc = CCLE.groups[[2]])

CorScoreCalc(CCLE.mito[unlist(CCLE.groups),],CCLE.seed.groups)
heatmap.2(cor(t(CCLE.mito[unlist(CCLE.groups),CCLE.seed.groups])),trace = "none")

# 11. Running on a HPC (example for Legion) ----

# What follows is an example script to be run on a HPC designed for the UCL Legion system
# Note this script is not designed to run as is written, and will need to be adapted to suit
# the needs of the user. Two scripts are needed, a batch script and an R script. The Batch script
# will repeatedly run the RScript with different inputs (different random seeds) the results
# of which are saved to a sub-directory in Scratch.

########################################
# Legion_Example_Batch_Script.sh
########################################

#!/bin/bash -l
# Batch script to run an array job on Legion with the upgraded
# software stack under SGE.
# 1. Force bash
#$ -S /bin/bash
# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:27:0
# 3. Request 1 gigabyte of RAM.
#$ -l mem=1G
# 4. Request 10 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=10G
# 5. Set up the job array.  In this instance we have requested 1000 tasks
# numbered 1 to 1000.
#$ -t 1-1000
# 6. Set the name of the job.
#$ -N [JobName]
# 7. Select the project that this job will run under.
# Find <your_project_id> by running the command "groups"
#$ -P [ProjectName]
# 8. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/[username]/Scratch/R_output

module unload compilers/intel/11.1/072
module unload mpi/qlogic/1.2.7/intel
module unload mkl/10.2.5/035


module load compilers/gnu/4.6.3
module load curl/7.21.3/gnu.4.6.3
module load atlas/3.8.3/gnu.4.6.3
module load gsl/1.15/gnu.4.6.3
module load fftw/3.3.1/double/gnu.4.6.3
module load hdf/5-1.8.7/gnu.4.6.3
module load netcdf/4.2.1.1/gnu.4.6.3
module load root/5.34.09/gnu.4.6.3
module load java/1.6.0_32
module load texlive/2012

module load r/3.0.1-atlas/gnu.4.6.3

Rscript /home/[username]/UCL_Legion_Example_Script.R $SGE_TASK_ID

######## End ##########


########################################
# UCL_Legion_Example_Script.R
########################################

# load necessary arguments
args <- commandArgs(TRUE)
.libpaths("/path/to/R/librarys")

require(MCBiclust)

# load gene expression data to analyse
# data includes:
# 1. Gene expression matrix, (gem)
# 2. Any gene set locations of interest, e.g. mitochondrial (mito.loc)
# 3. (optional) a list of pre-generated initial seeds
load("/path/to/gene/expression/data")

# set random seed
set.seed(as.numeric(args[1]))

# can select a known gene set
gem.mito <- gem[mito.loc,]

# or alternatively can select a random gene set, e.g
# gem.random <- gem[sample(seq(length = dim(gem)[1]),1000),]

# Run the FindSeed algorithm
gem.seed <- FindSeed(gem = gem.mito, seed.size = 10,
                     iterations = 10000,
                     initial.seed = random.seed.list[as.numeric(args[1])])

# Calculate the correlation vector
gem.gene.vec <- GeneVecFun(gem.mito, gem.seed,splits = 20)
gem.cor.vec <- CalcCorVector(gene.vec = gem.gene.vec,
                              gem = gem[,gem.seed])

# Calculate the gene set enrichment

gem.gene.names <- rownames(gem)

GSE.MW <- GOEnrichmentAnalysis(gene.names = gem.gene.names,
                               gene.values = gem.cor.vec,
                               sig.rate = 0.05)

# Save data to files
write.table(data.frame(Seed=gem.seed),file=paste("Seed",args[1],"gem_list.txt",sep="_"))
write.table(data.frame(CV=gem.cor.vec),file=paste("Seed",args[1],"cor_vec.txt",sep="_"))
write.table(GSE.MW,file=paste("Seed",args[1],"GSE_MW.txt",sep="_"))

######## End ##########

# Before logging on to legion any needed files must be copied to
# your home directory. This can be done in terminal with the scp
# command, e.g.

# scp UCL_Legion_Example_Script.R [username]@legion.rc.ucl.ac.uk:/home/[username]/

# First login on terminal with the following:
# ssh [username]@legion.rc.ucl.ac.uk 

# R and a text editior (vim) are already installed on legion.
# You will need to scp any necessary files to run an analysis, 
# as well as any data files needed and have additional R libraries installed somewhere within your
# user directory. All output files must be saved to the Scratch directory

# To run type:
# qsub Legion_Example_Batch_Script.sh

# Your job will be added to the queue and will be completed depending on 
# it's priority. Once complete all output files can be copied from the
# HPC using scp with the following command

# scp ucbprbe@login05.external.legion.ucl.ac.uk:/home/ucbprbe/Scratch/R_output/* .

# Addition information on using Legion can be found:
# https://wiki.rc.ucl.ac.uk/wiki/Legion

# For other HPC systems please adapt these instructions as necessary
