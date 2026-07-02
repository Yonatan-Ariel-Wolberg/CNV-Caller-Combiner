#!/usr/bin/env Rscript

# =============================================================================
# Detailed documentation: custom_ref_panel.R
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

library(ggplot2)
library(FNN)

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.pca <- args[1]
in.qcs <- args[2]
in.sex <- args[3]

## GET THE PCA DATA
pca.data <- read.table(in.pca, col.names=c("SAMPLE", "PC1", "PC2", "PC3", "PC4"), colClasses=c("character", rep("numeric", 4)))

## PLOT PCA
pdf("pca.pdf")
ggplot(pca.data, aes(x = PC1, y = PC2)) +
    geom_point()
dev.off()

## load("/home/phelelani/applications/clamms/data/example_qcs.Rdata")
qcs.data <- read.table(in.qcs, header=TRUE, sep='\t')
## selected_metrics <- c("AT_DROPOUT","GC_DROPOUT","ON_BAIT_VS_SELECTED","PCT_PF_UQ_READS","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_50X")
## qcs.data <- qcs.data[which(col.names %in% selected_metrics)]

# CREATE A SCALED COPY OF THE DATA FRAME 
the.data.scaled <- pca.data
for (i in 2:ncol(the.data.scaled)) {
    mini <- min(the.data.scaled[,i])
    maxi <- max(the.data.scaled[,i])
    the.data.scaled[,i] <- apply(the.data.scaled, 1, function(row) { 
		row[[i]] <- (as.numeric(row[[i]]) - mini) / (maxi - mini)
	} )
}

# GET K-NEAREST NEIGHBORS FOR EACH SAMPLE
k.param <- min(20, nrow(the.data.scaled) - 1)
knns <- get.knn(the.data.scaled[,c(seq(2,ncol(the.data.scaled)))], k=k.param, algorithm="kd_tree")

# GENERATE A SINGLE FILE FOR EACH SAMPLE LISTING ITS K-NEAREST NEIGHBOR SAMPLE IDS
for (i in 1:nrow(the.data.scaled)) {
    fname <- paste(the.data.scaled$SAMPLE[i], ".", k.param, "nns.txt", sep="")
    nn.sampleids <- the.data.scaled$SAMPLE[ knns$nn.index[i,] ]
    write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
}

# Read sex information
sex.data <- read.table(in.sex, col.names=c("SAMPLE", "SEX"), colClasses=c("character", "character"))

# GENERATE A SINGLE FILE FOR EACH SAMPLE LISTING ITS K-NEAREST NEIGHBOR SAMPLE IDS AND SEX FOR CLAMMS FIT_MODELS
for (i in 1:nrow(the.data.scaled)) {
    nn.sampleids <- the.data.scaled$SAMPLE[ knns$nn.index[i,] ]
    nn.paths <- paste(nn.sampleids, ".norm.cov.bed", sep="")
    nn.sex <- sex.data$SEX[match(nn.sampleids, sex.data$SAMPLE)]
    
    ref.panel.df <- data.frame(PATH=nn.paths, SEX=nn.sex)
    if (any(is.na(ref.panel.df$SEX))) {
        warning(paste("Some nearest neighbors for sample", the.data.scaled$SAMPLE[i], "were not found in sexinfo. Removing them."))
        ref.panel.df <- ref.panel.df[!is.na(ref.panel.df$SEX), ]
    }
    
    fname <- paste(the.data.scaled$SAMPLE[i], ".ref.panel.files.txt", sep="")
    write.table(ref.panel.df, fname, quote=F, sep="\t", row.names=F, col.names=F)
}

# TO CHECK HOW WELL EACH SAMPLE'S KNNS FIT, COMPUTE THE DISTANCE TO ITS KNN CLUSTER MEAN
the.data.scaled$DistanceToClusterMean <- sapply(1:nrow(the.data.scaled),  function(x) {
    this.knns <- knns$nn.index[x,];
    center <- colMeans(the.data.scaled[this.knns, 2:ncol(the.data.scaled)]);
    return(as.numeric(dist(rbind(as.numeric(the.data.scaled[x, 2:ncol(the.data.scaled)]), as.numeric(center)))))
})

## PLOT DISTANCE DISTRIBUTION
pdf("distance_distribution.pdf")
plot(ecdf(the.data.scaled$DistanceToClusterMean))
dev.off()
