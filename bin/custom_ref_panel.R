#!/usr/bin/env Rscript

# =============================================================================
# Detailed documentation: custom_ref_panel.R
# Purpose: This pipeline script calculates the custom reference panel files for
#   CLAMMS callers. It uses PCA coordinates (from SVD) and Picard metrics (QCs)
#   along with sample sex information to determine the K-Nearest Neighbors (KNN)
#   for each sample. These neighbors serve as the matched reference panel.
#
# Arguments (Positional):
#   - args[1]: Path to PCA coordinates file (e.g. `pca.coordinates.txt`)
#   - args[2]: Path to QC metrics file (e.g. `qcs_metrics`)
#   - args[3]: Path to Sex information file (e.g. `sex_info.txt`)
#
# Outputs:
#   - pca.pdf: PCA plot of PC1 vs PC2.
#   - <SAMPLE>.<K>nns.txt: A list of K-nearest neighbor sample IDs for each sample.
#   - <SAMPLE>.ref.panel.files.txt: Tab-separated path-to-coverage and sex mappings for the KNNs.
#   - distance_distribution.pdf: Cumulative distribution function (ECDF) plot of cluster mean distances.
# =============================================================================

library(ggplot2)
library(FNN)

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.pca <- args[1]
in.qcs <- args[2]
in.sex <- args[3]

## GET THE PCA DATA
# Read sample PCA coordinates (first column is sample name, columns 2-5 are PC1 to PC4)
pca.data <- read.table(in.pca, col.names=c("SAMPLE", "PC1", "PC2", "PC3", "PC4"), colClasses=c("character", rep("numeric", 4)))

## PLOT PCA
# Generate a diagnostic plot of PC1 vs PC2 to visualize cohort clustering
pdf("pca.pdf")
ggplot(pca.data, aes(x = PC1, y = PC2)) +
    geom_point()
dev.off()

## load("/home/phelelani/applications/clamms/data/example_qcs.Rdata")
# Read quality metrics
qcs.data <- read.table(in.qcs, header=TRUE, sep='\t')
## selected_metrics <- c("AT_DROPOUT","GC_DROPOUT","ON_BAIT_VS_SELECTED","PCT_PF_UQ_READS","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_50X")
## qcs.data <- qcs.data[which(col.names %in% selected_metrics)]

# CREATE A SCALED COPY OF THE DATA FRAME
# Scale each PC column to the range [0, 1] to prevent magnitude differences from dominating KNN distance metrics.
the.data.scaled <- pca.data
for (i in 2:ncol(the.data.scaled)) {
    mini <- min(the.data.scaled[,i])
    maxi <- max(the.data.scaled[,i])
    if (maxi > mini) {
        the.data.scaled[,i] <- (the.data.scaled[,i] - mini) / (maxi - mini)
    } else {
        the.data.scaled[,i] <- 0.0
    }
}

# GET K-NEAREST NEIGHBORS FOR EACH SAMPLE
# K is chosen as min(20, total_samples - 1) to build a cohort of up to 20 reference samples.
k.param <- max(1, min(20, nrow(the.data.scaled) - 1))
if (nrow(the.data.scaled) <= 1) {
    knns <- list(nn.index = matrix(1, nrow=1, ncol=1))
} else {
    knns <- get.knn(the.data.scaled[,c(seq(2,ncol(the.data.scaled)))], k=k.param, algorithm="kd_tree")
}

# GENERATE A SINGLE FILE FOR EACH SAMPLE LISTING ITS K-NEAREST NEIGHBOR SAMPLE IDS
for (i in 1:nrow(the.data.scaled)) {
    fname <- paste(the.data.scaled$SAMPLE[i], ".", k.param, "nns.txt", sep="")
    nn.sampleids <- the.data.scaled$SAMPLE[ knns$nn.index[i,] ]
    write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
}

# Read sex information
sex.data <- read.table(in.sex, col.names=c("SAMPLE", "SEX"), colClasses=c("character", "character"))

# GENERATE A SINGLE FILE FOR EACH SAMPLE LISTING ITS K-NEAREST NEIGHBOR SAMPLE IDS AND SEX FOR CLAMMS FIT_MODELS
# This maps each neighbor's path to its corresponding sex code (M/F/U).
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
# This serves as a quality metrics metric to detect outlier samples in the cohort.
the.data.scaled$DistanceToClusterMean <- sapply(1:nrow(the.data.scaled),  function(x) {
    this.knns <- knns$nn.index[x,];
    center <- colMeans(the.data.scaled[this.knns, 2:ncol(the.data.scaled)]);
    return(as.numeric(dist(rbind(as.numeric(the.data.scaled[x, 2:ncol(the.data.scaled)]), as.numeric(center)))))
})

## PLOT DISTANCE DISTRIBUTION
pdf("distance_distribution.pdf")
plot(ecdf(the.data.scaled$DistanceToClusterMean))
dev.off()
