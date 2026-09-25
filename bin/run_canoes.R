#!/usr/bin/env Rscript

# =============================================================================
# Detailed documentation: run_canoes.R
# Purpose: This pipeline orchestrator script runs the CANOES (Copy Number
#   Variation Detection in Clinical Exome Sequence Data) caller. It parses
#   exome target GC content, read counts, and sample IDs, checks and installs
#   required R dependencies, dynamically sets up parallel processing using
#   multiple CPU cores, invokes the core calling/genotyping functions,
#   and outputs formatted candidate CNV predictions.
#
# Arguments (Positional):
#   - args[1]: Path to GC content file (e.g. `gc.txt`)
#   - args[2]: Path to read counts matrix (e.g. `canoes_reads_new.txt`)
#   - args[3]: Path to sample names file (one sample per line)
#   - args[4]: (Optional) Number of CPU cores for parallelization. Default is 1.
#
# Outputs:
#   - <CHR>_CNVs_pass.csv: High-confidence CNV calls (few overall sample calls).
#   - <CHR>_CNVs_fail.csv: Low-confidence/outlier sample CNV calls.
#   - <CHR>_xcnvs.RData: Serialized R session image for post-processing/debugging.
#   - <CHR>_CNV_plots/: Folder containing visual coverage plots of CNV calls.
#   - <CHR>_CNV_genotype/: Folder containing genotype quality scores per sample.
# =============================================================================

## LOAD THE REQUIRED PACKAGES - INSTALL IF MISSING!
# Packages are checked and automatically installed from CRAN if missing.
if(!require(nnls)){
    install.packages("nnls", dependencies=TRUE)
    library("nnls")
}
if(!require(Hmisc)){
    install.packages("Hmisc", dependencies=TRUE)
    library("Hmisc")
}
if(!require(mgcv)){
    install.packages("mgcv", dependencies=TRUE)
    library("mgcv")
}
if(!require(plyr)){
    install.packages("plyr", dependencies=TRUE)
    library("plyr")
}

## GET USER INPUT
args = commandArgs(trailingOnly=TRUE)
in.gc <- args[1]
in.reads <- args[2]
in.samples <- args[3]
num_cores <- if (length(args) >= 4) as.integer(args[4]) else 1

## READ IN DATA
# Extract chromosome name from input filename by stripping the suffix.
chr <- gsub("_canoes_reads_new.txt", "", in.reads)
gc <- read.table(in.gc)$V2
canoes.reads <- read.table(in.reads)

## RENAME CANOES.READS COLUMNS
sample.names <- scan(in.samples, what="", sep="\n")
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)

## CREATE A VECTOR OF CONSECUTIVE TARGET IDS
target <- seq(1, nrow(canoes.reads))

## COMBINE THE DATA INTO ONE DATA FRAME
canoes.reads <- cbind(target, gc, canoes.reads)

## CALL CNVS IN THE SAMPLES
# Locate CANOES.R in the same directory as this script to source functions.
args0 <- commandArgs(FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
script_dir  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else "."
source(file.path(script_dir, "CANOES.R"))

library(parallel)

## CREATE A VECTOR TO HOLD THE RESULTS FOR EACH SAMPLE
xcnv.pass <- vector('list', length(sample.names))
xcnv.fail <- vector('list', length(sample.names))

## CALL CNVS IN EACH SAMPLE IN PARALLEL
# Call CNVs using mclapply. If a sample produces >= 50 calls, it is flagged as failing
# QC (possible noise outlier) and separated from the passing list.
cat("\nI am now finding CNVs in parallel using", num_cores, "cores\n")
results <- mclapply(1:length(sample.names), function(i) {
    sname <- sample.names[i]
    res <- list(pass = NULL, fail = NULL)
    tryCatch(
    {
        pass_calls <- CallCNVs(sname, canoes.reads)
        if (!is.null(pass_calls) && nrow(pass_calls) >= 50) {
            res$fail <- pass_calls
        } else {
            res$pass <- pass_calls
        }
    },
    error = function(e) {
        message(paste('An Error Occurred in sample:', sname))
        print(e)
    },
    warning = function(w) {
        message(paste('A Warning Occurred in sample:', sname))
        print(w)
    }
    )
    return(res)
}, mc.cores = num_cores)

for (i in 1:length(sample.names)){
    xcnv.pass[[i]] <- results[[i]]$pass
    xcnv.fail[[i]] <- results[[i]]$fail
}
cat("-- Done finding CNVs\n")

## COMBINE THE RESULTS INTO ONE DATA FRAME - WRITE TO CSV AND SAVE THE SESSION
xcnvs.pass <- do.call('rbind', xcnv.pass)
xcnvs.fail <- do.call('rbind', xcnv.fail)
write.table(xcnvs.pass, file=paste0(chr, "_CNVs_pass.csv"), sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)
write.table(xcnvs.fail, file=paste0(chr, "_CNVs_fail.csv"), sep = "\t", row.names=FALSE, col.names = TRUE, quote=FALSE)
save.image(file=paste0(chr, "_xcnvs.RData"))

## CREATE DIRECTORIES TO PUT GENOTYPE INFO AND PLOTS
dir_plots <- paste0(chr, "_CNV_plots")
dir_geno <- paste0(chr, "_CNV_genotype")
dir.create(dir_plots)
dir.create(dir_geno)

unique_samples <- unique(xcnvs.pass$SAMPLE)
if (length(unique_samples) > 0) {
    cat("\nI am now genotyping and plotting CNVs in parallel\n")
    mclapply(unique_samples, function(the.sample) {
        ## GENOTYPE ALL THE CNV CALLS MADE ABOVE PER SAMPLE
        # Calculate detailed quality metrics for the predicted intervals.
        the.name <- paste0(the.sample, "_", chr, "_genotype")
        gen_df <- GenotypeCNVs(xcnvs.pass, the.sample, canoes.reads)
        write.table(gen_df, file=paste0(dir_geno, "/", the.name, ".csv"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        ## PLOT ALL THE CNV CALLS PER SAMPLE    
        # Generate diagnostic coverages vs reference distribution plots.
        pdf(paste0(dir_plots, "/", the.sample, "_", chr, "_CNVs_plot.pdf"))
        the.sample.cnvs <- xcnvs.pass[xcnvs.pass$SAMPLE==the.sample,]
        if (the.sample %in% the.sample.cnvs$SAMPLE){
            for (i in 1:nrow(the.sample.cnvs)){
                PlotCNV(canoes.reads, the.sample.cnvs[i, "SAMPLE"], the.sample.cnvs[i, "TARGETS"])        
            }
        }
        dev.off()
    }, mc.cores = num_cores)
    cat("-- Done genotyping and plotting\n")
}
