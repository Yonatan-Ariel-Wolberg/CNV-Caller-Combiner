#!/usr/bin/env Rscript

# =============================================================================
# Detailed documentation: run_canoes.R
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================

## LOAD THE REQUIRED PACKAGES - INSTALL IF MISSING!
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
# Locate CANOES.R in the same directory as this script
args0 <- commandArgs(FALSE)
script_path <- sub("--file=", "", args0[grep("--file=", args0)])
script_dir  <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else "."
source(file.path(script_dir, "CANOES.R"))

library(parallel)

## CREATE A VECTOR TO HOLD THE RESULTS FOR EACH SAMPLE
xcnv.pass <- vector('list', length(sample.names))
xcnv.fail <- vector('list', length(sample.names))

## CALL CNVS IN EACH SAMPLE IN PARALLEL
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
        the.name <- paste0(the.sample, "_", chr, "_genotype")
        gen_df <- GenotypeCNVs(xcnvs.pass, the.sample, canoes.reads)
        write.table(gen_df, file=paste0(dir_geno, "/", the.name, ".csv"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        ## PLOT ALL THE CNV CALLS PER SAMPLE    
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
