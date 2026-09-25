#!/usr/bin/env Rscript

# =============================================================================
# Detailed documentation: run_exomedepth.R
# Purpose: This script executes the ExomeDepth CNV calling algorithm on a single
#   target sample using a pre-computed cohort read-count matrix. It leverages
#   beta-binomial modeling and automated reference sample selection.
#
# Arguments (Positional):
#   - args[1]: Path to count matrix TSV/TXT (header: chromosome start end name <sample1> <sample2> ...)
#   - args[2]: Target sample ID (must correspond to a column name in the matrix)
#   - args[3]: Path to output CSV file
#   - args[4]: (Optional) Minimum correlation threshold for reference selection. Default: 0.97
#   - args[5]: (Optional) Transition probability for the HMM. Default: 1e-4
#
# Outputs:
#   - <output_csv>: CSV file containing called CNVs (start.p, end.p, type, nexons,
#     start, end, chromosome, id, BF, reads.expected, reads.observed, reads.ratio, sample)
# =============================================================================

suppressPackageStartupMessages({
    if (!requireNamespace("ExomeDepth", quietly = TRUE)) {
        stop("ExomeDepth package is not installed in the R environment.")
    }
    library(ExomeDepth)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: run_exomedepth.R <count_matrix.txt> <sample_id> <output_csv> [min_correlation] [transition_prob]\n")
    quit(status = 1)
}

count_matrix_file <- args[1]
sample_id         <- args[2]
output_csv        <- args[3]
min_corr          <- if (length(args) >= 4) as.numeric(args[4]) else 0.97
trans_prob        <- if (length(args) >= 5) as.numeric(args[5]) else 1e-4

cat(paste0("[INFO] Loading count matrix from: ", count_matrix_file, "\n"))
counts <- read.table(count_matrix_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Ensure required coordinate columns exist
meta_cols <- c("chromosome", "start", "end", "name")
missing_meta <- setdiff(meta_cols, colnames(counts))
if (length(missing_meta) > 0) {
    stop(paste("Count matrix is missing required metadata columns:", paste(missing_meta, collapse = ", ")))
}

# Verify target sample exists in columns
if (!(sample_id %in% colnames(counts))) {
    stop(paste("Target sample '", sample_id, "' was not found in count matrix columns!"))
}

# Identify reference samples in the cohort
sample_cols <- setdiff(colnames(counts), meta_cols)
ref_sample_cols <- setdiff(sample_cols, sample_id)

if (length(ref_sample_cols) < 2) {
    stop(paste0("ExomeDepth requires at least 2 reference samples in the cohort. Found: ", length(ref_sample_cols)))
}

test_counts <- counts[[sample_id]]
ref_matrix  <- as.matrix(counts[, ref_sample_cols, drop = FALSE])

# Calculate exon bin lengths in kilobases for reference optimization
bin_lengths <- (counts$end - counts$start + 1) / 1000.0

cat(paste0("[INFO] Selecting optimal aggregate reference for sample: ", sample_id, " across ", length(ref_sample_cols), " potential reference exomes...\n"))

# Select reference set via ExomeDepth's step-wise forward model
ref_choice <- NULL
tryCatch({
    selected_ref <- select.reference.set(
        test.counts = test_counts,
        reference.counts = ref_matrix,
        bin.length = bin_lengths,
        n.bins.reduced = min(10000, nrow(counts))
    )
    ref_choice <- selected_ref$reference.choice
}, error = function(e) {
    cat(paste0("[WARNING] select.reference.set encountered an error: ", e$message, "\n"))
})

# Robust fallback if reference selection failed or returned empty set
if (is.null(ref_choice) || length(ref_choice) == 0) {
    cat("[WARNING] Fallback: Selecting top correlated reference samples using Pearson correlation.\n")
    corrs <- apply(ref_matrix, 2, function(x) cor(test_counts, x))
    top_n <- min(5, length(ref_sample_cols))
    ref_choice <- names(sort(corrs, decreasing = TRUE))[1:top_n]
}

cat(paste0("[INFO] Selected ", length(ref_choice), " reference samples: ", paste(ref_choice, collapse = ", "), "\n"))

# Sum counts across selected reference samples
ref_aggregated <- apply(
    X = ref_matrix[, ref_choice, drop = FALSE],
    MARGIN = 1,
    FUN = sum
)

# Instantiate ExomeDepth model object
cat("[INFO] Initializing ExomeDepth beta-binomial model...\n")
ed_obj <- new(
    'ExomeDepth',
    test = test_counts,
    reference = ref_aggregated,
    formula = 'cbind(test, reference) ~ 1'
)

# Call CNVs using Viterbi HMM decoding
cat("[INFO] Running CallCNVs across target intervals...\n")
cnvs <- CallCNVs(
    x = ed_obj,
    transition.probability = trans_prob,
    chromosome = counts$chromosome,
    start = counts$start,
    end = counts$end,
    name = counts$name
)

# Format and write results
if (!is.null(cnvs@CNV.calls) && nrow(cnvs@CNV.calls) > 0) {
    res <- cnvs@CNV.calls
    res$sample <- sample_id
    write.csv(res, file = output_csv, row.names = FALSE, quote = FALSE)
    cat(paste0("[INFO] Successfully called ", nrow(res), " CNVs for sample: ", sample_id, "\n"))
} else {
    empty_df <- data.frame(
        start.p = integer(),
        end.p = integer(),
        type = character(),
        nexons = integer(),
        start = integer(),
        end = integer(),
        chromosome = character(),
        id = character(),
        BF = numeric(),
        reads.expected = numeric(),
        reads.observed = numeric(),
        reads.ratio = numeric(),
        sample = character(),
        stringsAsFactors = FALSE
    )
    write.csv(empty_df, file = output_csv, row.names = FALSE, quote = FALSE)
    cat(paste0("[INFO] Zero CNVs identified for sample: ", sample_id, "\n"))
}
