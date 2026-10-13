# CNV & SV Caller Workflow Architecture & Specifications

This document provides exhaustive, publication-grade workflow diagrams, algorithmic foundations, container runtime mappings, and Nextflow DSL2 orchestration details for all Copy Number Variant (CNV) and Structural Variant (SV) callers integrated into the `CNV-Caller-Combiner` master ensemble pipeline.

---

## Table of Contents
1. [CANOES: Negative Binomial HMM Exome CNV Caller](#1-canoes-negative-binomial-hmm-exome-cnv-caller)
2. [CLAMMS: Mappability-Aware Mixture-Model Exome CNV Caller](#2-clamms-mappability-aware-mixture-model-exome-cnv-caller)
3. [XHMM: PCA SVD Matrix-Denoising Exome CNV Caller](#3-xhmm-pca-svd-matrix-denoising-exome-cnv-caller)
4. [CNVkit: On-Target & Antitarget Exome Segmentation Caller](#4-cnvkit-on-target--antitarget-exome-segmentation-caller)
5. [GATK-gCNV: Hierarchical Bayesian MCMC Exome CNV Caller](#5-gatk-gcnv-hierarchical-bayesian-mcmc-exome-cnv-caller)
6. [ExomeDepth: Beta-Binomial Exome CNV Caller](#6-exomedepth-beta-binomial-exome-cnv-caller)
7. [SvABA: Local De Novo Assembly Structural Variant Caller](#7-svaba-local-de-novo-assembly-structural-variant-caller)
8. [DRAGEN Germline Enrichment v4.4.6: Combined CNV & SV Caller](#8-dragen-germline-enrichment-v446-combined-cnv--sv-caller)
9. [Comparative Architectural & Statistical Summary Matrix](#9-comparative-architectural--statistical-summary-matrix)

---

## 1. CANOES: Negative Binomial HMM Exome CNV Caller

- **Primary Module**: [`modules/callers/canoes.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/canoes.nf)
- **Operational Guide**: [CANOES Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/canoes/SKILL.md)
- **Publication**: Backenroth et al., *Nucleic Acids Res.* (2014) [PMID: 24771342](https://pubmed.ncbi.nlm.nih.gov/24771342/)
- **Statistical Model**: Negative Binomial Generalized Linear Model with GC, $\text{GC}^2$, and target length covariate regression; Viterbi HMM decoding.

```mermaid
flowchart TD
    subgraph Inputs ["Input Artifacts & References"]
        PROBES["probes.bed (Exome Capture Intervals)"]
        FASTA["reference.fasta + .fai"]
        BAM_LIST["cohort_bams.list (Cohort BAM Paths)"]
        SAMPLE_LIST["sample_list.txt"]
    end

    subgraph GC_Prep ["Target Splitting & GC Modeling"]
        SPLIT_PROBES["Split Targets by Chromosome (chr1..chr22)"]
        GC_CALC["CALC_GC_CANOES\n(GATK AnnotateIntervals -imr OVERLAPPING_ONLY)"]
        GC_FILES["${chr}_gc.txt (GC Percentage per Target)"]
        
        PROBES --> SPLIT_PROBES
        SPLIT_PROBES --> GC_CALC
        FASTA --> GC_CALC
        GC_CALC --> GC_FILES
    end

    subgraph Batch_Counting ["BAM Batching & Depth Extraction"]
        BATCH["BATCH_BAMS\n(split -l canoes_batch_size [default: 100])"]
        MULTICOV["GEN_READ_COUNTS\n(bedtools multicov -q 20 per BAM Chunk)"]
        MERGE_COUNTS["MERGE_READ_COUNTS\n(Deterministic Horizontal Paste & Cut)"]
        READ_MAT["${chr}_canoes_reads.txt (Chromosome Count Matrix)"]
        
        BAM_LIST --> BATCH
        BATCH --> MULTICOV
        SPLIT_PROBES --> MULTICOV
        MULTICOV --> MERGE_COUNTS
        MERGE_COUNTS --> READ_MAT
    end

    subgraph HMM_Core ["CANOES Statistical & HMM Engine"]
        JOIN["Join Target GC Covariates + Read Counts per Chr"]
        RUN_CANOES["RUN_CANOES\n(run_canoes.R in phelelani/nf-exomecnv:canoes)\n- Strip 'chr' Prefixes\n- Fit NB Regression: log(mu) ~ alpha + beta*GC + beta*GC^2 + beta*log(Length)\n- Viterbi HMM State Path Decoding\n- Compute Phred Q_SOME"]
        PASS_CSV["${chr}_CNVs_pass.csv"]
        FAIL_CSV["${chr}_CNVs_fail.csv"]
        RDATA["${chr}_xcnvs.RData & CNV Plots"]
        COMBINE["COMBINE_CANOES_CNVS\n(Concatenate Autosomal Passes)"]
        SAMPLE_CSV["Sample_CNVs.csv"]
        
        GC_FILES --> JOIN
        READ_MAT --> JOIN
        SAMPLE_LIST --> RUN_CANOES
        JOIN --> RUN_CANOES
        RUN_CANOES --> PASS_CSV
        RUN_CANOES --> FAIL_CSV
        RUN_CANOES --> RDATA
        PASS_CSV --> COMBINE
        COMBINE --> SAMPLE_CSV
    end

    subgraph VCF_Post ["Standardization & Compression"]
        CONV["CONVERT_CANOES_TO_VCF\n(bin/canoes_csv_to_vcf.py via pysam)"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix -p vcf)\n- Set ##source=CANOES"]
        CANOES_VCF["${sample}_CANOES_output.vcf.gz + .tbi"]
        
        SAMPLE_CSV --> CONV
        FASTA --> CONV
        SAMPLE_LIST --> CONV
        CONV --> BGZIP
        BGZIP --> CANOES_VCF
    end

    subgraph Downstream ["Downstream Ensemble Integration"]
        NORM["NORM_CANOES_FULL\n(Harmonize SVTYPE, END, CIPOS, Q_SOME)"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST["CatBoost ML Matrix Extraction\n(Cohort Parquet Feature Database)"]
        
        CANOES_VCF --> NORM
        NORM --> CONSENSUS
        NORM --> CATBOOST
    end
```

---

## 2. CLAMMS: Mappability-Aware Mixture-Model Exome CNV Caller

- **Primary Module**: [`modules/callers/clamms.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/clamms.nf)
- **Operational Guide**: [CLAMMS Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/clamms/SKILL.md)
- **Publication**: Packer et al., *Bioinformatics* (2016) [PMID: 26382196](https://pubmed.ncbi.nlm.nih.gov/26382196/)
- **Statistical Model**: Dynamic sample-specific reference panel selection (PCA + Picard HS QC distance); Gaussian mixture models on mappability-adjusted coverage windows; sex-aware allosomal genotyping.

```mermaid
flowchart TD
    subgraph Inputs ["Input Artifacts & Metadata"]
        PROBES["probes.bed (Target Intervals)"]
        REF["reference.fasta"]
        MAP["mappability.bed (Empirical k-mer Mappability)"]
        SPEC["special_regions.bed (Non-PAR Pseudoautosomal Mask)"]
        BAMS["Cohort BAMs / BAIs"]
        SEX["sexinfo.txt (Karyotypic Sex Metadata)"]
    end

    subgraph Window_Prep ["Target Window Annotation"]
        GEN_WIN["GENERATE_WINDOWS\n(annotate_windows.sh)\n- Window Decomposition\n- Insert Size Calibration (200bp)\n- Mappability & Special Regions Annotation"]
        WINDOWS["windows.bed"]
        
        PROBES --> GEN_WIN
        REF --> GEN_WIN
        MAP --> GEN_WIN
        SPEC --> GEN_WIN
        GEN_WIN --> WINDOWS
    end

    subgraph Depth_Extraction ["Depth Collection & Normalization"]
        SAMTOOLS_DOC["SAMTOOLS_DOC\n(samtools bedcov -Q 30 across Windows)"]
        NORMALIZE["NORMALIZE_DOC\n(normalize_coverage)\nNormCov = RawDepth / (WindowLength * SampleMedianCoverage)"]
        NORM_COV["${sample}.norm.cov.bed"]
        
        WINDOWS --> SAMTOOLS_DOC
        BAMS --> SAMTOOLS_DOC
        SAMTOOLS_DOC --> NORMALIZE
        NORMALIZE --> NORM_COV
    end

    subgraph Picard_Covariates ["Picard Hybrid Selection QC"]
        PICARD_HS["GET_PICARD_QC_METRICS\n(Picard CollectHsMetrics:\nON/OFF Target, Fold-80 Penalty, AT/GC Dropout)"]
        PICARD_INS["GET_PICARD_MEAN_INSERT_SIZE\n(Picard CollectInsertSizeMetrics)"]
        COMB_QC["COMBINE_PICARD_QC_METRICS"]
        QC_TABLE["cohort_picard_qc.txt"]
        
        BAMS --> PICARD_HS
        BAMS --> PICARD_INS
        PICARD_HS --> COMB_QC
        PICARD_INS --> COMB_QC
        COMB_QC --> QC_TABLE
    end

    subgraph Ref_Panel ["Dynamic Reference Panel Selection"]
        PCA_CALC["CREATE_PCA_DATA\n(custom_ref_panel.sh - PCA on NormCov Matrix)"]
        REF_OPT["CREATE_CUSTOM_REF_PANEL\n(custom_ref_panel.R)\n- Minimize Covariate & PCA Distance\n- Align K Control Samples Matching Sex"]
        SAMPLE_PANELS["${sample}.ref_panel.txt"]
        
        NORM_COV --> PCA_CALC
        PCA_CALC --> REF_OPT
        QC_TABLE --> REF_OPT
        SEX --> REF_OPT
        REF_OPT --> SAMPLE_PANELS
    end

    subgraph Model_Calling ["Mixture Modeling & Genotyping"]
        TRAIN["TRAIN_MODELS\n(fit_models on Reference Windows)\nP(NormCov | CN=k) ~ Normal(mu_ik, sigma_ik^2)"]
        CALL["CALL_CNVS\n(call_cnv --sex ${sample_sex})\n- Assign Integer Copy States\n- Calculate Q_SOME, Q_EXACT_DEL, Q_EXACT_DUP"]
        COMBINE_BED["COMBINE_CLAMMS_CNVS\n(samples.cnv.bed)"]
        
        WINDOWS --> TRAIN
        SAMPLE_PANELS --> TRAIN
        NORM_COV --> TRAIN
        TRAIN --> CALL
        NORM_COV --> CALL
        SEX --> CALL
        CALL --> COMBINE_BED
    end

    subgraph VCF_Post ["Standardization & Downstream Routing"]
        CONV["CONVERT_CLAMMS_TO_VCF\n(bin/clamms_bed_to_vcf.py)"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix)\n- Set ##source=CLAMMS"]
        CLAMMS_VCF["${sample}_CLAMMS_output.vcf.gz + .tbi"]
        NORM_FULL["NORM_CLAMMS_FULL"]
        REFINEMENT["Breakpoint Refinement Engine\n(bin/refine_breakpoints.py)"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        
        COMBINE_BED --> CONV
        REF --> CONV
        CONV --> BGZIP
        BGZIP --> CLAMMS_VCF
        CLAMMS_VCF --> NORM_FULL
        NORM_FULL --> REFINEMENT
        NORM_FULL --> CONSENSUS
    end
```

---

## 3. XHMM: PCA SVD Matrix-Denoising Exome CNV Caller

- **Primary Module**: [`modules/callers/xhmm.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/xhmm.nf)
- **Operational Guide**: [XHMM Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/xhmm/SKILL.md)
- **Publication**: Fromer et al., *Am. J. Hum. Genet.* (2012) [PMID: 23040492](https://pubmed.ncbi.nlm.nih.gov/23040492/)
- **Statistical Model**: Singular Value Decomposition (SVD) Principal Component Analysis to subtract batch/plate noise; distance-weighted Viterbi HMM decoding on target z-scores.

```mermaid
flowchart TD
    subgraph Inputs ["Input Artifacts & Targets"]
        PROBES["probes.bed (Target Intervals)"]
        FASTA["reference.fasta"]
        BAMS["Cohort BAMs / BAIs"]
    end

    subgraph DOC_Collection ["Depth Collection & Matrix Aggregation"]
        GATK_DOC["GATK_DOC\n(GATK DepthOfCoverage across Targets)"]
        DOC_SUMMARIES["${sample}.sample_interval_summary"]
        COMBINE_DOC["COMBINE_DOC\n(xhmm --mergeGATKdepths)"]
        RAW_MATRIX["DATA.RD.txt (Raw Depth Matrix)"]
        
        PROBES --> GATK_DOC
        FASTA --> GATK_DOC
        BAMS --> GATK_DOC
        GATK_DOC --> DOC_SUMMARIES
        DOC_SUMMARIES --> COMBINE_DOC
        COMBINE_DOC --> RAW_MATRIX
    end

    subgraph GC_Target_Filtering ["Target QC & Extreme GC Elimination"]
        CALC_GC["CALC_GC_XHMM\n(GATK AnnotateIntervals)"]
        EXTREME_GC["extreme_gc_targets.txt (GC < 0.1 or > 0.9)"]
        FILTER_RD1["FILTER_SAMPLES\n(xhmm --matrix -r DATA.RD.txt)\n- Exclude Extreme GC Targets\n- Filter Low Target Depth (<10x)\n- Filter Low Sample Mean Depth (<25x)"]
        CENTERED_MAT["DATA.filtered_centered.RD.txt"]
        
        PROBES --> CALC_GC
        FASTA --> CALC_GC
        CALC_GC --> EXTREME_GC
        RAW_MATRIX --> FILTER_RD1
        EXTREME_GC --> FILTER_RD1
        FILTER_RD1 --> CENTERED_MAT
    end

    subgraph SVD_PCA ["SVD Denoising & Normalization"]
        PCA_SVD["RUN_PCA\n(xhmm --pca)\nDecompose M = U * Sigma * V^T"]
        PCA_MODEL["DATA.RD_PCA"]
        NORMALIZE["NORMALIZE\n(xhmm --normalize --pcaModel)\n- Subtract Leading Technical Components (>70% Var)\n- Compute Standardized Z-scores: Z_ij = (M_hat - mu_j) / sigma_j"]
        NORM_Z["DATA.pca_normalized.txt"]
        FILTER_RD2["FILTER_SAMPLES_ROUND2\n(xhmm --matrix)\n- Filter Targets with Extreme Residual Variance (SD > 30)"]
        FINAL_Z["DATA.same_filtered.pca_normalized.txt"]
        
        CENTERED_MAT --> PCA_SVD
        PCA_SVD --> PCA_MODEL
        CENTERED_MAT --> NORMALIZE
        PCA_MODEL --> NORMALIZE
        NORMALIZE --> NORM_Z
        NORM_Z --> FILTER_RD2
        FILTER_RD2 --> FINAL_Z
    end

    subgraph HMM_Calling ["Viterbi HMM Genotyping & VCF Conversion"]
        RUN_XHMM["RUN_XHMM\n(xhmm --discover & --genotype)\n- Distance-Penalized State Transition Matrix\n- Emit Phred Scores: SQ (Q_SOME), NQ/NDQ, EQ"]
        XCNV["DATA.xcnv"]
        CONV["CONVERT_XHMM_TO_VCF\n(bin/xhmm_xcnv_to_vcf.py)"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix)\n- Set ##source=XHMM"]
        XHMM_VCF["${sample}_XHMM_output.vcf.gz + .tbi"]
        
        FINAL_Z --> RUN_XHMM
        RUN_XHMM --> XCNV
        XCNV --> CONV
        FASTA --> CONV
        CONV --> BGZIP
        BGZIP --> XHMM_VCF
    end

    subgraph Downstream ["Ensemble Consensus & Machine Learning"]
        NORM_FULL["NORM_XHMM_FULL"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST["CatBoost Prioritization Database\n(SVD-Cleaned Z-Scores & Phred Quality)"]
        
        XHMM_VCF --> NORM_FULL
        NORM_FULL --> CONSENSUS
        NORM_FULL --> CATBOOST
    end
```

---

## 4. CNVkit: On-Target & Antitarget Exome Segmentation Caller

- **Primary Module**: [`modules/callers/cnvkit.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/cnvkit.nf)
- **Operational Guide**: [CNVkit Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/cnvkit/SKILL.md)
- **Publication**: Talevich et al., *PLoS Comput. Biol.* (2016) [PMID: 27100738](https://pubmed.ncbi.nlm.nih.gov/27100738/)
- **Statistical Model**: On-target and off-target (antitarget) binning; GC, length, and repeat bias regression; Circular Binary Segmentation (CBS); discrete copy number mapping.

```mermaid
flowchart TD
    subgraph Inputs ["Input References & Target Intervals"]
        FASTA["reference.fasta"]
        TARGETS["targets.bed (Exome Capture Targets)"]
        REFFLAT["refFlat.txt (Gene Annotations)"]
        BAMS["Sample & Cohort BAMs"]
    end

    subgraph Binning_Prep ["Antitarget Partitioning & Autobinning"]
        ACCESS["GENERATE_ACCESS\n(cnvkit.py access)\n- Exclude Telomeres, Centromeres, Uncalled N Bases"]
        ACCESS_BED["access.hg38.bed"]
        AUTOBIN["AUTOBIN\n(cnvkit.py autobin)\n- Dynamically Size On-Target Bins\n- Size Antitarget Bins (10-100kb Intervening Regions)"]
        TARGET_BED["*.target.bed"]
        ANTITARGET_BED["*.antitarget.bed"]
        
        FASTA --> ACCESS
        ACCESS --> ACCESS_BED
        ACCESS_BED --> AUTOBIN
        TARGETS --> AUTOBIN
        REFFLAT --> AUTOBIN
        BAMS --> AUTOBIN
        AUTOBIN --> TARGET_BED
        AUTOBIN --> ANTITARGET_BED
    end

    subgraph Depth_Collection ["Target & Antitarget Coverage Counting"]
        COV_SAMPLE["COVERAGE (Sample BAM)\n(cnvkit.py coverage -p 2)\nExtract .targetcoverage.cnn & .antitargetcoverage.cnn"]
        COV_POOLS["COVERAGE (Cohort Control BAMs)\n(cnvkit.py coverage -p 2)"]
        
        BAMS --> COV_SAMPLE
        BAMS --> COV_POOLS
        TARGET_BED --> COV_SAMPLE
        ANTITARGET_BED --> COV_SAMPLE
        TARGET_BED --> COV_POOLS
        ANTITARGET_BED --> COV_POOLS
    end

    subgraph Pooled_Ref ["Pooled Reference Construction"]
        MAKE_REF["CREATE_POOLED_REFERENCE\n(cnvkit.py reference)\n- GC Bias Correction (Local Polynomial Regression)\n- Target Density & Bin Size Adjustments\n- Repeat Mappability Normalization"]
        POOLED_REF["pooled_reference.cnn"]
        
        COV_POOLS --> MAKE_REF
        FASTA --> MAKE_REF
        MAKE_REF --> POOLED_REF
    end

    subgraph Fix_Segment_Call ["Correction, CBS Segmentation & Calling"]
        FIX["cnvkit.py fix\nLog2Ratio = log2(SampleDepth / RefDepth) - f(GC) - g(Size)"]
        CNR["${sample}.cnr (Normalized Bin Ratios)"]
        SEGMENT["cnvkit.py segment\nCircular Binary Segmentation (CBS) Change-Point Detection"]
        CNS["${sample}.cns (Segmented Intervals)"]
        CALL_INT["cnvkit.py call\nDiscrete Integer Copy State Mapping (0, 1, 2, 3...)"]
        VISUALS["cnvkit.py scatter & diagram\n(.scatter.pdf & .diagram.pdf)"]
        
        COV_SAMPLE --> FIX
        POOLED_REF --> FIX
        FIX --> CNR
        CNR --> SEGMENT
        SEGMENT --> CNS
        CNS --> CALL_INT
        CNR --> VISUALS
        CNS --> VISUALS
    end

    subgraph VCF_Post ["VCF Export & Ensemble Integration"]
        EXPORT["EXPORT_RESULTS\n(cnvkit.py export vcf)"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix)\n- Set ##source=CNVKIT"]
        CNVKIT_VCF["${sample}_CNVKIT_output.vcf.gz + .tbi"]
        NORM_FULL["NORM_CNVKIT_FULL\n(Map PROBES, FOLD_CHANGE, LOG2R)"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST["CatBoost Prioritization Database\n(Segment Probe Counts, Spread, Weight)"]
        
        CALL_INT --> EXPORT
        EXPORT --> BGZIP
        BGZIP --> CNVKIT_VCF
        CNVKIT_VCF --> NORM_FULL
        NORM_FULL --> CONSENSUS
        NORM_FULL --> CATBOOST
    end
```

---

## 5. GATK-gCNV: Hierarchical Bayesian MCMC Exome CNV Caller

- **Primary Module**: [`modules/callers/gatk-gcnv.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/gatk-gcnv.nf)
- **Operational Guide**: [GATK-gCNV Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/gatk-gcnv/SKILL.md)
- **Publication**: Babadi et al., *Nat. Genet.* (2023) [PMID: 37723223](https://pubmed.ncbi.nlm.nih.gov/37723223/)
- **Statistical Model**: Hierarchical Bayesian probabilistic graphical model; MCMC and variational inference modeling global depth, interval capture efficiency, latent bias projections, and copy ratio; contig ploidy priors; interval sharding.

```mermaid
flowchart TD
    subgraph Inputs ["Input Genomic References & Annotations"]
        FASTA["reference.fasta + .fai + .dict"]
        PROBES["targets.bed (Capture Targets)"]
        BAMS["Cohort BAMs / BAIs"]
        SEX["sexinfo.txt (Verified Karyotypic Sex)"]
    end

    subgraph Interval_Prep ["Interval Preprocessing & Priors"]
        PRIORS["GENERATE_PLOIDY_PRIORS\n(Autosome Diploid 0.97 Prior + Allosome Priors)"]
        PLOIDY_TSV["ploidy_priors.tsv"]
        PREP_INT["PREPROCESS_INTERVALS\n(GATK PreprocessIntervals --bin-length 0 --padding 250)"]
        PADDED_LIST["preprocessed.interval_list"]
        ANNOT_INT["ANNOTATE_INTERVALS\n(GATK AnnotateIntervals: GC Content & Complexity)"]
        ANNOT_TSV["annotated_intervals.tsv"]
        
        FASTA --> PRIORS
        PRIORS --> PLOIDY_TSV
        PROBES --> PREP_INT
        FASTA --> PREP_INT
        PREP_INT --> PADDED_LIST
        PADDED_LIST --> ANNOT_INT
        FASTA --> ANNOT_INT
        ANNOT_INT --> ANNOT_TSV
    end

    subgraph Count_Filter ["Read Count Collection & Target Filtering"]
        COLLECT["COLLECT_READ_COUNTS\n(GATK CollectReadCounts per Sample)"]
        COUNTS_HDF5["${sample}.counts.hdf5"]
        FILTER_INT["FILTER_INTERVALS\n(GATK FilterIntervals: Strip Low/Extreme Depth Intervals)"]
        FILTERED_LIST["cohort_filtered.interval_list"]
        
        PADDED_LIST --> COLLECT
        BAMS --> COLLECT
        COLLECT --> COUNTS_HDF5
        COUNTS_HDF5 --> FILTER_INT
        ANNOT_TSV --> FILTER_INT
        FILTER_INT --> FILTERED_LIST
    end

    subgraph Ploidy_Sex ["Contig Ploidy Baseline & Karyotype Alignment"]
        PLOIDY_DET["DETERMINE_PLOIDY_COHORT\n(GATK DetermineGermlineContigPloidy)"]
        CALL_SEX["ENFORCE_SEX_METADATA\n- Cross-Reference with sexinfo.txt\n- Adjust chrX / chrY Ploidy Baselines"]
        PLOIDY_CALLS["cohort-contig-ploidy-calls/"]
        
        FILTERED_LIST --> PLOIDY_DET
        COUNTS_HDF5 --> PLOIDY_DET
        PLOIDY_TSV --> PLOIDY_DET
        PLOIDY_DET --> CALL_SEX
        SEX --> CALL_SEX
        CALL_SEX --> PLOIDY_CALLS
    end

    subgraph Sharded_MCMC ["Interval Sharding & Bayesian MCMC Engine"]
        SCATTER["SCATTER_INTERVALS\n(IntervalListTools --SCATTER_COUNT 5000)"]
        SHARDS["shard_0000_list/ ... shard_NNNN_list/"]
        GCNV_MCMC["GERMLINE_CNV_CALLER_COHORT (per Shard)\n(GATK GermlineCNVCaller MCMC / Variational Inference)\n- Model: log(lambda) = log(mu_s) + log(f_t) + w_t^T * z_s + log(c_st)\n- Parameters: --cnv-coherence-length 1000.0\n--interval-psi-scale 1.0E-6 --sample-psi-scale 1.0E-6"]
        SHARD_CALLS["shard_calls/ per Cohort Shard"]
        
        FILTERED_LIST --> SCATTER
        SCATTER --> SHARDS
        SHARDS --> GCNV_MCMC
        COUNTS_HDF5 --> GCNV_MCMC
        ANNOT_TSV --> GCNV_MCMC
        PLOIDY_CALLS --> GCNV_MCMC
        GCNV_MCMC --> SHARD_CALLS
    end

    subgraph Post_Processing ["Shard Reassembly & Genotyping"]
        POSTPROCESS["POSTPROCESS_CALLS\n(GATK PostprocessGermlineCNVCalls)\n- Reassemble Sharded Intervals\n- Segment Genotyping & Boundary Quality (QSS, QSE)\n- Compute Phred QS Score & Continuous Copy Ratio (CR)"]
        DENOISED_RATIOS["${sample}.denoised_copy_ratios.tsv"]
        SEGMENTS_VCF["${sample}.genotyped_segments.vcf.gz"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix)\n- Set ##source=GCNV"]
        FINAL_VCF["${sample}.genotyped_segments.vcf.gz + .tbi"]
        
        SHARD_CALLS --> POSTPROCESS
        PLOIDY_CALLS --> POSTPROCESS
        POSTPROCESS --> DENOISED_RATIOS
        POSTPROCESS --> SEGMENTS_VCF
        SEGMENTS_VCF --> BGZIP
        BGZIP --> FINAL_VCF
    end

    subgraph Downstream ["Consensus & Machine Learning Prioritization"]
        NORM_FULL["NORM_GCNV_FULL\n(Map NP, QS, CR Annotations)"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST["CatBoost Prioritization Database\n(Quality Scores, Posterior Probs, Copy Ratios)"]
        
        FINAL_VCF --> NORM_FULL
        NORM_FULL --> CONSENSUS
        NORM_FULL --> CATBOOST
    end
```

---

## 6. ExomeDepth: Beta-Binomial Exome CNV Caller

- **Primary Module**: [`modules/callers/exomedepth.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/exomedepth.nf)
- **Operational Guide**: [ExomeDepth Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/exomedepth/SKILL.md)
- **Publication**: Plagnol et al., *Bioinformatics* (2012) [PMID: 22942019](https://pubmed.ncbi.nlm.nih.gov/22942019/)
- **Statistical Model**: Dynamic sample-specific reference panel selection (correlation $\ge 0.97$); over-dispersed Beta-Binomial read count distribution; 3-state HMM with transition penalty; Bayes Factor (BF) variant quantification.

```mermaid
flowchart TD
    subgraph Inputs ["Input Artifacts & Cohort Alignments"]
        PROBES["probes.bed (Exome Capture Targets)"]
        BAM_LIST["cohort_bams.list (Cohort BAM Filepaths)"]
        FASTA_FAI["reference.fasta.fai"]
        SAMPLE_LIST["sample_list.txt"]
    end

    subgraph Target_Prep ["Interval Standardization & Batching"]
        PREP_BED["PREPARE_EXOMEDEPTH_TARGETS\n(Standardize to 4-Column BED: chr, start, end, name)"]
        BED_4COL["exomedepth_targets.bed"]
        BATCH["BATCH_BAMS_EXOMEDEPTH\n(split -l exomedepth_batch_size [default: 100])"]
        BATCH_CHUNKS["batch_chunks_*.txt"]
        
        PROBES --> PREP_BED
        PREP_BED --> BED_4COL
        BAM_LIST --> BATCH
        BATCH --> BATCH_CHUNKS
    end

    subgraph Depth_Extraction ["Multicov Batch Counting & Matrix Assembly"]
        MULTICOV["GEN_EXOMEDEPTH_COUNTS\n(bedtools multicov -q 20 per BAM Chunk)"]
        CHUNK_COUNTS["batch_*_counts.txt"]
        MERGE_COUNTS["MERGE_EXOMEDEPTH_COUNTS\n(Assemble Unified Cohort Count Matrix)"]
        COUNT_MATRIX["exomedepth_cohort_counts.txt"]
        
        BED_4COL --> MULTICOV
        BATCH_CHUNKS --> MULTICOV
        MULTICOV --> CHUNK_COUNTS
        CHUNK_COUNTS --> MERGE_COUNTS
        MERGE_COUNTS --> COUNT_MATRIX
    end

    subgraph Statistical_Engine ["Dynamic Reference Selection & Beta-Binomial HMM"]
        RUN_ED["RUN_EXOMEDEPTH (Parallelized per Sample)\n(run_exomedepth.R in r-exomedepth:1.1.16)\n- Dynamic Reference Selection: Forward Selection (Correlation >= 0.97)\n- Over-dispersed Model: Y_i ~ BetaBinomial(N_i, p_i, Phi)\n- 3-State HMM (Del, Normal, Dup) with Transition Prob 1e-4\n- Bayes Factor: BF = log10( P(Counts | CNV) / P(Counts | Diploid) )"]
        SAMPLE_CSV["${sample}_CNVs.csv (Raw Call Records)"]
        
        COUNT_MATRIX --> RUN_ED
        SAMPLE_LIST --> RUN_ED
        RUN_ED --> SAMPLE_CSV
    end

    subgraph VCF_Post ["Standardization & Compression"]
        CONV["CONVERT_EXOMEDEPTH_TO_VCF\n(bin/exomedepth_csv_to_vcf.py)\n- Filter by exomedepth_bf_min\n- Format to Standard VCF"]
        BGZIP["BGZIP_SORT_INDEX_VCF\n(bcftools sort | bgzip | tabix)\n- Set ##source=EXOMEDEPTH"]
        ED_VCF["${sample}_EXOMEDEPTH_output.vcf.gz + .tbi"]
        
        SAMPLE_CSV --> CONV
        FASTA_FAI --> CONV
        CONV --> BGZIP
        BGZIP --> ED_VCF
    end

    subgraph Downstream ["Consensus Powerhouse & Machine Learning"]
        NORM_FULL["NORM_EXOMEDEPTH_FULL\n(Preserve BF, READS_RATIO, EXPECTED_COUNTS)"]
        CONSENSUS["Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST["CatBoost Prioritization Database\n(Single-Exon & Multi-Exon High-Sensitivity Matrix)"]
        
        ED_VCF --> NORM_FULL
        NORM_FULL --> CONSENSUS
        NORM_FULL --> CATBOOST
    end
```

---

## 7. SvABA: Local De Novo Assembly Structural Variant Caller

- **Primary Module**: [`modules/callers/svaba.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/svaba.nf)
- **Operational Guide**: [SvABA Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/svaba/SKILL.md)
- **Publication**: Wala et al., *Genome Res.* (2018) [PMID: 29535149](https://pubmed.ncbi.nlm.nih.gov/29535149/)
- **Algorithmic Model**: Target-restricted candidate read extraction (split reads, discordant pairs); localized de novo assembly via bidirected string graph (FM-index/BWT); contig realignment via BWA-MEM; base-pair resolution breakpoint calling.

```mermaid
flowchart TD
    subgraph Inputs ["Input Alignments & Target Definitions"]
        ALIGN_IN["Input Alignment: Sample BAM or CRAM + Index"]
        FASTA["reference.fasta"]
        TARGETS["target_bed (Exome Capture Targets)"]
    end

    subgraph Preflight ["BWA Index Pre-Flight & Alignment Slicing"]
        BWA_CHECK{"Check Reference BWA Index\n(${fasta}.bwt exists?)"}
        BUILD_BWA["bwa index ${fasta}\n(Generate Burrows-Wheeler Transform Index)"]
        BWA_INDEX["Reference BWA Indices (.bwt, .pac, .ann, .amb, .sa)"]
        
        IS_CRAM{"Is Input Format CRAM?"}
        SLICE_BAM["Target Slicing & Conversion\n(samtools view -b -T ${fasta} -L ${target_bed})\n- Prevent Decompression Bottlenecks\n- Generate Target-Padded BAM"]
        DIRECT_BAM["Direct BAM Path Pass-Through"]
        FINAL_ALIGN["Validated Local Alignment Stream"]
        
        FASTA --> BWA_CHECK
        BWA_CHECK -- No --> BUILD_BWA
        BWA_CHECK -- Yes --> BWA_INDEX
        BUILD_BWA --> BWA_INDEX
        
        ALIGN_IN --> IS_CRAM
        IS_CRAM -- Yes --> SLICE_BAM
        IS_CRAM -- No --> DIRECT_BAM
        TARGETS --> SLICE_BAM
        FASTA --> SLICE_BAM
        SLICE_BAM --> FINAL_ALIGN
        DIRECT_BAM --> FINAL_ALIGN
    end

    subgraph Assembly_Engine ["SvABA Local Assembly & Contig Realignment"]
        RUN_SVABA["RUN_SVABA_ASSEMBLY\n(svaba run -t ${ALIGN} -G ${fasta} -k ${target_bed} --germline -a ${sample_id} -p 4)\n1. Candidate Extraction: Discordant Pairs, Soft-clipped/Split Reads, Unmapped Mates\n2. String-Graph Assembly: Bidirected String Graph from FM-index\n3. Contig Realignment: Realign Assembled Contigs via Embedded BWA-MEM\n4. Breakpoint Calling: Single-Base Microhomology & Non-Templated Insertions"]
        RAW_SV["${sample}.svaba.raw.sv.vcf"]
        INDEL_VCF["${sample}.svaba.indel.vcf"]
        
        FINAL_ALIGN --> RUN_SVABA
        BWA_INDEX --> RUN_SVABA
        TARGETS --> RUN_SVABA
        RUN_SVABA --> RAW_SV
        RUN_SVABA --> INDEL_VCF
    end

    subgraph Standardization ["VCF Harmonization & Indexing"]
        BGZIP["BGZIP_SORT_INDEX_SV_TUPLE_VCF\n(bcftools sort | bgzip | tabix -p vcf)\n- Set ##source=SvABA"]
        SVABA_VCF["${sample}.svaba.sv.vcf.gz + .tbi"]
        
        RAW_SV --> BGZIP
        BGZIP --> SVABA_VCF
    end

    subgraph Downstream ["Consensus & Nucleotide-Resolution Refinement"]
        REFINEMENT["Breakpoint Refinement Engine\n(bin/refine_breakpoints.py)\n- Reconcile Exon-Bounded Depth Consensus with Assembly Breakends\n- Resolve Exact Single-Base Junctions"]
        CONSENSUS["Multi-Caller SV Consensus Engine\n(SURVIVOR & Truvari)"]
        
        SVABA_VCF --> REFINEMENT
        SVABA_VCF --> CONSENSUS
    end
```

---

## 8. DRAGEN Germline Enrichment v4.4.6: Combined CNV & SV Caller

- **Primary Module**: [`modules/callers/modules-icav2-dragen.nf`](file:///home/ywolberg/CNV-Caller-Combiner/modules/callers/modules-icav2-dragen.nf)
- **Operational Guide**: [DRAGEN Germline Enrichment v4.4.6 Skill](file:///home/ywolberg/CNV-Caller-Combiner/.agents/skills/dragen-germline-enrichment-v4-4-6/SKILL.md)
- **Foundational Publication**: Behera et al., *Nature Biotechnology* (2025) [DOI: 10.1038/s41587-024-02382-1](https://doi.org/10.1038/s41587-024-02382-1)
- **Architecture**: Unified hardware-accelerated FPGA pipeline integrating multigenome graph mapping, in-run Panel of Normals (PoN) read-depth CNV calling, and Manta split-read/string-graph SV assembly.

```mermaid
flowchart TD
    subgraph Inputs ["Input Artifacts & Cloud Payload Configuration"]
        READS["Input Reads: FASTQ / ORA or Alignments: BAM / CRAM"]
        REF_TAR["ref_tar (Multigenome Graph Reference + 64 Haplotypes + FPGA Hash Tables)"]
        TARGET_BED["target_bed (Capture Targets) & cnv_target_bed"]
        NORMAL_SAMPLES["Normal Cohort Controls (>= 5 Samples for In-Run PoN)"]
        PON_FILE["Optional Precomputed PoN: cnv_normals_file / cnv_combined_counts"]
        CONFIG_PARAMS["Pipeline Parameters:\nenable_map_align=true, enable_cnv=true,\ncnv_enable_in_run_pon=true, cnv_enable_gcbias_correction=true,\nenable_sv=true, output_format=CRAMv3.1"]
    end

    subgraph Orchestration ["ICAv2 Cloud Wrapper & FPGA Orchestration"]
        WRAPPER["CNV-Caller-Combiner Wrapper\n(modules-icav2-dragen.nf)\n- Package CRAMs/BAMs & Target BED into ICAv2 Payloads\n- Launch via icav2 projectpipelines start"]
        FPGA_INSTANCE["DRAGEN FPGA Instance Execution\n(Cloud dedicated fpga2-medium / fpga-medium)\n- Hardware-Accelerated Smith-Waterman & Graph Alignment"]
        
        READS --> WRAPPER
        TARGET_BED --> WRAPPER
        WRAPPER --> FPGA_INSTANCE
        REF_TAR --> FPGA_INSTANCE
        CONFIG_PARAMS --> FPGA_INSTANCE
    end

    subgraph Shared_Front ["Hardware Graph Alignment & QC"]
        GRAPH_ALIGN["Multigenome Graph Alignment\n(Eliminates Reference Bias in Low-Mappability Regions)"]
        DUP_MARK["Duplicate Marking & Coordinate Sorting"]
        BAM_CRAM_OUT["${sample}.cram + .crai (CRAM v3.1)"]
        
        FPGA_INSTANCE --> GRAPH_ALIGN
        GRAPH_ALIGN --> DUP_MARK
        DUP_MARK --> BAM_CRAM_OUT
    end

    subgraph CNV_Branch ["DRAGEN Read-Depth CNV Calling Engine"]
        subgraph PoN_Sub ["PoN Calibration (cnv-baseline-ica.nf)"]
            COUNT_EXTRACT["Target Read Binning\n(Collect *.target.counts.gz)"]
            GC_CORRECT["Systematic Bias Correction\n- Quadratic GC Polynomial Regression\n- Target Length & Flank Hybridization Adjustment\n- Mappability Weighting"]
            COMBINE_PON["DRAGEN_COMBINE_COUNTS\n(Multi-Sample Median Normalization across >= 5 Controls)"]
            PON_OUT["cohort_pon.combined.counts.txt.gz"]
            
            NORMAL_SAMPLES --> COUNT_EXTRACT
            COUNT_EXTRACT --> GC_CORRECT
            GC_CORRECT --> COMBINE_PON
            PON_FILE -.-> COMBINE_PON
            COMBINE_PON --> PON_OUT
        end

        subgraph Sample_CNV ["Sample Segmentation & Genotyping"]
            SAMPLE_COUNTS["Sample Target Count Extraction\n(*.target.counts.gc-corrected.gz)"]
            LOG2R_CALC["Log2Ratio Calculation:\nLog2R = log2(SampleNormalizedDepth / PoNMedianDepth)"]
            SEGMENTATION["Segmentation via CBS / HMM\n(Calibrated against Baseline Variance)"]
            CNV_GENOTYPE["Integer Copy State Calling:\n- Loss (CN=0), Het Del (CN=1), Neutral (CN=2)\n- Single Gain (CN=3), High Amplification (CN>=4)"]
            CNV_VCF["${sample}.cnv.vcf.gz + .tbi"]
            CNV_METRICS["${sample}.cnv_metrics.csv"]
            
            DUP_MARK --> SAMPLE_COUNTS
            SAMPLE_COUNTS --> LOG2R_CALC
            PON_OUT --> LOG2R_CALC
            LOG2R_CALC --> SEGMENTATION
            SEGMENTATION --> CNV_GENOTYPE
            CNV_GENOTYPE --> CNV_VCF
            CNV_GENOTYPE --> CNV_METRICS
        end
    end

    subgraph SV_Branch ["DRAGEN Structural Variant Discovery Engine"]
        EXTRACT_SV["Extract Candidate Variant Reads\n- Discordant Paired-End (PE) Clusters\n- Soft-Clipped & Split Reads (SR)"]
        MANTA_ASSEMBLY["String-Graph De Novo Assembly\n(Integrated Manta Assembly into Assembled Contigs)"]
        GRAPH_REMAP["Contig Remapping to Multigenome Graph Reference"]
        SV_CALL["Breakend Genotyping & Resolution:\nDEL, DUP, INV, INS, BND (Single-Nucleotide Precision)"]
        SV_VCF["${sample}.sv.vcf.gz + .tbi"]
        
        DUP_MARK --> EXTRACT_SV
        EXTRACT_SV --> MANTA_ASSEMBLY
        MANTA_ASSEMBLY --> GRAPH_REMAP
        GRAPH_REMAP --> SV_CALL
        SV_CALL --> SV_VCF
    end

    subgraph Downstream_Integration ["CNV-Caller-Combiner Master Downstream Integration"]
        NORM_CNV["FILTER_DRAGEN_QUALITY & NORM_DRAGEN_FULL\n(modules/normalise/normalise-quality.nf)"]
        CONSENSUS_CNV["CNV Consensus Merging Engine\n(SURVIVOR & Truvari)"]
        CATBOOST_ML["CatBoost ML Prioritization\n(Dual-Pathway Clinical Scorer)"]
        
        REFINEMENT_ENGINE["Breakpoint Refinement Engine\n(bin/refine_breakpoints.py)\nHybrid Reconciliation:\nRead-Depth CNV Segments + Nucleotide-Resolution SV Breakends"]
        CONSENSUS_SV["SV Consensus Merging Engine\n(Truvari SV Collapsing)"]
        
        CNV_VCF --> NORM_CNV
        NORM_CNV --> CONSENSUS_CNV
        NORM_CNV --> CATBOOST_ML
        
        NORM_CNV --> REFINEMENT_ENGINE
        SV_VCF --> REFINEMENT_ENGINE
        SV_VCF --> CONSENSUS_SV
    end
```

---

## 9. Comparative Architectural & Statistical Summary Matrix

| Tool | Category | Statistical Foundation | Primary Covariates & Normalization | Container Runtime | Output VCF Name Pattern | Downstream Role |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CANOES** | Exome CNV | Negative Binomial Generalized Linear Model + Viterbi HMM | GC, $\text{GC}^2$, target exon length, library size | `phelelani/nf-exomecnv:canoes` | `${sample}_CANOES_output.vcf.gz` | High specificity in medium cohorts ($N=100\text{--}500$); CatBoost feature contributor. |
| **CLAMMS** | Exome CNV | Dynamic PCA reference selection + Gaussian Mixture Models | Mappability tracks, Picard HS metrics (fold-80, AT/GC drop), sex | `phelelani/nf-exomecnv:clamms` | `${sample}_CLAMMS_output.vcf.gz` | Mappability-aware calling; allosomal sex-calibrated genotyping. |
| **XHMM** | Exome CNV | Singular Value Decomposition (SVD) PCA + Distance HMM | Target mean depth, extreme GC exclusion ($<0.1$ or $>0.9$), PCA variance subtraction | `phelelani/nf-exomecnv:xhmm` | `${sample}_XHMM_output.vcf.gz` | Removes technical batch and sequencing plate artifacts. |
| **CNVkit** | Exome & Genome CNV | Circular Binary Segmentation (CBS) on On/Off-target bins | Dynamic autobin sizing, refFlat annotations, pooled reference normalization | `biocontainers/cnvkit:0.9.12` | `${sample}_CNVKIT_output.vcf.gz` | Broad multi-exon and whole-gene CNVs; antitarget off-probe resolution. |
| **GATK-gCNV** | Germline Exome CNV | Hierarchical Bayesian Graphical Model (MCMC / Variational Inference) | GC waves, probe length, latent sample bias loadings, contig ploidy priors | `broadinstitute/gatk:4.6.2.0` | `${sample}.genotyped_segments.vcf.gz` | Gold-standard statistical authority; complex and allosomal loci. |
| **ExomeDepth** | Exome CNV | Over-dispersed Beta-Binomial Distribution + 3-state HMM | Step-wise sample correlation forward selection ($\ge 0.97$), Bayes Factor | `quay.io/biocontainers/r-exomedepth:1.1.16` | `${sample}_EXOMEDEPTH_output.vcf.gz` | Highest sensitivity for single-exon deletions and duplications. |
| **SvABA** | Structural Variant (SV) | Local de novo String-Graph Assembly (FM-index/BWT) + BWA-MEM | Target capture restriction (`-k`), candidate read filtering, assembly contigs | `quay.io/biocontainers/svaba:1.2.0` | `${sample}.svaba.sv.vcf.gz` | Primary provider of single-nucleotide breakpoint coordinates for refinement. |
| **DRAGEN v4.4.6** | Combined CNV & SV | Hardware FPGA Multigenome Graph + Median PoN CBS/HMM + Manta Assembly | Graph alignment, quadratic GC correction, target length, median PoN, PE/SR clusters | Dedicated ICAv2 Cloud FPGA (`fpga2-medium`) | `${sample}.cnv.vcf.gz` & `${sample}.sv.vcf.gz` | Hardware-accelerated rapid caller; dual read-depth and split-read engine. |
