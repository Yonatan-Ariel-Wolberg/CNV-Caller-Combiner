# CNV-Caller-Combiner: 9-Phase Workflow Architecture & Detailed Diagrams

This document provides exhaustive, publication-grade Mermaid.js workflow diagrams and operational specifications for all 9 sequential phases of the `CNV-Caller-Combiner` master exome pipeline.

---

## Phase 1: CRAM to BAM Conversion & Target-Padded Slicing

```mermaid
flowchart TD
    subgraph Inputs ["Input Alignments & Target Specifications"]
        CRAM["Input CRAM Files + .crai Index"]
        FASTA["Reference FASTA + .fai"]
        TARGET_BED["Exome Target Capture BED"]
    end

    subgraph Forensic_Engine ["40-Point Forensic Verification Engine"]
        PREFLIGHT["Section A: Pre-Flight Integrity & System Gates
- quickcheck, deterministic pipefail, BAM > 1MB, BAI > 0B
- mtime(index) >= mtime(alignment), 28-byte EOF block"]
        HEADER_PARITY["Section B: Header & Metadata Parity
- Coordinate sort (@HD SO:coordinate)
- Contig sequence dictionary & lengths (@SQ)
- Read groups (@RG), Program history (@PG), M5 checksum"]
        COUNT_PARITY["Section C: Alignment Counts & Flagstat Concordance
- Exact chromosome read count parity (idxstats: Delta = 0)
- 16-category flagstat diff, proper pairs %, duplicate counts"]
        BIOCHEM_PARITY["Section D: Biochemical & Quality Concordance
- Total mapped bases, read length distribution, Phred quality parity
- Insert size mean & SD, strand orientation, GC %
- CIGAR insertion, deletion, clip counts exact match"]
        SPEC_AUDIT["Section E: Structural Compliance Audit
- Unmapped formatting (MAPQ=0, CIGAR=*), QNAME ASCII integrity
- Base quality range (Q0-Q93), 1-based indexing (POS > 0)
- Core alignment stream MD5 checksum parity"]
    end

    subgraph Conversion_Modes ["Dual Conversion & Slicing Pathways"]
        subgraph Mode_Full ["Path A: Full CRAM-to-BAM Conversion"]
            CONV_FULL["CONVERT_CRAM_TO_BAM_VERIFIED
(bin/convert_cram_to_bam_verified.sh)
- samtools view -b -T ${fasta} ${cram}
- samtools index -b
- Execute 40 Forensic Validation Gates"]
            FULL_BAM["${sample}.bam + ${sample}.bam.bai
(Full-Genome Validated BAM)"]
            FULL_AUDIT["${sample}.forensic_audit.log
(40 Forensic Gate Certificates)"]
        end

        subgraph Mode_Slice ["Path B: Target-Padded Exome Slicing"]
            PAD_BED["Pad Target Intervals (+250 bp Flanks)
bedtools slop -b 250 -g ${fasta}.fai"]
            CONV_SLICE["CONVERT_CRAM_TO_TARGET_PADDED_BAM
(bin/convert_cram_to_target_padded_bam.sh)
- On-the-fly Slicing: samtools view -b -T ${fasta} -L ${padded_bed}
- Coordinate Indexing & 40 Verification Gates
- >80% Storage Reduction & Accelerated Caller I/O"]
            SLICED_BAM["${sample}.target_padded.bam + .bai
(Exome-Padded Sliced BAM)"]
            SLICE_AUDIT["${sample}.target_padded.forensic_audit.log"]
        end
    end

    subgraph Archival ["Lossless Archival Conversion (Optional)"]
        CONV_ARCHIVE["CONVERT_BAM_TO_CRAM_VERIFIED
(bin/convert_bam_to_cram_verified.sh)
- samtools view -C -T ${fasta}
- Lossless Re-compression & MD5 Checksum Verification"]
        ARCHIVE_CRAM["${sample}.archival.cram + .crai"]
    end

    subgraph Phase1_Outputs ["Phase 1 Certified Alignment Stream"]
        CERTIFIED_BAMS["Certified Coordinate-Sorted BAMs + BAIs
(Ready for Cohort QC & Downstream Callers)"]
        FLAGSTATS["Forensic Quality Audit Reports & Flagstat Logs"]
    end

    CRAM --> CONV_FULL
    FASTA --> CONV_FULL
    CONV_FULL --> PREFLIGHT
    PREFLIGHT --> HEADER_PARITY
    HEADER_PARITY --> COUNT_PARITY
    COUNT_PARITY --> BIOCHEM_PARITY
    BIOCHEM_PARITY --> SPEC_AUDIT
    SPEC_AUDIT --> FULL_BAM
    SPEC_AUDIT --> FULL_AUDIT

    CRAM --> CONV_SLICE
    FASTA --> PAD_BED
    TARGET_BED --> PAD_BED
    PAD_BED --> CONV_SLICE
    CONV_SLICE --> PREFLIGHT
    SPEC_AUDIT --> SLICED_BAM
    SPEC_AUDIT --> SLICE_AUDIT

    FULL_BAM -.-> CONV_ARCHIVE
    FASTA -.-> CONV_ARCHIVE
    CONV_ARCHIVE --> ARCHIVE_CRAM

    FULL_BAM --> CERTIFIED_BAMS
    SLICED_BAM --> CERTIFIED_BAMS
    FULL_AUDIT --> FLAGSTATS
    SLICE_AUDIT --> FLAGSTATS
```

---

## Phase 2: Cohort QC and Consensus Target BED Generation

```mermaid
flowchart TD
    subgraph Inputs ["Input Alignment Cohort & Vendor Capture Designs"]
        BAMS["Certified Cohort BAMs / CRAMs + Indices"]
        FASTA["Reference Genome FASTA + .fai"]
        VENDOR_BEDS["Vendor Target BED Files
(Agilent SureSelect, Twist Exome, Illumina TruSeq)"]
    end

    subgraph Target_Consensus ["Consensus Target BED Generation"]
        STANDARDIZE_BEDS["Standardize Target Intervals to 4-Column BED
(chr, start, end, target_id)"]
        BEDTOOLS_MERGE["BEDTOOLS_MERGE & INTERSECT
- Intersect Overlapping Exome Capture Kits
- Resolve Coordinate Gaps & Flanking Windows
- Strip Non-Standard Contigs"]
        CONSENSUS_BED["consensus_targets.bed / probes.bed
(Standardized Cohort Capture Target BED)"]
        
        VENDOR_BEDS --> STANDARDIZE_BEDS
        STANDARDIZE_BEDS --> BEDTOOLS_MERGE
        BEDTOOLS_MERGE --> CONSENSUS_BED
    end

    subgraph Sample_QC ["Target-Aware Sample Quality Assessment"]
        SAMTOOLS_DOC["SAMTOOLS_QC_BAM / SAMTOOLS_QC_CRAM
(samtools bedcov -Q 30 across Target BED)
- Dynamic 3-col vs 4-col Target Offset Handling (AWK $NF)
- Per-Target Breadth & Mean Depth Calculation"]
        QUALIMAP_RUN["QUALIMAP (qualimap bamqc)
- Mean Mapping Quality, Mean On-Target Coverage
- Insert Size Distribution, Duplicate Rate, GC Bias Plot"]
        SEX_INFER["INFER_SEX_AND_PLOIDY
- chrX & chrY Relative Depth Assessment
- Pseudoautosomal Region (PAR) Masking
- Karyotypic Sex Assignment (MALE / FEMALE)"]
        SEX_INFO["sexinfo.txt (Cohort Karyotypic Sex Metadata)"]
        CONTAM_EST["ESTIMATE_CONTAMINATION
- Cross-Sample Heterozygosity & Allele Swaps Audit"]
        
        BAMS --> SAMTOOLS_DOC
        CONSENSUS_BED --> SAMTOOLS_DOC
        BAMS --> QUALIMAP_RUN
        CONSENSUS_BED --> QUALIMAP_RUN
        SAMTOOLS_DOC --> SEX_INFER
        SEX_INFER --> SEX_INFO
        BAMS --> CONTAM_EST
    end

    subgraph QC_Aggregation ["Cohort Aggregation & Threshold Gating"]
        MULTIQC_RUN["MULTIQC
(Aggregate Qualimap, Samtools Flagstat, bedcov)"]
        QC_REPORT["multiqc_report.html & multiqc_data.json"]
        
        QC_GATE{"Cohort Quality Filtering Gates
1. Target Breadth >= 10x (Cutoff >= 75% Sanger / >= 90% Twist)
2. Mean Target Depth >= 30x
3. PCR Duplicate Rate <= 20%
4. Cohort Correlation >= 0.85 (Sanger) / >= 0.90"}
        
        PASS_SAMPLES["Passed Cohort Samples (Clean BAMs)"]
        FAIL_SAMPLES["Flagged Outlier Samples (High Dups / Low Depth)"]
        
        QUALIMAP_RUN --> MULTIQC_RUN
        SAMTOOLS_DOC --> MULTIQC_RUN
        CONTAM_EST --> MULTIQC_RUN
        MULTIQC_RUN --> QC_REPORT
        QC_REPORT --> QC_GATE
        BAMS --> QC_GATE
        QC_GATE -- Pass --> PASS_SAMPLES
        QC_GATE -- Fail --> FAIL_SAMPLES
    end

    subgraph Phase2_Outputs ["Phase 2 Certified Cohort Foundation"]
        COHORT_BAM_LIST["cohort_bams.list (Validated Cohort Alignment Manifest)"]
        VERIFIED_SEX["sexinfo.txt (Authentication for CLAMMS & GATK-gCNV)"]
        FINAL_PROBES["probes.bed (Unified Consensus Target Capture BED)"]
    end

    PASS_SAMPLES --> COHORT_BAM_LIST
    SEX_INFO --> VERIFIED_SEX
    CONSENSUS_BED --> FINAL_PROBES
```

---

## Phase 3: High-Level Parallel Calling, Quality Score Filtering and Normalisation

```mermaid
flowchart TD
    subgraph Inputs ["Phase 2 Certified Cohort Inputs"]
        COHORT_BAMS["cohort_bams.list (Passed Cohort BAMs)"]
        PROBES["probes.bed (Consensus Target BED)"]
        REF["reference.fasta"]
        SEX["sexinfo.txt (Karyotypic Sex Metadata)"]
    end

    subgraph HighLevel_Calling ["Parallel CNV & SV Caller Execution (High-Level)"]
        subgraph Local_CNV_Lane ["Local Multi-Caller Read-Depth CNV Lane"]
            CALL_CANOES["RUN_CANOES
(Negative Binomial HMM)"]
            CALL_CLAMMS["RUN_CLAMMS
(PCA Mixture Model)"]
            CALL_XHMM["RUN_XHMM
(SVD PCA Denoising)"]
            CALL_CNVKIT["RUN_CNVKIT
(On/Off-Target CBS)"]
            CALL_GCNV["RUN_GCNV
(Hierarchical Bayesian MCMC)"]
            CALL_EXOMEDEPTH["RUN_EXOMEDEPTH
(Beta-Binomial Forward Selection)"]
        end

        subgraph Local_SV_Lane ["Local Assembly SV Lane"]
            CALL_SVABA["RUN_SVABA
(Local String-Graph De Novo Assembly)"]
        end

        subgraph Cloud_Lane ["Cloud FPGA Unified Enrichment Lane (ICAv2)"]
            CALL_DRAGEN["RUN_DRAGEN_GERMLINE_ENRICHMENT
(Hardware FPGA Multigenome Graph Execution)
- Emits DRAGEN CNV Calls (*.cnv.vcf.gz)
- Emits DRAGEN SV Calls (*.sv.vcf.gz)"]
        end
    end

    subgraph Raw_Callsets ["Raw Candidate Variant Callsets"]
        RAW_CANOES["${sample}_CANOES_raw.vcf.gz"]
        RAW_CLAMMS["${sample}_CLAMMS_raw.vcf.gz"]
        RAW_XHMM["${sample}_XHMM_raw.vcf.gz"]
        RAW_CNVKIT["${sample}_CNVKIT_raw.vcf.gz"]
        RAW_GCNV["${sample}_GCNV_raw.vcf.gz"]
        RAW_EXOMEDEPTH["${sample}_EXOMEDEPTH_raw.vcf.gz"]
        RAW_SVABA["${sample}_SvABA_raw.sv.vcf.gz"]
        RAW_DRAGEN_CNV["${sample}_DRAGEN_raw.cnv.vcf.gz"]
        RAW_DRAGEN_SV["${sample}_DRAGEN_raw.sv.vcf.gz"]
    end

    COHORT_BAMS --> HighLevel_Calling
    PROBES --> HighLevel_Calling
    REF --> HighLevel_Calling
    SEX --> HighLevel_Calling

    CALL_CANOES --> RAW_CANOES
    CALL_CLAMMS --> RAW_CLAMMS
    CALL_XHMM --> RAW_XHMM
    CALL_CNVKIT --> RAW_CNVKIT
    CALL_GCNV --> RAW_GCNV
    CALL_EXOMEDEPTH --> RAW_EXOMEDEPTH
    CALL_SVABA --> RAW_SVABA
    CALL_DRAGEN --> RAW_DRAGEN_CNV
    CALL_DRAGEN --> RAW_DRAGEN_SV

    subgraph Quality_Filtering ["Stage 2: Strict Quality Score Filtering (quality-filters.nf)"]
        FILT_CANOES["FILTER_CANOES_QUALITY
(Q_SOME >= 80, Length >= 1000bp)"]
        FILT_CLAMMS["FILTER_CLAMMS_QUALITY
(Q_SOME >= 100, Q_EXACT_DEL >= 0.5, Q_EXACT_DUP >= -0.2)"]
        FILT_XHMM["FILTER_XHMM_QUALITY
(SQ >= 60, NDQ >= 60, EQ >= 60, Length >= 1000bp)"]
        FILT_CNVKIT["FILTER_CNVKIT_QUALITY
(Probes >= 5, Length >= 10000bp, GQ >= 1)"]
        FILT_GCNV["FILTER_GCNV_QUALITY
(QS_DEL >= 400, QS_DUP >= 400, Single_Interval_QS >= 1000)"]
        FILT_ED["FILTER_EXOMEDEPTH_QUALITY
(Bayes Factor >= 30.0, Length >= 1000bp)"]
        FILT_DRAGEN_CNV["FILTER_DRAGEN_QUALITY
(QUAL >= 10, FILTER == PASS)"]
        FILT_SVABA["FILTER_SVABA_QUALITY
(QUAL >= 30, Contig Realignment MapQ >= 20)"]
        FILT_DRAGEN_SV["FILTER_DRAGEN_SV_QUALITY
(QUAL >= 15, Read Support >= 5)"]

        RAW_CANOES --> FILT_CANOES
        RAW_CLAMMS --> FILT_CLAMMS
        RAW_XHMM --> FILT_XHMM
        RAW_CNVKIT --> FILT_CNVKIT
        RAW_GCNV --> FILT_GCNV
        RAW_EXOMEDEPTH --> FILT_ED
        RAW_DRAGEN_CNV --> FILT_DRAGEN_CNV
        RAW_SVABA --> FILT_SVABA
        RAW_DRAGEN_SV --> FILT_DRAGEN_SV
    end

    subgraph Normalisation ["Stage 3: Unified Quality Normalisation (normalise-quality.nf)"]
        NORM_CANOES["NORM_CANOES_FULL"]
        NORM_CLAMMS["NORM_CLAMMS_FULL"]
        NORM_XHMM["NORM_XHMM_FULL"]
        NORM_CNVKIT["NORM_CNVKIT_FULL"]
        NORM_GCNV["NORM_GCNV_FULL"]
        NORM_ED["NORM_EXOMEDEPTH_FULL"]
        NORM_DRAGEN_CNV["NORM_DRAGEN_FULL"]
        
        HARMONIZE["normalise_cnv_caller_quality_scores.py
- Transform Raw Tool Scores to Unified Phred Confidence: Q_NORM in [0, 999]
- Standardize SVTYPE (DEL, DUP), SVLEN, END Coordinates
- Standardize CIPOS & CIEND Confidence Intervals"]
        
        FILT_CANOES --> NORM_CANOES
        FILT_CLAMMS --> NORM_CLAMMS
        FILT_XHMM --> NORM_XHMM
        FILT_CNVKIT --> NORM_CNVKIT
        FILT_GCNV --> NORM_GCNV
        FILT_ED --> NORM_ED
        FILT_DRAGEN_CNV --> NORM_DRAGEN_CNV

        NORM_CANOES --> HARMONIZE
        NORM_CLAMMS --> HARMONIZE
        NORM_XHMM --> HARMONIZE
        NORM_CNVKIT --> HARMONIZE
        NORM_GCNV --> HARMONIZE
        NORM_ED --> HARMONIZE
        NORM_DRAGEN_CNV --> HARMONIZE
    end

    subgraph Phase3_Outputs ["Phase 3 Standardized Callsets"]
        NORM_CNV_VCFS["${sample}_${caller}_normalised.vcf.gz + .tbi
(Standardized Exome CNV Callsets)"]
        QUAL_SV_VCFS["${sample}_${sv_caller}_filtered.sv.vcf.gz + .tbi
(High-Confidence SV Breakends: SvABA & DRAGEN SV)"]
    end

    HARMONIZE --> NORM_CNV_VCFS
    FILT_SVABA --> QUAL_SV_VCFS
    FILT_DRAGEN_SV --> QUAL_SV_VCFS
```

---

## Phase 4: Intra-Sample Merging, Dual-Pathway Breakpoint Resolution, and Cohort Merging

```mermaid
flowchart TD
    subgraph Inputs ["Standardized Callsets from Phase 3"]
        NORM_CNVS["Normalised CNV Callsets per Sample
(CANOES, CLAMMS, XHMM, CNVkit, GATK-gCNV, ExomeDepth, DRAGEN CNV)"]
        QUAL_SVS["Filtered SV Assembly Callsets per Sample
(SvABA & DRAGEN SV)"]
        REF["reference.fasta"]
        PROBES["probes.bed"]
    end

    subgraph Step1_Intra_Merge ["Step 1: Intra-Sample Inter-Tool Merging (SURVIVOR or Truvari)"]
        MERGE_TOOL{"Select Ensemble Merger Engine
(SURVIVOR or Truvari)"}
        
        SURV_INTRA["SURVIVOR (Intra-Sample Mode)
- Reciprocal Overlap: >= 50%
- Boundary Proximity: <= 1000 bp
- Strandedness & SVTYPE Parity
- Agreement Support: >= 2 Callers"]
        TRUV_INTRA["TRUVARI (Intra-Sample Mode)
- pctsize >= 0.70, pctovl >= 0.50
- Intra-Sample Collapsing to Benchmark Union/Intersection"]
        
        INTRA_CNV["${sample}_intra_consensus.vcf.gz + .tbi
(Intra-Sample Exon-Bounded Read-Depth Consensus)"]
        
        NORM_CNVS --> MERGE_TOOL
        MERGE_TOOL -- SURVIVOR --> SURV_INTRA
        MERGE_TOOL -- Truvari --> TRUV_INTRA
        SURV_INTRA --> INTRA_CNV
        TRUV_INTRA --> INTRA_CNV
    end

    subgraph Step2_Breakpoint_Refine ["Step 2: Dual-Pathway Breakpoint Resolution Engine (bin/refine_breakpoints.py)"]
        ROUTE_PATH{"Select Execution Architecture:
Local vs Cloud-Based"}
        
        subgraph Local_Refine_Lane ["Local Architecture Pathway"]
            LOCAL_DEPTH["Local Consensus CNV Callset
(CANOES + CLAMMS + XHMM + CNVkit + GATK-gCNV + ExomeDepth)"]
            LOCAL_SV["Local SV Assembly Breakends
(SvABA Local String-Graph Contigs)"]
            REFINE_LOCAL["REFINE_BREAKPOINTS (Local Mode)
- Intersect Exon-Bounded Depth Intervals with SvABA Microhomology
- Reconcile Coarse Boundaries (+-100..500bp) to Exact Nucleotide Junctions
- Set CIPOS=(0,0) & Append REFINE_METHOD=SVABA_ASSEMBLY"]
            LOCAL_REFINED_VCF["${sample}_local_refined.vcf.gz"]
            
            LOCAL_DEPTH --> REFINE_LOCAL
            LOCAL_SV --> REFINE_LOCAL
            REFINE_LOCAL --> LOCAL_REFINED_VCF
        end

        subgraph Cloud_Refine_Lane ["Cloud-Based Architecture Pathway"]
            CLOUD_DEPTH["Cloud DRAGEN CNV Callset
(DRAGEN Germline Enrichment CNV Caller)"]
            CLOUD_SV["Cloud DRAGEN SV Callset
(DRAGEN Germline Enrichment SV Caller)"]
            REFINE_CLOUD["REFINE_BREAKPOINTS (Cloud Mode)
- Intersect DRAGEN Depth Segments with DRAGEN Manta Assembly Breakends
- Single-Base Junction Resolution at Breakpoint Junctions
- Append REFINE_METHOD=DRAGEN_MANTA"]
            CLOUD_REFINED_VCF["${sample}_cloud_refined.vcf.gz"]
            
            CLOUD_DEPTH --> REFINE_CLOUD
            CLOUD_SV --> REFINE_CLOUD
            REFINE_CLOUD --> CLOUD_REFINED_VCF
        end

        INTRA_CNV --> ROUTE_PATH
        QUAL_SVS --> ROUTE_PATH
        ROUTE_PATH -- Local Callers --> Local_Refine_Lane
        ROUTE_PATH -- Cloud ICAv2 Callers --> Cloud_Refine_Lane
        
        REFINED_SAMPLE_VCF["${sample}_refined_consensus.vcf.gz + .tbi
(Base-Pair Resolution Intra-Sample Consensus)"]
        LOCAL_REFINED_VCF --> REFINED_SAMPLE_VCF
        CLOUD_REFINED_VCF --> REFINED_SAMPLE_VCF
    end

    subgraph Step3_Cohort_Merge ["Step 3: Cohort Inter-Sample Merging (SURVIVOR or Truvari)"]
        COHORT_MERGE_TOOL{"Select Inter-Sample Merger
(SURVIVOR or Truvari)"}
        
        SURV_COHORT["MERGE_CONSENSUS_COHORT (SURVIVOR)
- Combine All Cohort Refined Sample VCFs
- Parameterize: max_distance=1000bp, min_support=1
- Preserve Caller Provenance & Tool Mask Strings"]
        TRUV_COHORT["MERGE_CONSENSUS_COHORT (Truvari)
- Truvari Collapse across Cohort
- Multi-Sample Site Clustering & Genotype Harmonization"]
        
        POST_PROCESS["Consensus Postprocessing (consensus-postprocess.nf)
- SANITIZE_COHORT_VCF
- BCFTOOLS_FLATTEN_COHORT (Generate Per-Sample Genotype Columns)"]
        
        COHORT_VCF["cohort_consensus_merged.vcf.gz + .tbi
(Unified Cohort Structural Variation Callset)"]
        
        REFINED_SAMPLE_VCF --> COHORT_MERGE_TOOL
        COHORT_MERGE_TOOL -- SURVIVOR --> SURV_COHORT
        COHORT_MERGE_TOOL -- Truvari --> TRUV_COHORT
        SURV_COHORT --> POST_PROCESS
        TRUV_COHORT --> POST_PROCESS
        POST_PROCESS --> COHORT_VCF
    end

    subgraph Phase4_Outputs ["Phase 4 Deliverables"]
        INTRA_DELIV["Sample-Level Refined Consensus Calls"]
        COHORT_DELIV["Cohort-Level Unified Consensus Callset
(Ready for Internal AF & Machine Learning)"]
    end

    REFINED_SAMPLE_VCF --> INTRA_DELIV
    COHORT_VCF --> COHORT_DELIV
```

---

## Phase 5: Internal Allele Frequency Filtering

```mermaid
flowchart TD
    subgraph Inputs ["Phase 4 Cohort Callsets"]
        COHORT_VCF["cohort_consensus_merged.vcf.gz + .tbi
(Cohort Unified Consensus Callset)"]
        PROBES["probes.bed (Target Capture BED)"]
    end

    subgraph Internal_AF_Engine ["Strict Cohort-Internal AF Engine (modules/filters/allele-frequency-filters.nf)"]
        CALC_AC_AN["Calculate Internal Allele Count (AC) & Allele Number (AN)
- Parse Sample Genotypes across N Cohort Individuals (AN = 2 * N)
- Eliminate External Database Ascertainment & Ancestry Bias
- AF_internal = AC / AN"]
        
        AF_METRICS["Internal Cohort Variant Frequency Distribution
- Common Polymorphisms (AF >= 1.0% or Carrier Count > 3)
- Rare Structural Variants (AF < 1.0% or Singletons/Doubletons)"]
        
        FILTER_RARE["FILTER_RARE_CNVS
(bin/filter_rare_cnvs.py & filter-rare-cnvs-ddduk-90.py)
- Isolate Clinically Relevant Rare CNVs (AF < 0.01)
- Segregate Common Structural Polymorphisms into Baseline"]
        
        CALC_AC_AN --> AF_METRICS
        AF_METRICS --> FILTER_RARE
    end

    subgraph Callset_Partitioning ["Callset Partitioning & Stratification"]
        RARE_COHORT_VCF["cohort_rare_cnvs.vcf.gz + .tbi
(Rare Structural Variants AF < 1%)"]
        COMMON_COHORT_VCF["cohort_common_polymorphisms.vcf.gz + .tbi
(Common Structural Variants AF >= 1%)"]
        AF_TABLE["cohort_allele_frequency_summary.tsv
(AC, AN, AF_internal, Carrier Sample IDs)"]
        
        FILTER_RARE --> RARE_COHORT_VCF
        FILTER_RARE --> COMMON_COHORT_VCF
        FILTER_RARE --> AF_TABLE
    end

    subgraph Sample_Decomposition ["Per-Sample Rare VCF Extraction"]
        FILTER_BY_REGIONS["FILTER_VCF_BY_RARE_REGIONS
(Decompose Cohort Rare VCF to Individual Patients)
bcftools view -s ${sample} -i 'GT!=0/0 & GT!=./.'"]
        SAMPLE_RARE_VCF["${sample}_rare_filtered.vcf.gz + .tbi
(Patient-Specific Rare CNV Callset)"]
        
        RARE_COHORT_VCF --> FILTER_BY_REGIONS
        FILTER_BY_REGIONS --> SAMPLE_RARE_VCF
    end

    subgraph Phase5_Outputs ["Phase 5 Deliverables"]
        RARE_DELIV["Clean Rare Structural Variant Callset
(Purged of Benign Population Polymorphisms)"]
        BENIGN_DELIV["Cohort Common Structural Baseline
(Informative for Population Architecture)"]
    end

    SAMPLE_RARE_VCF --> RARE_DELIV
    COMMON_COHORT_VCF --> BENIGN_DELIV
```

---

## Phase 6: Feature Extraction, CatBoost ML Training, Explainability, and In Silico Validation

```mermaid
flowchart TD
    subgraph Inputs ["Input Genomic Data & Rare Callsets"]
        RARE_VCFS["${sample}_rare_filtered.vcf.gz
(Phase 5 Rare CNV Calls)"]
        COHORT_BAMS["cohort_bams.list"]
        PROBES["probes.bed (Target Intervals)"]
        REF["reference.fasta"]
        TRUTH_LABELS["truth_labels.tsv (Orthogonal aCGH / MLPA Confirmations)"]
    end

    subgraph Feature_Extraction_Stage ["Step 1: Multidimensional Feature Extraction (feature_extraction.py)"]
        BIOPHYSICAL["1. Biophysical Genomic Features
- Segment Length (bp), Target GC Content (%), Flank GC
- Probe Count (NP), Probe Density (probes/kb)
- Sequence Complexity & Repeat Masking Density (LTR/LINE/SINE)"]
        CONCORDANCE["2. Ensemble Concordance Features
- Caller Agreement Count (N/7)
- Binary Tool Support Vector (CANOES, CLAMMS, XHMM, etc.)
- Reciprocal Overlap Ratios & Breakpoint Refinement Method"]
        L2R_METRICS["3. Probe-Level Depth Metrics (COMPUTE_TARGET_COVERAGE)
- Mean Log2 Depth Ratio (L2R), Median Absolute Deviation (MAD)
- L2R Variance, Denoised Copy Ratio Dispersion"]
        
        FEATURE_EXTRACT["FEATURE_EXTRACTION (modules-feature-extraction.nf)
Compile 20+ Genomic, Biophysical, and Statistical Features"]
        SAMPLE_FEATURES["${sample}_features.tsv"]
        
        RARE_VCFS --> FEATURE_EXTRACT
        COHORT_BAMS --> FEATURE_EXTRACT
        PROBES --> FEATURE_EXTRACT
        REF --> FEATURE_EXTRACT
        BIOPHYSICAL --> FEATURE_EXTRACT
        CONCORDANCE --> FEATURE_EXTRACT
        L2R_METRICS --> FEATURE_EXTRACT
        FEATURE_EXTRACT --> SAMPLE_FEATURES
    end

    subgraph Cohort_Database ["Step 2: Cohort Feature Database Generation"]
        BUILD_DB["BUILD_COHORT_FEATURE_DATABASE
(bin/build_cohort_feature_database.py)
- Aggregate All Per-Sample Feature Tables
- Handle Missingness & Log-Transform Skewed Covariates
- Serialize to Compressed Columnar Format"]
        PARQUET_DB["cohort_features.parquet / cohort_features.db
(Unified SQL / Parquet Cohort Feature Database)"]
        
        SAMPLE_FEATURES --> BUILD_DB
        BUILD_DB --> PARQUET_DB
    end

    subgraph CatBoost_Training ["Step 3: CatBoost Training & Stratified Cross-Validation (modules-train.nf)"]
        SPLIT_DATA["Stratified GroupKFold Cross-Validation (K=5)
- Grouped by Genomic Locus / Sample ID
- Eliminates Target Overlap Data Leakage"]
        TRAIN_MODEL["CatBoost Gradient-Boosted Decision Trees (train_predict_catboost.py)
- Hyperparameters: iterations=1000, depth=6, l2_leaf_reg=3
- Objective: Logloss & Focal Loss with Class-Weight Balancing
- Early Stopping on Evaluation Validation Set"]
        CATBOOST_MODEL["catboost_cnv_model.cbm
(Trained CatBoost Ensemble Classifier)"]
        
        PARQUET_DB --> SPLIT_DATA
        TRUTH_LABELS --> SPLIT_DATA
        SPLIT_DATA --> TRAIN_MODEL
        TRAIN_MODEL --> CATBOOST_MODEL
    end

    subgraph Model_Evaluation ["Step 4: Explainability & Discrimination Metrics Assessment"]
        SHAP_ANALYSIS["SHAP (SHapley Additive exPlanations) Values
- Global Feature Importance Ranking
- Local Individual Prediction Force Plots"]
        ROC_PR_METRICS["Discriminative Performance Metrics
- Receiver Operating Characteristic (AUC-ROC >= 0.95)
- Precision-Recall Curve (PR-AUC >= 0.90)"]
        RECLASS_METRICS["Clinical Value-Added Metrics
- Net Reclassification Improvement (NRI > 0)
- Integrated Discrimination Improvement (IDI > 0)
Quantifies Diagnostic Gain Over Standalone Callers"]
        
        TRAIN_MODEL --> SHAP_ANALYSIS
        TRAIN_MODEL --> ROC_PR_METRICS
        TRAIN_MODEL --> RECLASS_METRICS
    end

    subgraph InSilico_Validation ["Step 5: In Silico Validation & Cohort Prioritization (modules-predict.nf)"]
        RUN_PREDICT["RUN_CATBOOST_PREDICT
- Ingest Entire Cohort Feature Database
- Calculate Continuous Pathogenicity Probability: P_score in [0.0, 1.0]
- Apply Diagnostic Operating Cutoff: P_score >= 0.70"]
        SCORED_PARQUET["cohort_catboost_scored.parquet"]
        FILTER_ML["FILTER_VCF_BY_ML_SCORE
Append CATBOOST_PROB & Flag FILTER=PASS_PRIORITIZED"]
        ML_PRIORITIZED_VCF["${sample}_ml_prioritized.vcf.gz + .tbi
(ML-Validated High-Confidence Rare CNVs)"]
        
        CATBOOST_MODEL --> RUN_PREDICT
        PARQUET_DB --> RUN_PREDICT
        RUN_PREDICT --> SCORED_PARQUET
        SCORED_PARQUET --> FILTER_ML
        RARE_VCFS --> FILTER_ML
        FILTER_ML --> ML_PRIORITIZED_VCF
    end

    subgraph Phase6_Outputs ["Phase 6 Deliverables"]
        FEATURE_DELIV["Cohort Feature SQL/Parquet Database"]
        EVAL_DELIV["SHAP, AUC-ROC, PR-AUC, NRI, and IDI Evaluation Reports"]
        MODEL_DELIV["Validated ML-Prioritized Clinical Callsets"]
    end

    PARQUET_DB --> FEATURE_DELIV
    SHAP_ANALYSIS --> EVAL_DELIV
    ROC_PR_METRICS --> EVAL_DELIV
    RECLASS_METRICS --> EVAL_DELIV
    ML_PRIORITIZED_VCF --> MODEL_DELIV
```

---

## Phase 7: AnnotSV Annotation and ACMG 4/5 Clinical Prioritisation

```mermaid
flowchart TD
    subgraph Inputs ["Input Callsets & Genomic Builds"]
        ML_VCF["${sample}_ml_prioritized.vcf.gz + .tbi
(Phase 6 ML-Prioritized Callsets)"]
        GENOME_BUILD["Genome Build Identifier (GRCh38 / hg38)"]
    end

    subgraph Dynamic_DB ["Step 1: Dynamic AnnotSV Database Download & Compilation"]
        DOWNLOAD_TRIGGER{"Check Local AnnotSV Database Cache"}
        SETUP_SCRIPT["DOWNLOAD_ANNOTSV_DB
(bin/setup_reference_bundle.sh)
- Sync External Clinical Databases:
  ClinVar, ClinGen Dosage Sensitivity, Decipher, OMIM
  DGV, 1000 Genomes SV, gnomAD SV v2.1/v4.0
- Pre-compile Gene Overlap & Exon Coordinate Indices"]
        ANNOTSV_DB["AnnotSV Annotation Cache Directory
($ANNOTSV/share/AnnotSV/)"]
        
        GENOME_BUILD --> DOWNLOAD_TRIGGER
        DOWNLOAD_TRIGGER -- Cache Missing --> SETUP_SCRIPT
        DOWNLOAD_TRIGGER -- Cache Exists --> ANNOTSV_DB
        SETUP_SCRIPT --> ANNOTSV_DB
    end

    subgraph Variant_Annotation ["Step 2: Comprehensive Variant Annotation (RUN_ANNOTSV)"]
        RUN_ANNOTSV["RUN_ANNOTSV (modules/annotation/annotation.nf)
Container: phelelani/annotsv:3.4.2
- Dual-Level Record Formatting:
  1. Full Event Record (split = 0): Aggregate Overlap
  2. Transcripts Record (split = 1): Per-Gene Exonic Overlap
- Gene Overlaps, CDS Truncations, HI/TS Dosage Flags"]
        ANNOT_TSV["${sample}.annotated.tsv
(Comprehensive Multi-Database Annotation Table)"]
        ANNOT_VCF["${sample}.annotated.vcf.gz
(Annotated VCF with ANNOTSV Info Tags)"]
        
        ML_VCF --> RUN_ANNOTSV
        ANNOTSV_DB --> RUN_ANNOTSV
        RUN_ANNOTSV --> ANNOT_TSV
        RUN_ANNOTSV --> ANNOT_VCF
    end

    subgraph ACMG_Tiering ["Step 3: Automated ACMG/ClinGen 5-Class Pathogenicity Tiering"]
        ADD_ACMG["ADD_ACMG_ANNOTATIONS
(bin/add_acmg_annotations.py)
Automated 2019/2020 ACMG/ClinGen Guidelines:
- Section 1: Genomic Content (HI/TS, Gene Count)
- Section 2: Overlap with Established Benign Regions (DGV, gnomAD)
- Section 3: Evaluation of Gene Dosage Sensitivity
- Section 4: Detailed Case-Level & Inheritance Data
- Section 5: Phenotypic Specificity Matching"]
        
        ACMG_CLASSES["ACMG Pathogenicity Stratification:
Class 1: Benign
Class 2: Likely Benign
Class 3: Variant of Uncertain Significance (VUS)
Class 4: Likely Pathogenic
Class 5: Pathogenic"]
        
        ANNOT_TSV --> ADD_ACMG
        ANNOT_VCF --> ADD_ACMG
        ADD_ACMG --> ACMG_CLASSES
    end

    subgraph ACMG_Prioritisation ["Step 4: ACMG 4/5 Clinical Prioritisation (STRATIFY_VCF_BY_ACMG)"]
        STRATIFY["STRATIFY_VCF_BY_ACMG
(bin/stratify_vcf_by_acmg.py)
- Filter: ACMG_CLASS in [Class 4, Class 5]
- Filter for Single-Gene & Syndromic Actionable Events
- Partition Full Callset into Tiered VCFs"]
        
        TIER_4_5_VCF["${sample}_acmg_tier4_5_prioritized.vcf.gz + .tbi
(Actionable Pathogenic & Likely Pathogenic CNVs)"]
        ALL_TIERS_VCF["${sample}_acmg_all_classes.vcf.gz + .tbi
(Full Callset with Standardized ACMG Class Tags)"]
        CLINICAL_REPORT["${sample}_acmg_clinical_report.md & .tsv
(Diagnostic Clinical Summary Report)"]
        
        ACMG_CLASSES --> STRATIFY
        STRATIFY --> TIER_4_5_VCF
        STRATIFY --> ALL_TIERS_VCF
        STRATIFY --> CLINICAL_REPORT
    end

    subgraph Phase7_Outputs ["Phase 7 Deliverables"]
        ACTIONABLE_DELIV["Actionable Pathogenic Diagnostic Candidates (Class 4 & 5)"]
        ANNOTATED_DELIV["Fully Annotated Cohort Clinical Callsets"]
        REPORT_DELIV["Executive Diagnostic ACMG Summaries"]
    end

    TIER_4_5_VCF --> ACTIONABLE_DELIV
    ALL_TIERS_VCF --> ANNOTATED_DELIV
    CLINICAL_REPORT --> REPORT_DELIV
```

---

## Phase 8: Quality Score Metric Extraction and Analysis

```mermaid
flowchart TD
    subgraph Inputs ["Callsets Across All Pipeline Stages & Stratifications"]
        STAGE1_VCFS["Stage 1: Raw Caller Callsets"]
        STAGE2_VCFS["Stage 2: Quality-Filtered Callsets"]
        STAGE3_VCFS["Stage 3: Normalised Callsets"]
        STAGE4_VCFS["Stage 4: Refined Consensus Callsets"]
        STAGE5_VCFS["Stage 5: Rare-Filtered Callsets"]
        STAGE6_VCFS["Stage 6: CatBoost ML-Prioritized Callsets"]
        STAGE7_VCFS["Stage 7: ACMG Tiered Callsets (Class 1..5)"]
        PROBES["probes.bed (Target Intervals)"]
    end

    subgraph Metric_Extraction ["Step 1: Multi-Dimensional Quality Metric Extraction (quality-metric-evaluation.nf)"]
        ANALYZE_QC["ANALYZE_CNV_QUALITY_METRICS
(bin/analyze_cnv_quality_metrics.py)
- Variant Counts: Total Variants, DEL Count, DUP Count, DEL/DUP Ratio
- Length Metrics: Median Size, IQR, Total Footprint (Mb)
- Probe Coverage: Median Probes/Call, Single-Exon vs Multi-Exon Counts
- Confidence Metrics: Mean QUAL, Q_NORM Distribution
- Depth Metrics: Mean Target Depth, GC Percentage, Mappability"]
        
        SAMPLE_QC_CSV["${run_name}_cnv_vcf_metrics.csv
(Granular Per-Variant Metrics)"]
        SUMMARY_TXT["${run_name}_cnv_metrics_summary.txt
(Stage-Level Summary Statistics)"]
        
        STAGE1_VCFS --> ANALYZE_QC
        STAGE2_VCFS --> ANALYZE_QC
        STAGE3_VCFS --> ANALYZE_QC
        STAGE4_VCFS --> ANALYZE_QC
        STAGE5_VCFS --> ANALYZE_QC
        STAGE6_VCFS --> ANALYZE_QC
        STAGE7_VCFS --> ANALYZE_QC
        PROBES --> ANALYZE_QC
        ANALYZE_QC --> SAMPLE_QC_CSV
        ANALYZE_QC --> SUMMARY_TXT
    end

    subgraph Callset_Assessment ["Step 2: Callset Quality Auditing & Calibration Assessment"]
        ASSESS_CALLSETS["ASSESS_CALLSET_QUALITY_METRICS
- Verify Phred Score Calibration across Callers
- Single-Exon vs Multi-Exon Confidence Stratification
- Sample Outlier Detection (Extreme Call Burdens or Depth Shifts)"]
        AUDIT_REPORT["callset_quality_audit.tsv
(Cohort Quality Calibration Table)"]
        
        SAMPLE_QC_CSV --> ASSESS_CALLSETS
        ASSESS_CALLSETS --> AUDIT_REPORT
    end

    subgraph MultiStage_Compilation ["Step 3: Compilation of Multi-Stage Quality Report"]
        COMPILE_REPORT["COMPILE_ALL_STAGES_QUALITY_REPORT
- Synthesize Metrics across All 7 Callers & Consensus Stages
- Tabulate Quality Score Trajectory across Pipeline Funnel
- Generate Stage Summary Tables (.md, .tsv)"]
        ALL_STAGES_REPORT["all_stages_quality_report.md & .tsv
(Comprehensive Multi-Stage Quality Audit)"]
        
        SUMMARY_TXT --> COMPILE_REPORT
        AUDIT_REPORT --> COMPILE_REPORT
        COMPILE_REPORT --> ALL_STAGES_REPORT
    end

    subgraph Phase8_Outputs ["Phase 8 Quality Analytics Deliverables"]
        REPORT_OUT["Multi-Stage Comprehensive Quality Metric Report"]
        CSV_OUT["Granular Variant-Level Quality Metric Database"]
        AUDIT_OUT["Cohort Calibration & Outlier Diagnosis Audit"]
    end

    ALL_STAGES_REPORT --> REPORT_OUT
    SAMPLE_QC_CSV --> CSV_OUT
    AUDIT_REPORT --> AUDIT_OUT
```

---

## Phase 9: Performance Metric Evaluation, Figure Generation, and Attrition Funnel

```mermaid
flowchart TD
    subgraph Inputs ["Evaluated Callsets & Ground Truth Validation Benchmarks"]
        ALL_STAGE_VCFS["Callsets across All Pipeline Stages (Stages 1..10)
(Raw, Filtered, Normalised, Intra-Merged, Refined,
Cohort-Merged, Rare, ML-Prioritized, ACMG 4/5)"]
        TRUTH_BED["Gold-Standard Truthset BED
(Orthogonal aCGH, MLPA, WGS, or GIAB Benchmarks)"]
        EVAL_REPORTS["Multi-Stage Quality Reports from Phase 8"]
        PROBES["probes.bed (Target Intervals)"]
    end

    subgraph BED_Standardization ["Step 1: Callset Format Standardization (performance-metric-evaluation.nf)"]
        VCF2BED["VCF_TO_BED (bin/vcf_to_bed.py)
Convert Per-Sample VCFs to 5-Column BED Files
(chr, start, end, sample_id, svtype)"]
        COMBINE_BEDS["COMBINE_BEDS
Combine Per-Sample BEDs into Unified Stage Callset BEDs"]
        STAGE_BEDS["${stage}_unified_callset.bed"]
        
        ALL_STAGE_VCFS --> VCF2BED
        VCF2BED --> COMBINE_BEDS
        COMBINE_BEDS --> STAGE_BEDS
    end

    subgraph Performance_Benchmarking ["Step 2: Diagnostic Performance Evaluation (EVALUATE_COMPARATIVE)"]
        EVAL_COMP["EVALUATE_COMPARATIVE
(bin/evaluate_cnv_callset.py)
- Overlap Reciprocal Threshold: >= 50%
- Calculate Confusion Matrix Parameters:
  True Positives (TP), False Positives (FP)
  False Negatives (FN), True Negatives (TN)
- Compute Diagnostic Performance Indices:
  1. Sensitivity / Recall: R = TP / (TP + FN)
  2. Precision / PPV: P = TP / (TP + FP)
  3. F_beta Score (beta = 2, Clinically Weighted for Recall):
     F_2 = (5 * P * R) / (4 * P + R)
  4. Matthews Correlation Coefficient (MCC):
     MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
  5. False Positive Rate (FPR) & False Negative Rate (FNR)
  6. Specificity: TN / (TN + FP)"]
        
        PERF_METRICS["performance_metrics_summary.tsv
(Comparative Diagnostic Performance Table)"]
        ACMG_EVAL["EVALUATE_ACMG_CSV
(Performance Evaluation Stratified by ACMG Classes 1..5)"]
        
        STAGE_BEDS --> EVAL_COMP
        TRUTH_BED --> EVAL_COMP
        PROBES --> EVAL_COMP
        EVAL_COMP --> PERF_METRICS
        STAGE_BEDS --> ACMG_EVAL
        TRUTH_BED --> ACMG_EVAL
    end

    subgraph Attrition_Funnel_Engine ["Step 3: Variant Attrition Funnel Dynamics (GENERATE_ATTRITION_FUNNEL)"]
        FUNNEL_GEN["GENERATE_ATTRITION_FUNNEL
(bin/generate_attrition_funnel.py)
- Track Sequential Candidate Variant Counts & Retention %:
  Raw Calls -> QC Filtered -> Intra Consensus -> Break Refined
  -> Rare AF Filtered -> CatBoost ML -> ACMG 4/5 Prioritized
- Quantify False Positive Pruning vs True Variant Preservation"]
        
        FUNNEL_PLOT["attrition_funnel.png
(High-Resolution Diagnostic Funnel Diagram)"]
        FUNNEL_TABLE["attrition_funnel_table.md & .tsv
(Stage-by-Stage Attrition Statistics Table)"]
        
        EVAL_REPORTS --> FUNNEL_GEN
        PERF_METRICS --> FUNNEL_GEN
        FUNNEL_GEN --> FUNNEL_PLOT
        FUNNEL_GEN --> FUNNEL_TABLE
    end

    subgraph Visualization_Suite ["Step 4: Publication-Ready Figure Generation (PLOT_EVALUATION_METRICS)"]
        PLOT_METRICS["PLOT_EVALUATION_METRICS
(bin/plot_evaluation_metrics.py)
Generate Publication Figures:
1. Comparative ROC & Precision-Recall Curves
2. F2-Score & MCC Radar / Spider Charts
3. Caller Concordance UpSet Plots
4. Single-Exon vs Multi-Exon Sensitivity Curves
5. Clinical Utility Decision Curve Analysis (DCA)"]
        
        PUB_PLOTS["evaluation_figures/
- roc_curves.png, pr_curves.png
- radar_metrics.png, upset_concordance.png
- f2_score_comparison.png, dca_curves.png"]
        
        PERF_METRICS --> PLOT_METRICS
        ACMG_EVAL --> PLOT_METRICS
        PLOT_METRICS --> PUB_PLOTS
    end

    subgraph Phase9_Outputs ["Phase 9 Final Deliverables"]
        BENCHMARK_DELIV["Diagnostic Performance Benchmark Report (Recall, Precision, F2, MCC)"]
        FUNNEL_DELIV["Pipeline Variant Attrition Funnel Plot & Table"]
        FIGURES_DELIV["Publication-Grade Diagnostic Figures & Decision Curves"]
    end

    PERF_METRICS --> BENCHMARK_DELIV
    FUNNEL_PLOT --> FUNNEL_DELIV
    FUNNEL_TABLE --> FUNNEL_DELIV
    PUB_PLOTS --> FIGURES_DELIV
```

---

