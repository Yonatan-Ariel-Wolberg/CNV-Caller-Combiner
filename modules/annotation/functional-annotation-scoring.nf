#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE: functional-annotation-scoring.nf
// Purpose:
//   1. Annotates structural and copy number variants using AnnotSV purely for gene,
//      transcript, and functional impact features.
//   2. Prioritizes variants using a 100% FREQUENCY-AGNOSTIC biological scoring model.
//      No population frequency databases (e.g. gnomAD-SV, 1000 Genomes, DGV, ClinGen AF)
//      are used or required, eliminating database ascertainment bias in diverse cohorts.
//
// Designed to be simple, easy to read, maintain, and modify for all researchers.
// =====================================================================================

params.outdir          = './results'
params.annotsv_dir     = "${projectDir}/annotsv_db"
params.genome_build    = 'GRCh38'
params.min_score       = 0.0  // Set e.g. 40.0 to only keep Moderate & High priority calls
params.strvctvre_tabix = null // Optional path to precomputed StrVCTVRE tabix file

// =====================================================================================
// PROCESS 1: ANNOTATE_FUNCTIONAL_ANNOTSV
// =====================================================================================
// Annotates candidate CNV calls with overlapping genes, coding exons (CDS),
// gene loss-of-function intolerance (pLI, LOEUF), and OMIM phenotypes.
process ANNOTATE_FUNCTIONAL_ANNOTSV {
    tag "${sample_id}"
    label 'annotsv'
    publishDir "${params.outdir}/functional_scoring/1_annotated", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path("${sample_id}.annotsv.tsv"), emit: annotated_tsv

    script:
    def genome_build = params.get('genome_build', params.get('build', 'GRCh38'))
    """
    #!/bin/bash
    set -euo pipefail

    export ANNOTSV="\$(readlink -f ${params.annotsv_dir})"

    # Run AnnotSV in full annotation mode
    AnnotSV \\
        -SVinputFile "${vcf_file}" \\
        -genomeBuild "${genome_build}" \\
        -outputDir . \\
        -outputFile "${sample_id}.annotsv" \\
        -annotationMode full \\
        -overwrite 1 || true

    # Safety fallback: create an empty template header if input VCF had zero variants
    if [ ! -f "${sample_id}.annotsv.tsv" ]; then
        echo -e "Annotation_mode\\tSV_chrom\\tSV_start\\tSV_end\\tSV_type\\tGene_name\\tpLI\\tLOEUF\\tClinGen_HI\\tClinGen_TS\\tOMIM_morbid\\tLocation" > "${sample_id}.annotsv.tsv"
    fi
    """
}

// =====================================================================================
// PROCESS 2: ANNOTATE_STRVCTVRE
// =====================================================================================
// Annotates candidate CNVs with StrVCTVRE deleteriousness scores (Sharo et al. 2021)
// using precomputed tabix database or functional constraint fallback.
process ANNOTATE_STRVCTVRE {
    tag "${sample_id}"
    label 'python'
    publishDir "${params.outdir}/functional_scoring/2_strvctvre", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(annotated_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.strvctvre.tsv"), emit: strvctvre_tsv

    script:
    def strvctvre_arg = params.strvctvre_tabix ? "--strvctvre-tabix \"${params.strvctvre_tabix}\"" : ""
    """
    python3 ${projectDir}/bin/annotate_strvctvre.py \\
        --input "${annotated_tsv}" \\
        --output "${sample_id}.strvctvre.tsv" \\
        --genome-build "${params.genome_build}" \\
        ${strvctvre_arg}
    """
}

// =====================================================================================
// PROCESS 3: SCORE_FREQUENCY_AGNOSTIC
// =====================================================================================
// Scores variants (0 to 100 points) based entirely on intrinsic biological constraints:
//   - Gene Intolerance: High pLI (>= 0.90) or Low LOEUF (<= 0.35)
//   - Dosage Sensitivity: ClinGen Haploinsufficiency (loss) / Triplosensitivity (gain)
//   - Disease Relevance: OMIM Morbid Gene overlap
//   - Coding Disruption: Exon / CDS disruption vs intronic
//   - StrVCTVRE Deleteriousness: >= 0.58 pathogenic threshold / >= 0.80 high confidence
//   - Algorithm Concordance: Multi-caller support (CANOES, CLAMMS, XHMM, etc.)
//
// Outputs:
//   - Biological_Priority_Score (0 - 100)
//   - Priority_Tier (Tier 1: High Constraint, Tier 2: Moderate, Tier 3: Low/Tolerated)
//   - Biological_Criteria (Human-readable rationale for the score)
process SCORE_FREQUENCY_AGNOSTIC {
    tag "${sample_id}"
    label 'python'
    publishDir "${params.outdir}/functional_scoring/3_prioritized", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(annotated_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.prioritized.tsv"), emit: scored_tsv

    script:
    """
    python3 ${projectDir}/bin/score_frequency_agnostic.py \\
        --input "${annotated_tsv}" \\
        --output "${sample_id}.prioritized.tsv" \\
        --min-score ${params.min_score}
    """
}

// =====================================================================================
// WORKFLOW DEFINITION
// =====================================================================================
workflow FUNCTIONAL_ANNOTATION_SCORING {
    take:
        vcf_ch  // tuple(sample_id, vcf_file)

    main:
        // Step 1: Annotate with AnnotSV (genes, pLI, LOEUF, CDS, OMIM)
        ANNOTATE_FUNCTIONAL_ANNOTSV(vcf_ch)

        // Step 2: Annotate with StrVCTVRE deleteriousness scores
        ANNOTATE_STRVCTVRE(ANNOTATE_FUNCTIONAL_ANNOTSV.out.annotated_tsv)

        // Step 3: Score using frequency-agnostic biological constraints & StrVCTVRE
        SCORE_FREQUENCY_AGNOSTIC(ANNOTATE_STRVCTVRE.out.strvctvre_tsv)

    emit:
        annotated_tsv   = ANNOTATE_FUNCTIONAL_ANNOTSV.out.annotated_tsv
        strvctvre_tsv   = ANNOTATE_STRVCTVRE.out.strvctvre_tsv
        prioritized_tsv = SCORE_FREQUENCY_AGNOSTIC.out.scored_tsv
}
