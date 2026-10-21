#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// MODULE: phen2gene.nf
// Purpose:
//   Prioritizes patient CNVs using patient-specific Human Phenotype Ontology (HPO)
//   clinical manifestations via Phen2Gene (Zhao et al., Nucleic Acids Res 2020).
//
// Inputs:
//   - Two-column manifest: (sample_id, hpo_terms)
//   - Annotated CNV TSV (AnnotSV or consensus tables)
//
// Outputs:
//   - Patient-specific ranked gene tables (<sample_id>.gene_weights.tsv)
//   - HPO-prioritized candidate variant TSVs (<sample_id>.hpo_annotated.tsv)
// =====================================================================================

params.outdir = './results'
params.p2g_kb = null

// Process 1: Generate gene priority weights from an entire cohort manifest
process PHEN2GENE_FROM_MANIFEST {
    tag "cohort_manifest"
    label 'python'
    publishDir "${params.outdir}/phen2gene/weights", mode: 'copy', overwrite: true

    input:
    path hpo_manifest

    output:
    path "weights/*.gene_weights.tsv", emit: gene_weights

    script:
    def kb_arg = params.p2g_kb ? "--p2g-kb \"${params.p2g_kb}\"" : ""
    """
    python3 ${projectDir}/bin/run_phen2gene_wrapper.py \\
        --manifest "${hpo_manifest}" \\
        --output-dir weights \\
        ${kb_arg}
    """
}

// Process 2: Generate gene priority weights for a single sample tuple
process PHEN2GENE_SAMPLE {
    tag "${sample_id}"
    label 'python'
    publishDir "${params.outdir}/phen2gene/weights", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(hpo_terms)

    output:
    tuple val(sample_id), path("${sample_id}.gene_weights.tsv"), emit: gene_weights

    script:
    def kb_arg = params.p2g_kb ? "--p2g-kb \"${params.p2g_kb}\"" : ""
    """
    python3 ${projectDir}/bin/run_phen2gene_wrapper.py \\
        --sample-id "${sample_id}" \\
        --hpo-terms "${hpo_terms}" \\
        --output-dir . \\
        ${kb_arg}
    """
}

// Process 3: Annotate CNV candidates with patient-specific HPO gene match scores
process ANNOTATE_CNVS_WITH_HPO {
    tag "${sample_id}"
    label 'python'
    publishDir "${params.outdir}/phen2gene/annotated", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(cnv_tsv), path(gene_weights_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.hpo_annotated.tsv"), emit: hpo_annotated_tsv

    script:
    """
    python3 - << 'EOF'
import sys
import pandas as pd
import numpy as np

cnv_file = "${cnv_tsv}"
weights_file = "${gene_weights_tsv}"
out_file = "${sample_id}.hpo_annotated.tsv"

sep = "\\t" if cnv_file.endswith(".tsv") or ".tsv." in cnv_file else ","
df_cnv = pd.read_csv(cnv_file, sep=sep, low_memory=False)

if df_cnv.empty:
    df_cnv["Phen2Gene_Score"] = []
    df_cnv["Phen2Gene_Rank"] = []
    df_cnv["Matched_HPO_Terms"] = []
    df_cnv.to_csv(out_file, sep="\\t", index=False)
    sys.exit(0)

df_w = pd.read_csv(weights_file, sep="\\t") if open(weights_file).read().strip() else pd.DataFrame()
gene_map = {}
if not df_w.empty and "Gene_name" in df_w.columns:
    for _, r in df_w.iterrows():
        gene = str(r["Gene_name"]).strip().upper()
        gene_map[gene] = {
            "score": r.get("Phen2Gene_Score", np.nan),
            "rank": r.get("Rank", np.nan),
            "hpo": str(r.get("Matched_HPO_Terms", ""))
        }

gene_col = next((c for c in ["Gene_name", "gene", "GENE", "SYMBOL"] if c in df_cnv.columns), None)

p2g_scores, p2g_ranks, p2g_hpos = [], [], []
for _, row in df_cnv.iterrows():
    val_score = np.nan
    val_rank = np.nan
    val_hpo = ""
    if gene_col and pd.notna(row[gene_col]):
        genes = [g.strip().upper() for g in str(row[gene_col]).replace(";", ",").split(",") if g.strip()]
        for g in genes:
            if g in gene_map:
                s = gene_map[g]["score"]
                if np.isnan(val_score) or s > val_score:
                    val_score = s
                    val_rank = gene_map[g]["rank"]
                    val_hpo = gene_map[g]["hpo"]
    p2g_scores.append(val_score)
    p2g_ranks.append(val_rank)
    p2g_hpos.append(val_hpo)

df_cnv["Phen2Gene_Score"] = p2g_scores
df_cnv["Phen2Gene_Rank"] = p2g_ranks
df_cnv["Matched_HPO_Terms"] = p2g_hpos
df_cnv.to_csv(out_file, sep="\\t", index=False)
print(f"Annotated {len(df_cnv)} CNVs with Phen2Gene scores.")
EOF
    """
}

workflow PHEN2GENE_PRIORITIZATION {
    take:
        manifest_ch // path(manifest)
        cnv_tsv_ch  // tuple(sample_id, path(cnv_tsv))

    main:
        PHEN2GENE_FROM_MANIFEST(manifest_ch)
        
        // Flatten output TSVs and associate with sample_id
        weights_by_sample = PHEN2GENE_FROM_MANIFEST.out.gene_weights
            .flatten()
            .map { file -> tuple(file.name.replace('.gene_weights.tsv', ''), file) }

        // Join CNVs with matching sample gene weights
        joined_ch = cnv_tsv_ch.join(weights_by_sample)

        ANNOTATE_CNVS_WITH_HPO(joined_ch)

    emit:
        hpo_annotated = ANNOTATE_CNVS_WITH_HPO.out.hpo_annotated_tsv
}
