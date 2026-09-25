nextflow.enable.dsl=2

// =====================================================================================
// QUALITY SCORE & SIZE FILTERING THRESHOLDS
// =====================================================================================
// Parameter adjustments across the 3 benchmark profiles (Balanced, High Precision, Max Sensitivity).
// For comprehensive documentation and biological justifications, see:
// docs/PARAMETER_ADJUSTMENTS.md
//
// Key caller-specific threshold tunings:
//   - CANOES: Q_SOME (65 balanced / 80 high_prec / 50 max_sens), min_size (500 / 1000 / 0 bp)
//   - CLAMMS: Q_SOME (50 / 100 / 30), Q_EXACT_DEL (0.0 / 0.5 / -1.0), Q_EXACT_DUP (-0.5 / -0.2 / -1.0)
//   - XHMM:   SQ (50 / 60 / 30), NDQ & EQ (30 / 60 / 0), min_size (500 / 1000 / 0 bp)
//   - CNVkit: min_probes (3 / 5 / 2), min_size (5000 / 10000 / 1000 bp), max_size <= 5 Mb
//   - GATK:   QS_DEL (100 / 400 / 50), QS_DUP (50 / 400 / 20), single_interval_qs (100 / 1000 / 50)
//   - ExomeDepth: BF (10.0 / 30.0 / 5.0), min_size (500 / 1000 / 0 bp)
//   - DRAGEN: QUAL (30 / 50 / 20), min_bc (2 / 3 / 1), min_size (5000 / 10000 / 2000 bp)
// =====================================================================================

def get_filter_threshold(String key, def default_balanced) {
    def raw_profile = params.get('filter-profile', params.get('filter_profile', 'balanced'))
    def profile = (raw_profile == 'default') ? 'balanced' : raw_profile
    def profile_map = [
        'balanced': [
            canoes_q_some_min: 65, canoes_min_size: 500,
            clamms_q_some_min: 50, clamms_q_exact_del_min: 0.0, clamms_q_exact_dup_min: -0.5,
            xhmm_sq_min: 50, xhmm_ndq_min: 30, xhmm_eq_min: 30, xhmm_min_size: 500,
            cnvkit_min_size: 5000, cnvkit_min_probes: 3,
            gcnv_qs_del_min: 100, gcnv_qs_dup_min: 50, gcnv_single_interval_qs_min: 100,
            dragen_qual_min: 30, dragen_min_bc: 2, dragen_min_size: 5000,
            exomedepth_bf_min: 10.0, exomedepth_min_size: 500
        ],
        'high_precision': [
            canoes_q_some_min: 80, canoes_min_size: 1000,
            clamms_q_some_min: 100, clamms_q_exact_del_min: 0.5, clamms_q_exact_dup_min: -0.2,
            xhmm_sq_min: 60, xhmm_ndq_min: 60, xhmm_eq_min: 60, xhmm_min_size: 1000,
            cnvkit_min_size: 10000, cnvkit_min_probes: 5,
            gcnv_qs_del_min: 400, gcnv_qs_dup_min: 400, gcnv_single_interval_qs_min: 1000,
            dragen_qual_min: 50, dragen_min_bc: 3, dragen_min_size: 10000,
            exomedepth_bf_min: 30.0, exomedepth_min_size: 1000
        ],
        'max_sensitivity': [
            canoes_q_some_min: 50, canoes_min_size: 0,
            clamms_q_some_min: 30, clamms_q_exact_del_min: -1.0, clamms_q_exact_dup_min: -1.0,
            xhmm_sq_min: 30, xhmm_ndq_min: 0, xhmm_eq_min: 0, xhmm_min_size: 0,
            cnvkit_min_size: 1000, cnvkit_min_probes: 2,
            gcnv_qs_del_min: 50, gcnv_qs_dup_min: 20, gcnv_single_interval_qs_min: 50,
            dragen_qual_min: 20, dragen_min_bc: 1, dragen_min_size: 2000,
            exomedepth_bf_min: 5.0, exomedepth_min_size: 0
        ]
    ]
    profile_map['default'] = profile_map['balanced']
    def chosen_profile = profile_map.get(profile, profile_map['balanced'])
    def profile_default = chosen_profile.get(key, default_balanced)
    def val = params.get(key, null)
    return (val != null) ? val : profile_default
}

process FILTER_CANOES_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def q_some_min = get_filter_threshold('canoes_q_some_min', 65)
    def min_size   = get_filter_threshold('canoes_min_size', 500)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.canoes.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Filter CANOES calls by format field Q_SOME >= ${q_some_min}
        bcftools view -i "FMT/Q_SOME >= ${q_some_min}" reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                # Determine SVTYPE: prefer ALT column (\$5) if present, fallback to GT
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                else {
                    split(\$10, format_values, ":");
                    gt_val = format_values[1];
                    if (gt_val == "1") svtype="DEL";
                    else if (gt_val == "2") svtype="DUP";
                }
                
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                
                # Parse END coordinate and enforce minimum size filter (>= ${min_size} bp)
                end_val = "";
                n_fields = split(\$8, info_fields, ";");
                for (i = 1; i <= n_fields; i++) {
                    if (info_fields[i] ~ /^END=/) {
                        split(info_fields[i], end_part, "=");
                        end_val = end_part[2];
                        break;
                    }
                }
                if (end_val != "") {
                    len = end_val - \$2;
                    if (len < ${min_size}) next;
                }
                
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_CLAMMS_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def q_some_min = get_filter_threshold('clamms_q_some_min', 50)
    def q_del_min  = get_filter_threshold('clamms_q_exact_del_min', 0.0)
    def q_dup_min  = get_filter_threshold('clamms_q_exact_dup_min', -0.5)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.clamms.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Filter CLAMMS calls by format fields (configurable thresholds):
        # Q_SOME >= ${q_some_min}, Q_EXACT >= ${q_del_min} for DEL, Q_EXACT >= ${q_dup_min} for DUP
        bcftools view -i 'FMT/Q_SOME >= '"${q_some_min}"' && ((ALT="<DEL>" && FMT/Q_EXACT >= '"${q_del_min}"') || (ALT="<DUP>" && FMT/Q_EXACT >= '"${q_dup_min}"') || (ALT!="<DEL>" && ALT!="<DUP>" && FMT/Q_EXACT >= '"${q_dup_min}"'))' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                # Determine SVTYPE: ALWAYS check ALT column (\$5) first to prevent GT misclassification
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                else {
                    split(\$10, format_values, ":");
                    gt_val = format_values[1];
                    if (gt_val == "1") svtype="DEL";
                    else if (gt_val == "2") svtype="DUP";
                }
                
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_XHMM_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def sq_min   = get_filter_threshold('xhmm_sq_min', 50)
    def ndq_min  = get_filter_threshold('xhmm_ndq_min', 30)
    def eq_min   = get_filter_threshold('xhmm_eq_min', 30)
    def min_size = get_filter_threshold('xhmm_min_size', 500)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.xhmm.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Filter XHMM calls by quality metrics (configurable thresholds):
        # SQ >= ${sq_min}; EQ/NDQ optional >= ${ndq_min} / ${eq_min}
        bcftools view -i '((ALT="<DEL>" || ALT="DEL") && FMT/SQ[0:0] >= '"${sq_min}"' && (FMT/NDQ[0] >= '"${ndq_min}"' || FMT/NDQ[0]=".") && (FMT/EQ[0:0] >= '"${eq_min}"' || FMT/EQ[0:0]=".")) || ((ALT="<DUP>" || ALT="DUP") && FMT/SQ[0:1] >= '"${sq_min}"' && (FMT/NDQ[0] >= '"${ndq_min}"' || FMT/NDQ[0]=".") && (FMT/EQ[0:1] >= '"${eq_min}"' || FMT/EQ[0:1]=".")) || (ALT!="<DEL>" && ALT!="<DUP>" && (FMT/SQ[0:0] >= '"${sq_min}"' || FMT/SQ[0:1] >= '"${sq_min}"'))' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                # Determine SVTYPE: ALWAYS check ALT column (\$5) first
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                else {
                    split(\$10, format_vals, ":");
                    gt = format_vals[1];
                    if (gt ~ /^1/) svtype="DEL";
                    else if (gt ~ /^2/) svtype="DUP";
                }
                
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                
                # Minimum size filter (>= ${min_size} bp)
                end_val = "";
                n_fields = split(\$8, info_fields, ";");
                for (i = 1; i <= n_fields; i++) {
                    if (info_fields[i] ~ /^END=/) {
                        split(info_fields[i], end_part, "=");
                        end_val = end_part[2];
                        break;
                    }
                }
                if (end_val != "") {
                    len = end_val - \$2;
                    if (len < ${min_size}) next;
                }
                
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_CNVKIT_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def min_size   = get_filter_threshold('cnvkit_min_size', 5000)
    def min_probes = get_filter_threshold('cnvkit_min_probes', 3)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.cnvkit.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Format SVTYPE info tag to DEL/DUP for SURVIVOR compatibility and keep only PASS/unfiltered records
        # Apply configurable filters: Min length >= ${min_size} bp, min probes/bins >= ${min_probes}
        bcftools view -i 'FILTER="PASS" || FILTER="."' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL" || \$5 ~ /DEL/) svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP" || \$5 ~ /DUP/) svtype="DUP";
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                
                # Parse END coordinate, length, and probe count
                end_val = "";
                probes_val = 0;
                n_fields = split(\$8, info_fields, ";");
                for (i = 1; i <= n_fields; i++) {
                    if (info_fields[i] ~ /^END=/) {
                        split(info_fields[i], end_part, "=");
                        end_val = end_part[2];
                    } else if (info_fields[i] ~ /^PROBES=/) {
                        split(info_fields[i], probe_part, "=");
                        probes_val = probe_part[2];
                    }
                }
                if (end_val != "") {
                    len = end_val - \$2;
                    if (len < ${min_size}) next;
                }
                if (probes_val != "" && probes_val > 0 && probes_val < ${min_probes}) next;
                
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_GCNV_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)
    path probes_bed

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def qs_del_min   = get_filter_threshold('gcnv_qs_del_min', 100)
    def qs_dup_min   = get_filter_threshold('gcnv_qs_dup_min', 50)
    def single_int_qs = get_filter_threshold('gcnv_single_interval_qs_min', 100)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.gcnv.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Filter only actual deletions or duplications with high confidence
        # Supports both raw GCNV (FMT/QS) and normalized VCFs (QUAL)
        bcftools view -T "${probes_bed}" -i '((ALT="<DEL>" || ALT="DEL") && (FMT/QS>='"${qs_del_min}"' || QUAL>='"${qs_del_min}"')) || ((ALT="<DUP>" || ALT="DUP") && (FMT/QS>='"${qs_dup_min}"' || QUAL>='"${qs_dup_min}"'))' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                # Parse format keys to find NP index
                split(\$9, format_keys, ":");
                np_idx = 0;
                for (i=1; i<=length(format_keys); i++) {
                    if (format_keys[i] == "NP") {
                        np_idx = i;
                        break;
                    }
                }
                
                # Parse sample fields
                split(\$10, fields, ":");
                if (np_idx > 0) {
                    np = fields[np_idx];
                    # Single-interval calls (NP=1) require highest QS (QUAL >= ${single_int_qs})
                    if (np < 2 && \$6 < ${single_int_qs}) next;
                }
                
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                
                # Convert haploid GT to diploid to prevent SURVIVOR merge format field parse error
                gt = fields[1]
                if (gt !~ /[\\/|]/) {
                    if (gt == ".") {
                        fields[1] = "./."
                    } else {
                        fields[1] = gt "/" gt
                    }
                    s = fields[1]
                    for (i=2; i<=length(fields); i++) {
                        s = s ":" fields[i]
                    }
                    \$10 = s
                }
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_DRAGEN_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def qual_min = get_filter_threshold('dragen_qual_min', 30)
    def min_bc   = get_filter_threshold('dragen_min_bc', 2)
    def min_size = get_filter_threshold('dragen_min_size', 5000)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.dragen.qs-filtered.cnv.vcf.gz"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf

        bcftools view -i '((FILTER="PASS" || FILTER=".") && QUAL>='"${qual_min}"' && (FORMAT/BC>='"${min_bc}"' || FORMAT/BC="."))' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;
                
                # Length filter >= ${min_size} bp
                end_val = "";
                n_fields = split(\$8, info_fields, ";");
                for (i = 1; i <= n_fields; i++) {
                    if (info_fields[i] ~ /^END=/) {
                        split(info_fields[i], end_part, "=");
                        end_val = end_part[2];
                        break;
                    }
                }
                if (end_val != "") {
                    len = end_val - \$2;
                    if (len < ${min_size}) next;
                }
                
                print;
            }' > temp_filtered.vcf
            
        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

process FILTER_EXOMEDEPTH_QUALITY {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*.qs-filtered.cnv.vcf.gz'), emit: filtered_data
    path('*.qs-filtered.cnv.vcf.gz.tbi'), emit: filtered_indices, optional: true

    script:
    def bf_min   = get_filter_threshold('exomedepth_bf_min', 10.0)
    def min_size = get_filter_threshold('exomedepth_min_size', 500)
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.vcf.gz *.vcf; do
        [ -f "\$vcf" ] || continue
        [[ "\$vcf" =~ \\.(tbi|csi)\$ ]] && continue
        base=\$(basename "\$vcf" | sed -E 's/\\.(vcf|vcf\\.gz)\$//')
        base_clean=\$(echo "\$base" | sed -E 's/\\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)(\\..*)?\$//I' | sed -E 's/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*//I' | sed -E 's/(\\.(sorted|normalised|qs-filtered|filtered|cnv|annotated)|_output|_genotyped_segments)*\$//' | sed 's/\\.md\$//')
        filtered_gz="\${base_clean}.exomedepth.qs-filtered.cnv.vcf.gz"

        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf

        # Filter ExomeDepth calls by BF >= ${bf_min} and minimum size >= ${min_size} bp
        bcftools view -i '((FILTER="PASS" || FILTER=".") && (FMT/BF >= '"${bf_min}"' || INFO/BF >= '"${bf_min}"' || QUAL >= '"${bf_min}"' || FMT/BF="."))' reheadered.vcf \
            | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} {
                svtype="CNV";
                if (\$5 == "<DEL>" || \$5 == "DEL") svtype="DEL";
                else if (\$5 == "<DUP>" || \$5 == "DUP") svtype="DUP";
                if (\$8 ~ /SVTYPE=/) {
                    sub(/SVTYPE=[^;]+/, "SVTYPE=" svtype, \$8);
                } else {
                    \$8 = \$8 ";SVTYPE=" svtype;
                }
                if (\$8 !~ /END=/) next;

                end_val = "";
                n_fields = split(\$8, info_fields, ";");
                for (i = 1; i <= n_fields; i++) {
                    if (info_fields[i] ~ /^END=/) {
                        split(info_fields[i], end_part, "=");
                        end_val = end_part[2];
                        break;
                    }
                }
                if (end_val != "") {
                    len = end_val - \$2;
                    if (len < ${min_size}) next;
                }
                print;
            }' > temp_filtered.vcf

        bcftools sort temp_filtered.vcf -o "\${filtered_gz}" -O z
        tabix -f -p vcf "\${filtered_gz}" || true
        echo "\${filtered_gz}" >> vcf_list.txt
        rm -f temp_filtered.vcf reheadered.vcf new_sample.txt
    done
    """
}

