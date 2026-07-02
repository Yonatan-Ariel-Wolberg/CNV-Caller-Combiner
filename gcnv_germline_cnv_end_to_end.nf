#!/usr/bin/env nextflow

// =============================================================================
// Detailed documentation: gcnv_germline_cnv_end_to_end.nf
// Purpose: This script supports CNV calling and/or post-processing tasks in
//   this project, with parameters and I/O conventions aligned to repository workflows.
// Scope: End-to-end run path for GATK-gCNV.
// Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
// Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
// Operational notes: Keep parameter names synchronized with related workflow/process
//   definitions and update this file alongside upstream interface changes.
// =============================================================================
nextflow.enable.dsl=2

/*
 * -------------------------------------------------------------------------------------
 *  GATK-gCNV GERMLINE CNV END-TO-END ANALYSIS PIPELINE
 *  Consensus BED Generation -> GATK Caller -> Cohort Filtering & ACMG -> Evaluation & Plotting
 * -------------------------------------------------------------------------------------
 */

// 1. Reference & Input Paths
params.fasta                 = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.fasta_fai             = "${params.fasta}.fai"
params.dict                  = "/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.dict"
params.standard_bed          = "/dataG/ddd/data/resources/canoes/probes_sanger.bed"
params.control_interval_list = "/dataG/ddd/data/resources/canoes/probes_sanger.interval_list"
params.bam_dir               = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/bams"

// 2. Pipeline Modes & Options
params.run_modes             = "Control,Mosdepth" // Modes to run. Can be "Control", "Mosdepth", or "Control,Mosdepth"
params.run_calling           = true  // Execute CNV calling locally using GATK-gCNV
params.pass_percent          = 90    // Consensus threshold: regions covered in >= X% of samples

// 3. Pre-called VCF Fallback (used if params.run_calling = false)
params.precalled_control_vcfs  = "/home/ywolberg/nf-multi-exome-cnv/results/sbatch/ddd-uk/10-samples/gcnv/out_GCNV/vcfs/*.vcf.gz"
params.precalled_mosdepth_vcfs = "" // path to pre-called mosdepth GATK VCFs if run_calling = false
params.precalculated_consensus_bed = null

// 4. Cohort Filtering & ACMG Annotations
params.total_samples         = 10
params.max_af                = 0.05
params.outdir                = "./results_gcnv_e2e"
params.annotsv_db            = ""
params.annotsv_version       = "3.5.10"

// 5. Evaluator Paths
params.map_file              = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/DDDUK-MAP/EGAF-EGAN-ddduk.map"
params.truth_bed             = "/spaces/ywolberg/DECIPHERING_DD_DATA/DDD_UK_DATA/ddduk-truthset-cnvs/hg38-ddduk-ega-cnv-truthset.bed"

// Define dynamic repository and bin directory paths depending on script location
def repoDir = file("${baseDir}/main.nf").exists() ? "${baseDir}" : (file("${baseDir}/../main.nf").exists() ? "${baseDir}/.." : "${baseDir}")
def binDir = file("${repoDir}/bin").exists() ? "${repoDir}/bin" : "${repoDir}"


/*
 * =====================================================================================
 *  PHASE 1: MOSDEPTH TARGET BED GENERATION & GATK INTERVAL LIST CONVERSION
 * =====================================================================================
 */

process RUN_MOSDEPTH {
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0'
    cpus 4
    memory '24 GB'
    maxForks 10
    
    input:
    tuple val(sample), path(bam), path(bai)
    path fasta
    path fasta_fai
    path bed

    output:
    path "${sample}.regions.bed.gz", emit: usable_beds
    
    script:
    """
    mosdepth --threads ${task.cpus} \\
             --fasta ${fasta} \\
             --by ${bed} \\
             ${sample} ${bam}
    """
}

process BUILD_CONSENSUS {
    container 'quay.io/biocontainers/bedtools:2.27.1--h077b44d_9'
    cpus 2
    memory '8 GB'
    
    publishDir "${params.outdir}/consensus_bed", mode: 'copy'

    input:
    path usable_beds
    path bed

    output:
    path "ddduk_karyotypic_consensus.bed", emit: karyotypic_bed

    script:
    """
    TOTAL_SAMPLES=\$(ls *.regions.bed.gz | wc -l)
    MIN_PASS=\$(( (TOTAL_SAMPLES * ${params.pass_percent}) / 100 ))

    # Extract target interval coverages for each sample
    for f in *.regions.bed.gz; do
        gunzip -c "\$f" | cut -f4 > tmp_\$f
    done
    
    # Get coordinates from the first file
    first_file=\$(ls *.regions.bed.gz | head -n1)
    gunzip -c "\$first_file" | cut -f1-3 > coords.txt
    
    # Combine coordinates and coverages
    paste coords.txt tmp_*.regions.bed.gz > matrix.txt
    
    # Filter intervals based on sample coverage count
    awk -v min="\$MIN_PASS" -v n="\$TOTAL_SAMPLES" 'BEGIN{FS="\\t"; OFS="\\t"} {
        count = 0;
        for(i=4; i<=3+n; i++) {
            if(\$i >= 20.0) count++;
        }
        if(count >= min) {
            print \$1, \$2, \$3;
        }
    }' matrix.txt > ddduk_clean_targeted_consensus.bed
    
    # Clean up temp files and sort consensus bed
    rm -f tmp_* matrix.txt coords.txt
    sort -k1,1 -k2,2n ddduk_clean_targeted_consensus.bed | bedtools merge -i stdin -d 1500 > ddduk_karyotypic_consensus.bed
    """
}

process CONVERT_BED_TO_INTERVAL_LIST {
    container 'broadinstitute/gatk:4.6.2.0'
    label 'gatk'
    
    input:
    path bed
    path dict

    output:
    path "ddduk_karyotypic_consensus.interval_list", emit: interval_list

    script:
    """
    gatk BedToIntervalList \\
        -I ${bed} \\
        -O ddduk_karyotypic_consensus.interval_list \\
        -SD ${dict}
    """
}


/*
 * =====================================================================================
 *  PHASE 2: GATK-gCNV NESTED CALLING WORKFLOW
 * =====================================================================================
 */

process RUN_GCNV_CALLING {
    tag "${run_name}"
    maxForks 1
    cpus 4
    memory '16 GB'
    publishDir "${params.outdir}/${run_name}/raw_calls", mode: 'copy'

    input:
    tuple val(run_name), path(bed_file), path(interval_list)

    output:
    tuple val(run_name), path("out_vcfs/*.vcf.gz"), emit: vcfs

    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Resolve absolute paths of input files from task directory and copy to persistent folder
    mkdir -p ${repoDir}/.nested_gcnv_work/${run_name}/inputs
    cp -L ${interval_list} ${repoDir}/.nested_gcnv_work/${run_name}/inputs/targets.interval_list
    
    # Create and move to persistent workspace directory
    cd ${repoDir}/.nested_gcnv_work/${run_name}
    
    # Execute nested nextflow run for GATK-gCNV
    nextflow run ${repoDir}/main.nf \\
        -profile wits,small \\
        --workflow gcnv \\
        -params-file ${repoDir}/params/ddd-uk/params-gatk-gcnv-wits-ddd-uk.json \\
        --exome_targets "${repoDir}/.nested_gcnv_work/${run_name}/inputs/targets.interval_list" \\
        --outdir ./tmp_out \\
        --samples_path "${params.bam_dir}/*.{bam,bam.bai}" \\
        -resume
        
    # Move back to parent task directory and copy raw calls to output folder
    cd \$OLDPWD
    mkdir -p out_vcfs
    cp ${repoDir}/.nested_gcnv_work/${run_name}/tmp_out/out_GCNV/vcfs/*.vcf.gz out_vcfs/
    """
}


/*
 * =====================================================================================
 *  PHASE 3: COHORT FILTERING, MERGING & ACMG PRIORITIZATION
 * =====================================================================================
 */

process bcftoolsFilterVcfs {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/filtered_vcfs", mode: 'copy'

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path('vcf_list.txt'), path('*_filtered.cnv.vcf'), emit: filtered_data

    script:
    """
    #!/bin/bash
    set -euo pipefail
    > vcf_list.txt
    for vcf in *.normalised.vcf.gz; do
        base=\$(basename "\$vcf" .normalised.vcf.gz)
        base_clean=\$(echo "\$base" | sed 's/\\.md//g' | sed 's/_GCNV_genotyped_segments//g' | sed 's/\\.genotyped_segments//g')
        filtered_vcf="\${base_clean}_filtered.cnv.vcf"
        
        # Reheader VCF to ensure sample ID matches base_clean exactly
        echo "\${base_clean}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\$vcf" -o reheadered.vcf
        
        # Filter only actual deletions (CN<2, QUAL>=100) or duplications (CN>2, QUAL>=50)
        # and format the SVTYPE info tag to DEL/DUP, and require at least 3 target intervals (NP >= 3)
        bcftools view -i '(ALT="<DEL>" && QUAL>=100) || (ALT="<DUP>" && QUAL>=50)' reheadered.vcf \
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
                    if (np < 3) next;
                }
                
                svtype="CNV";
                if (\$5 == "<DEL>") svtype="DEL";
                else if (\$5 == "<DUP>") svtype="DUP";
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
            }' > "\$filtered_vcf"
            
        if grep -q -m1 -v '^#' "\$filtered_vcf"; then
            echo "\$filtered_vcf" >> vcf_list.txt
        else
            rm -f "\$filtered_vcf"
        fi
    done
    """
}

process SURVIVORMergeCohort {
    label 'survivor'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(vcf_list), path(filtered_vcfs)

    output:
    tuple val(run_name), path('cohort_merged.vcf'), emit: merged_vcf

    script:
    """
    SURVIVOR merge vcf_list.txt 1000 1 1 0 0 0 cohort_merged_raw.vcf
    awk '
    BEGIN { OFS="\\t" }
    /^#/ { print; next }
    {
        num_format_fields = split(\$9, tmp, ":")
        for (i = 10; i <= NF; i++) {
            n = split(\$i, colons, ":")
            if (n > num_format_fields) {
                val = colons[1]
                for (j = 2; j < num_format_fields; j++) {
                    val = val ":" colons[j]
                }
                last_val = colons[num_format_fields]
                for (j = num_format_fields + 1; j <= n; j++) {
                    last_val = last_val "_" colons[j]
                }
                \$i = val ":" last_val
                n = num_format_fields
                split(\$i, colons, ":")
            }
            for (j = 1; j <= n; j++) {
                if (colons[j] == "NA" || colons[j] == "NaN" || colons[j] == "NAN" || colons[j] == "") {
                    colons[j] = "."
                }
            }
            val = colons[1]
            for (j = 2; j <= n; j++) {
                val = val ":" colons[j]
            }
            \$i = val
        }
        print
    }
    ' cohort_merged_raw.vcf > cohort_merged.vcf
    """
}

process bcftoolsFlattenCohort {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf)

    output:
    tuple val(run_name), path('cohort_cnvs_by_sample.tsv'), emit: cohort_tsv

    script:
    """
    bcftools query -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\t%TY\\n]' cohort_merged.vcf \\
        | awk 'BEGIN{OFS="\\t"} \$7 != "NaN" && \$7 != "0/0" && \$7 != "0" && \$7 != "." && \$7 != "" {print \$1, \$2, \$3, \$4, \$5, \$6}' \\
        > cohort_cnvs_by_sample.tsv
    """
}

process filterRareCnvs {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(cohort_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), path('rare_regions.tsv'), emit: rare_data

    script:
    """
    # cache invalidate for small cohort variant filtering fix
    python3 "${binDir}/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --total-samples ${params.total_samples} \\
        --max-af ${params.max_af} \\
        --output-regions "rare_regions.tsv"
    """
}

process filterVcfByRareRegions {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf), path(rare_regions)

    output:
    tuple val(run_name), path('cohort_rare.vcf'), emit: rare_vcf

    script:
    """
    if [ ! -s "${rare_regions}" ]; then
        bcftools view -h "${merged_vcf}" -o cohort_rare.vcf
    else
        bcftools view -T "${rare_regions}" "${merged_vcf}" -o cohort_rare.vcf
    fi
    """
}

process DOWNLOAD_ANNOTSV_DB {
    publishDir "${repoDir}/annotsv_db", mode: 'copy'

    output:
    path "AnnotSV_annotations", emit: db_dir

    script:
    def version_tag = params.annotsv_version.startsWith('v') ? params.annotsv_version : "v${params.annotsv_version}"
    def version_clean = params.annotsv_version.replaceAll(/^v/, '')
    def parts = version_clean.split('\\.')
    def major_minor = parts.size() >= 2 ? "${parts[0]}.${parts[1]}" : version_clean
    """
    git clone --branch ${version_tag} --depth 1 https://github.com/lgmgeo/AnnotSV.git
    cd AnnotSV
    
    # Pre-download the correct human annotations file. Keep it under the name the Makefile expects
    # (Annotations_Human_${major_minor}.tar.gz) and copy it to the versioned clean name just in case.
    curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_${major_minor}.tar.gz
    cp Annotations_Human_${major_minor}.tar.gz Annotations_Human_${version_clean}.tar.gz 2>/dev/null || true
    
    make PREFIX=\$(pwd)/../AnnotSV_annotations install
    make PREFIX=\$(pwd)/../AnnotSV_annotations install-human-annotation
    cd ..
    rm -rf AnnotSV
    """
}

process runAnnotSV {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/cohort_merged", mode: 'copy'

    input:
    tuple val(run_name), path(rare_vcf)
    path db_dir

    output:
    tuple val(run_name), path('cohort_rare.annotated.tsv'), emit: annotated_tsv

    script:
    """
    echo "Running AnnotSV 3.5.10 on ${run_name}..."
    export ANNOTSV=\$(readlink -f ${db_dir})
    
    # Run AnnotSV
    AnnotSV \\
        -SVinputFile "${rare_vcf}" \\
        -outputDir . \\
        -outputFile cohort_rare.annotated \\
        -genomeBuild GRCh38 \\
        -overwrite 0 || true
        
    # If the output file was not created (e.g., empty input), touch it with standard headers
    if [ ! -f cohort_rare.annotated.tsv ]; then
        echo -e "Annotation_mode\\tSV_chrom\\tSV_start\\tSV_end\\tSV_type\\tACMG_class\\tAnnotSV_ranking_score\\tAnnotSV_ranking_criteria\\tGene_name" > cohort_rare.annotated.tsv
    fi
    """
}

process addACMGAnnotations {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/acmg_reports", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), emit: rare_cnvs_annotated
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs_*.csv'), emit: score_separated_csvs

    script:
    """
    python3 "${binDir}/add-acmg-annotations.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv"
    """
}

process ANALYZE_CNV_QUALITY_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}/metrics_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
    path probes_bed

    output:
    tuple val(run_name), path("${run_name}_cnv_metrics_summary.txt"), path("${run_name}_cnv_vcf_metrics.csv"), emit: metrics_reports

    script:
    """
    python3 "${binDir}/analyze_cnv_quality_metrics.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --vcf-dir . \\
        --probes-bed "${probes_bed}" \\
        --run-name "${run_name}" \\
        --output-summary "${run_name}_cnv_metrics_summary.txt" \\
        --output-csv "${run_name}_cnv_vcf_metrics.csv"
    """
}


/*
 * =====================================================================================
 *  PHASE 4: RIGOROUS EVALUATION & VISUALIZATION
 * =====================================================================================
 */

process EVALUATE_ACMG_CSV {
    tag "${run_name}-${csv.baseName}"
    publishDir "${params.outdir}/evaluate-results/${run_name}", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(csv), path(probes_bed)
    path map_file
    path truth_bed

    output:
    path "*_performance_report.txt", emit: report

    script:
    def base_name = csv.baseName
    """
    # cache invalidate for evaluation sample ID mapping fix
    python3 "${binDir}/evaluate_separated_helper.py" \\
        --input-csv "${csv}" \\
        --map-file "${map_file}" \\
        --truth-bed "${truth_bed}" \\
        --probes-bed "${probes_bed}" \\
        --workflow-type "${run_name}" \\
        --output "${run_name}_${base_name}_performance_report.txt"
    """
}

process PLOT_EVALUATION_METRICS {
    publishDir "${params.outdir}/evaluate-results", mode: 'copy', overwrite: true

    input:
    path reports
    
    output:
    path 'performance_metrics_comparison.png', emit: comparison_plot

    script:
    """
    # Structure files into subdirectories so parser detects Control vs Mosdepth
    mkdir -p Control Mosdepth
    for f in *_performance_report.txt; do
        if [[ "\$f" == *"Control"* ]]; then
            cp "\$f" Control/
        elif [[ "\$f" == *"Mosdepth"* ]]; then
            cp "\$f" Mosdepth/
        else
            cp "\$f" Control/
        fi
    done
    
    python3 "${binDir}/plot_evaluation_performance_dynamic.py" \\
        --input-dir . \\
        --output-png performance_metrics_comparison.png
    """
}


/*
 * =====================================================================================
 *  WORKFLOW ENTRY POINT
 * =====================================================================================
 */

workflow {

    // 1. Parse run modes to execute
    def modes = params.run_modes.split(',').collect { it.trim() }

    // Set up AnnotSV database channel (uses host pre-installed DB or downloads dynamically)
    if (params.annotsv_db) {
        annotsv_db_ch = Channel.value(file(params.annotsv_db))
    } else {
        def local_db = file("${repoDir}/annotsv_db/AnnotSV_annotations")
        if (local_db.exists()) {
            annotsv_db_ch = Channel.value(local_db)
        } else {
            DOWNLOAD_ANNOTSV_DB()
            annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
        }
    }
    
    // 2. Generate consensus BED dynamically if 'Mosdepth' run mode is active
    if (modes.contains('Mosdepth')) {
        if (params.precalculated_consensus_bed && file(params.precalculated_consensus_bed).exists()) {
            mosdepth_bed_ch = Channel.value(file(params.precalculated_consensus_bed))
        } else {
            cram_ch = Channel.fromPath("${params.bam_dir}/*.bam")
                .map { bam -> 
                    def bai = file("${bam}.bai")
                    if (!bai.exists()) bai = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
                    if (!bai.exists()) error "Missing BAM index file for ${bam}"
                    return [bam.baseName, bam, bai]
                }
            
            mosdepth_out = RUN_MOSDEPTH(
                cram_ch, 
                file(params.fasta), 
                file(params.fasta_fai), 
                file(params.standard_bed)
            )
            
            BUILD_CONSENSUS(mosdepth_out.usable_beds.collect(), file(params.standard_bed))
            mosdepth_bed_ch = BUILD_CONSENSUS.out.karyotypic_bed
        }
        CONVERT_BED_TO_INTERVAL_LIST(mosdepth_bed_ch, file(params.dict))
        mosdepth_interval_list_ch = CONVERT_BED_TO_INTERVAL_LIST.out.interval_list
    } else {
        mosdepth_bed_ch = Channel.empty()
        mosdepth_interval_list_ch = Channel.empty()
    }

    // 3. Construct target_beds_ch pairing run_name with target BED and GATK interval_list
    control_bed_ch = modes.contains('Control') ? Channel.of([ 'Control', file(params.standard_bed), file(params.control_interval_list) ]) : Channel.empty()
    mosdepth_bed_entry_ch = modes.contains('Mosdepth') ? mosdepth_bed_ch.combine(mosdepth_interval_list_ch).map { bed, il -> [ 'Mosdepth', bed, il ] } : Channel.empty()
    
    target_beds_ch = control_bed_ch.concat(mosdepth_bed_entry_ch)

    // 4. CNV Variant Calling Phase (Dynamic GATK Calling vs Precalled Fallback)
    if (params.run_calling) {
        // Filter modes: if results already exist for a mode, bypass caller execution to save time
        called_modes_ch = target_beds_ch.filter { run_name, bed, il ->
            if (file("${params.outdir}/${run_name}/raw_calls/out_vcfs/EGAF00001051389.md.genotyped_segments.vcf.gz").exists()) {
                log.info "Bypassing GATK CNV calling for mode ${run_name} because raw calls were found in results directory."
                return false
            }
            return true
        }

        bypassed_vcfs_ch = target_beds_ch.filter { run_name, bed, il ->
            if (file("${params.outdir}/${run_name}/raw_calls/out_vcfs/EGAF00001051389.md.genotyped_segments.vcf.gz").exists()) {
                return true
            }
            return false
        }.map { run_name, bed, il ->
            def files = file("${params.outdir}/${run_name}/raw_calls/out_vcfs/*.vcf.gz")
            return [ run_name, files ]
        }

        RUN_GCNV_CALLING(called_modes_ch)
        vcfs_ch = RUN_GCNV_CALLING.out.vcfs.concat(bypassed_vcfs_ch)
    } else {
        // Bypass calling and select pre-called GATK VCF sets
        vcfs_ch = Channel.fromList(modes)
            .map { run_name ->
                def glob_path = (run_name == 'Control') ? params.precalled_control_vcfs : params.precalled_mosdepth_vcfs
                if (!glob_path) error "Pre-called VCF path is not defined for mode ${run_name}."
                def files = file(glob_path)
                return [ run_name, files ]
            }
    }

    // 5. Cohort Analysis & ACMG Priorities (Runs in parallel across run modes)
    bcftoolsFilterVcfs(vcfs_ch)
    
    if (!params.skip_annotation) {
        SURVIVORMergeCohort(bcftoolsFilterVcfs.out.filtered_data)
        bcftoolsFlattenCohort(SURVIVORMergeCohort.out.merged_vcf)
        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
        
        // Join merged VCF with rare regions
        slice_joined = SURVIVORMergeCohort.out.merged_vcf.join(filterRareCnvs.out.rare_data)
        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
        
        runAnnotSV(filterVcfByRareRegions.out.rare_vcf, annotsv_db_ch)
        
        // Join unannotated rare CSV with AnnotSV TSV output
        acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })

        // 6. CNV Quality Metrics Assessment
        metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
            .join(runAnnotSV.out.annotated_tsv)
            .join(bcftoolsFilterVcfs.out.filtered_data.map { [ it[0], it[2] ] })
            
        metrics_joined = metrics_inputs.combine(target_beds_ch.map { [ it[0], it[1] ] }, by: 0)
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_joined.map { [ it[0], it[1], it[2], it[3] ] }, // [run_name, rare_cnvs_csv, annotated_tsv, filtered_vcfs]
            metrics_joined.map { it[4] } // bed_file
        )

        // 7. Score-separated Evaluation
        eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
        eval_joined = eval_inputs.combine(target_beds_ch.map { [ it[0], it[1] ] }, by: 0)
        
        EVALUATE_ACMG_CSV(
            eval_joined.map { [ it[0], it[1], it[2] ] }, // [ run_name, csv_file, bed_file ]
            file(params.map_file),
            file(params.truth_bed)
        )

        // 8. Dynamic Performance Comparison plotting
        PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
    }
}
