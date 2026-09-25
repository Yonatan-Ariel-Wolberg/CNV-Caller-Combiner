nextflow.enable.dsl = 2

// Default parameters (configured in nextflow.config, fallbacks here)
params.workflow       = params.get('workflow', 'qc')
params.input_dir      = params.get('input_dir', "$launchDir")
params.output_dir     = params.get('output_dir', "./qc_results")
params.ref_genome     = params.get('ref_genome', null)
params.target_bed     = params.get('target_bed', null) // Path to your exome kit's target/probe BED file (highly recommended)
params.mapq_threshold = params.get('mapq_threshold', 20)
params.threads        = params.get('threads', 4)

// Early Validation
def cram_files = file("${params.input_dir}/*.cram")
if (cram_files.size() > 0 && !params.ref_genome) {
    error "ERROR: CRAM files found in '${params.input_dir}', but no reference genome was provided. Please supply it using --ref_genome <path/to/fasta>"
}

/* * PROCESS 1: Handle BAM files (Samtools metrics)
 */
process SAMTOOLS_QC_BAM {
    tag "$bam_file.name"
    label 'samtools'
    publishDir "${params.output_dir}/samtools_metrics", mode: 'copy'

    input:
    path bam_file
    path target_bed

    output:
    path "*.flagstat"          , emit: flagstat
    path "*.stats"             , emit: stats
    path "*.idxstats"          , emit: idxstats
    path "*.coverage_summary"   , emit: coverage
    path "*.per_probe_depth.txt", emit: probe_depth, optional: true
    path "*.low_quality.bam"    , emit: low_qual_bam
    tuple path(bam_file), path("*.bai"), emit: bam_and_index 

    script:
    def prefix = bam_file.baseName
    def bed_arg = params.target_bed ? "-b ${target_bed}" : ""
    """
    set -euo pipefail

    # 1. Quickcheck: Verify BAM integrity
    if ! samtools quickcheck ${bam_file}; then
        echo "[ERROR] ${bam_file} failed the integrity check." >&2
        exit 1
    fi

    # Create standard BAI index on the fly (required for idxstats, bedcov, and QualiMap)
    samtools index -b -@ ${task.cpus} ${bam_file}

    # 2. Flagstat: Overview of mapping quality
    samtools flagstat -@ ${task.cpus} ${bam_file} > ${prefix}.flagstat

    # 3. Stats: Comprehensive quality metrics
    samtools stats -@ ${task.cpus} ${bam_file} > ${prefix}.stats

    # 4. Idxstats: Map reads per chromosome
    samtools idxstats ${bam_file} > ${prefix}.idxstats

    # 5. On-Target Coverage Calculation (Only looks at target regions if BED is provided)
    samtools depth ${bed_arg} ${bam_file} | awk '
    {
        sum += \$3;
        if(\$3 >= 1)  c1++;
        if(\$3 >= 10) c10++;
        if(\$3 >= 20) c20++;
        if(\$3 >= 30) c30++;
        total_bases++
    }
    END {
        if (total_bases > 0) {
            print "Sample: ${prefix}"
            print "Target-Region Average Depth: " sum/total_bases "x"
            print "Target Breadth of Coverage (>=1x): " (c1/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=10x): " (c10/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=20x): " (c20/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=30x): " (c30/total_bases)*100 "%"
        } else {
            print "No on-target coverage recorded. Did you supply the correct target BED?"
        }
    }' > ${prefix}.coverage_summary

    # 6. Per-Probe/Interval Depth Calculation (using samtools bedcov)
    if [ -f "${target_bed}" ]; then
        echo -e "Chr\\tStart\\tEnd\\tProbe_ID\\tTotal_Per_Base_Depth\\tAverage_Probe_Depth" > ${prefix}.per_probe_depth.txt
        samtools bedcov ${target_bed} ${bam_file} | awk -F'\\t' '
        {
            len = \$3 - \$2;
            probe_id = (NF >= 5) ? \$(NF-1) : \$1 ":" \$2 "-" \$3;
            tot_depth = \$NF;
            avg_depth = (len > 0) ? tot_depth / len : 0;
            print \$1 "\\t" \$2 "\\t" \$3 "\\t" probe_id "\\t" tot_depth "\\t" avg_depth
        }' >> ${prefix}.per_probe_depth.txt
    fi

    # 7. View: Extract low quality/unmapped reads
    samtools view -h -@ ${task.cpus} ${bam_file} | \\
    awk -v thresh="${params.mapq_threshold}" -F'\\t' '
    /^@/ { print; next }
    {
        if (\$5 < thresh || int(\$2/4)%2 == 1 || int(\$2/512)%2 == 1) {
            print \$0
        }
    }' | samtools view -bS -@ ${task.cpus} - > ${prefix}.low_quality.bam
    """
}

/* * PROCESS 2: Handle CRAM files (Samtools metrics)
 */
process SAMTOOLS_QC_CRAM {
    tag "$cram_file.name"
    label 'samtools'
    publishDir "${params.output_dir}/samtools_metrics", mode: 'copy'

    input:
    path cram_file
    path ref_fasta
    path target_bed

    output:
    path "*.flagstat"          , emit: flagstat
    path "*.stats"             , emit: stats
    path "*.idxstats"          , emit: idxstats
    path "*.coverage_summary"   , emit: coverage
    path "*.per_probe_depth.txt", emit: probe_depth, optional: true
    path "*.low_quality.bam"    , emit: low_qual_bam
    tuple path(cram_file), path("*.crai"), emit: cram_and_index

    script:
    def prefix = cram_file.baseName
    def bed_arg = params.target_bed ? "-b ${target_bed}" : ""
    """
    set -euo pipefail

    if ! samtools quickcheck ${cram_file}; then
        echo "[ERROR] ${cram_file} failed the integrity check." >&2
        exit 1
    fi

    samtools index -@ ${task.cpus} ${cram_file}

    samtools flagstat -@ ${task.cpus} ${cram_file} > ${prefix}.flagstat
    samtools stats -@ ${task.cpus} --reference ${ref_fasta} ${cram_file} > ${prefix}.stats
    samtools idxstats ${cram_file} > ${prefix}.idxstats

    samtools depth ${bed_arg} --reference ${ref_fasta} ${cram_file} | awk '
    {
        sum += \$3;
        if(\$3 >= 1)  c1++;
        if(\$3 >= 10) c10++;
        if(\$3 >= 20) c20++;
        if(\$3 >= 30) c30++;
        total_bases++
    }
    END {
        if (total_bases > 0) {
            print "Sample: ${prefix}"
            print "Target-Region Average Depth: " sum/total_bases "x"
            print "Target Breadth of Coverage (>=1x): " (c1/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=10x): " (c10/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=20x): " (c20/total_bases)*100 "%"
            print "Target Breadth of Coverage (>=30x): " (c30/total_bases)*100 "%"
        } else {
            print "No on-target coverage recorded. Did you supply the correct target BED?"
        }
    }' > ${prefix}.coverage_summary

    if [ -f "${target_bed}" ]; then
        echo -e "Chr\\tStart\\tEnd\\tProbe_ID\\tTotal_Per_Base_Depth\\tAverage_Probe_Depth" > ${prefix}.per_probe_depth.txt
        samtools bedcov --reference ${ref_fasta} ${target_bed} ${cram_file} | awk -F'\\t' '
        {
            len = \$3 - \$2;
            probe_id = (NF >= 5) ? \$(NF-1) : \$1 ":" \$2 "-" \$3;
            tot_depth = \$NF;
            avg_depth = (len > 0) ? tot_depth / len : 0;
            print \$1 "\\t" \$2 "\\t" \$3 "\\t" probe_id "\\t" tot_depth "\\t" avg_depth
        }' >> ${prefix}.per_probe_depth.txt
    fi

    samtools view -h -@ ${task.cpus} --reference ${ref_fasta} ${cram_file} | \\
    awk -v thresh="${params.mapq_threshold}" -F'\\t' '
    /^@/ { print; next }
    {
        if (\$5 < thresh || int(\$2/4)%2 == 1 || int(\$2/512)%2 == 1) {
            print \$0
        }
    }' | samtools view -bS -@ ${task.cpus} - > ${prefix}.low_quality.bam
    """
}

/* * PROCESS 3: QualiMap BAMQC (On-Target Metrics)
 */
process QUALIMAP {
    tag "$alignment_file.name"
    label 'qualimap'
    publishDir "${params.output_dir}/qualimap_results", mode: 'copy'

    input:
    tuple path(alignment_file), path(index_file)
    path target_bed

    output:
    path "${prefix}_qualimap", emit: qc_dir

    script:
    prefix = alignment_file.baseName
    // If BED is provided, we tell Qualimap to use it as a target design region
    def gff_arg = params.target_bed ? "-gff ${target_bed}" : ""
    """
    export JAVA_OPTS="-Xmx6G"

    qualimap bamqc \\
        -bam ${alignment_file} \\
        -outdir ${prefix}_qualimap \\
        -outformat HTML \\
        --java-mem-size=6G \\
        ${gff_arg}
    """
}

/* * PROCESS 4: MultiQC
 */
process MULTIQC {
    label 'multiqc'
    publishDir "${params.output_dir}/multiqc_report", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data

    script:
    """
    multiqc .
    """
}

/* * PROCESS 5: Infer Sex and Ploidy
 */
process INFER_SEX_AND_PLOIDY {
    label 'pysam'
    publishDir "${params.output_dir ?: params.outdir + '/qc_results'}/sex_and_ploidy", mode: 'copy', overwrite: true

    input:
    path idxstats_files

    output:
    path "cohort_sex_and_ploidy.tsv", emit: sex_table
    path "clamms_sexinfo.tsv", emit: clamms_sexinfo

    script:
    """
    python3 ${projectDir}/bin/infer_sex_and_ploidy.py \\
        --idxstats-files ${idxstats_files} \\
        --output-tsv cohort_sex_and_ploidy.tsv \\
        --clamms-sexinfo clamms_sexinfo.tsv
    """
}

/* * PROCESS 6: Reference-Free Contamination Estimator
 */
process ESTIMATE_CONTAMINATION {
    label 'pysam'
    publishDir "${params.output_dir ?: params.outdir + '/qc_results'}/contamination", mode: 'copy', overwrite: true

    input:
    path alignments
    path target_bed
    path ref_fasta

    output:
    path "cohort_contamination_estimates.tsv", emit: contamination_report

    script:
    def ref_arg = (ref_fasta && ref_fasta.name != 'NO_REF_GENOME' && ref_fasta.exists()) ? "--ref-fasta ${ref_fasta}" : ""
    """
    python3 ${projectDir}/bin/estimate_contamination.py \\
        --bam-files ${alignments} \\
        --target-bed ${target_bed} \\
        ${ref_arg} \\
        --output-tsv cohort_contamination_estimates.tsv
    """
}

/* * WORKFLOW
 */
workflow {
    input_ch = Channel.fromPath("${params.input_dir}/*.{bam,cram}")

    input_ch.branch {
        bams: it.name.endsWith('.bam')
        crams: it.name.endsWith('.cram')
    }.set { branched_files }

    // Channel for optional target bed file
    bed_ch = params.target_bed ? Channel.fromPath(params.target_bed).collect() : Channel.value([])
    ref_ch = params.ref_genome ? Channel.fromPath(params.ref_genome).collect() : Channel.value([])

    // Execute BAM pipeline
    SAMTOOLS_QC_BAM(branched_files.bams, bed_ch)
    
    // Execute CRAM pipeline
    if (params.ref_genome) {
        SAMTOOLS_QC_CRAM(branched_files.crams, ref_ch, bed_ch)
    }

    // QualiMap only supports BAM format (HTSJDK limitation; crashes on CRAM)
    QUALIMAP(SAMTOOLS_QC_BAM.out.bam_and_index, bed_ch)

    // Collect all raw reports for MultiQC
    bam_metrics = SAMTOOLS_QC_BAM.out.flagstat.mix(SAMTOOLS_QC_BAM.out.stats).mix(SAMTOOLS_QC_BAM.out.idxstats)
    cram_metrics = params.ref_genome ? SAMTOOLS_QC_CRAM.out.flagstat.mix(SAMTOOLS_QC_CRAM.out.stats).mix(SAMTOOLS_QC_CRAM.out.idxstats) : Channel.empty()

    all_qc_metrics = bam_metrics.mix(cram_metrics).mix(QUALIMAP.out.qc_dir.ifEmpty([])).collect()

    MULTIQC(all_qc_metrics)

    // Execute automated sex & ploidy inference from idxstats
    idxstats_ch = SAMTOOLS_QC_BAM.out.idxstats.mix(params.ref_genome ? SAMTOOLS_QC_CRAM.out.idxstats : Channel.empty()).collect()
    INFER_SEX_AND_PLOIDY(idxstats_ch)

    // Execute automated reference-free contamination estimator across all BAMs and CRAMs
    all_alignments = branched_files.bams.mix(branched_files.crams).collect()
    ref_contam_ch = params.ref_genome ? ref_ch : Channel.value(file("NO_REF_GENOME"))
    ESTIMATE_CONTAMINATION(all_alignments, bed_ch, ref_contam_ch)
}