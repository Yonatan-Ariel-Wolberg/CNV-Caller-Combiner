#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { NORMALISE_CNV_QUALITY_SCORES } from '../common/modules-common.nf'

// =============================================================================
// INDELIBLE MODULE
// =============================================================================

// REQUIRED FILES
params.ref            = ''
params.priors         = ''
params.indelible_conf = ''
ref               = file(params.ref            ?: 'NO_FILE', type: 'file')
priors            = file(params.priors         ?: 'NO_FILE', type: 'file')
indelible_conf    = file(params.indelible_conf ?: 'NO_FILE', type: 'file')
outdir            = file(params.outdir, type: 'dir')

// 1. THE FETCH COMMAND EXTRACTS THE READS FROM THE BAM FILE, IT TAKES 2 ARGUMENTS:
process RUN_FETCH {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.sc_reads"), path(bam), path(bai), emit: sc_reads
    
    script:
    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s
    indelible.py fetch --config ${indelible_conf} --i ${bam} --o ${bam}.sc_reads
    """
}

// 2. THE AGGREGATE MERGES INFORMATION ACROSS READS TOWARDS A POSITION-LEVEL VIEW OF THE DATA:
process RUN_AGGREGATE {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(sc_read), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.counts"), emit: counts
    
    script:
    """
    indelible.py aggregate --i ${sc_read} --b ${bam} --o ${bam}.counts --r ${ref} --config ${indelible_conf}
    """
}

// 3. THE SCORE COMMAND SCORES POSITIONS BASED ON THE READ INFORMATION AND SEQUENCE CONTEXT:
process RUN_SCORE {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(count)

    output:
    tuple val(sample), path("${count}.scored"), emit: scores
    path("${count}.scored"), emit: database_in
    
    script:
    """
    indelible.py score --i ${count} --o ${count}.scored --config ${indelible_conf}
    """
}

// 4. THE DATABASE COMMAND GENERATES THE ALLELE FREQUENCY AND BREAKPOINT DATABASE REQUIRED FOR THE NEXT STEP – ANNOTATE
process RUN_DATABASE {
    tag { "Indel_DB" }
    label 'indelible'
    cpus 6
    publishDir "${outdir}/out_INDELIBLE/database", mode: 'copy', overwrite: true

    input:
    path(score)

    output:
    path("InDelible_db.tsv"), emit: indel_database
    
    script:
    """
    ls *.scored > scores.txt
    indelible.py database --f scores.txt --o InDelible_db.tsv --r ${ref} --priors ${priors} --config ${indelible_conf} --tb ${task.cpus}
    """
}

// 5. THE ANNOTATE COMMAND ENRICHES THE RESULT WITH GENE/EXON ANNOTATIONS AND MERGES THE DATABASE RESULTS WITH THE POSITION FILE:
process RUN_ANNOTATE {
    tag { "${score}" }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/annotations", mode: 'copy', overwrite: true

    input:
    path(database)
    tuple val(sample), path(score)

    output:
    tuple val(sample), path("${score}.annotated"), emit: annotated
    
    script:
    """
    indelible.py annotate --i ${score} --o ${score}.annotated --d ${database} --config ${indelible_conf}
    """
}

// 6. ONE CAN THEN LOOK FOR DE NOVO MUTATION EVENTS USING THE DENOVO COMMAND:
// TRIO
process RUN_DENOVO_TRIO {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_trio", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo

    script:
    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// MOM
process RUN_DENOVO_MOM {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_mom", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_mom

    script:
    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// DAD
process RUN_DENOVO_DAD {
    label 'indelible'
    tag { sample }

    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_dad", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_dad

    script:
    """
    indelible.py denovo --c ${annotation} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// Filter the results from Indelible
process FILTER_INDELIBLE {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/filtered", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(annotated)
    
    output:
    tuple val(sample), path("${sample}.annotated.filtered.tsv"), emit: filtered_cnvs
    
    script:
    """
    awk '{ if ((\$39 < 2) && (\$40 < 2)) { print } }' ${annotated} > ${sample}.annotated.filtered.tsv
    """
}

// Process to convert INDELIBLE filtered TSV to VCF using the Python script
process CONVERT_INDELIBLE_TO_VCF {
    tag "${sample_id}"
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/vcfs", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(filtered_tsv)

    output:
    path("${sample_id}_INDELIBLE_output.vcf"), emit: vcfs

    script:
    def fai_arg = params.get('fai', null) ? "--fai_file ${file(params.fai)}" : ""
    """
    indelible_tsv_to_vcf.py \\
        --input_file ${filtered_tsv} \\
        --sample_id ${sample_id} \\
        --output_dir . \\
        ${fai_arg}
    """
}

// Process to compress, sort and index each INDELIBLE VCF.
// TOOL=INDELIBLE is already written into the INFO field by CONVERT_INDELIBLE_TO_VCF,
// so no additional bcftools-annotate step is needed here.  A second attempt to add
// the same ##INFO header would produce a non-standard duplicate header line.
process BGZIP_SORT_INDEX_VCF {
    tag "${vcf_file.simpleName}"
    label 'bcftools'
    publishDir "${outdir}/out_INDELIBLE/vcfs", mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path("*.sorted.vcf.gz"),     emit: sorted_vcf
    path("*.sorted.vcf.gz.tbi"), emit: sorted_vcf_index

    script:
    def sample_name = vcf_file.simpleName
    def sorted_gz   = "${sample_name}.sorted.vcf.gz"
    """
    bcftools sort ${vcf_file} -o ${sorted_gz} -O z
    tabix -p vcf ${sorted_gz}
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow INDELIBLE {
    take:
    crams_ch        // channel: tuple(sample, bam, bai)
    cram_trios_ch   // channel: tuple(sample, child_bam, child_bai, mom_bam, mom_bai, dad_bam, dad_bai)
    cram_mom_ch     // channel: tuple(sample, child_bam, child_bai, mom_bam, mom_bai)
    cram_dad_ch     // channel: tuple(sample, child_bam, child_bai, dad_bam, dad_bai)

    main:
    // Step 1: Extract split/clipped reads from BAM files
    RUN_FETCH(crams_ch)

    // Step 2: Aggregate read information to a position-level view
    RUN_AGGREGATE(RUN_FETCH.out.sc_reads)

    // Step 3: Score positions based on read information and sequence context
    RUN_SCORE(RUN_AGGREGATE.out.counts)

    // Step 4: Generate the allele frequency and breakpoint database
    RUN_DATABASE(RUN_SCORE.out.database_in.collect())

    // Step 5: Annotate positions with gene/exon information
    // Use .first() to convert the single-item database queue channel into a
    // value channel that broadcasts to every scored-file in the scores channel.
    // Without .first(), Nextflow's default zip behaviour would pair the single
    // database with only the FIRST scored file, leaving all other samples
    // un-annotated.
    RUN_ANNOTATE(RUN_DATABASE.out.indel_database.first(), RUN_SCORE.out.scores)

    // Step 6: Identify de novo mutations in complete trios
    RUN_DENOVO_TRIO(cram_trios_ch.join(RUN_ANNOTATE.out.annotated))

    // Step 7: Identify de novo mutations with mother only
    RUN_DENOVO_MOM(cram_mom_ch.join(RUN_ANNOTATE.out.annotated))

    // Step 8: Identify de novo mutations with father only
    RUN_DENOVO_DAD(cram_dad_ch.join(RUN_ANNOTATE.out.annotated))

    // Step 9: Filter annotated results by confidence thresholds
    FILTER_INDELIBLE(RUN_ANNOTATE.out.annotated)

    // Step 10: Convert filtered TSV to VCF
    CONVERT_INDELIBLE_TO_VCF(FILTER_INDELIBLE.out.filtered_cnvs)

    // Step 11: bgzip-compress, coordinate-sort, and tabix-index each per-sample VCF
    BGZIP_SORT_INDEX_VCF(CONVERT_INDELIBLE_TO_VCF.out.vcfs)

    // Step 12: Normalise quality scores to a common scale
    NORMALISE_CNV_QUALITY_SCORES(BGZIP_SORT_INDEX_VCF.out.sorted_vcf.flatten(), 'INDELIBLE', 'out_INDELIBLE')

    emit:
    filtered_cnvs        = FILTER_INDELIBLE.out.filtered_cnvs
    indelible_denovo     = RUN_DENOVO_TRIO.out.indelible_denovo
    indelible_denovo_mom = RUN_DENOVO_MOM.out.indelible_denovo_mom
    indelible_denovo_dad = RUN_DENOVO_DAD.out.indelible_denovo_dad
    sorted_vcf           = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index     = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
    normalised_vcf       = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf
    normalised_vcf_index = NORMALISE_CNV_QUALITY_SCORES.out.normalised_vcf_index
}
