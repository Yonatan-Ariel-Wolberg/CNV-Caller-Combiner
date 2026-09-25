#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// PROCESS: RUN_MOSDEPTH (Base-level regions coverage extraction)
// =====================================================================================
process RUN_MOSDEPTH {
    label 'mosdepth'
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

// =====================================================================================
// PROCESS: BUILD_CONSENSUS (Cohort-wide unfragmented target coverage filtration)
// =====================================================================================
process BUILD_CONSENSUS {
    label 'bedtools'
    cpus 2
    memory '8 GB'
    
    publishDir "${params.outdir}/consensus_bed", mode: 'copy'

    input:
    path usable_beds
    path bed

    output:
    path "ddduk_karyotypic_consensus.bed", emit: karyotypic_bed

    script:
    def sort_flag = params.workflow == 'dragen' ? '-k1,1V -k2,2n' : '-k1,1 -k2,2n'
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
    sort ${sort_flag} ddduk_clean_targeted_consensus.bed | bedtools merge -i stdin -d 1500 > ddduk_karyotypic_consensus.bed
    """
}

// =====================================================================================
// PROCESS: CONVERT_BED_TO_INTERVAL_LIST (GATK conversion helper)
// =====================================================================================
process CONVERT_BED_TO_INTERVAL_LIST {
    label 'gatk'
    publishDir "${params.outdir}/consensus_bed", mode: 'copy'
    
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

// =====================================================================================
// SUB-WORKFLOW: MOSDEPTH_TARGET_GENERATION
// =====================================================================================
workflow MOSDEPTH_TARGET_GENERATION {
    take:
        bam_ch      // Channel of [ sampleId, bam, bai ]
        fasta       // Reference Fasta
        fasta_fai   // Fasta Index
        standard_bed // Sanger design target BED file

    main:
        RUN_MOSDEPTH(bam_ch, fasta, fasta_fai, standard_bed)
        BUILD_CONSENSUS(RUN_MOSDEPTH.out.usable_beds.collect(), standard_bed)

    emit:
        karyotypic_bed = BUILD_CONSENSUS.out.karyotypic_bed
}
