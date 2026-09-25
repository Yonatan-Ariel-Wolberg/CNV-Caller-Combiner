#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BGZIP_SORT_INDEX_VCF; BGZIP_SORT_INDEX_SV_VCF } from '../common/bgzip-sort-index.nf'

// =====================================================================================
// PROCESSES FOR ICAv2 DRAGEN GERMLINE ENRICHMENT
// =====================================================================================

outdir = file(params.outdir, type: 'dir')

// Inspired by a workflow by Reagan Cannell @ https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow

process UPLOAD_CRAM_FILES {
    debug true
    tag "$sampleId"
    label 'icav2_dragen'
    maxForks params.maxUploadForks
    cpus 1
    time '4h'

    errorStrategy { task.exitStatus == 100 ? 'terminate' : 'retry' }
    maxRetries params.uploadRetries

    input:
    tuple val(sampleId), file(cramPair)

    output:
    path "data_upload.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def uploadPath = params.icaUploadPath
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    // Support both CRAM/CRAI and BAM/BAI pairs; pair[0] is the primary file, pair[1] the index
    def (cram_file, crai_file) = cramPair
    """
    #!/bin/bash
    set -euo pipefail

    sleep \$((RANDOM % 15))
    touch data_upload.txt

    clean_json() { sed -n '/^{/,\$p'; }

    get_or_upload_file() {
        local fpath=\$1
        local folder=\$2
        local pid=\$3
        local fname=\$(basename "\$fpath")

        local list_raw=\$(icav2 projectdata list --project-id "\$pid" --parent-folder "\$folder" --file-name "\$fname" --match-mode EXACT -o json 2>/dev/null || true)
        local existing_id=\$(echo "\$list_raw" | clean_json | jq -r --arg n "\$fname" '.items[]? | select(.name == \$n) | .id' | head -n 1 || true)
        local existing_status=\$(echo "\$list_raw" | clean_json | jq -r --arg n "\$fname" '.items[]? | select(.name == \$n) | .status' | head -n 1 || true)

        if [ -n "\$existing_id" ] && [ "\$existing_id" != "null" ]; then
            if [ "\$existing_status" == "PARTIAL" ]; then
                icav2 projectdata delete "\$existing_id" --project-id "\$pid" >/dev/null 2>&1 || true
            else
                local chk=\$(icav2 projectdata get "\$existing_id" --project-id "\$pid" -o json 2>/dev/null || true)
                if ! echo "\$chk" | grep -qE "NotFound|not found|404"; then
                    echo "\$existing_id" | tr -d '[:space:]'
                    return 0
                fi
            fi
        fi

        local real_fpath=\$(realpath "\$fpath")
        local upload_raw=\$(timeout 3600 icav2 projectdata upload "\$real_fpath" "\$folder" --project-id "\$pid" -o json 2>&1 || true)

        local new_id=\$(echo "\$upload_raw" | clean_json | jq -r '.id // empty' 2>/dev/null || true)
        if [ -z "\$new_id" ] || [ "\$new_id" == "null" ]; then
            new_id=\$(echo "\$upload_raw" | grep -o 'fil.[a-zA-Z0-9]*' | head -1 || true)
        fi

        echo "\$new_id" | tr -d '[:space:]'
    }

    CRAM_ID=\$(get_or_upload_file "${cram_file}" "${uploadPath}" "${projectId}")
    CRAI_ID=\$(get_or_upload_file "${crai_file}" "${uploadPath}" "${projectId}")

    if [ -z "\$CRAM_ID" ] || [ -z "\$CRAI_ID" ]; then exit 1; fi

    echo "sampleId:${sampleId}" >> data_upload.txt
    echo "${cramCode}:\$CRAM_ID" >> data_upload.txt
    echo "${craiCode}:\$CRAI_ID" >> data_upload.txt
    """
}

process GET_STATIC_FILES {
    debug true
    tag "Validating Static Files"
    label 'icav2_dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_static.txt", emit: dataFile

    script:
    def refCode = params.referenceAnalysisDataCode
    def refId = params.referenceFileId
    def bedCode = params.targetBedAnalysisDataCode
    def bedId = params.targetBedFileId
    def cramRefCode = params.cramReferenceAnalysisDataCode
    def cramRefId = params.cramReferenceFileId
    """
    #!/bin/bash
    cat ${dataFile} > data_static.txt
    echo "${refCode}:${refId}" >> data_static.txt
    echo "${bedCode}:${bedId}" >> data_static.txt
    echo "${cramRefCode}:${cramRefId}" >> data_static.txt
    """
}

process CHECK_FILE_STATUS {
    debug true
    tag "Checking Folder Status"
    label 'icav2_dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_checked.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def uploadPath = params.icaUploadPath
    def interval = params.fileStatusCheckInterval
    def limit = params.fileStatusCheckLimit
    """
    #!/bin/bash
    set -euo pipefail
    cat ${dataFile} > data_checked.txt

    check_count=0
    while true; do
        ((check_count+=1))
        echo "Attempt \${check_count}/${limit}..."

        # Simplified Check: Ask ICA to list ANY files in the folder that are NOT 'AVAILABLE'
        # This uses 1 API call per loop instead of N, completely preventing rate limits.
        not_ready=\$(icav2 projectdata list --project-id "${projectId}" --parent-folder "${uploadPath}" -o json 2>/dev/null | jq -r '.items[]? | select((.details.status // .status) != "AVAILABLE") | .name' || echo "API_ERROR")

        if [ "\$not_ready" == "API_ERROR" ]; then
            echo "Warning: API hiccup. Retrying..."
        elif [ -z "\$not_ready" ]; then
            echo "✓ All files in folder are AVAILABLE!"
            break
        else
            echo "Files still syncing. Waiting..."
        fi

        if [ \${check_count} -ge ${limit} ]; then
            echo "Timeout reached."
            exit 1
        fi
        sleep ${interval}
    done
    """
}

process START_ANALYSIS_BATCH {
    debug true
    tag "Launch Batch Analysis"
    label 'icav2_dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_analysis.txt", emit: dataFile

    script:
    def projectId = params.projectId
    def pipelineId = params.pipelineId
    def storageSize = params.storageSize
    def userRefPrefix = params.userReference
    def samplesPerNode = params.samplesPerNode ?: 8

    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    def refCode = params.referenceAnalysisDataCode
    def bedCode = params.targetBedAnalysisDataCode
    def cramRefCode = params.cramReferenceAnalysisDataCode
    """
    #!/bin/bash
    set -euo pipefail
    cat ${dataFile} > data_analysis.txt

    # Validate minimum 5 samples required for in-run panel of normals
    sample_count=\$(grep -c "^${cramCode}:" data_analysis.txt || true)
    if [ "\$sample_count" -lt 5 ]; then
        echo "Error: minimum 5 samples required for in-run panel of normals. Found \$sample_count." >&2
        exit 1
    fi
    echo "Sample count: \$sample_count (minimum 5 satisfied for in-run panel of normals)"

    ref_id=\$(grep '^${refCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    bed_id=\$(grep '^${bedCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    cram_ref_id=\$(grep '^${cramRefCode}:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')

    cram_ids_list=\$(grep '^${cramCode}:' ${dataFile} | cut -f2- -d: | tr -d '\\r' | sort -u | paste -sd "," -)
    crai_ids_list=\$(grep '^${craiCode}:' ${dataFile} | cut -f2- -d: | tr -d '\\r' | sort -u | paste -sd "," -)

    user_ref="${userRefPrefix}-batch"

    analysis_response_raw=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --project-id ${projectId} \
        --user-reference "\${user_ref}" \
        --storage-size ${storageSize} \
        --input ${refCode}:\${ref_id} \
        --input ${bedCode}:\${bed_id} \
        --input ${cramRefCode}:\${cram_ref_id} \
        --input ${cramCode}:\${cram_ids_list} \
        --input ${craiCode}:\${crai_ids_list} \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:true \
        --parameters output_format:CRAMv3.1 \
        --parameters enable_duplicate_marking:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:GVCF \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:true \
        --parameters cnv_enable_in_run_pon:true \
        --parameters cnv_enable_gcbias_correction:true \
        --parameters enable_sv:true \
        --parameters qc_coverage_report:cov_report \
        --parameters qc_coverage_filters:"'mapq<20,bq<20'" \
        --parameters qc_detect_contamination:true \
        --parameters enable_variant_annotation:true \
        --parameters samples_per_node:${samplesPerNode} 2>&1 || true)

    analysis_id=\$(echo "\$analysis_response_raw" | sed -n '/^{/,\$p' | jq -r ".id" 2>/dev/null || echo "")

    if [ -z "\${analysis_id}" ] || [ "\${analysis_id}" == "null" ]; then exit 1; fi

    echo "analysisId:\${analysis_id}" >> data_analysis.txt
    echo "analysisRef:\${user_ref}" >> data_analysis.txt
    echo "batchMode:true" >> data_analysis.txt
    """
}

process CHECK_ANALYSIS_STATUS {
    debug true
    tag "Monitor Analysis"
    label 'icav2_dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_status.txt", emit: dataFile

    script:
    def interval = params.analysisStatusCheckInterval
    def limit = params.analysisStatusCheckLimit
    def projectId = params.projectId
    """
    #!/bin/bash
    set -euo pipefail
    cat ${dataFile} > data_status.txt
    analysis_id=\$(grep '^analysisId:' ${dataFile} | cut -f2- -d: | tr -d '[:space:]')
    check_count=0

    while true; do
        ((check_count+=1))
        status_json=\$(icav2 projectanalyses get "\${analysis_id}" --project-id "${projectId}" -o json)
        status=\$(echo "\${status_json}" | jq -r ".status" | tr -d '[:space:]')

        if [[ "\${status}" == "SUCCEEDED" ]]; then
            echo "analysisStatus:SUCCEEDED" >> data_status.txt; break
        elif [[ "\${status}" == "FAILED" || "\${status}" == "FAILED_FINAL" || "\${status}" == "ABORTED" ]]; then
            echo "analysisStatus:FAILED" >> data_status.txt; break
        fi

        if [ \${check_count} -gt ${limit} ]; then
            echo "analysisStatus:TIMEOUT" >> data_status.txt; break
        fi
        sleep ${interval}
    done
    """
}

process DOWNLOAD_ANALYSIS_OUTPUT {
    debug true
    tag "Download Output"
    label 'icav2_dragen'
    cpus 1

    input:
    path(dataFile)

    output:
    path "data_final.txt", emit: dataFile

    script:
    def downloadPath = params.localDownloadPath ?: "${outdir}/out_DRAGEN"
    def projectId = params.projectId
    """
    #!/bin/bash
    set -euo pipefail
    cat ${dataFile} > data_final.txt
    status=\$(grep '^analysisStatus:' ${dataFile} | cut -f2- -d:)
    analysis_id=\$(grep '^analysisId:' ${dataFile} | cut -f2- -d:)

    if [ "\$status" == "SUCCEEDED" ]; then
        output_json=\$(icav2 projectanalyses output "\${analysis_id}" --project-id "${projectId}" -o json)
        folder_id=\$(echo "\${output_json}" | jq -r ".items[0].data[0].dataId")

        if [ -n "\${folder_id}" ] && [ "\${folder_id}" != "null" ]; then
            icav2 projectdata download "\${folder_id}" "${downloadPath}" --project-id "${projectId}"
            echo "outputFolderId:\${folder_id}" >> data_final.txt
            echo "deleteData:true" >> data_final.txt
        else
            echo "deleteData:false" >> data_final.txt
        fi
    else
        echo "deleteData:false" >> data_final.txt
    fi
    """
}

process DELETE_DATA {
    debug true
    tag "Cleanup"
    label 'icav2_dragen'
    cpus 1
    errorStrategy 'ignore'

    input:
    path(dataFile)

    output:
    stdout

    script:
    def projectId = params.projectId
    def cramCode = params.cramAnalysisDataCode
    def craiCode = params.cramIndexAnalysisDataCode
    """
    #!/bin/bash
    delete_flag=\$(grep '^deleteData:' ${dataFile} | cut -f2- -d:)
    if [ "\$delete_flag" == "true" ]; then
        for id in \$(grep '^${cramCode}:' ${dataFile} | cut -f2- -d: | tr -d '\\r'); do
            icav2 projectdata delete "\$id" --project-id "${projectId}" || true
        done
        for id in \$(grep '^${craiCode}:' ${dataFile} | cut -f2- -d: | tr -d '\\r'); do
            icav2 projectdata delete "\$id" --project-id "${projectId}" || true
        done
        folder_id=\$(grep '^outputFolderId:' ${dataFile} | cut -f2- -d: | tr -d '\\r')
        icav2 projectdata delete "\${folder_id}" --project-id "${projectId}" || true
    fi
    """
}

process ADD_DRAGEN_TOOL_ANNOTATION {
    debug true
    tag "Add TOOL=DRAGEN annotation"
    label 'bcftools'
    publishDir "${outdir}/out_DRAGEN/annotated_vcfs", mode: 'copy', overwrite: true
    cpus 1

    input:
    path(dataFile)

    output:
    path "*_DRAGEN.annotated.vcf.gz", emit: annotated_vcfs, optional: true

    script:
    def downloadPath = params.localDownloadPath ?: "${outdir}/out_DRAGEN"
    """
    #!/bin/bash
    set -euo pipefail

    echo "Annotating downloaded DRAGEN VCF files with TOOL=DRAGEN in INFO field..."

    printf '##INFO=<ID=TOOL,Number=1,Type=String,Description="Calling tool">\\n' > extra_header.txt

    found=0
    while IFS= read -r vcf; do
        found=1
        sample_name=\$(basename "\$vcf" | sed 's/\\.cnv\\.vcf\\.gz\$//' | sed 's/\\.vcf\\.gz\$//' | sed 's/\\.vcf\$//')

        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' "\$vcf" | \\
            awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, "DRAGEN"}' | \\
            sort -k1,1 -k2,2n | \\
            bgzip -c > "\${sample_name}_tool_annot.bed.gz"
        tabix -p bed "\${sample_name}_tool_annot.bed.gz"

        bcftools annotate \\
            -a "\${sample_name}_tool_annot.bed.gz" \\
            -c CHROM,FROM,TO,INFO/TOOL \\
            -h extra_header.txt \\
            "\$vcf" \\
            -O z -o "\${sample_name}_DRAGEN.annotated.vcf.gz"

        rm -f "\${sample_name}_tool_annot.bed.gz" "\${sample_name}_tool_annot.bed.gz.tbi"
    done < <(find ${downloadPath} -type f -name "*.cnv.vcf.gz" 2>/dev/null || true)

    if [ \$found -eq 0 ]; then
        echo "No VCF files found in ${downloadPath}. Skipping annotation."
    else
        echo "TOOL=DRAGEN annotation complete."
    fi
    """
}

// =====================================================================================
// SUB-WORKFLOW TO CHAIN THE PROCESSES TOGETHER
// =====================================================================================

workflow DRAGEN {
    take:
    cram_ch  // channel: tuple(sampleId, [cram_file, crai_file])

    main:
    // Validate minimum 5 samples required for in-run panel of normals
    cram_ch
        .collect()
        .subscribe { list ->
            if (list.size() < 5) {
                error("Error: minimum 5 samples required for in-run panel of normals. Found ${list.size()}.")
            }
        }

    validated_cram_ch = cram_ch

    // Step 1: Upload CRAM and CRAI files to ICA for each sample
    UPLOAD_CRAM_FILES(validated_cram_ch)

    // Step 2: Combine all per-sample upload data into a single file
    combined_upload_ch = UPLOAD_CRAM_FILES.out.dataFile
        .collectFile(name: 'combined_data.txt')

    // Step 3: Append static reference file IDs to the combined data file
    GET_STATIC_FILES(combined_upload_ch)

    // Step 4: Verify all uploaded files are available in ICA
    CHECK_FILE_STATUS(GET_STATIC_FILES.out.dataFile)

    // Step 5: Launch the DRAGEN batch analysis on ICA
    START_ANALYSIS_BATCH(CHECK_FILE_STATUS.out.dataFile)

    // Step 6: Monitor the analysis until completion
    CHECK_ANALYSIS_STATUS(START_ANALYSIS_BATCH.out.dataFile)

    // Step 7: Download the analysis output from ICA
    DOWNLOAD_ANALYSIS_OUTPUT(CHECK_ANALYSIS_STATUS.out.dataFile)

    // Step 8: Clean up temporary CRAM files and output folder from ICA
    DELETE_DATA(DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile)

    // Step 9: Annotate downloaded VCF files with TOOL=DRAGEN in INFO field
    ADD_DRAGEN_TOOL_ANNOTATION(DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile)

    // Step 10: Sort and index each annotated VCF via bgzip-sort-index.nf
    BGZIP_SORT_INDEX_VCF(ADD_DRAGEN_TOOL_ANNOTATION.out.annotated_vcfs.flatten(), 'DRAGEN', 'out_DRAGEN')

    emit:
    result         = DOWNLOAD_ANALYSIS_OUTPUT.out.dataFile
    annotated_vcfs = ADD_DRAGEN_TOOL_ANNOTATION.out.annotated_vcfs
    sorted_vcf     = BGZIP_SORT_INDEX_VCF.out.sorted_vcf
    sorted_vcf_index = BGZIP_SORT_INDEX_VCF.out.sorted_vcf_index
}
