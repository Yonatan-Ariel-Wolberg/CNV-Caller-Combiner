#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.transform.Field

// =====================================================================================
// MODULE INCLUDES
// =====================================================================================
include { INDELIBLE } from './modules/callers/modules-indelible.nf'
include { CANOES } from './modules/callers/modules-canoes.nf'
include { XHMM } from './modules/callers/modules-xhmm.nf'
include { CLAMMS } from './modules/callers/modules-clamms.nf'
include { DRAGEN } from './modules/callers/modules-icav2-dragen.nf'
include { CNVKIT } from './modules/callers/modules-cnvkit.nf'
include { GATK_GCNV } from './modules/callers/modules-gatk-gcnv.nf'
include { SURVIVOR } from './modules/sv-mergers/modules-survivor.nf'
include { TRUVARI } from './modules/sv-mergers/modules-truvari.nf'
include { FEATURE_EXTRACTION } from './modules/ml/modules-feature-extraction.nf'
include { NORMALISE } from './modules/normalise/modules-normalise.nf'
include { TRAIN } from './modules/ml/modules-train.nf'
include { EVALUATE } from './modules/evaluate/modules-evaluate.nf'
include { SAMTOOLS_QC_BAM; SAMTOOLS_QC_CRAM; QUALIMAP; MULTIQC } from './target-aware-exome-sequencing-quality-assess.nf'

// =====================================================================================
// GLOBAL SETUP
// =====================================================================================
workflow_mode = params.workflow
params.test_size = params.get('test_size', -1)

@Field def CALLER_DIR_PARAMS = [
    'canoes_dir',
    'clamms_dir',
    'xhmm_dir',
    'cnvkit_dir',
    'gcnv_dir',
    'dragen_dir',
    'indelible_dir',
]

// Valid caller names accepted by normalise_cnv_caller_quality_scores.py.
// Note: GATK gCNV is referred to as 'GATK' (not 'GCNV') by the normalise script.
@Field def VALID_NORMALISE_CALLERS = ['CANOES', 'CLAMMS', 'XHMM', 'GATK', 'CNVKIT', 'DRAGEN', 'INDELIBLE']
@Field def BYTES_PER_GIB = 1024D * 1024D * 1024D
@Field def DEFAULT_MIN_FREE_GB = 5
@Field def MAX_VCF_SCHEMA_VALIDATION_FILES = 5

@Field def REQUIRED_PARAMS_BY_WORKFLOW = [
    // Required params are aligned to the workflow-specific params/*.json templates.
    'indelible': ['outdir', 'crams', 'ref', 'priors', 'indelible_conf', 'fai'],
    'canoes': ['outdir', 'samplesheet_bams', 'ref', 'fai', 'probes', 'canoes_batch_size'],
    'xhmm': ['outdir', 'samplesheet_bams', 'ref', 'probes', 'xhmm_conf', 'xhmm_batch_size'],
    'clamms': ['outdir', 'samplesheet_bams', 'ref', 'fai', 'probes', 'interval_list', 'mappability', 'special_reg', 'sexinfo'],
    'dragen': [
        'outdir', 'projectId', 'pipelineId', 'userReference', 'storageSize',
        'referenceAnalysisDataCode', 'targetBedAnalysisDataCode', 'cramAnalysisDataCode',
        'cramIndexAnalysisDataCode', 'cramReferenceAnalysisDataCode', 'referenceFileId',
        'cramReferenceFileId', 'targetBedFileId', 'localDownloadPath',
        'cramFilePairsUploadPath', 'icaUploadPath', 'maxUploadForks', 'uploadRetries',
        'fileStatusCheckInterval', 'fileStatusCheckLimit', 'analysisStatusCheckInterval',
        'analysisStatusCheckLimit'
    ],
    'cnvkit': ['outdir', 'bams', 'fasta', 'targets', 'refflat', 'test_size', 'test_list'],
    'gcnv': ['outdir', 'samples_path', 'fasta', 'fai', 'dict', 'exome_targets', 'bin_length', 'padding', 'is_wgs', 'scatter_count'],
    'survivor': ['outdir', 'canoes_dir', 'clamms_dir', 'xhmm_dir', 'cnvkit_dir', 'gcnv_dir', 'dragen_dir', 'indelible_dir'],
    'truvari': ['outdir', 'canoes_dir', 'clamms_dir', 'xhmm_dir', 'cnvkit_dir', 'gcnv_dir', 'dragen_dir', 'indelible_dir'],
    'survivor_with_features': [
        'outdir', 'canoes_dir', 'clamms_dir', 'xhmm_dir', 'cnvkit_dir', 'gcnv_dir',
        'dragen_dir', 'indelible_dir', 'canoes_norm_dir', 'clamms_norm_dir',
        'xhmm_norm_dir', 'cnvkit_norm_dir', 'gcnv_norm_dir', 'dragen_norm_dir',
        'indelible_norm_dir', 'merger_mode'
    ],
    'truvari_with_features': [
        'outdir', 'canoes_dir', 'clamms_dir', 'xhmm_dir', 'cnvkit_dir', 'gcnv_dir',
        'dragen_dir', 'indelible_dir', 'canoes_norm_dir', 'clamms_norm_dir',
        'xhmm_norm_dir', 'cnvkit_norm_dir', 'gcnv_norm_dir', 'dragen_norm_dir',
        'indelible_norm_dir', 'merger_mode'
    ],
    'normalise': ['outdir', 'vcf_dir', 'caller'],
    'feature_extraction': ['outdir', 'merged_vcf_dir', 'merger_mode'],
    'train': ['outdir', 'features_dir', 'truth_labels'],
    'evaluate': ['outdir', 'vcf_dir', 'caller', 'truth_bed', 'probes_bed'],
    'full': ['truth_labels'],
    'qc': ['input_dir'],
]

def is_param_set(param_name) {
    def value = params.get(param_name, null)
    if (value == null) return false
    if (value instanceof CharSequence) return value.toString().trim().length() > 0
    return true
}

def validate_required_params(workflow_name) {
    if (REQUIRED_PARAMS_BY_WORKFLOW.containsKey(workflow_name)) {
        def missing = REQUIRED_PARAMS_BY_WORKFLOW[workflow_name].findAll { !is_param_set(it) }
        if (!missing.isEmpty()) {
            exit 1, "Error: Missing required parameter(s) for --workflow ${workflow_name}: ${missing.collect { '--' + it }.join(', ')}"
        }
    }

    if (['survivor', 'truvari', 'survivor_with_features', 'truvari_with_features'].contains(workflow_name)) {
        def configured_caller_dirs = CALLER_DIR_PARAMS.findAll { is_param_set(it) }
        if (configured_caller_dirs.size() < 2) {
            exit 1, "Error: --workflow ${workflow_name} requires at least TWO caller VCF directories. Configure at least 2 of: ${CALLER_DIR_PARAMS.collect { '--' + it }.join(', ')}"
        }
    }

    if (workflow_name == 'full') {
        def configured_caller_count = 0

        if (['samplesheet_bams', 'fai'].every { is_param_set(it) }) {
            // One config block enables 3 callers: CANOES, CLAMMS, XHMM
            configured_caller_count += 3
        }
        if (['bams', 'fasta', 'targets', 'refflat'].every { is_param_set(it) }) {
            configured_caller_count += 1
        }
        if (['samples_path', 'fasta', 'fai', 'dict', 'exome_targets'].every { is_param_set(it) }) {
            configured_caller_count += 1
        }
        if (is_param_set('cramFilePairsUploadPath')) {
            // DRAGEN
            configured_caller_count += 1
        }
        if (is_param_set('crams')) {
            // INDELIBLE
            configured_caller_count += 1
        }

        if (configured_caller_count < 2) {
            exit 1, "Error: --workflow full requires at least 2 configured CNV callers out of 7 (CANOES, CLAMMS, XHMM, CNVKIT, GCNV, DRAGEN, INDELIBLE). Configure at least one of: (--samplesheet_bams + --fai) [enables 3 callers], (--bams + --fasta + --targets + --refflat), (--samples_path + --fasta + --fai + --dict + --exome_targets), (--cramFilePairsUploadPath), (--crams)."
        }
    }

    if (workflow_name == 'normalise' && is_param_set('caller')) {
        def caller_val = params.caller.toString().trim().toUpperCase()
        if (!VALID_NORMALISE_CALLERS.contains(caller_val)) {
            exit 1, "Error: --caller '${params.caller}' is not a supported caller for --workflow normalise. Valid values are: ${VALID_NORMALISE_CALLERS.collect { "'${it}'" }.join(', ')}"
        }
    }

    if (workflow_name == 'gcnv') {
        if (is_param_set('bin_length')) {
            def bin_length_val = params.bin_length.toString().toInteger()
            if (bin_length_val < 0) {
                exit 1, "Error: --bin_length must be a non-negative integer for --workflow gcnv. Got: ${params.bin_length}"
            }
        }
        if (is_param_set('padding')) {
            def padding_val = params.padding.toString().toInteger()
            if (padding_val < 0) {
                exit 1, "Error: --padding must be a non-negative integer for --workflow gcnv. Got: ${params.padding}"
            }
        }
        if (is_param_set('scatter_count')) {
            def scatter_count_val = params.scatter_count.toString().toInteger()
            if (scatter_count_val <= 0) {
                exit 1, "Error: --scatter_count must be a positive integer for --workflow gcnv. Got: ${params.scatter_count}"
            }
        }
    }

    if (workflow_name == 'canoes' && is_param_set('canoes_batch_size')) {
        def canoes_batch_val = params.canoes_batch_size.toString().toInteger()
        if (canoes_batch_val <= 0) {
            exit 1, "Error: --canoes_batch_size must be a positive integer for --workflow canoes. Got: ${params.canoes_batch_size}"
        }
    }

    if (workflow_name == 'xhmm' && is_param_set('xhmm_batch_size')) {
        def xhmm_batch_val = params.xhmm_batch_size.toString().toInteger()
        if (xhmm_batch_val <= 0) {
            exit 1, "Error: --xhmm_batch_size must be a positive integer for --workflow xhmm. Got: ${params.xhmm_batch_size}"
        }
    }

    if (workflow_name == 'cnvkit' && is_param_set('test_size')) {
        def test_size_val = params.test_size.toString().toInteger()
        if (test_size_val != -1 && test_size_val <= 0) {
            exit 1, "Error: --test_size must be -1 (disabled) or a positive integer for --workflow cnvkit. Got: ${params.test_size}"
        }
    }

    if (['feature_extraction', 'survivor_with_features', 'truvari_with_features'].contains(workflow_name) && is_param_set('merger_mode')) {
        def valid_merger_modes = ['survivor', 'truvari']
        def merger_mode_val = params.merger_mode.toString().trim().toLowerCase()
        if (!valid_merger_modes.contains(merger_mode_val)) {
            exit 1, "Error: --merger_mode '${params.merger_mode}' is not valid for --workflow ${workflow_name}. Valid values are: 'survivor', 'truvari'"
        }
    }
}

def to_java_file(def path_value) {
    if (path_value == null) {
        return null
    }
    return (path_value instanceof File) ? path_value : new File(path_value.toString())
}

def resolve_input_path(def raw_path, File base_dir = null) {
    if (raw_path == null) {
        return null
    }
    def p = raw_path.toString().trim()
    if (!p) {
        return null
    }
    def f = new File(p)
    if (f.isAbsolute() || base_dir == null) {
        return f
    }
    return new File(base_dir, p)
}

def fail_validation(String message) {
    exit 1, message
}

def require_readable_file(def path_value, String param_name) {
    if (!path_value) {
        fail_validation("Error: missing required path for ${param_name}")
    }
    def p = to_java_file(path_value)
    if (!p.exists() || !p.isFile() || !p.canRead()) {
        fail_validation("Error: unreadable file for ${param_name}: ${path_value}")
    }
    return p
}

def require_readable_dir(def path_value, String param_name) {
    if (!path_value) {
        fail_validation("Error: missing required directory for ${param_name}")
    }
    def p = to_java_file(path_value)
    if (!p.exists() || !p.isDirectory() || !p.canRead()) {
        fail_validation("Error: unreadable directory for ${param_name}: ${path_value}")
    }
    return p
}

def require_writable_dir(def path_value, String label) {
    if (!path_value) {
        fail_validation("Error: output path not writable (${label}): ${path_value}")
    }
    def p = to_java_file(path_value)
    if (!p.exists() && !p.mkdirs()) {
        fail_validation("Error: output path not writable (${label}): ${path_value}")
    }
    if (!p.exists() || !p.isDirectory() || !p.canWrite()) {
        fail_validation("Error: output path not writable (${label}): ${path_value}")
    }
    return p
}

def validate_truth_labels_schema(def truth_labels_path) {
    def truth_labels_file = require_readable_file(truth_labels_path, 'truth_labels')
    def header_line = null
    truth_labels_file.withReader { reader ->
        header_line = reader.readLine()
    }

    if (header_line == null || header_line.trim().isEmpty()) {
        fail_validation("Error: malformed truth_labels '${truth_labels_path}': file is empty or missing header")
    }

    def required_cols = ['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label']
    def header_cols = header_line.split('\t', -1).collect { it.trim() }
    def missing = required_cols.findAll { !header_cols.contains(it) }
    if (!missing.isEmpty()) {
        fail_validation("Error: malformed truth_labels '${truth_labels_path}': missing required columns ${missing.join(', ')}")
    }
}

def validate_samplesheet_and_bam_indices(def samplesheet_path) {
    def samplesheet_file = require_readable_file(samplesheet_path, 'samplesheet_bams')
    def header_line = null
    samplesheet_file.withReader { reader ->
        header_line = reader.readLine()
    }
    if (header_line == null || header_line.trim().isEmpty()) {
        fail_validation("Error: malformed samplesheet '${samplesheet_path}': file is empty or missing header")
    }

    def header_cols = header_line.split('\t', -1).collect { it.trim() }
    if (!header_cols.contains('SampleID') || !header_cols.contains('BAM')) {
        fail_validation("Error: malformed samplesheet '${samplesheet_path}': expected tab-separated columns 'SampleID' and 'BAM'")
    }

    def sample_idx = header_cols.indexOf('SampleID')
    def bam_idx = header_cols.indexOf('BAM')
    def seen_sample_ids = [] as Set
    def row_count = 0
    samplesheet_file.eachLine { line, line_number ->
        if (line_number == 1 || line.trim().isEmpty()) {
            return
        }
        def fields = line.split('\t', -1)
        if (fields.size() <= Math.max(sample_idx, bam_idx)) {
            fail_validation("Error: malformed samplesheet '${samplesheet_path}': line ${line_number} has missing columns")
        }

        def sample_id = fields[sample_idx].trim()
        def bam_raw = fields[bam_idx].trim()
        if (!sample_id || !bam_raw) {
            fail_validation("Error: malformed samplesheet '${samplesheet_path}': line ${line_number} has empty SampleID/BAM value")
        }
        if (seen_sample_ids.contains(sample_id)) {
            fail_validation("Error: malformed samplesheet '${samplesheet_path}': duplicate SampleID '${sample_id}'")
        }
        seen_sample_ids.add(sample_id)

        def bam_path = resolve_input_path(bam_raw, samplesheet_file.parentFile)
        if (bam_path == null || !bam_path.exists() || !bam_path.isFile() || !bam_path.canRead()) {
            fail_validation("Error: malformed samplesheet '${samplesheet_path}': unreadable BAM path '${bam_raw}' for sample '${sample_id}'")
        }

        def index_candidates = []
        if (bam_path.name.toLowerCase().endsWith('.bam')) {
            index_candidates << new File(bam_path.toString() + '.bai')
            index_candidates << new File(bam_path.parentFile, bam_path.name.replaceAll(/(?i)\.bam$/, '.bai'))
        } else if (bam_path.name.toLowerCase().endsWith('.cram')) {
            index_candidates << new File(bam_path.toString() + '.crai')
            index_candidates << new File(bam_path.parentFile, bam_path.name.replaceAll(/(?i)\.cram$/, '.crai'))
        }
        if (!index_candidates.isEmpty()) {
            def has_readable_index = index_candidates.any { it.exists() && it.isFile() && it.canRead() }
            if (!has_readable_index) {
                fail_validation("Error: missing BAM index for sample '${sample_id}' in samplesheet '${samplesheet_path}'")
            }
        }
        row_count++
    }

    if (row_count == 0) {
        fail_validation("Error: malformed samplesheet '${samplesheet_path}': no data rows found")
    }
}

def validate_full_reference_assets_consistency() {
    ['fasta', 'fai', 'dict', 'targets', 'refflat', 'exome_targets'].each { param_name ->
        if (is_param_set(param_name)) {
            require_readable_file(params[param_name], param_name)
        }
    }

    if (is_param_set('samples_path') && (!is_param_set('fasta') || !is_param_set('fai') || !is_param_set('dict') || !is_param_set('exome_targets'))) {
        fail_validation("Error: inconsistent reference assets for --workflow full: gCNV input requires --fasta, --fai, --dict, and --exome_targets together with --samples_path")
    }

    if (is_param_set('bams') && (!is_param_set('fasta') || !is_param_set('targets') || !is_param_set('refflat'))) {
        fail_validation("Error: inconsistent reference assets for --workflow full: CNVkit input requires --fasta, --targets, and --refflat together with --bams")
    }
}

def validate_icav2_runtime_assets() {
    if (is_param_set('icav2_container')) {
        def container_file = require_readable_file(params.icav2_container, 'icav2_container')
        if (!container_file.name.toLowerCase().endsWith('.sif')) {
            fail_validation("Error: ICAv2 asset unreadable or invalid: --icav2_container must point to a .sif file (${params.icav2_container})")
        }
    }

    if (is_param_set('ica_credentials_dir')) {
        require_readable_dir(params.ica_credentials_dir, 'ica_credentials_dir')
    }

    if (is_param_set('localDownloadPath')) {
        require_writable_dir(params.localDownloadPath, 'localDownloadPath')
    }
}

def open_maybe_gzip_reader(File f) {
    if (f.name.toLowerCase().endsWith('.gz')) {
        FileInputStream fis = new FileInputStream(f)
        try {
            def gis = new java.util.zip.GZIPInputStream(fis)
            return new BufferedReader(new InputStreamReader(gis))
        } catch (Throwable t) {
            fis.close()
            throw t
        }
    }
    return new BufferedReader(new FileReader(f))
}

def validate_single_vcf_schema(File vcf_file, String context_label) {
    def saw_header = false
    def valid_header = false
    def reader = open_maybe_gzip_reader(vcf_file)
    try {
        String line
        while ((line = reader.readLine()) != null) {
            if (line.startsWith('#CHROM')) {
                saw_header = true
                def fields = line.split('\t', -1)
                valid_header = (fields.size() >= 8 &&
                    fields[0] == '#CHROM' &&
                    fields[1] == 'POS' &&
                    fields[2] == 'ID' &&
                    fields[3] == 'REF' &&
                    fields[4] == 'ALT' &&
                    fields[5] == 'QUAL' &&
                    fields[6] == 'FILTER' &&
                    fields[7] == 'INFO')
                break
            }
        }
    } finally {
        reader.close()
    }

    if (!saw_header || !valid_header) {
        fail_validation("Error: VCF schema incompatibility in ${context_label}: ${vcf_file}")
    }
}

def validate_vcf_schema_preconditions() {
    if (['normalise', 'evaluate'].contains(workflow_mode) && is_param_set('vcf_dir')) {
        def vcf_dir = require_readable_dir(params.vcf_dir, 'vcf_dir')
        def vcfs = vcf_dir.listFiles()?.findAll {
            it.isFile() && (it.name.toLowerCase().endsWith('.vcf') || it.name.toLowerCase().endsWith('.vcf.gz'))
        } ?: []
        vcfs.take(MAX_VCF_SCHEMA_VALIDATION_FILES).each { vcf_file ->
            validate_single_vcf_schema(vcf_file, "--workflow ${workflow_mode} --vcf_dir")
        }
    }

    if (['feature_extraction', 'survivor_with_features', 'truvari_with_features'].contains(workflow_mode) && is_param_set('merged_vcf_dir')) {
        def merged_vcf_dir = require_readable_dir(params.merged_vcf_dir, 'merged_vcf_dir')
        def suffix = (params.get('merger_mode', 'survivor').toString().trim().toLowerCase() == 'truvari') ? '_truvari_merged' : '_survivor_union'
        def merged_vcfs = merged_vcf_dir.listFiles()?.findAll {
            it.isFile() &&
            it.name.contains(suffix) &&
            (it.name.toLowerCase().endsWith('.vcf') || it.name.toLowerCase().endsWith('.vcf.gz'))
        } ?: []
        merged_vcfs.take(MAX_VCF_SCHEMA_VALIDATION_FILES).each { vcf_file ->
            validate_single_vcf_schema(vcf_file, "--workflow ${workflow_mode} --merged_vcf_dir")
        }
    }
}

def validate_io_permissions_and_disk(def output_dir_path, Integer min_free_gb = DEFAULT_MIN_FREE_GB) {
    if (output_dir_path) {
        def output_dir = require_writable_dir(output_dir_path, 'outdir')
        def free_gb = (output_dir.usableSpace / BYTES_PER_GIB)
        if (free_gb < min_free_gb) {
            fail_validation("Error: insufficient disk space in output path '${output_dir_path}'. Required at least ${min_free_gb} GB, found ${String.format('%.2f', free_gb)} GB.")
        }
    }
}

def validate_runtime_inputs(String workflow_name) {
    validate_io_permissions_and_disk(params.get('outdir', null))

    if (['canoes', 'xhmm', 'clamms', 'full'].contains(workflow_name) && is_param_set('samplesheet_bams')) {
        validate_samplesheet_and_bam_indices(params.samplesheet_bams)
    }
    if (!params.skip_training && ['train', 'full'].contains(workflow_name) && is_param_set('truth_labels')) {
        validate_truth_labels_schema(params.truth_labels)
    }
    if (workflow_name == 'full') {
        validate_full_reference_assets_consistency()
    }
    if (['dragen', 'full'].contains(workflow_name)) {
        validate_icav2_runtime_assets()
    }
    if (['normalise', 'evaluate', 'feature_extraction', 'survivor_with_features', 'truvari_with_features'].contains(workflow_name)) {
        validate_vcf_schema_preconditions()
    }
    if (workflow_name == 'qc') {
        validate_io_permissions_and_disk(params.get('output_dir', null))
        require_readable_dir(params.input_dir, 'input_dir')
    }
}

def is_validation_only_mode() {
    // `-stub-run` should validate config/inputs without launching process execution
    // (and therefore without container resolution or image pulls).
    workflow.stubRun
}

def extract_sample_id_from_vcf(f) {
    // Extract bare sample ID from per-caller VCF names used across modules,
    // including DRAGEN's ${sample_id}.cnv.vcf(.gz) convention.
    // Ordering matters: strip caller suffix and VCF extensions first, then trim
    // downstream ".sorted" / ".normalised" suffixes from intermediate filenames.
    return f.name
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
        .replaceAll(/(?i)\.cnv\.vcf(\.gz)?$/, '')
        .replaceAll(/(?i)\.vcf(\.gz)?$/, '')
        .replaceAll(/(?i)\.sorted$/, '')
        .replaceAll(/(?i)\.normalised$/, '')
}

def infer_caller_from_vcf(f) {
    def n = f.name.toUpperCase()
    if (n.contains('_CANOES'))    return 'canoes'
    if (n.contains('_CLAMMS'))    return 'clamms'
    if (n.contains('_XHMM'))      return 'xhmm'
    if (n.contains('_CNVKIT'))    return 'cnvkit'
    if (n.contains('_GCNV'))      return 'gatk_gcnv'
    if (n.contains('_DRAGEN') || n.contains('.CNV.')) return 'dragen'
    if (n.contains('_INDELIBLE')) return 'indelible'
    error "Could not infer caller from VCF filename: ${f.name}"
}

def build_tool_vcfs_str(vcfs) {
    // Build feature_extraction --tool_vcfs argument format:
    // caller1=/path/a.vcf.gz,caller2=/path/b.vcf.gz,...
    return vcfs
        .collect { vcf -> "${infer_caller_from_vcf(vcf)}=${vcf}" }
        .join(',')
}

// =====================================================================================
// HELPER FUNCTION: GATHER VCFS FOR CONSENSUS MODULES
// =====================================================================================
def gather_vcfs() {
    def ch = Channel.empty()
    def dir_count = 0
    
    CALLER_DIR_PARAMS.each { caller_dir ->
        if (params.get(caller_dir, false)) {
            ch = ch.mix(Channel.fromPath(params[caller_dir] + "/*.vcf*").map { f -> [extract_sample_id_from_vcf(f), f] })
            dir_count++
        }
    }
    
    if (dir_count < 2) {
        exit 1, "Error: You must provide VCF directories for at least TWO different callers (e.g., --canoes_dir and --clamms_dir) to run consensus modules."
    }
    
    return ch
}

def group_caller_vcfs(vcf_ch) {
    return vcf_ch
        .groupTuple()
        .filter { sample_id, vcfs -> vcfs.size() >= 2 }
        .ifEmpty { error "No samples with VCFs from 2 or more callers were found. Check that the caller VCF directories contain files with matching sample IDs across at least 2 callers." }
}

// =====================================================================================
// SUB-WORKFLOW DEFINITIONS
// =====================================================================================

workflow RUN_INDELIBLE {
    take: 
        crams
        cram_trios
        cram_mom
        cram_dad
    main:
        INDELIBLE(crams, cram_trios, cram_mom, cram_dad)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                INDELIBLE.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'INDELIBLE'
            )
        }
}

workflow RUN_CANOES {
    take: 
        bams
        chroms
        fai
    main:
        bam_list = bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        CANOES(bam_list, fai, Channel.from(chroms))
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CANOES.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CANOES'
            )
        }
}

workflow RUN_XHMM {
    take: 
        bams
    main:
        XHMM(bams)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                XHMM.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'XHMM'
            )
        }
}

workflow RUN_CLAMMS {
    take: 
        bams
        fai
    main:
        sample_list = bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
        CLAMMS(bams, fai, sample_list)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CLAMMS.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CLAMMS'
            )
        }
}

workflow RUN_DRAGEN {
    take: 
        cramPairs
    main:
        DRAGEN(cramPairs)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                DRAGEN.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'DRAGEN'
            )
        }
}

workflow RUN_CNVKIT {
    take: 
        bams
        fasta
        targets
        refflat
    main:
        // Guard against empty BAM input so workflow fails fast with a clear error.
        bams_nonempty = bams.ifEmpty { error "RUN_CNVKIT received no BAMs. Check --bams input/glob." }
        CNVKIT(bams_nonempty, fasta, targets, refflat, bams_nonempty.first().map { it[1] })
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                CNVKIT.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'CNVKIT'
            )
        }
}

workflow RUN_GCNV {
    take: 
        bams
        fasta
        fai
        dict
        targets
    main:
        GATK_GCNV(bams, fasta, fai, dict, targets)
        if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
            EVALUATE(
                GATK_GCNV.out.normalised_vcf.flatten(),
                file(params.truth_bed),
                file(params.probes_bed),
                'GCNV'
            )
        }
}

workflow RUN_EVALUATE {
    take:
        vcf_ch
        truth_bed
        probes_bed
        caller_name
    main:
        EVALUATE(vcf_ch, truth_bed, probes_bed, caller_name)
}

workflow RUN_SURVIVOR {
    take: 
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        SURVIVOR(grouped_vcfs)
    emit:
        union_vcf = SURVIVOR.out.union_vcf
}

workflow RUN_TRUVARI {
    take: 
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        TRUVARI(grouped_vcfs)
    emit:
        merged_vcf = TRUVARI.out.merged_vcf
}

workflow RUN_FEATURE_EXTRACTION {
    take:
        // Channel of tuples:
        // [ sample_id, merged_vcf, collapsed_vcf, tool_vcfs_str, merger_mode,
        //   bam_file, reference_fasta, bed_file, mappability_file, indelible_counts ]
        // collapsed_vcf should be [] when not available (SURVIVOR mode or absent).
        feature_inputs_ch
    main:
        FEATURE_EXTRACTION(feature_inputs_ch)
}

// Helper closure: build the feature_inputs tuple for a single merged VCF.
// collapsed_vcf_f: a file path or [] (empty list)
def build_feature_inputs(sample_id, merged_vcf, collapsed_vcf_f) {
    def parts = []
    if (params.get('canoes_norm_dir',    false)) parts << "canoes=${params.canoes_norm_dir}/${sample_id}_CANOES.normalised.vcf.gz"
    if (params.get('clamms_norm_dir',    false)) parts << "clamms=${params.clamms_norm_dir}/${sample_id}_CLAMMS.normalised.vcf.gz"
    if (params.get('xhmm_norm_dir',      false)) parts << "xhmm=${params.xhmm_norm_dir}/${sample_id}_XHMM.normalised.vcf.gz"
    if (params.get('cnvkit_norm_dir',    false)) parts << "cnvkit=${params.cnvkit_norm_dir}/${sample_id}_CNVKIT.normalised.vcf.gz"
    if (params.get('gcnv_norm_dir',      false)) parts << "gatk_gcnv=${params.gcnv_norm_dir}/${sample_id}_GCNV.normalised.vcf.gz"
    if (params.get('dragen_norm_dir',    false)) parts << "dragen=${params.dragen_norm_dir}/${sample_id}_DRAGEN.normalised.vcf.gz"
    if (params.get('indelible_norm_dir', false)) parts << "indelible=${params.indelible_norm_dir}/${sample_id}_INDELIBLE.normalised.vcf.gz"
    def tool_vcfs_str = parts.join(',')

    def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
    def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
    def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
    def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []
    def indelible_f = params.get('indelible_counts', false) ? file(params.indelible_counts) : []
    def mode        = params.get('merger_mode', 'survivor')

    return [sample_id, merged_vcf, collapsed_vcf_f, tool_vcfs_str, mode,
            bam_f, fasta_f, bed_f, map_f, indelible_f]
}

workflow RUN_NORMALISE {
    take:
        // Channel of tuples: [ sample_name, vcf_path, caller_name ]
        vcf_ch
    main:
        NORMALISE(vcf_ch)
}

workflow RUN_TRAIN {
    take:
        // Channel of individual *_features.tsv file paths (one per sample).
        features_tsv_ch
        // Single-element channel containing the path to the truth-labels TSV.
        // Required columns: sample_id, chrom, start, end, cnv_type, truth_label
        truth_labels_ch
        // Optional single-element channel containing probes BED.
        // Used for probe-overlap fallback when feature/truth CNV coordinates differ.
        probes_bed_ch
    main:
        TRAIN(features_tsv_ch, truth_labels_ch, probes_bed_ch)
}

workflow RUN_SURVIVOR_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        SURVIVOR(grouped_vcfs)
        feature_inputs_ch = SURVIVOR.out.union_vcf
            .map { sample_id, merged_vcf -> build_feature_inputs(sample_id, merged_vcf, []) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_TRUVARI_WITH_FEATURES {
    take:
        vcf_ch
    main:
        grouped_vcfs = group_caller_vcfs(vcf_ch)
        TRUVARI(grouped_vcfs)
        feature_inputs_ch = TRUVARI.out.merged_vcf
            .join(TRUVARI.out.collapsed_vcf)
            .map { sample_id, merged_vcf, collapsed_vcf -> build_feature_inputs(sample_id, merged_vcf, collapsed_vcf) }
        FEATURE_EXTRACTION(feature_inputs_ch)
}

workflow RUN_QC {
    take:
        input_dir
        output_dir
        ref_genome
        target_bed
        mapq_threshold
        threads
    main:
        input_ch = Channel.fromPath("${input_dir}/*.{bam,cram}")

        input_ch.branch {
            bams: it.name.endsWith('.bam')
            crams: it.name.endsWith('.cram')
        }.set { branched_files }

        // Channel for optional target bed file
        bed_ch = target_bed ? Channel.fromPath(target_bed).collect() : Channel.value([])

        // Execute BAM pipeline
        SAMTOOLS_QC_BAM(branched_files.bams, bed_ch)
        
        // Execute CRAM pipeline
        if (ref_genome) {
            ref_ch = Channel.fromPath(ref_genome).collect()
            SAMTOOLS_QC_CRAM(branched_files.crams, ref_ch, bed_ch)
        }

        // QualiMap
        bam_pairs = SAMTOOLS_QC_BAM.out.bam_and_index
        cram_pairs = ref_genome ? SAMTOOLS_QC_CRAM.out.cram_and_index : Channel.empty()
        qualimap_input = bam_pairs.mix(cram_pairs)

        QUALIMAP(qualimap_input, bed_ch)

        // Collect all raw reports for MultiQC
        bam_metrics = SAMTOOLS_QC_BAM.out.flagstat.mix(SAMTOOLS_QC_BAM.out.stats).mix(SAMTOOLS_QC_BAM.out.idxstats)
        cram_metrics = ref_genome ? SAMTOOLS_QC_CRAM.out.flagstat.mix(SAMTOOLS_QC_CRAM.out.stats).mix(SAMTOOLS_QC_CRAM.out.idxstats) : Channel.empty()

        all_qc_metrics = bam_metrics.mix(cram_metrics).mix(QUALIMAP.out.qc_dir).collect()

        MULTIQC(all_qc_metrics)
}

// =====================================================================================
// WORKFLOW SWITCH - FULLY INDEPENDENT EXECUTION BLOCKS
// =====================================================================================

workflow {
    validate_required_params(workflow_mode)
    validate_runtime_inputs(workflow_mode)
    if (is_validation_only_mode()) {
        log.info "Validation-only execution mode detected (stub-run); exiting before process execution."
        return
    }
    switch (workflow_mode) {

        case['qc']:
            RUN_QC(
                params.input_dir,
                params.output_dir,
                params.ref_genome,
                params.target_bed,
                params.mapq_threshold,
                params.threads
            )
            break
        
        case['indelible']:
            Channel.fromFilePairs([params.crams + '/*{.cram,.cram.crai}'])
                .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }
                .filter { it -> it[1] =~ '_01_1' }
                .ifEmpty { error "indelible workflow: no proband CRAMs found in '${params.crams}'. Check that files follow the *_01_1.cram naming convention." }
                .set { ch_crams }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: 6)
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] ] }
                .set { ch_cram_trios }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                .filter { it -> it[1].size() == 4 && it[1][3] =~ '02_2.cram' }
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                .set { ch_cram_mom }

            Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                .filter { it -> it[1].size() == 4 && it[1][3] =~ '03_3.cram' }
                .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                .set { ch_cram_dad }

            RUN_INDELIBLE(ch_crams, ch_cram_trios, ch_cram_mom, ch_cram_dad)
            break

        case['canoes']:
            chroms = (1..22).toList().collect { 'chr' + "${it}" }
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "canoes workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_CANOES(ch_bams, chroms, file(params.fai))
            break

        case['xhmm']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "xhmm workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_XHMM(ch_bams)
            break

        case['clamms']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "clamms workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_CLAMMS(ch_bams, file(params.fai))
            break

        case['dragen']:
            // Support both CRAM/CRAI and BAM/BAI file pairs while ignoring orphan index files.
            Channel.fromPath(params.cramFilePairsUploadPath, checkIfExists: true)
                .filter { f -> f.name.endsWith('.cram') || f.name.endsWith('.bam') }
                .map { primary ->
                    def sampleId = primary.name.replaceAll(/\.(cram|bam)$/, '')
                    def idxCandidates = primary.name.endsWith('.cram')
                        ? [file("${primary}.crai"), file(primary.toString().replaceAll(/\.cram$/, '.crai'))]
                        : [file("${primary}.bai"), file("${primary}.bam.bai"), file(primary.toString().replaceAll(/\.bam$/, '.bai'))]
                    def indexFile = idxCandidates.find { it.exists() }
                    [sampleId, [primary, indexFile]]
                }
                .filter { sampleId, pair -> pair[1] != null }
                .ifEmpty { error "dragen workflow: no CRAM or BAM file pairs found in '${params.cramFilePairsUploadPath}'. Check that the path contains indexed CRAM/BAM files." }
                .set { ch_cramPairs }
            RUN_DRAGEN(ch_cramPairs)
            break

        case['cnvkit']:
            def glob_pattern = params.bams.contains('{') ? params.bams : params.bams.replaceAll(/\.bam$/, '.{bam,bai,bam.bai}')
            Channel.fromPath(glob_pattern)
                .map { it -> 
                    def sample_id = it.name.replaceAll(/\.bam\.bai$/, '').replaceAll(/\.bam$/, '').replaceAll(/\.bai$/, '').replaceAll(/\.md$/, '')
                    [sample_id, it] 
                }
                .groupTuple(size: 2)
                .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                .filter { id, bam, bai -> 
                    def tl = params.get('test_list', 'all')
                    (tl != "all" && tl != "") ? tl.tokenize(',').contains(id) : true 
                }
                .set { ch_bams_filtered }
            
            ch_bams_final = (params.get('test_size', -1).toInteger() > 0) ? ch_bams_filtered.take(params.get('test_size', -1).toInteger()) : ch_bams_filtered
            RUN_CNVKIT(ch_bams_final, file(params.fasta), file(params.targets), file(params.refflat))
            break

        case['gcnv']:
            Channel.fromPath(params.samples_path)
                .filter { it.name.endsWith('.bam') || it.name.endsWith('.cram') }
                .map { it -> 
                    def index_fallback = [file("${it}.bai"), file("${it}.bam.bai"), file(it.toString().replaceAll(/\.bam$/, '.bai'))]
                        .find { f -> f.exists() }
                    def index_file = index_fallback ? index_fallback : (it.name.endsWith('.bam') ? file("${it}.bai") : file("${it}.crai"))
                    return [ it.baseName, it, index_file ] 
                }
                .ifEmpty { error "gcnv workflow: no sample files found matching '${params.samples_path}'. Check that the glob matches existing BAM or CRAM files." }
                .set { ch_bams }
            RUN_GCNV(ch_bams, file(params.fasta), file(params.fai), file(params.dict), file(params.exome_targets))
            break

        case['survivor']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR(ch_vcfs)
            
            if (params.get('annotsv_db', false)) {
                annotsv_db_ch = Channel.fromPath(params.annotsv_db)
            } else {
                DOWNLOAD_ANNOTSV_DB()
                annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
            }
            
            def probes_bed_file = params.get('probes_bed', '/dataG/ddd/data/resources/canoes/probes_sanger.bed')
            RUN_CONSENSUS_ANNOTATION('Consensus_Survivor', RUN_SURVIVOR.out.union_vcf.map { sample_id, vcf -> vcf }, annotsv_db_ch, file(probes_bed_file))
            break

        case['truvari']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI(ch_vcfs)
            
            if (params.get('annotsv_db', false)) {
                annotsv_db_ch = Channel.fromPath(params.annotsv_db)
            } else {
                DOWNLOAD_ANNOTSV_DB()
                annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
            }
            
            def probes_bed_file = params.get('probes_bed', '/dataG/ddd/data/resources/canoes/probes_sanger.bed')
            RUN_CONSENSUS_ANNOTATION('Consensus_Truvari', RUN_TRUVARI.out.merged_vcf.map { sample_id, vcf -> vcf }, annotsv_db_ch, file(probes_bed_file))
            break

        case['survivor_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR_WITH_FEATURES(ch_vcfs)
            break

        case['truvari_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI_WITH_FEATURES(ch_vcfs)
            break

        case['normalise']:
            // Normalise quality scores for pre-existing caller VCFs without running
            // the full caller pipeline.
            // Required: --vcf_dir  (directory containing raw .vcf / .vcf.gz files)
            //           --caller   (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN | INDELIBLE)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f ->
                    def sample_name = f.name
                        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
                        .replaceAll(/\.(sorted\.)?vcf(\.gz)?$/, '')
                    [sample_name, f, params.caller.toUpperCase()]
                }
                .ifEmpty { error "normalise workflow: no VCF files found in '${params.vcf_dir}'. Check that the directory contains .vcf or .vcf.gz files." }
                .set { ch_vcfs }
            RUN_NORMALISE(ch_vcfs)
            break

        case['feature_extraction']:
            // Run feature extraction on an already-produced merged VCF.
            // Required: --merged_vcf_dir  (directory containing <sample_id>_survivor_union.vcf
            //                              or <sample_id>_truvari_merged.vcf files)
            // Optional annotation params (any subset may be omitted):
            //   --canoes_norm_dir, --clamms_norm_dir, --xhmm_norm_dir, --cnvkit_norm_dir,
            //   --gcnv_norm_dir, --dragen_norm_dir, --indelible_norm_dir
            //   --bam_file, --reference_fasta, --bed_file, --mappability_file,
            //   --indelible_counts, --merger_mode (default: survivor)
            if (params.get('merger_mode', 'survivor') == 'truvari') {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_truvari_merged.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_truvari.*/, '').replaceAll(/(?i)\.vcf(\.gz)?$/, '')
                        def collapsed_f = file("${f.parent}/${sample_id}_truvari_collapsed.vcf")
                        def collapsed_vcf_f = collapsed_f.exists() ? collapsed_f : []
                        build_feature_inputs(sample_id, f, collapsed_vcf_f)
                    }
                    .ifEmpty { error "feature_extraction workflow (truvari mode): no merged VCF files found in '${params.merged_vcf_dir}'. Expected files matching *_truvari_merged.vcf or *_truvari_merged.vcf.gz." }
                    .set { ch_feature_inputs }
            } else {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_survivor_union.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_survivor.*/, '').replaceAll(/(?i)\.vcf(\.gz)?$/, '')
                        build_feature_inputs(sample_id, f, [])
                    }
                    .ifEmpty { error "feature_extraction workflow (survivor mode): no merged VCF files found in '${params.merged_vcf_dir}'. Expected files matching *_survivor_union.vcf or *_survivor_union.vcf.gz." }
                    .set { ch_feature_inputs }
            }
            RUN_FEATURE_EXTRACTION(ch_feature_inputs)
            break

        case['train']:
            if (params.skip_training) {
                log.info "Skipping model training as requested (skip_training = true)."
                break
            }
            // Train an XGBoost classifier on CNV feature matrices produced by
            // the feature_extraction workflow.
            // Required: --features_dir   (directory containing *_features.tsv files
            //                             for positive and negative instances)
            //           --truth_labels   (csv containing EGAN_ID, CNV_ID, and True/False label)
            // Optional: --probes_bed     (bed containing probes used to filter target regions)
            Channel.value(file(params.truth_labels))
                .set { ch_truth }
            def probesBedParam = params.get('probes_bed', false)
            Channel.value(probesBedParam ? file(probesBedParam) : [])
                .set { ch_probes }
            TRAIN(Channel.fromPath(params.features_dir + '/*_features.tsv').collect(), ch_truth, ch_probes)
            break

        case['evaluate']:
            // Evaluate an existing cohort-wide CNV callset VCF directory against the gold truthset
            // Required: --vcf_dir    (directory containing sample VCFs to evaluate)
            //           --caller     (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN | INDELIBLE)
            //           --truth_bed  (gold standard truth set BED file)
            //           --probes_bed (target Sanger probes BED file)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f -> [extract_sample_id_from_vcf(f), f] }
                .ifEmpty { error "evaluate workflow: no VCF files found in '${params.vcf_dir}'. Check that the directory contains .vcf or .vcf.gz files." }
                .set { ch_vcfs }
            RUN_EVALUATE(ch_vcfs, file(params.truth_bed), file(params.probes_bed), params.caller.toUpperCase())
            break

        case['full']:
            def caller_vcf_channels = []
            def mergerMode = params.get('merger_mode', 'survivor')

            if (params.get('samplesheet_bams', false) && params.get('fai', false)) {
                Channel.fromPath(params.samplesheet_bams)
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                    .set { ch_samplesheet_bams }
                chroms = (1..22).toList().collect { 'chr' + "${it}" }
                bam_list = ch_samplesheet_bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
                CANOES(bam_list, file(params.fai), Channel.from(chroms))
                sample_list = ch_samplesheet_bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
                CLAMMS(ch_samplesheet_bams, file(params.fai), sample_list)
                XHMM(ch_samplesheet_bams)
                caller_vcf_channels << CANOES.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
                caller_vcf_channels << CLAMMS.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
                caller_vcf_channels << XHMM.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (
                params.get('bams', false) &&
                params.get('fasta', false) &&
                params.get('targets', false) &&
                params.get('refflat', false)
            ) {
                def glob_pattern_full = params.bams.contains('{') ? params.bams : params.bams.replaceAll(/\.bam$/, '.{bam,bai,bam.bai}')
                Channel.fromPath(glob_pattern_full)
                    .map { it -> 
                        def sample_id = it.name.replaceAll(/\.bam\.bai$/, '').replaceAll(/\.bam$/, '').replaceAll(/\.bai$/, '').replaceAll(/\.md$/, '')
                        [sample_id, it] 
                    }
                    .groupTuple(size: 2)
                    .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                    .set { ch_bams_cnvkit }
                def ch_bams_cnvkit_nonempty = ch_bams_cnvkit.ifEmpty { error "full workflow CNVkit step received no BAMs. Check the --bams parameter value matches existing BAM files." }
                CNVKIT(
                    ch_bams_cnvkit_nonempty,
                    file(params.fasta),
                    file(params.targets),
                    file(params.refflat),
                    ch_bams_cnvkit_nonempty.first().map { it[1] }
                )
                caller_vcf_channels << CNVKIT.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (
                params.get('samples_path', false) &&
                params.get('fasta', false) &&
                params.get('fai', false) &&
                params.get('dict', false) &&
                params.get('exome_targets', false)
            ) {
                Channel.fromPath(params.samples_path)
                    .filter { it.name.endsWith('.bam') || it.name.endsWith('.cram') }
                    .map { it ->
                        def index_fallback = [file("${it}.bai"), file("${it}.bam.bai"), file(it.toString().replaceAll(/\.bam$/, '.bai'))]
                            .find { f -> f.exists() }
                        def index_file = index_fallback ? index_fallback : (it.name.endsWith('.bam') ? file("${it}.bai") : file("${it}.crai"))
                        return [ it.baseName, it, index_file ]
                    }
                    .set { ch_bams_gcnv }
                GATK_GCNV(
                    ch_bams_gcnv,
                    file(params.fasta),
                    file(params.fai),
                    file(params.dict),
                    file(params.exome_targets)
                )
                caller_vcf_channels << GATK_GCNV.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (params.get('cramFilePairsUploadPath', false)) {
                Channel.fromPath(params.cramFilePairsUploadPath, checkIfExists: true)
                    .filter { f -> f.name.endsWith('.cram') || f.name.endsWith('.bam') }
                    .map { primary ->
                        def sampleId = primary.name.replaceAll(/\.(cram|bam)$/, '')
                        def idxCandidates = primary.name.endsWith('.cram')
                            ? [file("${primary}.crai"), file(primary.toString().replaceAll(/\.cram$/, '.crai'))]
                            : [file("${primary}.bai"), file("${primary}.bam.bai"), file(primary.toString().replaceAll(/\.bam$/, '.bai'))]
                        def indexFile = idxCandidates.find { it.exists() }
                        [sampleId, [primary, indexFile]]
                    }
                    .filter { sampleId, pair -> pair[1] != null }
                    .ifEmpty { error "full workflow (DRAGEN path): no CRAM or BAM file pairs found in '${params.cramFilePairsUploadPath}'. Check that the path contains indexed CRAM/BAM files." }
                    .set { ch_cramPairs_full_dragen }
                DRAGEN(ch_cramPairs_full_dragen)
                caller_vcf_channels << DRAGEN.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (params.get('crams', false)) {
                Channel.fromFilePairs([params.crams + '/*{.cram,.cram.crai}'])
                    .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }
                    .filter { it -> it[1] =~ '_01_1' }
                    .ifEmpty { error "full workflow (INDELIBLE path): no proband CRAMs found in '${params.crams}'. Check that files follow the *_01_1.cram naming convention." }
                    .set { ch_crams_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: 6)
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] ] }
                    .set { ch_cram_trios_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                    .filter { it -> it[1].size() == 4 && it[1][3] =~ '02_2.cram' }
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                    .set { ch_cram_mom_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                    .filter { it -> it[1].size() == 4 && it[1][3] =~ '03_3.cram' }
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                    .set { ch_cram_dad_full }

                INDELIBLE(ch_crams_full, ch_cram_trios_full, ch_cram_mom_full, ch_cram_dad_full)
                caller_vcf_channels << INDELIBLE.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (caller_vcf_channels.size() < 2) {
                exit 1, "Error: full workflow requires at least two CNV caller inputs for consensus and model training. Please configure parameters for at least 2 callers."
            }

            def ch_vcfs = caller_vcf_channels[0]
            caller_vcf_channels.drop(1).each { ch_part ->
                ch_vcfs = ch_vcfs.mix(ch_part)
            }

            grouped_vcfs = group_caller_vcfs(ch_vcfs)

            def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
            def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
            def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
            def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []
            def indelible_f = params.get('indelible_counts', false) ? file(params.indelible_counts) : []

            if (mergerMode == 'truvari') {
                TRUVARI(grouped_vcfs)
                feature_inputs_ch = TRUVARI.out.merged_vcf
                    .join(TRUVARI.out.collapsed_vcf)
                    .join(grouped_vcfs)
                    .map { sample_id, merged_vcf, collapsed_vcf, vcfs ->
                        [sample_id, merged_vcf, collapsed_vcf, build_tool_vcfs_str(vcfs), 'truvari',
                         bam_f, fasta_f, bed_f, map_f, indelible_f]
                    }
                FEATURE_EXTRACTION(feature_inputs_ch)
            } else {
                SURVIVOR(grouped_vcfs)
                feature_inputs_ch = SURVIVOR.out.union_vcf
                    .join(grouped_vcfs)
                    .map { sample_id, merged_vcf, vcfs ->
                        [sample_id, merged_vcf, [], build_tool_vcfs_str(vcfs), 'survivor',
                         bam_f, fasta_f, bed_f, map_f, indelible_f]
                    }
                FEATURE_EXTRACTION(feature_inputs_ch)
            }

            if (params.skip_training) {
                log.info "Skipping model training as requested (skip_training = true)."
            } else {
                Channel.value(file(params.truth_labels))
                    .set { ch_truth_full }
                def probesBedParamFull = params.get('probes_bed', false)
                Channel.value(probesBedParamFull ? file(probesBedParamFull) : [])
                    .set { ch_probes_full }
                TRAIN(FEATURE_EXTRACTION.out.features_tsv, ch_truth_full, ch_probes_full)
            }
            break

        default:
            exit 1, """
OOOPS!! SEEMS LIKE WE HAVE A WORKFLOW ERROR!

No workflow 'mode' given!
Please use one of the following options for workflows:
    --workflow indelible    
    --workflow canoes       
    --workflow xhmm       
    --workflow clamms      
    --workflow dragen      
    --workflow cnvkit      
    --workflow gcnv       
    --workflow survivor     
    --workflow truvari
    --workflow survivor_with_features
    --workflow truvari_with_features
    --workflow feature_extraction
    --workflow normalise
    --workflow train
    --workflow evaluate
    --workflow full
    --workflow qc
"""
            break
    }
}


// =====================================================================================
// POST-CONSENSUS COHORT MERGING, ANNOTATION & EVALUATION
// =====================================================================================

process DOWNLOAD_ANNOTSV_DB {
    publishDir "${projectDir}/annotsv_db", mode: 'copy'

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
    
    # Pre-download the correct human annotations file
    curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_${major_minor}.tar.gz
    cp Annotations_Human_${major_minor}.tar.gz Annotations_Human_${version_clean}.tar.gz 2>/dev/null || true
    
    make PREFIX=\$({readlink -f \$(pwd)/../AnnotSV_annotations}) install
    make PREFIX=\$({readlink -f \$(pwd)/../AnnotSV_annotations}) install-human-annotation
    cd ..
    rm -rf AnnotSV
    """
}

process MERGE_CONSENSUS_COHORT {
    label 'survivor'
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    path vcfs

    output:
    path 'cohort_merged.vcf', emit: merged_vcf

    script:
    """
    # Create the list of VCFs to merge
    > vcf_list.txt
    for vcf in ${vcfs}; do
        echo "\$vcf" >> vcf_list.txt
    done

    # Merge cohort-wide using SURVIVOR
    SURVIVOR merge vcf_list.txt 1000 1 1 0 0 30 cohort_merged.vcf
    """
}

process bcftoolsFlattenCohort {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf)

    output:
    tuple val(run_name), path('cohort_cnvs_by_sample.tsv'), emit: cohort_tsv

    script:
    """
    bcftools query -f '[%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%SAMPLE\\t%GT\\t%TY\\n]' ${merged_vcf} \\
        | awk 'BEGIN{OFS="\\t"} \$7 != "NaN" && \$7 != "0/0" && \$7 != "0" && \$7 != "." && \$7 != "" {print \$1, \$2, \$3, \$4, \$5, \$6}' \\
        > cohort_cnvs_by_sample.tsv
    """
}

process filterRareCnvs {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(cohort_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), path('rare_regions.tsv'), emit: rare_data

    script:
    """
    python3 "${projectDir}/bin/filter-rare-cnvs-ddduk-90.py" \\
        --input-tsv "${cohort_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv" \\
        --total-samples ${params.get('total_samples', 90)} \\
        --max-af ${params.get('max_af', 0.05)} \\
        --output-regions "rare_regions.tsv"
    """
}

process filterVcfByRareRegions {
    label 'bcftools'
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(merged_vcf), path(rare_regions)

    output:
    tuple val(run_name), path('cohort_rare.vcf'), emit: rare_vcf

    script:
    """
    # Sort and compress regions
    sort -k1,1 -k2,2n ${rare_regions} | bgzip -c > regions.bed.gz
    tabix -p bed regions.bed.gz

    # Filter cohort VCF by rare regions
    bcftools view -R regions.bed.gz ${merged_vcf} -o cohort_rare.vcf
    """
}

process runAnnotSV {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_cohort", mode: 'copy'

    input:
    tuple val(run_name), path(rare_vcf)
    path db_dir

    output:
    tuple val(run_name), path('cohort_rare.annotated.tsv'), emit: annotated_tsv

    script:
    """
    export ANNOTSV=\$(readlink -f ${db_dir})
    
    # Run AnnotSV
    AnnotSV \\
        -SVinputFile "${rare_vcf}" \\
        -outputDir . \\
        -outputFile cohort_rare.annotated \\
        -genomeBuild GRCh38 \\
        -overwrite 0 || true
        
    # Fallback headers if file not created
    if [ ! -f cohort_rare.annotated.tsv ]; then
        echo -e "Annotation_mode\\tSV_chrom\\tSV_start\\tSV_end\\tSV_type\\tACMG_class\\tAnnotSV_ranking_score\\tAnnotSV_ranking_criteria\\tGene_name" > cohort_rare.annotated.tsv
    fi
    """
}

process addACMGAnnotations {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_acmg_reports", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv)

    output:
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs.csv'), emit: rare_cnvs_annotated
    tuple val(run_name), path('DDD_Cohort_Rare_CNVs_*.csv'), emit: score_separated_csvs

    script:
    """
    python3 "${projectDir}/bin/add-acmg-annotations.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --output-csv "DDD_Cohort_Rare_CNVs.csv"
    """
}

process ANALYZE_CNV_QUALITY_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/consensus_metrics", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(rare_cnvs_csv), path(annotated_tsv), path(filtered_vcfs)
    path probes_bed

    output:
    tuple val(run_name), path("${run_name}_cnv_metrics_summary.txt"), path("${run_name}_cnv_vcf_metrics.csv"), emit: metrics_reports

    script:
    """
    python3 "${projectDir}/bin/analyze_cnv_quality_metrics.py" \\
        --input-csv "${rare_cnvs_csv}" \\
        --annotsv-tsv "${annotated_tsv}" \\
        --vcf-dir . \\
        --probes-bed "${probes_bed}" \\
        --run-name "${run_name}" \\
        --output-summary "${run_name}_cnv_metrics_summary.txt" \\
        --output-csv "${run_name}_cnv_vcf_metrics.csv"
    """
}

process EVALUATE_ACMG_CSV {
    tag "${run_name}-${csv.baseName}"
    publishDir "${params.outdir}/consensus_evaluate/${run_name}", mode: 'copy', overwrite: true

    input:
    tuple val(run_name), path(csv), path(probes_bed)
    path map_file
    path truth_bed

    output:
    path "*_performance_report.txt", emit: report

    script:
    def base_name = csv.baseName
    """
    python3 "${projectDir}/bin/evaluate_separated_helper.py" \\
        --input-csv "${csv}" \\
        --map-file "${map_file}" \\
        --truth-bed "${truth_bed}" \\
        --probes-bed "${probes_bed}" \\
        --workflow-type "${run_name}" \\
        --output "${run_name}_${base_name}_performance_report.txt"
    """
}

process PLOT_EVALUATION_METRICS {
    publishDir "${params.outdir}/consensus_evaluate", mode: 'copy', overwrite: true

    input:
    path reports
    
    output:
    path 'performance_metrics_comparison.png', emit: comparison_plot

    script:
    """
    # Structure files dynamically into subdirectories based on their run name prefix
    for f in *_performance_report.txt; do
        folder=\$(echo "\$f" | cut -d'_' -f1-2)
        mkdir -p "\$folder"
        cp "\$f" "\$folder"/
    done
    
    python3 "${projectDir}/bin/plot_evaluation_performance_dynamic.py" \\
        --input-dir . \\
        --output-png performance_metrics_comparison.png
    """
}

workflow RUN_CONSENSUS_ANNOTATION {
    take:
        run_name
        sample_vcfs_ch
        annotsv_db_ch
        probes_bed
    main:
        // Cohort-wide merge
        MERGE_CONSENSUS_COHORT(sample_vcfs_ch.collect())
        
        cohort_ch = Channel.of(run_name).combine(MERGE_CONSENSUS_COHORT.out.merged_vcf)
        
        bcftoolsFlattenCohort(cohort_ch)
        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
        
        slice_joined = cohort_ch.join(filterRareCnvs.out.rare_data)
        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
        
        runAnnotSV(filterVcfByRareRegions.out.rare_vcf, annotsv_db_ch)
        
        acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })

        // Quality assessment
        metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
            .join(runAnnotSV.out.annotated_tsv)
            .combine(sample_vcfs_ch.collect())
            
        metrics_joined = metrics_inputs.combine(Channel.of(probes_bed))
        
        ANALYZE_CNV_QUALITY_METRICS(
            metrics_joined.map { [ it[0], it[1], it[2], it[3] ] },
            metrics_joined.map { it[4] }
        )

        // ACMG Evaluators
        eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
        eval_joined = eval_inputs.combine(Channel.of(probes_bed))
        
        EVALUATE_ACMG_CSV(
            eval_joined.map { [ it[0], it[1], it[2] ] },
            file(params.map_file),
            file(params.truth_bed)
        )

        PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
}

        case['survivor_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR_WITH_FEATURES(ch_vcfs)
            break

        case['truvari_with_features']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI_WITH_FEATURES(ch_vcfs)
            break

        case['normalise']:
            // Normalise quality scores for pre-existing caller VCFs without running
            // the full caller pipeline.
            // Required: --vcf_dir  (directory containing raw .vcf / .vcf.gz files)
            //           --caller   (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN | INDELIBLE)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f ->
                    def sample_name = f.name
                        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
                        .replaceAll(/\.(sorted\.)?vcf(\.gz)?$/, '')
                    [sample_name, f, params.caller.toUpperCase()]
                }
                .ifEmpty { error "normalise workflow: no VCF files found in '${params.vcf_dir}'. Check that the directory contains .vcf or .vcf.gz files." }
                .set { ch_vcfs }
            RUN_NORMALISE(ch_vcfs)
            break

        case['feature_extraction']:
            // Run feature extraction on an already-produced merged VCF.
            // Required: --merged_vcf_dir  (directory containing <sample_id>_survivor_union.vcf
            //                              or <sample_id>_truvari_merged.vcf files)
            // Optional annotation params (any subset may be omitted):
            //   --canoes_norm_dir, --clamms_norm_dir, --xhmm_norm_dir, --cnvkit_norm_dir,
            //   --gcnv_norm_dir, --dragen_norm_dir, --indelible_norm_dir
            //   --bam_file, --reference_fasta, --bed_file, --mappability_file,
            //   --indelible_counts, --merger_mode (default: survivor)
            if (params.get('merger_mode', 'survivor') == 'truvari') {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_truvari_merged.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_truvari.*/, '').replaceAll(/(?i)\.vcf(\.gz)?$/, '')
                        def collapsed_f = file("${f.parent}/${sample_id}_truvari_collapsed.vcf")
                        def collapsed_vcf_f = collapsed_f.exists() ? collapsed_f : []
                        build_feature_inputs(sample_id, f, collapsed_vcf_f)
                    }
                    .ifEmpty { error "feature_extraction workflow (truvari mode): no merged VCF files found in '${params.merged_vcf_dir}'. Expected files matching *_truvari_merged.vcf or *_truvari_merged.vcf.gz." }
                    .set { ch_feature_inputs }
            } else {
                Channel.fromPath(params.merged_vcf_dir + '/**/*_survivor_union.vcf*')
                    .map { f ->
                        def sample_id = f.name.replaceAll(/_survivor.*/, '').replaceAll(/(?i)\.vcf(\.gz)?$/, '')
                        build_feature_inputs(sample_id, f, [])
                    }
                    .ifEmpty { error "feature_extraction workflow (survivor mode): no merged VCF files found in '${params.merged_vcf_dir}'. Expected files matching *_survivor_union.vcf or *_survivor_union.vcf.gz." }
                    .set { ch_feature_inputs }
            }
            RUN_FEATURE_EXTRACTION(ch_feature_inputs)
            break

        case['train']:
            if (params.skip_training) {
                log.info "Skipping model training as requested (skip_training = true)."
                break
            }
            // Train an XGBoost classifier on CNV feature matrices produced by
            // the feature_extraction workflow.
            // Required: --features_dir   (directory containing *_features.tsv files
            //                             produced by the feature_extraction workflow)
            //           --truth_labels   (TSV file with columns:
            //                             sample_id, chrom, start, end, cnv_type, truth_label
            //                             where truth_label=1 means true CNV)
            // Optional: --probes_bed     (capture-target BED used to match rows by
            //                             shared probes when coordinates differ)
            Channel.fromPath(params.features_dir + '/**/*_features.tsv')
                .ifEmpty { error "train workflow: no feature TSV files found in '${params.features_dir}'. Expected *_features.tsv files produced by the feature_extraction workflow." }
                .set { ch_features }
            Channel.value(file(params.truth_labels))
                .set { ch_truth }
            def probesBedParam = params.get('probes_bed', false)
            Channel.value(probesBedParam ? file(probesBedParam) : [])
                .set { ch_probes }
            RUN_TRAIN(ch_features, ch_truth, ch_probes)
            break

        case['evaluate']:
            // Evaluate CNV caller performance (precision and sensitivity) using
            // pre-existing per-sample VCFs produced by one of the 7 callers.
            // VCFs are converted to BED, merged into a unified call set, and
            // compared against the truth set at the probe (capture target) level.
            //
            // Required: --vcf_dir    (directory containing per-sample .vcf / .vcf.gz files)
            //           --caller     (caller label used in output names, e.g. CANOES)
            //           --truth_bed  (truth set BED: CHR, START, STOP, CNV_TYPE, SAMPLE_ID)
            //           --probes_bed (capture target BED: CHR, START, STOP[, ...])
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .ifEmpty { error "evaluate workflow: no VCF files found in '${params.vcf_dir}'. Check that the directory contains .vcf or .vcf.gz files." }
                .set { ch_vcfs }
            RUN_EVALUATE(
                ch_vcfs,
                file(params.truth_bed),
                file(params.probes_bed),
                params.caller
            )
            break

        case['full']:
            def caller_vcf_channels = []
            def mergerMode = params.get('merger_mode', 'survivor')

            if (params.get('samplesheet_bams', false) && params.get('fai', false)) {
                Channel.fromPath(params.samplesheet_bams)
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                    .set { ch_samplesheet_bams }
                chroms = (1..22).toList().collect { 'chr' + "${it}" }
                bam_list = ch_samplesheet_bams.collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
                CANOES(bam_list, file(params.fai), Channel.from(chroms))
                sample_list = ch_samplesheet_bams.map { it -> it[0] + '\n' }.collectFile(name: 'sample_list.txt')
                CLAMMS(ch_samplesheet_bams, file(params.fai), sample_list)
                XHMM(ch_samplesheet_bams)
                caller_vcf_channels << CANOES.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
                caller_vcf_channels << CLAMMS.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
                caller_vcf_channels << XHMM.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (
                params.get('bams', false) &&
                params.get('fasta', false) &&
                params.get('targets', false) &&
                params.get('refflat', false)
            ) {
                def glob_pattern_full = params.bams.contains('{') ? params.bams : params.bams.replaceAll(/\.bam$/, '.{bam,bai,bam.bai}')
                Channel.fromPath(glob_pattern_full)
                    .map { it -> 
                        def sample_id = it.name.replaceAll(/\.bam\.bai$/, '').replaceAll(/\.bam$/, '').replaceAll(/\.bai$/, '').replaceAll(/\.md$/, '')
                        [sample_id, it] 
                    }
                    .groupTuple(size: 2)
                    .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                    .set { ch_bams_cnvkit }
                def ch_bams_cnvkit_nonempty = ch_bams_cnvkit.ifEmpty { error "full workflow CNVkit step received no BAMs. Check the --bams parameter value matches existing BAM files." }
                CNVKIT(
                    ch_bams_cnvkit_nonempty,
                    file(params.fasta),
                    file(params.targets),
                    file(params.refflat),
                    ch_bams_cnvkit_nonempty.first().map { it[1] }
                )
                caller_vcf_channels << CNVKIT.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (
                params.get('samples_path', false) &&
                params.get('fasta', false) &&
                params.get('fai', false) &&
                params.get('dict', false) &&
                params.get('exome_targets', false)
            ) {
                Channel.fromPath(params.samples_path)
                    .filter { it.name.endsWith('.bam') || it.name.endsWith('.cram') }
                    .map { it ->
                        def index_fallback = [file("${it}.bai"), file("${it}.bam.bai"), file(it.toString().replaceAll(/\.bam$/, '.bai'))]
                            .find { f -> f.exists() }
                        def index_file = index_fallback ? index_fallback : (it.name.endsWith('.bam') ? file("${it}.bai") : file("${it}.crai"))
                        return [ it.baseName, it, index_file ]
                    }
                    .set { ch_bams_gcnv }
                GATK_GCNV(
                    ch_bams_gcnv,
                    file(params.fasta),
                    file(params.fai),
                    file(params.dict),
                    file(params.exome_targets)
                )
                caller_vcf_channels << GATK_GCNV.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (params.get('cramFilePairsUploadPath', false)) {
                Channel.fromPath(params.cramFilePairsUploadPath, checkIfExists: true)
                    .filter { f -> f.name.endsWith('.cram') || f.name.endsWith('.bam') }
                    .map { primary ->
                        def sampleId = primary.name.replaceAll(/\.(cram|bam)$/, '')
                        def idxCandidates = primary.name.endsWith('.cram')
                            ? [file("${primary}.crai"), file(primary.toString().replaceAll(/\.cram$/, '.crai'))]
                            : [file("${primary}.bai"), file("${primary}.bam.bai"), file(primary.toString().replaceAll(/\.bam$/, '.bai'))]
                        def indexFile = idxCandidates.find { it.exists() }
                        [sampleId, [primary, indexFile]]
                    }
                    .filter { sampleId, pair -> pair[1] != null }
                    .ifEmpty { error "full workflow (DRAGEN path): no CRAM or BAM file pairs found in '${params.cramFilePairsUploadPath}'. Check that the path contains indexed CRAM/BAM files." }
                    .set { ch_cramPairs_full_dragen }
                DRAGEN(ch_cramPairs_full_dragen)
                caller_vcf_channels << DRAGEN.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (params.get('crams', false)) {
                Channel.fromFilePairs([params.crams + '/*{.cram,.cram.crai}'])
                    .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }
                    .filter { it -> it[1] =~ '_01_1' }
                    .ifEmpty { error "full workflow (INDELIBLE path): no proband CRAMs found in '${params.crams}'. Check that files follow the *_01_1.cram naming convention." }
                    .set { ch_crams_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: 6)
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] ] }
                    .set { ch_cram_trios_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                    .filter { it -> it[1].size() == 4 && it[1][3] =~ '02_2.cram' }
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                    .set { ch_cram_mom_full }

                Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
                    .filter { it -> it[1].size() == 4 && it[1][3] =~ '03_3.cram' }
                    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
                    .set { ch_cram_dad_full }

                INDELIBLE(ch_crams_full, ch_cram_trios_full, ch_cram_mom_full, ch_cram_dad_full)
                caller_vcf_channels << INDELIBLE.out.normalised_vcf.flatten().map { f -> [extract_sample_id_from_vcf(f), f] }
            }

            if (caller_vcf_channels.size() < 2) {
                exit 1, "Error: full workflow requires at least two CNV caller inputs for consensus and model training. Please configure parameters for at least 2 callers."
            }

            def ch_vcfs = caller_vcf_channels[0]
            caller_vcf_channels.drop(1).each { ch_part ->
                ch_vcfs = ch_vcfs.mix(ch_part)
            }

            grouped_vcfs = group_caller_vcfs(ch_vcfs)

            def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
            def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
            def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
            def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []
            def indelible_f = params.get('indelible_counts', false) ? file(params.indelible_counts) : []

            if (mergerMode == 'truvari') {
                TRUVARI(grouped_vcfs)
                feature_inputs_ch = TRUVARI.out.merged_vcf
                    .join(TRUVARI.out.collapsed_vcf)
                    .join(grouped_vcfs)
                    .map { sample_id, merged_vcf, collapsed_vcf, vcfs ->
                        [sample_id, merged_vcf, collapsed_vcf, build_tool_vcfs_str(vcfs), 'truvari',
                         bam_f, fasta_f, bed_f, map_f, indelible_f]
                    }
                FEATURE_EXTRACTION(feature_inputs_ch)
            } else {
                SURVIVOR(grouped_vcfs)
                feature_inputs_ch = SURVIVOR.out.union_vcf
                    .join(grouped_vcfs)
                    .map { sample_id, merged_vcf, vcfs ->
                        [sample_id, merged_vcf, [], build_tool_vcfs_str(vcfs), 'survivor',
                         bam_f, fasta_f, bed_f, map_f, indelible_f]
                    }
                FEATURE_EXTRACTION(feature_inputs_ch)
            }

            if (params.skip_training) {
                log.info "Skipping model training as requested (skip_training = true)."
            } else {
                Channel.value(file(params.truth_labels))
                    .set { ch_truth_full }
                def probesBedParamFull = params.get('probes_bed', false)
                Channel.value(probesBedParamFull ? file(probesBedParamFull) : [])
                    .set { ch_probes_full }
                TRAIN(FEATURE_EXTRACTION.out.features_tsv, ch_truth_full, ch_probes_full)
            }
            break

        default:
            exit 1, """
OOOPS!! SEEMS LIKE WE HAVE A WORKFLOW ERROR!

No workflow 'mode' given!
Please use one of the following options for workflows:
    --workflow indelible    
    --workflow canoes       
    --workflow xhmm       
    --workflow clamms      
    --workflow dragen      
    --workflow cnvkit      
    --workflow gcnv       
    --workflow survivor     
    --workflow truvari
    --workflow survivor_with_features
    --workflow truvari_with_features
    --workflow feature_extraction
    --workflow normalise
    --workflow train
    --workflow evaluate
    --workflow full
    --workflow qc
"""
            break
    }
}
