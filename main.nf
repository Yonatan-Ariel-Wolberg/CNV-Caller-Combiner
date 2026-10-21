#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    RUN_CANOES; RUN_XHMM; RUN_CLAMMS; RUN_DRAGEN; RUN_CNVKIT; RUN_GCNV; RUN_EXOMEDEPTH; RUN_EVALUATE; 
    RUN_SURVIVOR; RUN_TRUVARI; RUN_FEATURE_EXTRACTION; 
    RUN_SURVIVOR_WITH_FEATURES; RUN_TRUVARI_WITH_FEATURES; RUN_QC; RUN_CONSENSUS_ANNOTATION
} from './modules/workflows/workflows.nf'
include { NORMALISE } from './modules/normalise/modules-normalise.nf'

workflow RUN_NORMALISE {
    take:
        vcf_ch
    main:
        NORMALISE(vcf_ch)
}

workflow RUN_TRAIN {
    take:
        features_tsv_ch
        truth_labels_ch
        probes_bed_ch
    main:
        TRAIN(features_tsv_ch, truth_labels_ch, probes_bed_ch)
}



import groovy.transform.Field

// =====================================================================================
// MODULE INCLUDES
// =====================================================================================
include { SAMTOOLS_QC_BAM; SAMTOOLS_QC_CRAM; QUALIMAP; MULTIQC; INFER_SEX_AND_PLOIDY; ESTIMATE_CONTAMINATION } from './modules/sample-qc/target-aware-exome-sequencing-quality-assess.nf'
include { 
    FILTER_CANOES_QUALITY; FILTER_CLAMMS_QUALITY; FILTER_XHMM_QUALITY; 
    FILTER_CNVKIT_QUALITY; FILTER_GCNV_QUALITY; FILTER_DRAGEN_QUALITY;
    FILTER_EXOMEDEPTH_QUALITY 
} from './modules/filters/quality-filters.nf'
include { 
    FILTER_CANOES_QUALITY as FILTER_CANOES_FULL;
    FILTER_CLAMMS_QUALITY as FILTER_CLAMMS_FULL;
    FILTER_XHMM_QUALITY as FILTER_XHMM_FULL;
    FILTER_CNVKIT_QUALITY as FILTER_CNVKIT_FULL;
    FILTER_GCNV_QUALITY as FILTER_GCNV_FULL;
    FILTER_DRAGEN_QUALITY as FILTER_DRAGEN_FULL;
    FILTER_EXOMEDEPTH_QUALITY as FILTER_EXOMEDEPTH_FULL 
} from './modules/filters/quality-filters.nf'
include {
    NORMALISE_CNV_QUALITY_SCORES as NORM_CANOES_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_CLAMMS_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_XHMM_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_CNVKIT_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_GCNV_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_DRAGEN_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_EXOMEDEPTH_FULL;
    NORMALISE_CNV_QUALITY_SCORES as NORM_CANOES_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_CLAMMS_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_XHMM_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_CNVKIT_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_GCNV_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_DRAGEN_GATHER;
    NORMALISE_CNV_QUALITY_SCORES as NORM_EXOMEDEPTH_GATHER
} from './modules/normalise/normalise-quality.nf'
include { BGZIP_SORT_INDEX_SV_VCF } from './modules/common/bgzip-sort-index.nf'
include { EVALUATE_COMPARATIVE; EVALUATE_ACMG_CSV; PLOT_EVALUATION_METRICS } from './modules/evaluation/performance-metric-evaluation.nf'

include { CANOES } from './modules/callers/modules-canoes.nf'
include { XHMM } from './modules/callers/modules-xhmm.nf'
include { CLAMMS } from './modules/callers/modules-clamms.nf'
include { DRAGEN } from './modules/callers/modules-icav2-dragen.nf'
include { CNVKIT } from './modules/callers/modules-cnvkit.nf'
include { GATK_GCNV } from './modules/callers/modules-gatk-gcnv.nf'
include { EXOMEDEPTH } from './modules/callers/modules-exomedepth.nf'
include { SURVIVOR } from './modules/sv-mergers/modules-survivor.nf'
include { TRUVARI } from './modules/sv-mergers/modules-truvari.nf'
include { SVABA } from './modules/callers/modules-svaba.nf'
include { REFINE_CONSENSUS_BREAKPOINTS } from './modules/refinement/breakpoint-refinement.nf'
include { FEATURE_EXTRACTION } from './modules/ml/modules-feature-extraction.nf'
include { COMPUTE_TARGET_COVERAGE } from './modules/ml/modules-feature-extraction.nf'
include { TRAIN } from './modules/ml/modules-train.nf'
include { RUN_CATBOOST_PREDICT } from './modules/ml/modules-predict.nf'
include { MERGE_CONSENSUS_COHORT; SPLIT_COHORT_VCF; BCFTOOLS_FLATTEN_COHORT; SANITIZE_COHORT_VCF; bcftoolsFlattenCohort; sanitizeCohortVcf } from './modules/consensus/consensus-postprocess.nf'
include { BUILD_COHORT_FEATURE_DATABASE } from './modules/ml/modules-cohort-database.nf'
include { DOWNLOAD_ANNOTSV_DB; RUN_ANNOTSV; ADD_ACMG_ANNOTATIONS; STRATIFY_VCF_BY_ACMG; runAnnotSV; addACMGAnnotations; stratifyVcfByAcmg } from './modules/annotation/annotation.nf'
include { FILTER_RARE_CNVS; FILTER_VCF_BY_RARE_REGIONS; filterRareCnvs; filterVcfByRareRegions } from './modules/filters/allele-frequency-filters.nf'
include { ANALYZE_CNV_QUALITY_METRICS; GENERATE_ATTRITION_FUNNEL; ASSESS_CALLSET_QUALITY_METRICS; COMPILE_ALL_STAGES_QUALITY_REPORT } from './modules/evaluation/quality-metric-evaluation.nf'
include { CONVERT_CRAM_TO_BAM_VERIFIED; CONVERT_BAM_TO_CRAM_VERIFIED; RUN_CRAM_TO_BAM_CONVERSION; RUN_BAM_TO_CRAM_CONVERSION; CONVERT_CRAM_TO_TARGET_PADDED_BAM; RUN_CRAM_TO_TARGET_PADDED_BAM_CONVERSION } from './modules/conversions/modules-conversions.nf'
include { GENERATE_IGV_REPORT } from './modules/visualization/modules-igv.nf'


// =====================================================================================
// GLOBAL SETUP
// =====================================================================================
workflow_mode = params.workflow
params.test_size = params.get('test_size', -1)

// Harmonize probe / target BED aliases across all callers and evaluation workflows
def resolved_probes = params.probes_bed ?: params.probes ?: params.targets ?: null
if (resolved_probes) {
    params.probes_bed = resolved_probes
    params.probes     = resolved_probes
    params.targets    = resolved_probes
}

@Field def CALLER_DIR_PARAMS = [
    'canoes_dir',
    'clamms_dir',
    'xhmm_dir',
    'cnvkit_dir',
    'gcnv_dir',
    'dragen_dir',
    'exomedepth_dir',
]

// Valid caller names accepted by normalise_cnv_caller_quality_scores.py.
// Note: GATK gCNV is referred to as 'GATK' (not 'GCNV') by the normalise script.
@Field def VALID_NORMALISE_CALLERS = ['CANOES', 'CLAMMS', 'XHMM', 'GATK', 'CNVKIT', 'DRAGEN', 'EXOMEDEPTH']
@Field def BYTES_PER_GIB = 1024D * 1024D * 1024D
@Field def DEFAULT_MIN_FREE_GB = 5
@Field def MAX_VCF_SCHEMA_VALIDATION_FILES = 5

@Field def REQUIRED_PARAMS_BY_WORKFLOW = [
    // Required params are aligned to the workflow-specific params/*.json templates.
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
    'exomedepth': ['outdir', 'samplesheet_bams', 'ref', 'fai', 'probes'],
    'survivor': ['outdir', 'canoes_dir', 'clamms_dir'],
    'truvari': ['outdir', 'canoes_dir', 'clamms_dir'],
    'survivor_with_features': ['outdir', 'merger_mode', 'canoes_dir', 'clamms_dir'],
    'truvari_with_features': ['outdir', 'merger_mode', 'canoes_dir', 'clamms_dir'],
    'normalise': ['outdir', 'vcf_dir', 'caller'],
    'feature_extraction': ['outdir', 'merged_vcf_dir', 'merger_mode'],
    'train': ['outdir', 'features_dir', 'truth_labels'],
    'evaluate': ['outdir', 'vcf_dir', 'caller', 'truth_bed', 'probes_bed'],
    'full': ['truth_labels'],
    'local': ['truth_labels'],
    'full_with_dragen': ['truth_labels'],
    'qc': ['input_dir'],
]

def is_param_set(param_name) {
    def value = params.get(param_name, null)
    if (value == null) return false
    if (value instanceof CharSequence) return value.toString().trim().length() > 0
    return true
}

def is_caller_enabled(caller_name, workflow_name) {
    def caller = caller_name.toString().toLowerCase().trim()
    def wf = workflow_name ? workflow_name.toString().trim() : ''
    def skip_key = "skip_${caller}".toString()
    def run_key = "run_${caller}".toString()
    // Explicit skip flag (e.g. --skip_canoes, --skip_gcnv_training)
    if (params.get(skip_key, false) || (caller == 'gatk' && params.get('skip_gcnv_training', false))) {
        return false
    }
    // Explicit caller list passed via --callers (e.g. --callers canoes,xhmm,exomedepth)
    if (params.get('callers', false)) {
        def requested = params.callers.toString().toLowerCase().tokenize(',').collect { it.trim() }
        if (requested.contains('all')) return true
        if (requested.contains('local')) return caller != 'dragen'
        if (caller == 'gatk' && (requested.contains('gatk') || requested.contains('gcnv') || requested.contains('gatk_gcnv') || requested.contains('gatk-gcnv'))) return true
        return requested.contains(caller)
    }
    // Explicit run flag (e.g. --run_canoes true, --run_exomedepth true)
    if (params.containsKey(run_key)) {
        return params[run_key].toString().toBoolean()
    }
    if (caller == 'gatk' && (params.containsKey('run_gcnv') || params.containsKey('run_gatk_gcnv'))) {
        return (params.get('run_gcnv', false) || params.get('run_gatk_gcnv', false)).toString().toBoolean()
    }
    // Workflow-based defaults
    if (caller == 'dragen') {
        return wf == 'full_with_dragen' || (wf == 'full' && params.get('cramFilePairsUploadPath', false))
    }
    if (wf == 'local') {
        return caller != 'dragen'  // All 6 local callers enabled
    }
    if (wf in ['full', 'full_with_dragen']) {
        return true
    }
    return false
}

def validate_required_params(workflow_name) {
    if (REQUIRED_PARAMS_BY_WORKFLOW.containsKey(workflow_name)) {
        def required_keys = ['survivor', 'truvari', 'survivor_with_features', 'truvari_with_features'].contains(workflow_name) ?
            REQUIRED_PARAMS_BY_WORKFLOW[workflow_name].findAll { !it.endsWith('_dir') } :
            REQUIRED_PARAMS_BY_WORKFLOW[workflow_name]
        if (['full', 'local', 'full_with_dragen'].contains(workflow_name) && params.get('skip_ml_training', false)) {
            required_keys = required_keys.findAll { it != 'truth_labels' }
        }
        def missing = required_keys.findAll { !is_param_set(it) }
        if (!missing.isEmpty()) {
            exit 1, "Error: Missing required parameter(s) for --workflow ${workflow_name}: ${missing.collect { '--' + it }.join(', ')}"
        }
    }

    if (['survivor', 'truvari', 'survivor_with_features', 'truvari_with_features'].contains(workflow_name)) {
        def configured_caller_dirs = CALLER_DIR_PARAMS.findAll { is_param_set(it) }
        def min_callers = params.get('min_callers', 1) as int
        if (min_callers > 1 && configured_caller_dirs.size() < 2) {
            exit 1, "Error: --workflow ${workflow_name} requires at least TWO caller VCF directories. Configure at least 2 of: ${CALLER_DIR_PARAMS.collect { '--' + it }.join(', ')}"
        } else if (configured_caller_dirs.size() < min_callers) {
            exit 1, "Error: --workflow ${workflow_name} requires at least ${min_callers} caller VCF directory. Configure at least ${min_callers} of: ${CALLER_DIR_PARAMS.collect { '--' + it }.join(', ')}"
        }
    }

    if (['full', 'local', 'full_with_dragen'].contains(workflow_name)) {
        def configured_caller_count = 0
        def has_samplesheet = is_param_set('samplesheet_bams') || is_param_set('samplesheet_crams')

        // CANOES
        if (is_caller_enabled('canoes', workflow_name) && has_samplesheet && (is_param_set('fai') || is_param_set('ref') || is_param_set('fasta'))) {
            configured_caller_count += 1
        }
        // CLAMMS
        if (is_caller_enabled('clamms', workflow_name) && has_samplesheet && (is_param_set('fai') || is_param_set('ref') || is_param_set('fasta'))) {
            configured_caller_count += 1
        }
        // XHMM
        if (is_caller_enabled('xhmm', workflow_name) && has_samplesheet) {
            configured_caller_count += 1
        }
        // CNVKIT
        def has_cnvkit_bams = is_param_set('bams') || has_samplesheet
        def has_cnvkit_fasta = is_param_set('fasta') || is_param_set('ref') || is_param_set('reference_fasta')
        def has_cnvkit_targets = is_param_set('targets') || is_param_set('exome_targets') || is_param_set('probes') || is_param_set('probes_bed')
        if (is_caller_enabled('cnvkit', workflow_name) && has_cnvkit_bams && has_cnvkit_fasta && has_cnvkit_targets && is_param_set('refflat')) {
            configured_caller_count += 1
        }
        // GATK
        def has_gatk_bams = is_param_set('samples_path') || has_samplesheet
        def has_gatk_fasta = is_param_set('fasta') || is_param_set('ref') || is_param_set('reference_fasta')
        def has_gatk_targets = is_param_set('exome_targets') || is_param_set('targets') || is_param_set('probes') || is_param_set('probes_bed')
        if (is_caller_enabled('gatk', workflow_name) && has_gatk_bams && has_gatk_fasta && is_param_set('dict') && has_gatk_targets) {
            configured_caller_count += 1
        }
        // DRAGEN
        if (is_caller_enabled('dragen', workflow_name) && workflow_name != 'local' && is_param_set('cramFilePairsUploadPath')) {
            configured_caller_count += 1
        }
        // EXOMEDEPTH
        def has_ed_probes = is_param_set('probes') || is_param_set('targets') || is_param_set('exome_targets') || is_param_set('probes_bed')
        if (is_caller_enabled('exomedepth', workflow_name) && has_samplesheet && has_ed_probes) {
            configured_caller_count += 1
        }

        // configured_caller_count += 3
        def default_min = (workflow_name == 'dragen' || (params.get('enable_cnv', true) == false && is_param_set('cramFilePairsUploadPath'))) ? 0 : 1
        def min_callers = params.get('min_callers', default_min) as int
        if (min_callers > 1 && configured_caller_count < 2) {
            exit 1, "Error: --workflow full requires at least 2 configured CNV callers out of 7."
        } else if (configured_caller_count < min_callers) {
            exit 1, "Error: --workflow ${workflow_name} requires at least ${min_callers} configured and enabled CNV caller(s)."
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
        def container_str = params.icav2_container.toString()
        if (!container_str.startsWith('docker://') && !container_str.startsWith('shub://') && !container_str.startsWith('oras://')) {
            def container_file = require_readable_file(params.icav2_container, 'icav2_container')
            if (!container_file.name.toLowerCase().endsWith('.sif')) {
                fail_validation("Error: ICAv2 asset unreadable or invalid: --icav2_container must point to a .sif file or remote container URI (${params.icav2_container})")
            }
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
            (it.name.contains(suffix) || it.name.contains('.intra-merge.')) &&
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
    if (!params.skip_ml_training && ['train', 'full'].contains(workflow_name) && is_param_set('truth_labels')) {
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
    // Extract bare sample ID from per-caller VCF names across all stages:
    // Raw: {sample}.<caller>.cnv.vcf.gz
    // Filtered: {sample}.<caller>.qs-filtered.cnv.vcf.gz
    // Normalised: {sample}.<caller>.normalised.qs-filtered.cnv.vcf.gz
    // Normalised (bypass filter): {sample}.<caller>.normalised.cnv.vcf.gz
    // Legacy formats: {sample}_<CALLER>_output.vcf, {sample}.sorted.vcf.gz, etc.
    return f.name
        .replaceAll(/(?i)\.(canoes|clamms|xhmm|cnvkit|gcnv|dragen|exomedepth)\..*/, '')
        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|EXOMEDEPTH).*/, '')
        .replaceAll(/(?i)\.intra-merge.*$/, '')
        .replaceAll(/(?i)_survivor.*$/, '')
        .replaceAll(/(?i)_truvari.*$/, '')
        .replaceAll(/(?i)\.normalised.*$/, '')
        .replaceAll(/(?i)\.quality-filtered.*$/, '')
        .replaceAll(/(?i)\.qs-filtered.*$/, '')
        .replaceAll(/(?i)\.filtered.*$/, '')
        .replaceAll(/(?i)\.sorted.*$/, '')
        .replaceAll(/(?i)\.cnv.*$/, '')
        .replaceAll(/(?i)\.vcf(\.gz)?$/, '')
        .replaceAll(/(?i)\.md$/, '')
}

def infer_caller_from_vcf(f) {
    def n = f.name.toUpperCase()
    def p = f.toString().toUpperCase()
    if (n.contains('.CANOES.') || n.contains('_CANOES') || p.contains('/OUT_CANOES/'))        return 'canoes'
    if (n.contains('.CLAMMS.') || n.contains('_CLAMMS') || p.contains('/OUT_CLAMMS/'))        return 'clamms'
    if (n.contains('.XHMM.')   || n.contains('_XHMM')   || p.contains('/OUT_XHMM/'))          return 'xhmm'
    if (n.contains('.CNVKIT.') || n.contains('_CNVKIT') || p.contains('/OUT_CNVKIT/'))        return 'cnvkit'
    if (n.contains('.GCNV.')   || n.contains('_GCNV')   || p.contains('/OUT_GCNV/'))          return 'gatk_gcnv'
    if (n.contains('.DRAGEN.') || n.contains('_DRAGEN') || p.contains('/OUT_DRAGEN/') || n.contains('.CNV.')) return 'dragen'
    if (n.contains('.EXOMEDEPTH.') || n.contains('_EXOMEDEPTH') || p.contains('/OUT_EXOMEDEPTH/')) return 'exomedepth'
    error "Could not infer caller from VCF filename: ${f.name} (path: ${f.toString()})"
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
workflow GATHER_AND_FILTER_VCFS {
    main:
        def ch = Channel.empty()
        def ch_comp = Channel.empty()
        
        def skip_norm = params.get('skip_normalisation', false) || params.get('skip_quality_normalisation', false)
        if (params.get('canoes_dir', false)) {
            canoes_in = Channel.of('CANOES').combine(Channel.fromPath(params.canoes_dir + "/*.vcf*").collect())
            FILTER_CANOES_QUALITY(canoes_in)
            def vcf_filtered = FILTER_CANOES_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_CANOES_GATHER(vcf_filtered, 'CANOES', 'out_CANOES')
                vcf_norm = NORM_CANOES_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["CANOES_filtered", f] })
        }
        if (params.get('clamms_dir', false)) {
            clamms_in = Channel.of('CLAMMS').combine(Channel.fromPath(params.clamms_dir + "/*.vcf*").collect())
            FILTER_CLAMMS_QUALITY(clamms_in)
            def vcf_filtered = FILTER_CLAMMS_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_CLAMMS_GATHER(vcf_filtered, 'CLAMMS', 'out_CLAMMS')
                vcf_norm = NORM_CLAMMS_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["CLAMMS_filtered", f] })
        }
        if (params.get('xhmm_dir', false)) {
            xhmm_in = Channel.of('XHMM').combine(Channel.fromPath(params.xhmm_dir + "/*.vcf*").collect())
            FILTER_XHMM_QUALITY(xhmm_in)
            def vcf_filtered = FILTER_XHMM_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_XHMM_GATHER(vcf_filtered, 'XHMM', 'out_XHMM')
                vcf_norm = NORM_XHMM_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["XHMM_filtered", f] })
        }
        if (params.get('cnvkit_dir', false)) {
            cnvkit_in = Channel.of('CNVKIT').combine(Channel.fromPath(params.cnvkit_dir + "/*.vcf*").collect())
            FILTER_CNVKIT_QUALITY(cnvkit_in)
            def vcf_filtered = FILTER_CNVKIT_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_CNVKIT_GATHER(vcf_filtered, 'CNVKIT', 'out_CNVKIT')
                vcf_norm = NORM_CNVKIT_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["CNVKIT_filtered", f] })
        }
        if (params.get('gcnv_dir', false)) {
            gcnv_in = Channel.of('GCNV').combine(Channel.fromPath(params.gcnv_dir + "/*.vcf*").collect())
            FILTER_GCNV_QUALITY(gcnv_in, file(params.probes_bed))
            def vcf_filtered = FILTER_GCNV_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_GCNV_GATHER(vcf_filtered, 'GATK', 'out_GCNV')
                vcf_norm = NORM_GCNV_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["GCNV_filtered", f] })
        }
        if (params.get('dragen_dir', false)) {
            dragen_in = Channel.of('DRAGEN').combine(Channel.fromPath(params.dragen_dir + "/*.vcf*").collect())
            FILTER_DRAGEN_QUALITY(dragen_in)
            def vcf_filtered = FILTER_DRAGEN_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_DRAGEN_GATHER(vcf_filtered, 'DRAGEN', 'out_DRAGEN')
                vcf_norm = NORM_DRAGEN_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["DRAGEN_filtered", f] })
        }
        if (params.get('exomedepth_dir', false)) {
            exomedepth_in = Channel.of('EXOMEDEPTH').combine(Channel.fromPath(params.exomedepth_dir + "/*.vcf*").collect())
            FILTER_EXOMEDEPTH_QUALITY(exomedepth_in)
            def vcf_filtered = FILTER_EXOMEDEPTH_QUALITY.out.filtered_data.map { it[2] }.flatten()
            def vcf_norm = vcf_filtered
            if (!skip_norm) {
                NORM_EXOMEDEPTH_GATHER(vcf_filtered, 'EXOMEDEPTH', 'out_EXOMEDEPTH')
                vcf_norm = NORM_EXOMEDEPTH_GATHER.out.normalised_vcf
            }
            ch = ch.mix(vcf_norm.map { f -> [extract_sample_id_from_vcf(f), f] })
            ch_comp = ch_comp.mix(vcf_filtered.map { f -> ["EXOMEDEPTH_filtered", f] })
        }
        
    emit:
        vcfs = ch
        labeled_vcfs = ch_comp
}

def group_caller_vcfs(vcf_ch) {
    def min_c = params.get('min_callers_per_sample', 1)
    return vcf_ch
        .groupTuple()
        .filter { sample_id, vcfs -> min_c <= 1 ? (vcfs.size() >= 1) : (vcfs.size() >= 2) }
        .ifEmpty { error "No samples with caller VCFs were found. Check that the caller VCF directories contain files with matching sample IDs." }
}

def gather_vcfs() {
    def dir_count = 0
    CALLER_DIR_PARAMS.each { caller_dir ->
        if (params.get(caller_dir, false)) dir_count++
    }
    if (dir_count < 2) {
        exit 1, "Error: You must provide VCF directories for at least TWO different callers (e.g., --canoes_dir and --clamms_dir) to run consensus modules."
    }
    GATHER_AND_FILTER_VCFS()
    return GATHER_AND_FILTER_VCFS.out.vcfs
}

// =====================================================================================
// SUB-WORKFLOW DEFINITIONS
// =====================================================================================























// Helper closure: build the feature_inputs tuple for a single merged VCF.
// collapsed_vcf_f: a file path or [] (empty list)
def build_feature_inputs(sample_id, merged_vcf, collapsed_vcf_f) {
    def parts = []
    def find_norm_vcf = { dir, caller ->
        def candidates = [
            file("${dir}/${sample_id}.${caller}.normalised.qs-filtered.cnv.vcf.gz"),
            file("${dir}/${sample_id}.${caller}.normalised.cnv.vcf.gz"),
            file("${dir}/${sample_id}_${caller.toUpperCase()}.normalised.vcf.gz"),
            file("${dir}/${sample_id}_DRAGEN.normalised.vcf.gz"),
            file("${dir}/${sample_id}.${caller}.normalised.vcf.gz"),
            file("${dir}/${sample_id}.${caller}.cnv.vcf.gz")
        ]
        def match = candidates.find { it.exists() }
        return match ? match.toString() : candidates[0].toString()
    }
    if (params.get('canoes_norm_dir',    false)) parts << "canoes=${find_norm_vcf(params.canoes_norm_dir, 'canoes')}"
    if (params.get('clamms_norm_dir',    false)) parts << "clamms=${find_norm_vcf(params.clamms_norm_dir, 'clamms')}"
    if (params.get('xhmm_norm_dir',      false)) parts << "xhmm=${find_norm_vcf(params.xhmm_norm_dir, 'xhmm')}"
    if (params.get('cnvkit_norm_dir',    false)) parts << "cnvkit=${find_norm_vcf(params.cnvkit_norm_dir, 'cnvkit')}"
    if (params.get('gcnv_norm_dir',      false)) parts << "gatk_gcnv=${find_norm_vcf(params.gcnv_norm_dir, 'gcnv')}"
    if (params.get('dragen_norm_dir',    false)) parts << "dragen=${find_norm_vcf(params.dragen_norm_dir, 'dragen')}"
    if (params.get('exomedepth_norm_dir',false)) parts << "exomedepth=${find_norm_vcf(params.exomedepth_norm_dir, 'exomedepth')}"
    def tool_vcfs_str = parts.join(',')

    def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
    def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
    def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
    def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []
    def mode        = params.get('merger_mode', 'survivor')

    return [sample_id, merged_vcf, collapsed_vcf_f, tool_vcfs_str, mode,
            bam_f, fasta_f, bed_f, map_f]
}











def print_pipeline_help() {
    log.info """
===============================================================================
CNV-Caller-Combiner: Multi-Caller Germline CNV Pipeline & Consensus Framework
===============================================================================

Usage:
    nextflow run cnv-caller-combiner.nf -profile <profile> --workflow <mode> [options]

Workflows:
    --workflow canoes                 Run CANOES CNV caller
    --workflow xhmm                   Run XHMM CNV caller
    --workflow clamms                 Run CLAMMS CNV caller
    --workflow dragen                 Run Illumina DRAGEN CNV pipeline
    --workflow cnvkit                 Run CNVkit pipeline
    --workflow gcnv                   Run GATK gCNV pipeline
    --workflow exomedepth             Run ExomeDepth caller
    --workflow survivor               Run SURVIVOR intra-sample consensus merging
    --workflow truvari                Run Truvari intra-sample consensus merging
    --workflow survivor_with_features Run SURVIVOR merging with feature extraction
    --workflow truvari_with_features  Run Truvari merging with feature extraction
    --workflow feature_extraction     Extract ML features from merged VCFs
    --workflow normalise              Harmonize and scale caller quality scores
    --workflow train                  Train CatBoost / XGBoost classifiers
    --workflow evaluate               Evaluate callset against truth set BED
    --workflow full                   Run end-to-end consensus & ML pipeline
    --workflow local                  Run end-to-end local caller combiner
    --workflow full_with_dragen       Run end-to-end including DRAGEN
    --workflow qc                     Run target-aware quality control
    --workflow cram_to_bam            Convert CRAM files to BAM format
    --workflow bam_to_cram            Convert BAM files to CRAM format

Options:
    --outdir <dir>                    Output directory [defaults to ./results]
    --ml_pathway <A|B>                ML pathway: 'A' (AF -> ML) or 'B' (ML -> AF) [defaults to A]
    --skip_ml_training                Skip CatBoost/XGBoost training and prediction
    --skip_cohort_merge               Skip inter-sample cohort merging
    --skip_af_filtering               Skip internal allele frequency filtering
    --skip_annotation                 Skip AnnotSV variant annotation
    --help                            Print this help message
"""
}

// =====================================================================================
// WORKFLOW SWITCH - FULLY INDEPENDENT EXECUTION BLOCKS
// =====================================================================================

workflow RUN_BENCHMARK_PIPELINE {
            def caller_vcf_channels = []
            def mergerMode = params.get('merger_mode', 'survivor')
            def raw_vcfs_comp = Channel.empty()
            def filtered_vcfs_comp = Channel.empty()
            def normalised_vcfs_comp = Channel.empty()

            def ch_samplesheet_bams = Channel.empty()
            def ch_unsliced_bams = Channel.empty()
            def ref_fa_common = file(params.get('ref', params.get('reference_fasta', params.get('fasta', ''))))
            def ref_fai_common = file(params.get('fai', "${ref_fa_common}.fai"))

            if (params.get('samplesheet_bams', false)) {
                def ss_file = file(params.samplesheet_bams)
                def first_row = ss_file.exists() ? ss_file.readLines().find { it.trim() && !it.startsWith('SampleID') } : null
                if (first_row && (first_row.contains('.cram') || first_row.endsWith('.cram'))) {
                    // CRAM paths supplied in samplesheet_bams: automatically convert on the fly!
                    def ch_input_crams = Channel.fromPath(params.samplesheet_bams)
                        .splitCsv(header: true, sep: '\t')
                        .map { row -> [ "${row.SampleID}".replaceAll(/\.md$/, ''), file("${row.BAM}"), file("${row.BAM}.crai") ] }
                    CONVERT_CRAM_TO_BAM_VERIFIED(ch_input_crams, ref_fa_common, ref_fai_common)
                    ch_samplesheet_bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
                    ch_unsliced_bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
                } else {
                    ch_samplesheet_bams = Channel.fromPath(params.samplesheet_bams)
                        .splitCsv(header: true, sep: '\t')
                        .map { row ->
                            def sample_id = "${row.SampleID}".replaceAll(/\.md$/, '')
                            def bam_file = file("${row.BAM}")
                            def bai_path = row.BAI ?: (file("${row.BAM}.bai").exists() ? "${row.BAM}.bai" : "${row.BAM}".replaceAll(/\.bam$/, '.bai'))
                            def bai_file = file(bai_path)
                            [ sample_id, bam_file, bai_file ]
                        }
                    ch_unsliced_bams = ch_samplesheet_bams
                }
            } else if (params.get('samplesheet_crams', false)) {
                // Automated CRAM-to-BAM on-the-fly conversion with 30 forensic verification gates
                def ch_input_crams = Channel.fromPath(params.samplesheet_crams)
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> [ "${row.SampleID}".replaceAll(/\.md$/, ''), file("${row.CRAM}"), file(row.CRAI ?: "${row.CRAM}.crai") ] }
                def target_bed_for_slicing = params.get('target_bed', params.get('probes', params.get('exome_targets', params.get('targets', false))))
                if (target_bed_for_slicing && !params.get('disable_target_slicing', false)) {
                    CONVERT_CRAM_TO_TARGET_PADDED_BAM(ch_input_crams, ref_fa_common, ref_fai_common, file(target_bed_for_slicing))
                    ch_samplesheet_bams = CONVERT_CRAM_TO_TARGET_PADDED_BAM.out.bam_and_index
                    // For CNVkit, ensure unsliced alignment is provided so off-target (antitarget) depth calculations are preserved
                    CONVERT_CRAM_TO_BAM_VERIFIED(ch_input_crams, ref_fa_common, ref_fai_common)
                    ch_unsliced_bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
                } else {
                    CONVERT_CRAM_TO_BAM_VERIFIED(ch_input_crams, ref_fa_common, ref_fai_common)
                    ch_samplesheet_bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
                    ch_unsliced_bams = CONVERT_CRAM_TO_BAM_VERIFIED.out.bam_and_index
                }
            }

            def should_run_canoes = is_caller_enabled('canoes', workflow_mode)
            if (should_run_canoes && (params.get('fai', false) || params.get('ref', false) || params.get('fasta', false)) && (params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false))) {
                def chroms = (1..22).toList().collect { 'chr' + "${it}" }
                def canoes_bam_list = ch_samplesheet_bams.collectFile() { item -> [ 'canoes_bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
                CANOES(canoes_bam_list, file(params.get('fai', "${ref_fa_common}.fai")), Channel.fromList(chroms))
                def canoes_raw = CANOES.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(canoes_raw.map { f -> ["CANOES_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    FILTER_CANOES_FULL(Channel.of('CANOES_full').combine(canoes_raw.collect()))
                    def canoes_filtered = FILTER_CANOES_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(canoes_filtered.map { f -> ["CANOES_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_CANOES_FULL(canoes_filtered, 'CANOES', 'out_CANOES')
                        caller_vcf_channels << NORM_CANOES_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_CANOES_FULL.out.normalised_vcf.map { f -> ["CANOES", f] })
                    } else {
                        caller_vcf_channels << canoes_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(canoes_filtered.map { f -> ["CANOES", f] })
                    }
                } else {
                    caller_vcf_channels << canoes_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(canoes_raw.map { f -> ["CANOES", f] })
                }
            }

            def should_run_clamms = is_caller_enabled('clamms', workflow_mode)
            if (should_run_clamms && (params.get('fai', false) || params.get('ref', false) || params.get('fasta', false)) && (params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false))) {
                def clamms_sample_list = ch_samplesheet_bams.map { it -> it[0] + '\n' }.collectFile(name: 'clamms_sample_list.txt')
                CLAMMS(ch_samplesheet_bams, file(params.get('fai', "${ref_fa_common}.fai")), clamms_sample_list)
                def clamms_raw = CLAMMS.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(clamms_raw.map { f -> ["CLAMMS_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    FILTER_CLAMMS_FULL(Channel.of('CLAMMS_full').combine(clamms_raw.collect()))
                    def clamms_filtered = FILTER_CLAMMS_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(clamms_filtered.map { f -> ["CLAMMS_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_CLAMMS_FULL(clamms_filtered, 'CLAMMS', 'out_CLAMMS')
                        caller_vcf_channels << NORM_CLAMMS_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_CLAMMS_FULL.out.normalised_vcf.map { f -> ["CLAMMS", f] })
                    } else {
                        caller_vcf_channels << clamms_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(clamms_filtered.map { f -> ["CLAMMS", f] })
                    }
                } else {
                    caller_vcf_channels << clamms_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(clamms_raw.map { f -> ["CLAMMS", f] })
                }
            }

            def should_run_xhmm = is_caller_enabled('xhmm', workflow_mode)
            if (should_run_xhmm && (params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false))) {
                XHMM(ch_samplesheet_bams)
                def xhmm_raw = XHMM.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(xhmm_raw.map { f -> ["XHMM_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    FILTER_XHMM_FULL(Channel.of('XHMM_full').combine(xhmm_raw.collect()))
                    def xhmm_filtered = FILTER_XHMM_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(xhmm_filtered.map { f -> ["XHMM_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_XHMM_FULL(xhmm_filtered, 'XHMM', 'out_XHMM')
                        caller_vcf_channels << NORM_XHMM_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_XHMM_FULL.out.normalised_vcf.map { f -> ["XHMM", f] })
                    } else {
                        caller_vcf_channels << xhmm_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(xhmm_filtered.map { f -> ["XHMM", f] })
                    }
                } else {
                    caller_vcf_channels << xhmm_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(xhmm_raw.map { f -> ["XHMM", f] })
                }
            }

            def should_run_cnvkit = is_caller_enabled('cnvkit', workflow_mode)
            def cnvkit_targets_path = params.get('targets', params.get('exome_targets', params.get('probes', params.get('probes_bed', ''))))
            def cnvkit_fasta_path = params.get('fasta', params.get('ref', params.get('reference_fasta', '')))
            def has_cnvkit_inputs = should_run_cnvkit && (params.get('bams', false) || params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false)) && cnvkit_fasta_path && cnvkit_targets_path && params.get('refflat', false)

            if (has_cnvkit_inputs) {
                def ch_bams_cnvkit_resolved
                if (params.get('bams', false)) {
                    def glob_pattern_full = params.bams.contains('{') ? params.bams : params.bams.replaceAll(/\.bam$/, '.{bam,bai,bam.bai}')
                    ch_bams_cnvkit_resolved = Channel.fromPath(glob_pattern_full)
                        .map { it -> 
                            def sample_id = it.name.replaceAll(/\.bam\.bai$/, '').replaceAll(/\.bam$/, '').replaceAll(/\.bai$/, '').replaceAll(/\.md$/, '')
                            [sample_id, it] 
                        }
                        .groupTuple(size: 2)
                        .map { id, files -> [id, files.find { it.extension == 'bam' }, files.find { it.extension == 'bai' }] }
                        .ifEmpty { error "full workflow CNVkit step received no BAMs. Check the --bams parameter value matches existing BAM files." }
                } else {
                    // CNVkit requires unsliced alignments so that antitarget (off-target) coverage is preserved
                    ch_bams_cnvkit_resolved = ch_unsliced_bams
                }
                CNVKIT(
                    ch_bams_cnvkit_resolved,
                    file(cnvkit_fasta_path),
                    file(cnvkit_targets_path),
                    file(params.refflat),
                    ch_bams_cnvkit_resolved.first().map { it[1] }
                )
                def cnvkit_raw = CNVKIT.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(cnvkit_raw.map { f -> ["CNVKIT_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    FILTER_CNVKIT_FULL(Channel.of('CNVKIT_full').combine(cnvkit_raw.collect()))
                    def cnvkit_filtered = FILTER_CNVKIT_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(cnvkit_filtered.map { f -> ["CNVKIT_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_CNVKIT_FULL(cnvkit_filtered, 'CNVKIT', 'out_CNVKIT')
                        caller_vcf_channels << NORM_CNVKIT_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_CNVKIT_FULL.out.normalised_vcf.map { f -> ["CNVKIT", f] })
                    } else {
                        caller_vcf_channels << cnvkit_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(cnvkit_filtered.map { f -> ["CNVKIT", f] })
                    }
                } else {
                    caller_vcf_channels << cnvkit_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(cnvkit_raw.map { f -> ["CNVKIT", f] })
                }
            }

            def should_run_gatk = is_caller_enabled('gatk', workflow_mode)
            def gatk_targets_path = params.get('exome_targets', params.get('targets', params.get('probes', params.get('probes_bed', ''))))
            def gatk_fasta_path = params.get('fasta', params.get('ref', params.get('reference_fasta', '')))
            def gatk_fai_path = params.get('fai', "${gatk_fasta_path}.fai")
            def has_gatk_inputs = should_run_gatk && !params.skip_gcnv_training && (params.get('samples_path', false) || params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false)) && gatk_fasta_path && gatk_fai_path && params.get('dict', false) && gatk_targets_path

            if (has_gatk_inputs) {
                def ch_bams_gcnv_resolved
                if (params.get('samples_path', false)) {
                    ch_bams_gcnv_resolved = Channel.fromPath(params.samples_path)
                        .filter { it.name.endsWith('.bam') || it.name.endsWith('.cram') }
                        .map { it ->
                            def index_fallback = [file("${it}.bai"), file("${it}.bam.bai"), file(it.toString().replaceAll(/\.bam$/, '.bai'))]
                                .find { f -> f.exists() }
                            def index_file = index_fallback ? index_fallback : (it.name.endsWith('.bam') ? file("${it}.bai") : file("${it}.crai"))
                            return [ it.baseName.replaceAll(/\.md$/, ''), it, index_file ]
                        }
                } else {
                    ch_bams_gcnv_resolved = ch_samplesheet_bams
                }
                GATK_GCNV(
                    ch_bams_gcnv_resolved,
                    file(gatk_fasta_path),
                    file(gatk_fai_path),
                    file(params.dict),
                    file(gatk_targets_path)
                )
                def gcnv_raw = GATK_GCNV.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(gcnv_raw.map { f -> ["GCNV_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    def gcnv_probes_path = params.probes_bed ?: params.probes ?: params.targets ?: gatk_targets_path
                    def gcnv_probes = file(gcnv_probes_path)
                    FILTER_GCNV_FULL(Channel.of('GCNV_full').combine(gcnv_raw.collect()), gcnv_probes)
                    def gcnv_filtered = FILTER_GCNV_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(gcnv_filtered.map { f -> ["GCNV_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_GCNV_FULL(gcnv_filtered, 'GATK', 'out_GCNV')
                        caller_vcf_channels << NORM_GCNV_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_GCNV_FULL.out.normalised_vcf.map { f -> ["GCNV", f] })
                    } else {
                        caller_vcf_channels << gcnv_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(gcnv_filtered.map { f -> ["GCNV", f] })
                    }
                } else {
                    caller_vcf_channels << gcnv_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(gcnv_raw.map { f -> ["GCNV", f] })
                }
            }

            def should_run_dragen = is_caller_enabled('dragen', workflow_mode)
            if (should_run_dragen && workflow_mode != 'local') {
                if (params.get('cramFilePairsUploadPath', false)) {
                    Channel.fromPath(params.cramFilePairsUploadPath, checkIfExists: true)
                        .filter { f -> f.name.endsWith('.cram') || f.name.endsWith('.bam') }
                        .map { primary ->
                            def sampleId = primary.name.replaceAll(/\.(cram|bam)$/, '').replaceAll(/\.md$/, '')
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
                    def dragen_raw = DRAGEN.out.sorted_vcf.flatten()
                    raw_vcfs_comp = raw_vcfs_comp.mix(dragen_raw.map { f -> ["DRAGEN_raw", f] })
                    if (!params.get('skip_quality_filtering', false)) {
                        FILTER_DRAGEN_FULL(Channel.of('DRAGEN_full').combine(dragen_raw.collect()))
                        def dragen_filtered = FILTER_DRAGEN_FULL.out.filtered_data.map { it[2] }.flatten()
                        filtered_vcfs_comp = filtered_vcfs_comp.mix(dragen_filtered.map { f -> ["DRAGEN_filtered", f] })
                        if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                            NORM_DRAGEN_FULL(dragen_filtered, 'DRAGEN', 'out_DRAGEN')
                            caller_vcf_channels << NORM_DRAGEN_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                            normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_DRAGEN_FULL.out.normalised_vcf.map { f -> ["DRAGEN", f] })
                        } else {
                            caller_vcf_channels << dragen_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                            normalised_vcfs_comp = normalised_vcfs_comp.mix(dragen_filtered.map { f -> ["DRAGEN", f] })
                        }
                    } else {
                        caller_vcf_channels << dragen_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(dragen_raw.map { f -> ["DRAGEN", f] })
                    }
                }
            }

            def should_run_exomedepth = is_caller_enabled('exomedepth', workflow_mode)
            def exomedepth_probes_path = params.get('probes', params.get('targets', params.get('exome_targets', params.get('probes_bed', ''))))

            if (should_run_exomedepth && (params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false)) && exomedepth_probes_path) {
                def exomedepth_probes = file(exomedepth_probes_path)
                def exomedepth_bam_list = ch_samplesheet_bams.collectFile() { item -> [ 'exomedepth_bam_list.txt', "${item[1]}" + '\n' ] }
                EXOMEDEPTH(exomedepth_bam_list, exomedepth_probes, file(params.get('fai', "${ref_fa_common}.fai")))
                def exomedepth_raw = EXOMEDEPTH.out.sorted_vcf.flatten()
                raw_vcfs_comp = raw_vcfs_comp.mix(exomedepth_raw.map { f -> ["EXOMEDEPTH_raw", f] })
                if (!params.get('skip_quality_filtering', false)) {
                    FILTER_EXOMEDEPTH_FULL(Channel.of('EXOMEDEPTH_full').combine(exomedepth_raw.collect()))
                    def exomedepth_filtered = FILTER_EXOMEDEPTH_FULL.out.filtered_data.map { it[2] }.flatten()
                    filtered_vcfs_comp = filtered_vcfs_comp.mix(exomedepth_filtered.map { f -> ["EXOMEDEPTH_filtered", f] })
                    if (!params.get('skip_normalisation', false) && !params.get('skip_quality_normalisation', false)) {
                        NORM_EXOMEDEPTH_FULL(exomedepth_filtered, 'EXOMEDEPTH', 'out_EXOMEDEPTH')
                        caller_vcf_channels << NORM_EXOMEDEPTH_FULL.out.normalised_vcf.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(NORM_EXOMEDEPTH_FULL.out.normalised_vcf.map { f -> ["EXOMEDEPTH", f] })
                    } else {
                        caller_vcf_channels << exomedepth_filtered.map { f -> [extract_sample_id_from_vcf(f), f] }
                        normalised_vcfs_comp = normalised_vcfs_comp.mix(exomedepth_filtered.map { f -> ["EXOMEDEPTH", f] })
                    }
                } else {
                    caller_vcf_channels << exomedepth_raw.map { f -> [extract_sample_id_from_vcf(f), f] }
                    normalised_vcfs_comp = normalised_vcfs_comp.mix(exomedepth_raw.map { f -> ["EXOMEDEPTH", f] })
                }
            }



            if (caller_vcf_channels.size() < 1) {
                if (workflow_mode == 'dragen' || params.get('cramFilePairsUploadPath', false) || params.get('dragen_sv_dir', false)) {
                    log.info "DRAGEN execution active without CNV caller output or 0 CNV callers configured. Bypassing consensus merging and evaluation."
                } else {
                    exit 1, "Error: workflow requires at least one CNV caller input. Please configure parameters for at least 1 caller."
                }
            }

            def is_pathway_b = !(params.get('ml_pathway', 'B') ?: 'B').toString().trim().equalsIgnoreCase('A')
            def ch_intrasample_merged = Channel.empty()
            def ch_intrasample_tbi = Channel.empty()
            def ch_refined_intrasample = Channel.empty()
            def ch_annotated_vcf = Channel.empty()
            def ch_acmg_vcfs = Channel.empty()

            if (caller_vcf_channels.size() >= 1) {
                def ch_vcfs = caller_vcf_channels[0]
                caller_vcf_channels.drop(1).each { ch_part ->
                    ch_vcfs = ch_vcfs.mix(ch_part)
                }

                grouped_vcfs = group_caller_vcfs(ch_vcfs)

                def bam_f       = params.get('bam_file',         false) ? file(params.bam_file)         : []
                def fasta_f     = params.get('reference_fasta',  false) ? file(params.reference_fasta)  : []
                def bed_f       = params.get('bed_file',         false) ? file(params.bed_file)         : []
                def map_f       = params.get('mappability_file', false) ? file(params.mappability_file) : []

                // -----------------------------------------------------------------------------
                // CALLER-INDEPENDENT TARGET COVERAGE & LOG2 CALCULATION
                // -----------------------------------------------------------------------------
                def ch_coverage_cnn = Channel.empty()
                def ch_cov_targets = bed_f ?: (params.get('probes_bed', params.get('probes', params.get('target_bed', false))) ? file(params.get('probes_bed', params.get('probes', params.get('target_bed')))) : [])

                if (params.get('run_cnvkit', false)) {
                    ch_coverage_cnn = CNVKIT.out.target_cov
                } else if ((params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false)) && ch_cov_targets && file(ch_cov_targets.toString()).exists()) {
                    COMPUTE_TARGET_COVERAGE(ch_samplesheet_bams, ch_cov_targets)
                    ch_coverage_cnn = COMPUTE_TARGET_COVERAGE.out.coverage_cnn
                }

                // -----------------------------------------------------------------------------
                // STAGE 5: INTRA-SAMPLE, INTER-TOOL MERGING (SURVIVOR or Truvari)
                // -----------------------------------------------------------------------------
                ch_intrasample_merged = Channel.empty()
                ch_intrasample_tbi = Channel.empty()
                if (mergerMode == 'truvari') {
                    TRUVARI(grouped_vcfs)
                    ch_intrasample_merged = TRUVARI.out.merged_vcf
                    ch_intrasample_tbi = TRUVARI.out.merged_tbi
                    feature_inputs_ch = TRUVARI.out.merged_vcf
                        .join(TRUVARI.out.collapsed_vcf)
                        .join(grouped_vcfs)
                        .join(ch_coverage_cnn, remainder: true)
                        .filter { item ->
                            item.size() >= 4 && item[1] != null && item[1] != []
                        }
                        .map { item ->
                            def sample_id     = item[0]
                            def merged_vcf    = item[1]
                            def collapsed_vcf = item[2]
                            def vcfs          = item[3]
                            def cov_file      = (item.size() > 4 && item[4] != null) ? item[4] : []
                            [sample_id, merged_vcf, collapsed_vcf, build_tool_vcfs_str(vcfs), 'truvari',
                             bam_f, fasta_f, bed_f, map_f, cov_file]
                        }
                } else {
                    SURVIVOR(grouped_vcfs)
                    ch_intrasample_merged = SURVIVOR.out.union_vcf
                    ch_intrasample_tbi = SURVIVOR.out.union_tbi
                    feature_inputs_ch = SURVIVOR.out.union_vcf
                        .join(grouped_vcfs)
                        .join(ch_coverage_cnn, remainder: true)
                        .filter { item ->
                            item.size() >= 3 && item[1] != null && item[1] != []
                        }
                        .map { item ->
                            def sample_id  = item[0]
                            def merged_vcf = item[1]
                            def vcfs       = item[2]
                            def cov_file   = (item.size() > 3 && item[3] != null) ? item[3] : []
                            [sample_id, merged_vcf, [], build_tool_vcfs_str(vcfs), 'survivor',
                             bam_f, fasta_f, bed_f, map_f, cov_file]
                        }
                }

                // Optional SV Calling via SvABA if BAMs and target BED provided, or ingestion of DRAGEN SV calls
                def ch_sv_vcfs = Channel.empty()
                if (params.get('run_svaba', false) && (params.get('samplesheet_bams', false) || params.get('samplesheet_crams', false)) && params.get('fasta', false)) {
                    SVABA(ch_samplesheet_bams, file(params.fasta), file(params.fai), file(params.get('exome_targets', params.get('probes_bed', ''))))
                    ch_sv_vcfs = SVABA.out.sv_vcf
                } else if (params.get('dragen_sv_dir', false)) {
                    def dragen_sv_raw = Channel.fromPath("${params.dragen_sv_dir}/*.sv.vcf*")
                        .filter { !it.name.endsWith('.tbi') && !it.name.endsWith('.csi') }
                    BGZIP_SORT_INDEX_SV_VCF(dragen_sv_raw, 'DRAGEN', 'out_DRAGEN/sv_vcfs')
                    ch_sv_vcfs = BGZIP_SORT_INDEX_SV_VCF.out.sv_vcf_tuple
                }

                // -----------------------------------------------------------------------------
                // STAGE 6: POST-CONSENSUS BREAKPOINT REFINEMENT
                // -----------------------------------------------------------------------------
                ch_refined_intrasample = ch_intrasample_merged
                if (params.get('run_refinement', true) && (params.get('run_svaba', false) || params.get('dragen_sv_dir', false))) {
                    def ch_to_refine = ch_intrasample_merged.join(ch_intrasample_tbi).map { id, vcf, tbi -> [id, vcf, tbi] }
                    REFINE_CONSENSUS_BREAKPOINTS(ch_to_refine.join(ch_sv_vcfs))
                    ch_refined_intrasample = REFINE_CONSENSUS_BREAKPOINTS.out.refined_vcf.map { id, vcf, tbi -> [id, vcf] }
                }

                // -----------------------------------------------------------------------------
                // OPTIONAL STAGE: HEADLESS IGV HTML VIEWER GENERATION
                // -----------------------------------------------------------------------------
                if (params.get('generate_igv_reports', false) || params.get('run_igv', false)) {
                    def ch_igv_bams = ch_samplesheet_bams.map { id, bam, bai -> [id, bam, bai] }
                    def ch_igv_inputs = ch_refined_intrasample.join(ch_igv_bams).map { id, vcf, bam, bai -> [id, vcf, bam, bai] }
                    GENERATE_IGV_REPORT(ch_igv_inputs, ref_fa_common, ref_fai_common)
                }

                // -----------------------------------------------------------------------------
                // STAGE 7: INTER-SAMPLE COHORT MERGING (SURVIVOR / Truvari)
                // -----------------------------------------------------------------------------
                def ch_cohort_merged = Channel.empty()
                if (!params.get('skip_cohort_merge', false)) {
                    MERGE_CONSENSUS_COHORT(ch_refined_intrasample.map { it[1] }.collect())
                    ch_cohort_merged = Channel.of('full_cohort').combine(MERGE_CONSENSUS_COHORT.out.merged_vcf)
                }

                FEATURE_EXTRACTION(feature_inputs_ch)

                if (!params.get('skip_cohort_merge', false)) {
                    def ch_ref_fasta_db = params.get('fasta', params.get('reference_fasta', false)) ? file(params.get('fasta', params.get('reference_fasta'))) : file("${projectDir}/NO_FASTA")
                    def ch_target_bed_db = params.get('probes_bed', params.get('target_bed', false)) ? file(params.get('probes_bed', params.get('target_bed'))) : file("${projectDir}/NO_BED")
                    BUILD_COHORT_FEATURE_DATABASE(
                        ch_cohort_merged,
                        FEATURE_EXTRACTION.out.features_tsv.map { it[1] }.collect(),
                        ch_ref_fasta_db,
                        ch_target_bed_db
                    )
                }

                is_pathway_b = !(params.get('ml_pathway', 'B') ?: 'B').toString().trim().equalsIgnoreCase('A')
                def ch_af_filtered = Channel.empty()
                def ch_catboost_predicted = Channel.empty()

                if (!is_pathway_b) {
                    // -------------------------------------------------------------------------
                    // PATHWAY A (Default): AF Filtering -> CatBoost In Silico Prediction
                    // -------------------------------------------------------------------------
                    if (!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) {
                        bcftoolsFlattenCohort(ch_cohort_merged)
                        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
                        def slice_joined = ch_cohort_merged.join(filterRareCnvs.out.rare_data)
                        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
                        ch_af_filtered = filterVcfByRareRegions.out.rare_vcf
                    }

                    if (params.skip_ml_training) {
                        log.info "Skipping machine learning classifier training as requested (skip_ml_training = true)."
                    } else {
                        Channel.value(file(params.truth_labels))
                            .set { ch_truth_full }
                        def probesBedParamFull = params.get('probes_bed', false)
                        Channel.value(probesBedParamFull ? file(probesBedParamFull) : [])
                            .set { ch_probes_full }

                        // Run CatBoost In Silico Prediction on AF-filtered VCF (or cohort merged if AF skipped)
                        def vcfs_for_ml = (!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) ?
                            ch_af_filtered.map { it[1] }.collect() :
                            ch_cohort_merged.map { it[1] }.collect()

                        RUN_CATBOOST_PREDICT(
                            FEATURE_EXTRACTION.out.features_tsv.map { it[1] }.collect(),
                            vcfs_for_ml,
                            ch_truth_full,
                            ch_probes_full
                        )
                        ch_catboost_predicted = RUN_CATBOOST_PREDICT.out.predicted_vcfs
                    }
                } else {
                    // -------------------------------------------------------------------------
                    // PATHWAY B: CatBoost In Silico Prediction -> AF Filtering
                    // -------------------------------------------------------------------------
                    if (params.skip_ml_training) {
                        log.info "Skipping machine learning classifier training as requested (skip_ml_training = true)."
                    } else {
                        Channel.value(file(params.truth_labels))
                            .set { ch_truth_full }
                        def probesBedParamFull = params.get('probes_bed', false)
                        Channel.value(probesBedParamFull ? file(probesBedParamFull) : [])
                            .set { ch_probes_full }

                        // Run CatBoost In Silico Prediction directly on cohort merged VCF
                        def cohort_vcfs_for_ml = ch_cohort_merged.map { it[1] }.collect()
                        RUN_CATBOOST_PREDICT(
                            FEATURE_EXTRACTION.out.features_tsv.map { it[1] }.collect(),
                            cohort_vcfs_for_ml,
                            ch_truth_full,
                            ch_probes_full
                        )
                        ch_catboost_predicted = RUN_CATBOOST_PREDICT.out.predicted_vcfs
                    }

                    if (!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) {
                        def vcf_to_filter = (!params.skip_ml_training && params.get('truth_labels', false)) ?
                            RUN_CATBOOST_PREDICT.out.predicted_vcfs.flatten().filter { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.csi') }.map { f -> ['full_cohort', f] } :
                            ch_cohort_merged
                        bcftoolsFlattenCohort(vcf_to_filter)
                        filterRareCnvs(bcftoolsFlattenCohort.out.cohort_tsv)
                        def slice_joined = vcf_to_filter.join(filterRareCnvs.out.rare_data)
                        filterVcfByRareRegions(slice_joined.map { [ it[0], it[1], it[3] ] })
                        ch_af_filtered = filterVcfByRareRegions.out.rare_vcf
                    }
                }

                // -----------------------------------------------------------------------------
                // STAGE 12 & 13: ANNOTSV VARIANT ANNOTATION & ACMG STRATIFICATION
                // -----------------------------------------------------------------------------
                ch_annotated_vcf = Channel.empty()
                ch_acmg_vcfs = Channel.empty()
                if (!params.get('skip_annotation', false) && params.get('annotsv_db', false)) {
                    def annotsv_db_path = file(params.annotsv_db)
                    def vcf_for_annotsv = is_pathway_b ?
                        ((!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) ?
                            ch_af_filtered :
                            ((!params.skip_ml_training && params.get('truth_labels', false)) ?
                                RUN_CATBOOST_PREDICT.out.predicted_vcfs.flatten().filter { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.csi') }.map { f -> ['full_cohort', f] } :
                                ch_cohort_merged)) :
                        ((!params.skip_ml_training && params.get('truth_labels', false)) ?
                            RUN_CATBOOST_PREDICT.out.predicted_vcfs.flatten().filter { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.csi') }.map { f -> ['full_cohort', f] } :
                            ((!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) ?
                                ch_af_filtered :
                                ch_cohort_merged))
                    runAnnotSV(vcf_for_annotsv, annotsv_db_path)
                    if (!params.get('skip_af_filtering', false)) {
                        def acmg_joined = filterRareCnvs.out.rare_data.join(runAnnotSV.out.annotated_tsv)
                        addACMGAnnotations(acmg_joined.map { [ it[0], it[1], it[3] ] })
                        def vcf_acmg_joined = vcf_for_annotsv.join(runAnnotSV.out.annotated_tsv)
                        stratifyVcfByAcmg(vcf_acmg_joined)
                        ch_annotated_vcf = stratifyVcfByAcmg.out.annotated_vcf
                        ch_acmg_vcfs = stratifyVcfByAcmg.out.acmg_vcfs

                        // Quality assessment
                        def resolved_eval_probes = file(params.probes_bed ?: params.probes ?: params.targets)
                        def metrics_inputs = addACMGAnnotations.out.rare_cnvs_annotated
                            .join(runAnnotSV.out.annotated_tsv)
                            .combine(ch_refined_intrasample.map { it[1] }.collect().map { [ it ] })
                            .map { run_id, rare_csv, annot_tsv, vcfs ->
                                [ run_id, rare_csv, annot_tsv, vcfs ]
                            }
                        ANALYZE_CNV_QUALITY_METRICS(
                            metrics_inputs,
                            resolved_eval_probes
                        )

                        // ACMG Evaluators
                        if (params.get('map_file', false) && params.get('truth_bed', false)) {
                            def eval_inputs = addACMGAnnotations.out.score_separated_csvs.transpose()
                            def eval_joined = eval_inputs.combine(Channel.of(resolved_eval_probes))
                            EVALUATE_ACMG_CSV(
                                eval_joined.map { [ it[0], it[1], it[2] ] },
                                file(params.map_file),
                                file(params.truth_bed)
                            )
                            PLOT_EVALUATION_METRICS(EVALUATE_ACMG_CSV.out.report.collect())
                        }
                    }
                }
            }

            // -----------------------------------------------------------------------------
            // MULTI-STAGE EVALUATION & QUALITY ASSESSMENT ACROSS ALL 10 CALLSETS
            // -----------------------------------------------------------------------------
            def qc_str = params.get('run_qc', false) ? "QC" : "NoQC"
            def is_mosdepth = params.get('probes_bed', '').contains('mosdepth') || params.get('standard_bed', '').contains('mosdepth')
            def bed_str = is_mosdepth ? "Mosdepth" : "Control"
            def merger_str = mergerMode.toUpperCase()
            def label_suffix = "_${bed_str}_${qc_str}_${merger_str}"

            def all_evals = Channel.empty()
            def all_quality_evals = Channel.empty()

            // Stage 1: Raw CNV Callset
            all_evals = all_evals.mix(raw_vcfs_comp.map { label, f -> ["1_raw_${label}${label_suffix}", f] })
            all_quality_evals = all_quality_evals.mix(raw_vcfs_comp.map { label, f -> ["1_raw_${label}", f] })

            // Stage 2: Quality-Filtered CNV Callset
            all_evals = all_evals.mix(filtered_vcfs_comp.map { label, f -> ["2_quality_filtered_${label}${label_suffix}", f] })
            all_quality_evals = all_quality_evals.mix(filtered_vcfs_comp.map { label, f -> ["2_quality_filtered_${label}", f] })

            // Stage 3: Normalised CNV Callset
            all_evals = all_evals.mix(normalised_vcfs_comp.map { label, f -> ["3_normalised_${label}${label_suffix}", f] })
            all_quality_evals = all_quality_evals.mix(normalised_vcfs_comp.map { label, f -> ["3_normalised_${label}", f] })

            if (caller_vcf_channels.size() >= 1) {
                // Stage 4: Intra-Sample Inter-Tool Merged CNV Callset
                all_evals = all_evals.mix(ch_intrasample_merged.map { sample_id, f -> ["4_merged_intrasample${label_suffix}", f] })
                all_quality_evals = all_quality_evals.mix(ch_intrasample_merged.map { sample_id, f -> ["4_merged_intrasample", f] })

                // Stage 5: Refined-Breakpoint CNV Callset
                if (params.get('run_refinement', true) && (params.get('run_svaba', false) || params.get('dragen_sv_dir', false))) {
                    all_evals = all_evals.mix(ch_refined_intrasample.map { sample_id, f -> ["5_refined_intrasample${label_suffix}", f] })
                    all_quality_evals = all_quality_evals.mix(ch_refined_intrasample.map { sample_id, f -> ["5_refined_intrasample", f] })
                }

                // Stage 6: Inter-Sample Cohort Merged CNV Callset
                if (!params.get('skip_cohort_merge', false)) {
                    all_evals = all_evals.mix(MERGE_CONSENSUS_COHORT.out.merged_vcf.map { f -> ["6_cohort_merged${label_suffix}", f] })
                    all_quality_evals = all_quality_evals.mix(MERGE_CONSENSUS_COHORT.out.merged_vcf.map { f -> ["6_cohort_merged", f] })
                }

                if (is_pathway_b) {
                    // Pathway B: CatBoost (Analytical Validity) precedes Allele Frequency Filtering (Clinical Utility)
                    // Stage 7: in silico Validated CNV Callset (CatBoost predictions on cohort-merged consensus)
                    if (!params.skip_ml_training && params.get('truth_labels', false)) {
                        def ch_catboost_vcfs = RUN_CATBOOST_PREDICT.out.predicted_vcfs.flatten().filter { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.csi') }
                        all_evals = all_evals.mix(ch_catboost_vcfs.map { f -> ["7_insilico_validated${label_suffix}", f] })
                        all_quality_evals = all_quality_evals.mix(ch_catboost_vcfs.map { f -> ["7_insilico_validated", f] })
                    }

                    // Stage 8: Internal Allele Frequency-Filtered CNV Callset (Rare Disease Diagnostic Track)
                    if (!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) {
                        all_evals = all_evals.mix(filterVcfByRareRegions.out.rare_vcf.map { run_name, f -> ["8_af_filtered${label_suffix}", f] })
                        all_quality_evals = all_quality_evals.mix(filterVcfByRareRegions.out.rare_vcf.map { run_name, f -> ["8_af_filtered", f] })
                    }
                } else {
                    // Pathway A: Allele Frequency Filtering precedes CatBoost
                    // Stage 7: Internal Allele Frequency-Filtered CNV Callset
                    if (!params.get('skip_cohort_merge', false) && !params.get('skip_af_filtering', false)) {
                        all_evals = all_evals.mix(filterVcfByRareRegions.out.rare_vcf.map { run_name, f -> ["7_af_filtered${label_suffix}", f] })
                        all_quality_evals = all_quality_evals.mix(filterVcfByRareRegions.out.rare_vcf.map { run_name, f -> ["7_af_filtered", f] })
                    }

                    // Stage 8: in silico Validated CNV Callset (CatBoost predictions)
                    if (!params.skip_ml_training && params.get('truth_labels', false)) {
                        def ch_catboost_vcfs = RUN_CATBOOST_PREDICT.out.predicted_vcfs.flatten().filter { f -> !f.name.endsWith('.tbi') && !f.name.endsWith('.csi') }
                        all_evals = all_evals.mix(ch_catboost_vcfs.map { f -> ["8_insilico_validated${label_suffix}", f] })
                        all_quality_evals = all_quality_evals.mix(ch_catboost_vcfs.map { f -> ["8_insilico_validated", f] })
                    }
                }

                // Stage 9: Variant Annotated CNV Callset
                all_evals = all_evals.mix(ch_annotated_vcf.map { run_name, f -> ["9_variant_annotated${label_suffix}", f] })
                all_quality_evals = all_quality_evals.mix(ch_annotated_vcf.map { run_name, f -> ["9_variant_annotated", f] })

                // Stage 10: ACMG_{4_5, 5, 4, 3, 2, 1, NA} CNV Callsets
                def ch_acmg_eval = ch_acmg_vcfs.map { run_name, files -> files }.flatten().map { f ->
                    def score = f.name.replaceAll(/^ACMG/, '').replaceAll(/\..*/, '')
                    ["10_ACMG_${score}${label_suffix}", f]
                }
                def ch_acmg_qual = ch_acmg_vcfs.map { run_name, files -> files }.flatten().map { f ->
                    def score = f.name.replaceAll(/^ACMG/, '').replaceAll(/\..*/, '')
                    ["10_ACMG_${score}", f]
                }
                all_evals = all_evals.mix(ch_acmg_eval)
                all_quality_evals = all_quality_evals.mix(ch_acmg_qual)
            }

            // Quality Metric Extraction across all 10 callsets
            def probes_f_qual = params.get('probes_bed', params.get('probes', params.get('target_bed', false))) ? file(params.get('probes_bed', params.get('probes', params.get('target_bed')))) : file("${projectDir}/NO_BED")
            def grouped_quality_inputs = all_quality_evals.groupTuple()
            ASSESS_CALLSET_QUALITY_METRICS(grouped_quality_inputs, probes_f_qual)
            COMPILE_ALL_STAGES_QUALITY_REPORT(ASSESS_CALLSET_QUALITY_METRICS.out.summary_file.collect())

            // Performance Metric Evaluation (vs gold truth BED)
            def comp_probes_path = params.probes_bed ?: params.probes ?: params.targets
            if (params.get('truth_bed', false) && comp_probes_path) {
                def truth_f = file(params.truth_bed)
                def probes_f = file(comp_probes_path)
                EVALUATE_COMPARATIVE(all_evals, truth_f, probes_f)
                def comp_reports_ch = EVALUATE_COMPARATIVE.out.performance_reports.map { it instanceof List ? it[-1] : it }
                GENERATE_ATTRITION_FUNNEL(comp_reports_ch.collect())
            }
}

workflow {
    if (params.get('help', false)) {
        print_pipeline_help()
        return
    }
    validate_required_params(workflow_mode)
    validate_runtime_inputs(workflow_mode)
    if (is_validation_only_mode()) {
        log.info "Validation-only execution mode detected (stub-run); exiting before process execution."
        return
    }
    switch (params.workflow ?: workflow_mode) {

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

        case ['cram_to_bam', 'cram_to_target_bam']:
            def ref_fa_c2b = file(params.get('ref', params.get('reference_fasta', params.get('fasta', ''))))
            def ref_fai_c2b = file("${ref_fa_c2b}.fai")
            if (params.get('samplesheet_crams', false)) {
                Channel.fromPath(params.samplesheet_crams)
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> [ "${row.SampleID}".replaceAll(/\\.md\$/, ''), file("${row.CRAM}"), file(row.CRAI ?: "${row.CRAM}.crai") ] }
                    .set { ch_cram_pairs }
            } else if (params.get('input_dir', false)) {
                Channel.fromFilePairs("${params.input_dir}/*.{cram,cram.crai}", size: 2)
                    .map { id, files -> [ id.replaceAll(/\\.md\$/, ''), files[0].name.endsWith('.cram') ? files[0] : files[1], files[0].name.endsWith('.crai') ? files[0] : files[1] ] }
                    .set { ch_cram_pairs }
            } else {
                exit 1, "Please provide --samplesheet_crams or --input_dir for cram_to_bam conversion."
            }
            if (workflow_mode == 'cram_to_target_bam' || params.get('target_bed', params.get('probes', false))) {
                def target_bed_c2b = file(params.get('target_bed', params.get('probes', '')))
                RUN_CRAM_TO_TARGET_PADDED_BAM_CONVERSION(ch_cram_pairs, ref_fa_c2b, ref_fai_c2b, target_bed_c2b)
            } else {
                RUN_CRAM_TO_BAM_CONVERSION(ch_cram_pairs, ref_fa_c2b, ref_fai_c2b)
            }
            break

        case['bam_to_cram']:
            def ref_fa_b2c = file(params.get('ref', params.get('reference_fasta', params.get('fasta', ''))))
            def ref_fai_b2c = file("${ref_fa_b2c}.fai")
            if (params.get('samplesheet_bams', false)) {
                Channel.fromPath(params.samplesheet_bams)
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> [ "${row.SampleID}".replaceAll(/\\.md\$/, ''), file("${row.BAM}"), file("${row.BAM}".replaceAll(/\\.bam\$/, '.bam.bai')) ] }
                    .set { ch_bam_pairs }
            } else if (params.get('input_dir', false)) {
                Channel.fromFilePairs("${params.input_dir}/*.{bam,bam.bai}", size: 2)
                    .map { id, files -> [ id.replaceAll(/\\.md\$/, ''), files[0].name.endsWith('.bam') ? files[0] : files[1], files[0].name.endsWith('.bai') ? files[0] : files[1] ] }
                    .set { ch_bam_pairs }
            } else {
                exit 1, "Please provide --samplesheet_bams or --input_dir for bam_to_cram conversion."
            }
            RUN_BAM_TO_CRAM_CONVERSION(ch_bam_pairs, ref_fa_b2c, ref_fai_b2c)
            break
        


        case['canoes']:
            def canoes_chroms = (1..22).toList().collect { 'chr' + "${it}" }
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}".replaceAll(/\.md$/, ''), "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "canoes workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_CANOES(ch_bams, canoes_chroms, file(params.fai))
            break

        case['xhmm']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}".replaceAll(/\.md$/, ''), "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "xhmm workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_XHMM(ch_bams)
            break

        case['clamms']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}".replaceAll(/\.md$/, ''), "${row.BAM}", "${row.BAM}".replaceAll(/\.bam$/, '.bam.bai') ] }
                .ifEmpty { error "clamms workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            RUN_CLAMMS(ch_bams, file(params.fai))
            break

        case['dragen']:
            // Support both CRAM/CRAI and BAM/BAI file pairs while ignoring orphan index files.
            Channel.fromPath(params.cramFilePairsUploadPath, checkIfExists: true)
                .filter { f -> f.name.endsWith('.cram') || f.name.endsWith('.bam') }
                .map { primary ->
                    def sampleId = primary.name.replaceAll(/\\.(cram|bam)\$/, '').replaceAll(/\\.md\$/, '')
                    def idxCandidates = primary.name.endsWith('.cram')
                        ? [file("${primary}.crai"), file(primary.toString().replaceAll(/\\.cram\$/, '.crai'))]
                        : [file("${primary}.bai"), file("${primary}.bam.bai"), file(primary.toString().replaceAll(/\\.bam\$/, '.bai'))]
                    def indexFile = idxCandidates.find { it.exists() }
                    [sampleId, [primary, indexFile]]
                }
                .filter { sampleId, pair -> pair[1] != null }
                .ifEmpty { error "dragen workflow: no CRAM or BAM file pairs found in '${params.cramFilePairsUploadPath}'. Check that the path contains indexed CRAM/BAM files." }
                .set { ch_cramPairs }
            RUN_DRAGEN(ch_cramPairs)
            break

        case['cnvkit']:
            def glob_pattern = params.bams.contains('{') ? params.bams : params.bams.replaceAll(/\\.bam\$/, '.{bam,bai,bam.bai}')
            Channel.fromPath(glob_pattern)
                .map { it -> 
                    def sample_id = it.name.replaceAll(/\\.bam\\.bai\$/, '').replaceAll(/\\.bam\$/, '').replaceAll(/\\.bai\$/, '').replaceAll(/\\.md\$/, '')
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
                    def index_fallback = [file("${it}.bai"), file("${it}.bam.bai"), file(it.toString().replaceAll(/\\.bam\$/, '.bai'))]
                        .find { f -> f.exists() }
                    def index_file = index_fallback ? index_fallback : (it.name.endsWith('.bam') ? file("${it}.bai") : file("${it}.crai"))
                    return [ it.baseName.replaceAll(/\\.md\$/, ''), it, index_file ] 
                }
                .ifEmpty { error "gcnv workflow: no sample files found matching '${params.samples_path}'. Check that the glob matches existing BAM or CRAM files." }
                .set { ch_bams }
            RUN_GCNV(ch_bams, file(params.fasta), file(params.fai), file(params.dict), file(params.exome_targets))
            break

        case['exomedepth']:
            def exomedepth_probes_path = params.get('probes', params.get('targets', params.get('exome_targets', params.get('probes_bed', ''))))
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> "${row.BAM}\n" }
                .collectFile(name: 'exomedepth_bam_list.txt')
                .ifEmpty { error "exomedepth workflow: no samples found in samplesheet '${params.samplesheet_bams}'. Check that the file is tab-separated with SampleID and BAM columns." }
                .set { ch_bams }
            EXOMEDEPTH(ch_bams, file(exomedepth_probes_path), file(params.fai))
            break

        case['survivor']:
            def ch_vcfs = gather_vcfs()
            RUN_SURVIVOR(ch_vcfs)
            
            def qc_str = params.get('run_qc', false) ? "QC" : "NoQC"
            def is_mosdepth = params.get('probes_bed', '').contains('mosdepth') || params.get('standard_bed', '').contains('mosdepth')
            def bed_str = is_mosdepth ? "Mosdepth" : "Control"
            def run_name = "${bed_str}_${qc_str}_SURVIVOR"

            if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
                def truth_f = file(params.truth_bed)
                def probes_f = file(params.probes_bed)
                def filtered_comp = GATHER_AND_FILTER_VCFS.out.labeled_vcfs
                    .map { label, f -> ["${label}_${bed_str}_${qc_str}_SURVIVOR", f] }
                def merged_comp = RUN_SURVIVOR.out.union_vcf
                    .map { sample_id, f -> ["merged_intrasample_${bed_str}_${qc_str}_SURVIVOR", f] }
                EVALUATE_COMPARATIVE(filtered_comp.mix(merged_comp), truth_f, probes_f)
            }

            if (params.get('annotsv_db', false)) {
                annotsv_db_ch = Channel.fromPath(params.annotsv_db)
            } else {
                DOWNLOAD_ANNOTSV_DB()
                annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
            }
            
            def probes_bed_file = params.get('probes_bed', '/dataG/ddd/data/resources/canoes/probes_sanger.bed')
            RUN_CONSENSUS_ANNOTATION(run_name, RUN_SURVIVOR.out.union_vcf.map { sample_id, vcf -> vcf }, annotsv_db_ch, file(probes_bed_file))
            break

        case['truvari']:
            def ch_vcfs = gather_vcfs()
            RUN_TRUVARI(ch_vcfs)
            
            def qc_str_tr = params.get('run_qc', false) ? "QC" : "NoQC"
            def is_mosdepth_tr = params.get('probes_bed', '').contains('mosdepth') || params.get('standard_bed', '').contains('mosdepth')
            def bed_str_tr = is_mosdepth_tr ? "Mosdepth" : "Control"
            def run_name_tr = "${bed_str_tr}_${qc_str_tr}_TRUVARI"

            if (params.get('truth_bed', false) && params.get('probes_bed', false)) {
                def truth_f = file(params.truth_bed)
                def probes_f = file(params.probes_bed)
                def filtered_comp = GATHER_AND_FILTER_VCFS.out.labeled_vcfs
                    .map { label, f -> ["${label}_${bed_str_tr}_${qc_str_tr}_TRUVARI", f] }
                def merged_comp = RUN_TRUVARI.out.merged_vcf
                    .map { sample_id, f -> ["merged_intrasample_${bed_str_tr}_${qc_str_tr}_TRUVARI", f] }
                EVALUATE_COMPARATIVE(filtered_comp.mix(merged_comp), truth_f, probes_f)
            }

            if (params.get('annotsv_db', false)) {
                annotsv_db_ch = Channel.fromPath(params.annotsv_db)
            } else {
                DOWNLOAD_ANNOTSV_DB()
                annotsv_db_ch = DOWNLOAD_ANNOTSV_DB.out.db_dir
            }
            
            def probes_bed_file = params.get('probes_bed', '/dataG/ddd/data/resources/canoes/probes_sanger.bed')
            RUN_CONSENSUS_ANNOTATION(run_name_tr, RUN_TRUVARI.out.merged_vcf.map { sample_id, vcf -> vcf }, annotsv_db_ch, file(probes_bed_file))
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
            //           --caller   (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f ->
                    def sample_name = f.name
                        .replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN).*/, '')
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
            //   --gcnv_norm_dir, --dragen_norm_dir
            //   --bam_file, --reference_fasta, --bed_file, --mappability_file,
            //   -- --merger_mode (defaults to survivor)
            if (params.get('merger_mode', 'survivor') == 'truvari') {
                Channel.fromPath(params.merged_vcf_dir + '/**/*{_truvari_merged.vcf*,.intra-merge.cnv.vcf*}')
                    .map { f ->
                        def sample_id = extract_sample_id_from_vcf(f)
                        def collapsed_f = file("${f.parent}/${sample_id}_truvari_collapsed.vcf")
                        def collapsed_vcf_f = collapsed_f.exists() ? collapsed_f : []
                        build_feature_inputs(sample_id, f, collapsed_vcf_f)
                    }
                    .set { feature_inputs_ch }
                FEATURE_EXTRACTION(feature_inputs_ch)
            } else {
                Channel.fromPath(params.merged_vcf_dir + '/**/*{_survivor_union.vcf*,.intra-merge.cnv.vcf*}')
                    .map { f ->
                        def sample_id = extract_sample_id_from_vcf(f)
                        build_feature_inputs(sample_id, f, [])
                    }
                    .set { feature_inputs_ch }
                FEATURE_EXTRACTION(feature_inputs_ch)
            }
            break

        case['train']:
            // Train the ML model (XGBoost/CatBoost) using features extraction outputs.
            // Required: --features_dir   (directory containing *_features.tsv files
            //                             for positive and negative instances)
            //           --truth_labels   (csv containing EGAN_ID, CNV_ID, and True/False label)
            // Optional: --probes_bed     (bed containing probes used to filter target regions)
            Channel.value(file(params.truth_labels))
                .set { ch_truth }
            def probesBedParam = params.get('probes_bed', false)
            Channel.value(probesBedParam ? file(probesBedParam) : [])
                .set { ch_probes }
            RUN_TRAIN(Channel.fromPath(params.features_dir + '/*_features.tsv').collect(), ch_truth, ch_probes)
            break

        case['evaluate']:
            // Evaluate an existing cohort-wide CNV callset VCF directory against the gold truthset
            // Required: --vcf_dir    (directory containing sample VCFs to evaluate)
            //           --caller     (CANOES | CLAMMS | XHMM | GATK | CNVKIT | DRAGEN)
            //           --truth_bed  (gold standard truth set BED file)
            //           --probes_bed (target Sanger probes BED file)
            Channel.fromPath(params.vcf_dir + '/*.vcf*')
                .map { f -> [extract_sample_id_from_vcf(f), f] }
                .ifEmpty { error "evaluate workflow: no VCF files found in '${params.vcf_dir}'. Check that the directory contains .vcf or .vcf.gz files." }
                .set { ch_vcfs }
            RUN_EVALUATE(ch_vcfs, file(params.truth_bed), file(params.probes_bed), params.caller.toUpperCase())
            break

        case['full']:
        case['local']:
        case['full_with_dragen']:
            RUN_BENCHMARK_PIPELINE()
            break

        default:
            exit 1, """
OOOPS!! SEEMS LIKE WE HAVE A WORKFLOW ERROR!

No workflow 'mode' given!
Please use one of the following options for workflows:
       
    --workflow canoes       
    --workflow xhmm       
    --workflow clamms      
    --workflow dragen      
    --workflow cnvkit      
    --workflow gcnv       
    --workflow exomedepth
    --workflow survivor     
    --workflow truvari
    --workflow survivor_with_features
    --workflow truvari_with_features
    --workflow feature_extraction
    --workflow normalise
    --workflow train
    --workflow evaluate
    --workflow full
    --workflow local
    --workflow full_with_dragen
    --workflow qc
    --workflow cram_to_bam
    --workflow bam_to_cram
"""
            break
    }
}


// =====================================================================================
// POST-CONSENSUS COHORT MERGING, ANNOTATION & EVALUATION
// =====================================================================================






















