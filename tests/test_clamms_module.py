#!/usr/bin/env python3
"""
Tests for the CLAMMS Nextflow module (modules/callers/modules-clamms.nf).

Validates that the module correctly implements the CLAMMS pipeline for exome
sequencing CNV calling from a cohort (no matched normals), producing per-sample
VCF files. Based on the reference implementation in phelelani/nf-exomecnv and
the CLAMMS tool documentation (rgcgithub/clamms).

Checks:
  1. GENERATE_WINDOWS: probes sorted, INSERT_SIZE exported, correct argument
     order for annotate_windows.sh via $CLAMMS_DIR.
  2. SAMTOOLS_DOC: -Q 30 quality filter and awk mean-depth computation present.
  3. NORMALIZE_DOC: $CLAMMS_DIR prefix, sed chr-strip, output is .norm.cov.bed.
  4. CREATE_PCA_DATA: calls custom_ref_panel.sh (not a placeholder), outputs
     pca.coordinates.txt.
  5. CREATE_CUSTOM_REF_PANEL: calls custom_ref_panel.R with correct argument
     order, label is 'R', output glob is *.ref.panel.files.txt.
  6. COMBINE_PICARD_QC_METRICS: label is 'picard', no erroneous mv command,
     output named qcs_metrics.
  7. TRAIN_MODELS: chr-strip step present, correct $CLAMMS_DIR/fit_models
     argument order (ref_panel first, then windows).
  8. CALL_CNVS: correct $CLAMMS_DIR/call_cnv argument order (norm_cov first,
     then models), --sex flag present, output is .cnv.bed.
  9. FILTER_CLAMMS_CNVS: awk uses $10>0 (not $10>=0), output is a tuple of
     both samples.cnv.bed and samples.cnv.filtered.bed.
 10. Workflow: ref_panel extension .ref.panel.files.txt, filtered_cnvs mapped
     to extract the filtered BED before CONVERT_CLAMMS_TO_VCF.
"""

import os
import re

import pytest


NF_MODULE = "modules/callers/modules-clamms.nf"
REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_module(repo_root="."):
    path = os.path.join(repo_root, NF_MODULE)
    with open(path) as fh:
        return fh.read()


def _extract_process(module_text, process_name):
    """Return the body of a named process block, or None if not found."""
    match = re.search(
        rf"process {re.escape(process_name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def _extract_script(process_body):
    """Return the heredoc script block content from a process body."""
    match = re.search(r'"""\s*(.+?)\s*"""', process_body, re.DOTALL)
    return match.group(1) if match else ""


def _extract_doc_section(text, heading):
    """Return the text of a named section up to (but not including) the next heading."""
    section_start = text.find(heading)
    if section_start == -1:
        return None
    heading_match = re.search(r'\n#{1,6}(\s|\n|$)', text[section_start + 1:])
    section_end = (
        section_start + 1 + heading_match.start()
        if heading_match else len(text)
    )
    return text[section_start:section_end]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def module_text():
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    return _read_module(repo_root)


@pytest.fixture(scope="module")
def generate_windows_body(module_text):
    body = _extract_process(module_text, "GENERATE_WINDOWS")
    assert body is not None, "GENERATE_WINDOWS process not found"
    return body


@pytest.fixture(scope="module")
def samtools_doc_body(module_text):
    body = _extract_process(module_text, "SAMTOOLS_DOC")
    assert body is not None, "SAMTOOLS_DOC process not found"
    return body


@pytest.fixture(scope="module")
def normalize_doc_body(module_text):
    body = _extract_process(module_text, "NORMALIZE_DOC")
    assert body is not None, "NORMALIZE_DOC process not found"
    return body


@pytest.fixture(scope="module")
def create_pca_data_body(module_text):
    body = _extract_process(module_text, "CREATE_PCA_DATA")
    assert body is not None, "CREATE_PCA_DATA process not found"
    return body


@pytest.fixture(scope="module")
def create_custom_ref_panel_body(module_text):
    body = _extract_process(module_text, "CREATE_CUSTOM_REF_PANEL")
    assert body is not None, "CREATE_CUSTOM_REF_PANEL process not found"
    return body


@pytest.fixture(scope="module")
def combine_picard_qc_body(module_text):
    body = _extract_process(module_text, "COMBINE_PICARD_QC_METRICS")
    assert body is not None, "COMBINE_PICARD_QC_METRICS process not found"
    return body


@pytest.fixture(scope="module")
def train_models_body(module_text):
    body = _extract_process(module_text, "TRAIN_MODELS")
    assert body is not None, "TRAIN_MODELS process not found"
    return body


@pytest.fixture(scope="module")
def call_cnvs_body(module_text):
    body = _extract_process(module_text, "CALL_CNVS")
    assert body is not None, "CALL_CNVS process not found"
    return body


@pytest.fixture(scope="module")
def filter_clamms_cnvs_body(module_text):
    body = _extract_process(module_text, "FILTER_CLAMMS_CNVS")
    assert body is not None, "FILTER_CLAMMS_CNVS process not found"
    return body


# ===========================================================================
# 1. GENERATE_WINDOWS
# ===========================================================================

class TestGenerateWindows:
    """GENERATE_WINDOWS must sort probes, set INSERT_SIZE, and call
    annotate_windows.sh via $CLAMMS_DIR with the correct argument order."""

    def test_probes_sorted_before_annotation(self, generate_windows_body):
        """Probes BED must be sorted with sort -k1,1 -k2,2n before annotating."""
        script = _extract_script(generate_windows_body)
        assert "sort -k1,1 -k2,2n" in script, (
            "Probes must be sorted with 'sort -k1,1 -k2,2n' before calling "
            "annotate_windows.sh (CLAMMS requirement)."
        )

    def test_mappability_sorted_before_annotation(self, generate_windows_body):
        """Mappability BED must be sorted with sort -k1,1 -k2,2n before annotating.
        CLAMMS annotate_windows.sh uses bedtools with sorted input; passing an
        unsorted mappability file causes incorrect window annotations or errors."""
        script = _extract_script(generate_windows_body)
        # There must be a sort step that produces mappability_sorted.bed
        assert "mappability_sorted.bed" in script, (
            "Mappability BED must be sorted to 'mappability_sorted.bed' before "
            "being passed to annotate_windows.sh (bedtools -sorted requirement)."
        )
        # The sorted file (not the raw param) must be passed to annotate_windows.sh
        annot_match = re.search(r"annotate_windows\.sh\s+\S+\s+\S+\s+(\S+)", script)
        assert annot_match, "annotate_windows.sh call not found in script"
        mappability_arg = annot_match.group(1)
        assert "mappability_sorted" in mappability_arg, (
            "annotate_windows.sh must receive mappability_sorted.bed (the "
            "pre-sorted copy), not the raw ${mappability} parameter."
        )

    def test_insert_size_exported(self, generate_windows_body):
        """INSERT_SIZE environment variable must be exported/set."""
        script = _extract_script(generate_windows_body)
        assert "INSERT_SIZE" in script, (
            "INSERT_SIZE must be set (e.g., export INSERT_SIZE=200) for "
            "annotate_windows.sh to compute correct GC content windows."
        )

    def test_uses_clamms_dir_prefix(self, generate_windows_body):
        """annotate_windows.sh must be called via $CLAMMS_DIR."""
        script = _extract_script(generate_windows_body)
        assert "CLAMMS_DIR" in script and "annotate_windows.sh" in script, (
            "annotate_windows.sh must be called with the $CLAMMS_DIR prefix."
        )

    def test_correct_argument_order(self, generate_windows_body):
        """annotate_windows.sh must receive targets first, then ref, mappability,
        INSERT_SIZE, and special_regions — NOT ref first."""
        script = _extract_script(generate_windows_body)
        # The sorted targets file should appear before the ref in the command
        annot_match = re.search(r"annotate_windows\.sh\s+(\S+)\s+(\S+)", script)
        assert annot_match, "annotate_windows.sh call not found in script"
        first_arg = annot_match.group(1)
        # First argument should be the sorted targets bed, not the ref genome
        assert "targets_sorted" in first_arg or ".bed" in first_arg, (
            "First argument to annotate_windows.sh must be the sorted targets BED "
            "(not the reference genome FASTA)."
        )

    def test_output_redirected_to_windows_bed(self, generate_windows_body):
        """Output must be redirected to windows.bed."""
        script = _extract_script(generate_windows_body)
        assert "> windows.bed" in script, (
            "annotate_windows.sh output must be redirected to windows.bed."
        )


# ===========================================================================
# 2. SAMTOOLS_DOC
# ===========================================================================

class TestSamtoolsDoc:
    """SAMTOOLS_DOC must filter by mapping quality and compute mean depth."""

    def test_quality_filter_q30(self, samtools_doc_body):
        """samtools bedcov must use -Q 30 to filter low-quality reads."""
        script = _extract_script(samtools_doc_body)
        assert "-Q 30" in script, (
            "samtools bedcov must include '-Q 30' to filter reads with "
            "mapping quality below 30 (CLAMMS documentation requirement)."
        )

    def test_awk_mean_depth_computation(self, samtools_doc_body):
        """An awk pipeline must compute mean depth as $NF/($3-$2)."""
        script = _extract_script(samtools_doc_body)
        assert "awk" in script, (
            "An awk command must follow samtools bedcov to compute mean depth "
            "as $NF/($3-$2) per window."
        )
        assert r"$NF" in script or r"\$NF" in script, (
            "awk must use $NF (last column, total coverage) to compute mean depth."
        )

    def test_output_format_four_columns(self, samtools_doc_body):
        """Output must be reformatted to 4 columns: chrom, start, end, mean_depth."""
        script = _extract_script(samtools_doc_body)
        assert "printf" in script or r'printf' in script, (
            "awk printf must be used to format the 4-column output "
            "(chrom, start, end, mean_depth)."
        )


# ===========================================================================
# 3. NORMALIZE_DOC
# ===========================================================================

class TestNormalizeDoc:
    """NORMALIZE_DOC must use $CLAMMS_DIR prefix, strip chr, and output .norm.cov.bed."""

    def test_uses_clamms_dir_prefix(self, normalize_doc_body):
        """normalize_coverage must be called via $CLAMMS_DIR."""
        script = _extract_script(normalize_doc_body)
        assert "CLAMMS_DIR" in script and "normalize_coverage" in script, (
            "normalize_coverage must be called with the $CLAMMS_DIR prefix."
        )

    def test_sed_strips_chr_prefix(self, normalize_doc_body):
        """A sed command must strip the 'chr' prefix from chromosome names."""
        script = _extract_script(normalize_doc_body)
        assert "sed" in script and "chr" in script, (
            "sed 's/^chr//g' (or equivalent) must be applied to strip the 'chr' "
            "prefix so chromosome names are consistent with CLAMMS expectations."
        )

    def test_output_extension_norm_cov_bed(self, normalize_doc_body):
        """Output file must use .norm.cov.bed extension (not .norm.coverage.bed)."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", normalize_doc_body, re.DOTALL
        )
        assert output_section, "output section not found in NORMALIZE_DOC"
        out_text = output_section.group(1)
        assert ".norm.cov.bed" in out_text, (
            "Normalized coverage output must use the .norm.cov.bed extension "
            "to match what custom_ref_panel.sh expects (ls *.norm.cov.bed)."
        )
        assert ".norm.coverage.bed" not in out_text, (
            "Output must not use .norm.coverage.bed; use .norm.cov.bed instead."
        )


# ===========================================================================
# 4. CREATE_PCA_DATA
# ===========================================================================

class TestCreatePcaData:
    """CREATE_PCA_DATA must call the real custom_ref_panel.sh script, not a placeholder."""

    def test_calls_custom_ref_panel_sh(self, create_pca_data_body):
        """Script must call custom_ref_panel.sh to generate PCA coordinates."""
        script = _extract_script(create_pca_data_body)
        assert "custom_ref_panel.sh" in script, (
            "CREATE_PCA_DATA must call custom_ref_panel.sh to compute PCA "
            "coordinates from normalized coverage data."
        )

    def test_no_placeholder_content(self, create_pca_data_body):
        """Script must not contain placeholder echo or dummy commands."""
        script = _extract_script(create_pca_data_body)
        # A real implementation should NOT have echo+placeholder lines
        has_placeholder = bool(re.search(r'echo\s+.*[Pp]laceholder', script))
        assert not has_placeholder, (
            "CREATE_PCA_DATA must not contain placeholder content; "
            "it must call the real custom_ref_panel.sh script."
        )

    def test_output_is_pca_coordinates_txt(self, create_pca_data_body):
        """Output must be pca.coordinates.txt (matching custom_ref_panel.sh output)."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", create_pca_data_body, re.DOTALL
        )
        assert output_section, "output section not found in CREATE_PCA_DATA"
        out_text = output_section.group(1)
        assert "pca.coordinates.txt" in out_text, (
            "CREATE_PCA_DATA output must be 'pca.coordinates.txt', which is what "
            "custom_ref_panel.sh produces."
        )


# ===========================================================================
# 5. CREATE_CUSTOM_REF_PANEL
# ===========================================================================

class TestCreateCustomRefPanel:
    """CREATE_CUSTOM_REF_PANEL must call the R script, have label 'R',
    and output *.ref.panel.files.txt."""

    def test_calls_custom_ref_panel_r(self, create_custom_ref_panel_body):
        """Script must call custom_ref_panel.R (not a shell script)."""
        script = _extract_script(create_custom_ref_panel_body)
        assert "custom_ref_panel.R" in script, (
            "CREATE_CUSTOM_REF_PANEL must call custom_ref_panel.R to create "
            "sample-specific reference panels using k-nearest neighbors."
        )

    def test_label_is_r(self, create_custom_ref_panel_body):
        """Process label must be 'R' (not 'clamms_bedtools')."""
        assert "label 'R'" in create_custom_ref_panel_body, (
            "CREATE_CUSTOM_REF_PANEL must use label 'R' since it runs an R script."
        )

    def test_output_extension_ref_panel_files_txt(self, create_custom_ref_panel_body):
        """Output ref_panel glob must use *.ref.panel.files.txt extension."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", create_custom_ref_panel_body, re.DOTALL
        )
        assert output_section, "output section not found in CREATE_CUSTOM_REF_PANEL"
        out_text = output_section.group(1)
        assert "*.ref.panel.files.txt" in out_text, (
            "CREATE_CUSTOM_REF_PANEL output glob must be '*.ref.panel.files.txt' "
            "to match what custom_ref_panel.R generates."
        )

    def test_r_script_receives_sexinfo(self, create_custom_ref_panel_body):
        """custom_ref_panel.R must receive the sexinfo file as an argument."""
        script = _extract_script(create_custom_ref_panel_body)
        assert "sexinfo" in script, (
            "custom_ref_panel.R must receive the sexinfo file so it can "
            "create sex-aware reference panels."
        )


# ===========================================================================
# 6. COMBINE_PICARD_QC_METRICS
# ===========================================================================

class TestCombinePicardQcMetrics:
    """COMBINE_PICARD_QC_METRICS must have label 'picard' and not rename the output."""

    def test_label_is_picard(self, combine_picard_qc_body):
        """Process label must be 'picard'."""
        assert "label 'picard'" in combine_picard_qc_body, (
            "COMBINE_PICARD_QC_METRICS must use label 'picard' since it runs "
            "Picard-related QC combination using combine_picard_qc_metrics.sh."
        )

    def test_no_mv_command(self, combine_picard_qc_body):
        """Script must not rename qcs_metrics with mv."""
        script = _extract_script(combine_picard_qc_body)
        assert "mv qcs_metrics" not in script, (
            "COMBINE_PICARD_QC_METRICS must not rename 'qcs_metrics' with mv; "
            "combine_picard_qc_metrics.sh already produces the file as 'qcs_metrics'."
        )

    def test_output_named_qcs_metrics(self, combine_picard_qc_body):
        """Output file must be named qcs_metrics."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", combine_picard_qc_body, re.DOTALL
        )
        assert output_section, "output section not found in COMBINE_PICARD_QC_METRICS"
        out_text = output_section.group(1)
        assert "qcs_metrics" in out_text, (
            "COMBINE_PICARD_QC_METRICS output must be named 'qcs_metrics' to "
            "match the filename produced by combine_picard_qc_metrics.sh."
        )


# ===========================================================================
# 7. TRAIN_MODELS
# ===========================================================================

class TestTrainModels:
    """TRAIN_MODELS must strip chr prefix from windows and pass args in correct order."""

    def test_chr_strip_step(self, train_models_body):
        """A sed command must strip the 'chr' prefix from windows before fitting."""
        script = _extract_script(train_models_body)
        assert "sed" in script and "chr" in script, (
            "TRAIN_MODELS must strip the 'chr' prefix from windows.bed with "
            "sed 's/chr//g' since CLAMMS uses non-prefixed chromosome names."
        )

    def test_uses_clamms_dir_prefix(self, train_models_body):
        """fit_models must be called via $CLAMMS_DIR."""
        script = _extract_script(train_models_body)
        assert "CLAMMS_DIR" in script and "fit_models" in script, (
            "fit_models must be called with the $CLAMMS_DIR prefix."
        )

    def test_ref_panel_before_windows_in_fit_models(self, train_models_body):
        """fit_models argument order: ref_panel file first, then windows.bed."""
        script = _extract_script(train_models_body)
        fit_match = re.search(r"fit_models\s+(\S+)\s+(\S+)", script)
        assert fit_match, "fit_models command not found in TRAIN_MODELS script"
        first_arg = fit_match.group(1)
        second_arg = fit_match.group(2)
        # The ref panel file should be first (contains .txt or ref_panel_file)
        # The windows file should be second (ends with .bed or contains windows)
        assert "windows" in second_arg or ".bed" in second_arg, (
            "fit_models second argument must be the windows BED file "
            "(correct order: fit_models <ref_panel> <windows.bed>)."
        )


# ===========================================================================
# 8. CALL_CNVS
# ===========================================================================

class TestCallCnvs:
    """CALL_CNVS must pass args in the correct order, include --sex, and
    output .cnv.bed files."""

    def test_uses_clamms_dir_prefix(self, call_cnvs_body):
        """call_cnv must be called via $CLAMMS_DIR."""
        script = _extract_script(call_cnvs_body)
        assert "CLAMMS_DIR" in script and "call_cnv" in script, (
            "call_cnv must be called with the $CLAMMS_DIR prefix."
        )

    def test_sex_flag_present(self, call_cnvs_body):
        """call_cnv must include the --sex flag for sex-aware CNV calling."""
        script = _extract_script(call_cnvs_body)
        assert "--sex" in script, (
            "call_cnv must include '--sex $sex' for sex-aware CNV calling. "
            "Without this, sex chromosome CNV calls will be incorrect."
        )

    def test_sex_looked_up_from_sexinfo(self, call_cnvs_body):
        """Sex value must be looked up from the sexinfo file."""
        script = _extract_script(call_cnvs_body)
        assert "sexinfo" in script, (
            "CALL_CNVS must look up the sample sex from the sexinfo file "
            "(e.g., grep ${sample_id} ${sexinfo} | cut -f 2)."
        )

    def test_norm_cov_before_models_in_call_cnv(self, call_cnvs_body):
        """call_cnv argument order: norm_cov first, then models."""
        script = _extract_script(call_cnvs_body)
        # Find the call_cnv line
        call_match = re.search(r"call_cnv\s+(\S+)\s+(\S+)", script)
        assert call_match, "call_cnv command not found in CALL_CNVS script"
        first_arg = call_match.group(1)
        second_arg = call_match.group(2)
        # norm_cov should come before models
        assert "norm_cov" in first_arg or "cov" in first_arg, (
            "call_cnv first argument must be the normalized coverage file "
            "(correct order: call_cnv <norm_cov> <models> --sex <sex>)."
        )
        assert "model" in second_arg or "models" in second_arg, (
            "call_cnv second argument must be the models BED file."
        )

    def test_output_extension_cnv_bed(self, call_cnvs_body):
        """Output file must use .cnv.bed extension."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", call_cnvs_body, re.DOTALL
        )
        assert output_section, "output section not found in CALL_CNVS"
        out_text = output_section.group(1)
        assert ".cnv.bed" in out_text, (
            "CALL_CNVS output must use the .cnv.bed extension "
            "(matching CLAMMS convention and expected by FILTER_CLAMMS_CNVS)."
        )
        assert "_CLAMMS_calls.bed" not in out_text, (
            "Output must not use _CLAMMS_calls.bed; use .cnv.bed instead."
        )


# ===========================================================================
# 9. FILTER_CLAMMS_CNVS
# ===========================================================================

class TestFilterClammsCnvs:
    """FILTER_CLAMMS_CNVS must use the correct filter threshold and emit a tuple."""

    def test_awk_filter_q_exact_strictly_positive(self, filter_clamms_cnvs_body):
        """awk filter for Q_EXACT (column 10) must be > 0, not >= 0."""
        script = _extract_script(filter_clamms_cnvs_body)
        # Match both escaped (\$10) and unescaped ($10) forms
        assert re.search(r'\\?\$10\s*>\s*0\b', script), (
            "FILTER_CLAMMS_CNVS awk must filter with Q_EXACT > 0 (column 10 > 0), "
            "not >= 0. Calls with Q_EXACT = 0 are of questionable quality."
        )

    def test_no_awk_filter_q_exact_gte_zero(self, filter_clamms_cnvs_body):
        """awk must NOT use >= 0 for Q_EXACT (that accepts Q_EXACT=0 calls)."""
        script = _extract_script(filter_clamms_cnvs_body)
        assert not re.search(r'\\?\$10\s*>=\s*0', script), (
            "FILTER_CLAMMS_CNVS must not use '$10 >= 0'; use '$10 > 0' to "
            "exclude borderline-quality calls."
        )

    def test_output_is_tuple_of_both_beds(self, filter_clamms_cnvs_body):
        """Output must be a tuple containing both samples.cnv.bed and
        samples.cnv.filtered.bed."""
        output_section = re.search(
            r"output:\s*(.+?)(?=\s*script:)", filter_clamms_cnvs_body, re.DOTALL
        )
        assert output_section, "output section not found in FILTER_CLAMMS_CNVS"
        out_text = output_section.group(1)
        assert "samples.cnv.bed" in out_text, (
            "FILTER_CLAMMS_CNVS output must include samples.cnv.bed "
            "(all CNV calls before filtering)."
        )
        assert "samples.cnv.filtered.bed" in out_text, (
            "FILTER_CLAMMS_CNVS output must include samples.cnv.filtered.bed "
            "(quality-filtered CNV calls)."
        )
        assert "tuple" in out_text, (
            "FILTER_CLAMMS_CNVS output must be a tuple emitting both BED files."
        )

    def test_input_glob_uses_cnv_bed(self, filter_clamms_cnvs_body):
        """cat glob must match *.cnv.bed (not *_CLAMMS_calls.bed)."""
        script = _extract_script(filter_clamms_cnvs_body)
        assert "*.cnv.bed" in script, (
            "FILTER_CLAMMS_CNVS script must use '*.cnv.bed' glob to concatenate "
            "per-sample CNV calls (matching CALL_CNVS .cnv.bed output)."
        )


# ===========================================================================
# 10. Workflow orchestration
# ===========================================================================

class TestWorkflow:
    """CLAMMS workflow must use correct file extensions and extract the right
    element from the filtered_cnvs tuple for VCF conversion."""

    def test_ref_panel_uses_new_extension(self, module_text):
        """Workflow must parse .ref.panel.files.txt (not .ref.panel.files)."""
        wf_match = re.search(
            r"workflow CLAMMS\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "CLAMMS workflow not found"
        wf_body = wf_match.group(1)
        assert ".ref.panel.files.txt" in wf_body, (
            "CLAMMS workflow must use '.ref.panel.files.txt' extension when "
            "parsing ref panel filenames (matching CREATE_CUSTOM_REF_PANEL output)."
        )
        # The old extension (without .txt) must not appear in the replace call
        assert ".ref.panel.files'" not in wf_body, (
            "Old '.ref.panel.files' extension (without .txt) must not be used; "
            "update to '.ref.panel.files.txt'."
        )

    def test_filtered_cnvs_mapped_to_filtered_bed(self, module_text):
        """Workflow must extract the filtered BED (second element) from
        filtered_cnvs tuple before passing to CONVERT_CLAMMS_TO_VCF."""
        wf_match = re.search(
            r"workflow CLAMMS\s*\{(.+)",
            module_text,
            re.DOTALL,
        )
        assert wf_match, "CLAMMS workflow not found"
        wf_body = wf_match.group(1)
        # The workflow should map filtered_cnvs to get the second element
        assert re.search(
            r"filtered_cnvs\.map\s*\{",
            wf_body
        ), (
            "CLAMMS workflow must map filtered_cnvs to extract the filtered BED "
            "file (second tuple element) before passing to CONVERT_CLAMMS_TO_VCF."
        )


# ===========================================================================
# 11. README and docs sexinfo documentation
# ===========================================================================

class TestReadmeClammsSexinfoDocumentation:
    """README should document the required sexinfo file format for CLAMMS."""

    README_PATH = os.path.join(REPO_ROOT, 'README.md')

    @pytest.fixture(scope='class')
    def readme_text(self):
        with open(self.README_PATH) as fh:
            return fh.read()

    def test_readme_mentions_sexinfo_requirements(self, readme_text):
        assert 'Sexinfo file requirements' in readme_text, (
            "README.md must contain a 'Sexinfo file requirements' section."
        )

    def test_readme_sexinfo_documents_sample_id_column(self, readme_text):
        section = _extract_doc_section(readme_text, 'Sexinfo file requirements')
        assert section is not None
        assert 'sample_id' in section, (
            "Sexinfo requirements section must mention the 'sample_id' column."
        )

    def test_readme_sexinfo_documents_sex_values(self, readme_text):
        section = _extract_doc_section(readme_text, 'Sexinfo file requirements')
        assert section is not None
        assert re.search(r'\bM\b', section) and re.search(r'\bF\b', section), (
            "Sexinfo requirements section must document both 'M' (male) and "
            "'F' (female) as valid sex values."
        )

    def test_readme_sexinfo_documents_non_exact_id_match_requirement(self, readme_text):
        section = _extract_doc_section(readme_text, 'Sexinfo file requirements')
        assert section is not None
        assert 'does not need to contain exactly' in section and 'set of IDs' in section, (
            "Sexinfo requirements must clarify that sexinfo does not need to be an "
            "exact ID-set match to samplesheet.tsv."
        )
        assert 'extra IDs are ignored' in section and re.search(
            r'must\s+include\s+all sample IDs', section
        ), (
            "Sexinfo requirements must clarify that extra IDs are ignored but all "
            "samplesheet IDs are required."
        )


class TestClammsDocSexinfoDocumentation:
    """docs/modules/clamms.md should document the required sexinfo file format."""

    CLAMMS_DOC_PATH = os.path.join(REPO_ROOT, 'docs', 'modules', 'clamms.md')

    @pytest.fixture(scope='class')
    def doc_text(self):
        with open(self.CLAMMS_DOC_PATH) as fh:
            return fh.read()

    def test_clamms_doc_mentions_sexinfo_requirements(self, doc_text):
        assert 'Sexinfo file requirements' in doc_text, (
            "docs/modules/clamms.md must contain a 'Sexinfo file requirements' section."
        )

    def test_clamms_doc_sexinfo_documents_sample_id_column(self, doc_text):
        section = _extract_doc_section(doc_text, 'Sexinfo file requirements')
        assert section is not None
        assert 'sample_id' in section, (
            "clamms.md sexinfo section must mention the 'sample_id' column."
        )

    def test_clamms_doc_sexinfo_documents_sex_values(self, doc_text):
        section = _extract_doc_section(doc_text, 'Sexinfo file requirements')
        assert section is not None
        assert re.search(r'\bM\b', section) and re.search(r'\bF\b', section), (
            "clamms.md sexinfo section must document both 'M' (male) and "
            "'F' (female) as valid sex values."
        )

    def test_clamms_doc_sexinfo_documents_non_exact_id_match_requirement(self, doc_text):
        section = _extract_doc_section(doc_text, 'Sexinfo file requirements')
        assert section is not None
        assert 'does not need to contain exactly' in section and 'set of IDs' in section, (
            "clamms.md sexinfo section must clarify that sexinfo does not need to "
            "be an exact ID-set match to samplesheet.tsv."
        )
        assert 'extra IDs are ignored' in section and re.search(
            r'must\s+include\s+all sample IDs', section
        ), (
            "clamms.md sexinfo section must clarify that extra IDs are ignored but "
            "all samplesheet IDs are required."
        )
