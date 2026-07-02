#!/usr/bin/env python3
"""
Tests for the evaluate module and supporting scripts.

Validates:
  1. bin/vcf_to_bed.py correctly converts a single-sample VCF to a 5-column BED
     file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID).
  2. bin/evaluate_caller_performance.py accepts the new 5-column BED format
     (CHR, START, STOP, CNV_TYPE, SAMPLE_ID) for both truth and call sets and
     correctly computes precision / sensitivity metrics.
  3. modules/evaluate/modules-evaluate.nf declares the VCF_TO_BED, COMBINE_BEDS, and
     EVALUATE_CALLER processes together with an EVALUATE workflow that chains them.
  4. main.nf includes the EVALUATE workflow, defines the RUN_EVALUATE sub-workflow,
     integrates optional evaluation after each of the 7 callers, and exposes the
     standalone 'evaluate' workflow mode.
"""

import io
import os
import re
import sys

import pytest

REPO_ROOT = os.path.join(os.path.dirname(__file__), '..')
BIN_DIR   = os.path.join(REPO_ROOT, 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))

NF_EVALUATE = os.path.join(REPO_ROOT, 'modules', 'evaluate', 'modules-evaluate.nf')
NF_MAIN     = os.path.join(REPO_ROOT, 'main.nf')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read(path):
    with open(path) as fh:
        return fh.read()


def _extract_process(text, name):
    match = re.search(
        rf'process {re.escape(name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)',
        text, re.DOTALL)
    return match.group(1) if match else None


def _extract_workflow(text, name):
    match = re.search(
        rf'workflow {re.escape(name)}\s*\{{(.+)',
        text, re.DOTALL)
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def evaluate_nf():
    return _read(NF_EVALUATE)


@pytest.fixture(scope='module')
def main_nf():
    return _read(NF_MAIN)


@pytest.fixture(scope='module')
def vcf_to_bed_text():
    path = os.path.join(BIN_DIR, 'vcf_to_bed.py')
    return _read(path)


@pytest.fixture(scope='module')
def eval_script_text():
    path = os.path.join(BIN_DIR, 'evaluate_caller_performance.py')
    return _read(path)


# ---------------------------------------------------------------------------
# vcf_to_bed.py – unit tests
# ---------------------------------------------------------------------------

class TestVcfToBed:
    """Functional tests for bin/vcf_to_bed.py."""

    @pytest.fixture()
    def del_vcf(self, tmp_path):
        """Minimal VCF with one DEL record."""
        vcf = tmp_path / 'SAMPLE1_CANOES.sorted.vcf'
        vcf.write_text(
            "##fileformat=VCFv4.1\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "chr1\t1001\t.\tN\t<DEL>\t.\tPASS\tEND=2000;SVTYPE=CNV\tGT\t0/1\n"
        )
        return str(vcf)

    @pytest.fixture()
    def dup_vcf(self, tmp_path):
        """Minimal VCF with one DUP record."""
        vcf = tmp_path / 'SAMPLE2_CLAMMS.sorted.vcf'
        vcf.write_text(
            "##fileformat=VCFv4.1\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE2\n"
            "chr2\t5001\t.\tN\t<DUP>\t.\tPASS\tEND=8000;SVTYPE=CNV\tGT\t0/1\n"
        )
        return str(vcf)

    def test_del_bed_columns(self, del_vcf, tmp_path):
        """BED file has 5 tab-separated columns for a DEL record."""
        import vcf_to_bed as mod
        out = str(tmp_path / 'out.bed')
        mod.vcf_to_bed(del_vcf, out)
        with open(out) as fh:
            lines = [l.rstrip('\n') for l in fh if l.strip()]
        assert len(lines) == 1
        cols = lines[0].split('\t')
        assert len(cols) == 5

    def test_del_bed_values(self, del_vcf, tmp_path):
        """CHR, START (0-based), STOP, CNV_TYPE, SAMPLE_ID are correct for DEL."""
        import vcf_to_bed as mod
        out = str(tmp_path / 'out.bed')
        mod.vcf_to_bed(del_vcf, out)
        with open(out) as fh:
            cols = fh.read().strip().split('\t')
        assert cols[0] == 'chr1'
        assert cols[1] == '1000'   # POS 1001 → 0-based 1000
        assert cols[2] == '2000'   # END from INFO
        assert cols[3] == 'DEL'
        assert cols[4] == 'SAMPLE1'

    def test_dup_cnv_type(self, dup_vcf, tmp_path):
        """CNV_TYPE is DUP for a DUP ALT record."""
        import vcf_to_bed as mod
        out = str(tmp_path / 'out.bed')
        mod.vcf_to_bed(dup_vcf, out)
        with open(out) as fh:
            cols = fh.read().strip().split('\t')
        assert cols[3] == 'DUP'
        assert cols[4] == 'SAMPLE2'

    def test_sample_id_override(self, del_vcf, tmp_path):
        """--sample_id argument overrides the sample name from the VCF header."""
        import vcf_to_bed as mod
        out = str(tmp_path / 'out.bed')
        mod.vcf_to_bed(del_vcf, out, sample_id_override='MYSAMPLE')
        with open(out) as fh:
            cols = fh.read().strip().split('\t')
        assert cols[4] == 'MYSAMPLE'

    def test_script_has_shebang(self, vcf_to_bed_text):
        assert vcf_to_bed_text.startswith('#!/usr/bin/env python3')

    def test_script_uses_argparse(self, vcf_to_bed_text):
        assert 'argparse' in vcf_to_bed_text

    def test_script_has_vcf_arg(self, vcf_to_bed_text):
        assert '--vcf' in vcf_to_bed_text

    def test_script_has_output_arg(self, vcf_to_bed_text):
        assert '--output' in vcf_to_bed_text


# ---------------------------------------------------------------------------
# evaluate_caller_performance.py – unit tests
# ---------------------------------------------------------------------------

class TestEvaluateCallerPerformance:
    """Functional tests for bin/evaluate_caller_performance.py."""

    @pytest.fixture()
    def probes_bed(self, tmp_path):
        """Three probe regions."""
        bed = tmp_path / 'probes.bed'
        bed.write_text(
            "chr1\t900\t1100\n"   # overlaps truth CNV
            "chr1\t1500\t1700\n"  # overlaps truth CNV
            "chr2\t0\t500\n"      # no CNV
        )
        return str(bed)

    @pytest.fixture()
    def truth_bed(self, tmp_path):
        """One DEL covering probes 1 & 2 for SAMPLE1."""
        bed = tmp_path / 'truth.bed'
        bed.write_text("chr1\t800\t2000\tDEL\tSAMPLE1\n")
        return str(bed)

    @pytest.fixture()
    def callset_bed_full(self, tmp_path):
        """Call set covering both truth probes (TP=2, FP=0, FN=0, TN=1)."""
        bed = tmp_path / 'calls_full.bed'
        bed.write_text("chr1\t800\t2000\tDEL\tSAMPLE1\n")
        return str(bed)

    @pytest.fixture()
    def callset_bed_empty(self, tmp_path):
        """Empty call set (TP=0, FN=2, FP=0, TN=1)."""
        bed = tmp_path / 'calls_empty.bed'
        bed.write_text("")
        return str(bed)

    def test_perfect_precision_and_sensitivity(
            self, truth_bed, callset_bed_full, probes_bed):
        """Sensitivity and precision are 1.0 when all truth probes are called."""
        import evaluate_caller_performance as mod
        probes   = mod.load_bed_file(probes_bed)
        truth    = mod.load_cnv_file(truth_bed)
        callset  = mod.load_cnv_file(callset_bed_full)
        cats     = mod.categorize_probes(probes, truth, callset)
        sens, prec, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert sens == 1.0
        assert prec == 1.0
        assert tn == 1

    def test_zero_sensitivity_empty_callset(
            self, truth_bed, callset_bed_empty, probes_bed):
        """Sensitivity is 0 when the call set is empty."""
        import evaluate_caller_performance as mod
        probes   = mod.load_bed_file(probes_bed)
        truth    = mod.load_cnv_file(truth_bed)
        callset  = mod.load_cnv_file(callset_bed_empty)
        cats     = mod.categorize_probes(probes, truth, callset)
        sens, prec, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert sens == 0.0

    def test_load_cnv_file_has_cnv_type_column(self, truth_bed):
        """load_cnv_file returns a DataFrame with a 'cnv_type' column."""
        import evaluate_caller_performance as mod
        df = mod.load_cnv_file(truth_bed)
        assert 'cnv_type' in df.columns
        assert 'sample' in df.columns

    def test_script_has_shebang(self, eval_script_text):
        assert eval_script_text.startswith('#!/usr/bin/env python3')

    def test_script_uses_argparse(self, eval_script_text):
        assert 'argparse' in eval_script_text

    def test_script_has_truth_bed_arg(self, eval_script_text):
        assert '--truth_bed' in eval_script_text

    def test_script_has_callset_bed_arg(self, eval_script_text):
        assert '--callset_bed' in eval_script_text

    def test_script_has_probes_bed_arg(self, eval_script_text):
        assert '--probes_bed' in eval_script_text

    def test_script_has_output_arg(self, eval_script_text):
        assert '--output' in eval_script_text

    def test_main_writes_output_file(
            self, truth_bed, callset_bed_full, probes_bed, tmp_path):
        """main() writes all metrics including new ones to the output file."""
        import evaluate_caller_performance as mod
        out = str(tmp_path / 'metrics.txt')
        sys.argv = ['evaluate_caller_performance.py',
                    '--truth_bed',   truth_bed,
                    '--callset_bed', callset_bed_full,
                    '--probes_bed',  probes_bed,
                    '--output',      out]
        mod.main()
        content = open(out).read()
        assert 'Sensitivity' in content
        assert 'Precision'   in content
        assert 'True Negatives' in content
        assert 'Specificity' in content
        assert 'F_beta' in content
        assert 'Matthews Correlation Coefficient' in content
        assert 'Confusion Matrix' in content

    def test_perfect_call_specificity(
            self, truth_bed, callset_bed_full, probes_bed):
        """Specificity is 1.0 when no FPs exist (the one non-truth probe is TN)."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_full)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert spec == 1.0

    def test_perfect_call_mcc(
            self, truth_bed, callset_bed_full, probes_bed):
        """MCC is 1.0 for a perfect classifier."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_full)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert mcc == 1.0

    def test_perfect_call_f_beta(
            self, truth_bed, callset_bed_full, probes_bed):
        """F_beta is 1.0 for a perfect classifier."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_full)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert f_beta == 1.0

    def test_empty_callset_specificity(
            self, truth_bed, callset_bed_empty, probes_bed):
        """Specificity is 1.0 when call set is empty (no FPs)."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_empty)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert spec == 1.0

    def test_empty_callset_f_beta_zero(
            self, truth_bed, callset_bed_empty, probes_bed):
        """F_beta is 0.0 when sensitivity and precision are both 0."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_empty)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert f_beta == 0.0

    def test_confusion_matrix_counts(
            self, truth_bed, callset_bed_full, probes_bed):
        """TP=2, FP=0, FN=0, TN=1 for the full-call fixture."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_full)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta, mcc = mod.compute_metrics(cats)
        assert tp == 2
        assert fp == 0
        assert fn == 0
        assert tn == 1

    def test_beta_parameter(
            self, truth_bed, callset_bed_full, probes_bed):
        """compute_metrics accepts a custom beta value."""
        import evaluate_caller_performance as mod
        probes  = mod.load_bed_file(probes_bed)
        truth   = mod.load_cnv_file(truth_bed)
        callset = mod.load_cnv_file(callset_bed_full)
        cats    = mod.categorize_probes(probes, truth, callset)
        _, _, tp, tn, fp, fn, spec, f_beta1, mcc = mod.compute_metrics(cats, beta=1)
        _, _, tp, tn, fp, fn, spec, f_beta2, mcc = mod.compute_metrics(cats, beta=2)
        assert f_beta1 == 1.0
        assert f_beta2 == 1.0


# ---------------------------------------------------------------------------
# modules/evaluate/modules-evaluate.nf – structural tests
# ---------------------------------------------------------------------------

class TestEvaluateModule:
    """Verify the Nextflow module structure."""

    def test_module_has_vcf_to_bed_process(self, evaluate_nf):
        assert 'process VCF_TO_BED' in evaluate_nf

    def test_module_has_combine_beds_process(self, evaluate_nf):
        assert 'process COMBINE_BEDS' in evaluate_nf

    def test_module_has_evaluate_caller_process(self, evaluate_nf):
        assert 'process EVALUATE_CALLER' in evaluate_nf

    def test_module_has_evaluate_workflow(self, evaluate_nf):
        assert 'workflow EVALUATE' in evaluate_nf

    def test_vcf_to_bed_calls_script(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'VCF_TO_BED')
        assert body is not None
        assert 'vcf_to_bed.py' in body

    def test_combine_beds_cats_files(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'COMBINE_BEDS')
        assert body is not None
        assert 'cat' in body

    def test_evaluate_caller_calls_script(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'EVALUATE_CALLER')
        assert body is not None
        assert 'evaluate_caller_performance.py' in body

    def test_evaluate_caller_passes_truth_bed(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'EVALUATE_CALLER')
        assert '--truth_bed' in body

    def test_evaluate_caller_passes_callset_bed(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'EVALUATE_CALLER')
        assert '--callset_bed' in body

    def test_evaluate_caller_passes_probes_bed(self, evaluate_nf):
        body = _extract_process(evaluate_nf, 'EVALUATE_CALLER')
        assert '--probes_bed' in body

    def test_evaluate_workflow_calls_vcf_to_bed(self, evaluate_nf):
        body = _extract_workflow(evaluate_nf, 'EVALUATE')
        assert body is not None
        assert 'VCF_TO_BED' in body

    def test_evaluate_workflow_calls_combine_beds(self, evaluate_nf):
        body = _extract_workflow(evaluate_nf, 'EVALUATE')
        assert 'COMBINE_BEDS' in body

    def test_evaluate_workflow_calls_evaluate_caller(self, evaluate_nf):
        body = _extract_workflow(evaluate_nf, 'EVALUATE')
        assert 'EVALUATE_CALLER' in body

    def test_evaluate_workflow_collects_beds(self, evaluate_nf):
        """COMBINE_BEDS must receive all per-sample BEDs via .collect()."""
        body = _extract_workflow(evaluate_nf, 'EVALUATE')
        assert 'collect()' in body

    def test_evaluate_workflow_emits_report(self, evaluate_nf):
        body = _extract_workflow(evaluate_nf, 'EVALUATE')
        assert 'performance_report' in body

    def test_processes_use_pysam_label(self, evaluate_nf):
        """All three evaluate processes must use label 'pysam'."""
        for proc in ('VCF_TO_BED', 'COMBINE_BEDS', 'EVALUATE_CALLER'):
            body = _extract_process(evaluate_nf, proc)
            assert body is not None, f"{proc} not found"
            assert "label 'pysam'" in body, f"{proc} missing label 'pysam'"


# ---------------------------------------------------------------------------
# main.nf – integration tests
# ---------------------------------------------------------------------------

class TestMainNfIntegration:
    """Verify main.nf includes EVALUATE and integrates it correctly."""

    def test_includes_evaluate_workflow(self, main_nf):
        assert "include { EVALUATE } from './modules/evaluate/modules-evaluate.nf'" in main_nf

    def test_defines_run_evaluate_sub_workflow(self, main_nf):
        assert 'workflow RUN_EVALUATE' in main_nf

    def test_has_evaluate_case(self, main_nf):
        assert "case['evaluate']" in main_nf

    def test_evaluate_case_uses_vcf_dir(self, main_nf):
        assert 'params.vcf_dir' in main_nf

    def test_evaluate_case_uses_truth_bed(self, main_nf):
        assert 'params.truth_bed' in main_nf

    def test_evaluate_case_uses_probes_bed(self, main_nf):
        assert 'params.probes_bed' in main_nf

    def test_evaluate_in_default_error_message(self, main_nf):
        assert '--workflow evaluate' in main_nf

    @pytest.mark.parametrize('caller', [
        'CANOES', 'XHMM', 'CLAMMS', 'CNVKIT', 'GCNV', 'DRAGEN', 'INDELIBLE'
    ])
    def test_evaluate_integrated_after_caller(self, caller, main_nf):
        """Each of the 7 callers must call EVALUATE conditionally."""
        # Find the RUN_<CALLER> sub-workflow block
        match = re.search(
            rf'workflow RUN_{caller}\s*\{{(.+?)(?=\nworkflow |\Z)',
            main_nf, re.DOTALL)
        assert match is not None, f'RUN_{caller} workflow not found'
        body = match.group(1)
        assert 'EVALUATE' in body, \
            f'RUN_{caller} does not call EVALUATE'
        assert 'truth_bed' in body, \
            f'RUN_{caller} does not check for truth_bed param'
