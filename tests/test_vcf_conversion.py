#!/usr/bin/env python3
"""
Tests for bin/canoes_csv_to_vcf.py and bin/clamms_bed_to_vcf.py.

Verifies that each script:
  1. Runs without errors on representative input data.
  2. Produces a valid VCF file with the mandatory header lines.
  3. Includes TOOL=<TOOLNAME> in the INFO field of every data record.
  4. Outputs individual per-sample VCFs (one sample column each).
"""

import os
import sys
import pytest

# Allow importing the scripts directly
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))
import clamms_bed_to_vcf as clamms_converter


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def fai_file(tmp_path):
    """A minimal FAI file covering chromosomes used in tests."""
    fai = tmp_path / "reference.fa.fai"
    fai.write_text(
        "chr1\t248956422\t112\t70\t71\n"
        "chr2\t242193529\t252513167\t70\t71\n"
        "chr3\t198295559\t498166716\t70\t71\n"
    )
    return str(fai)


@pytest.fixture()
def sample_file(tmp_path):
    """A sample-list file with two sample IDs."""
    sf = tmp_path / "samples.txt"
    sf.write_text("SAMPLE1\nSAMPLE2\n")
    return str(sf)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _read_vcf_records(vcf_path):
    """Return (header_lines, data_lines) from a VCF file."""
    header, data = [], []
    with open(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
            else:
                data.append(line)
    return header, data


# ===========================================================================
# CANOES
# ===========================================================================

@pytest.fixture()
def canoes_csv(tmp_path):
    """Minimal CANOES TSV with one call per sample.

    Uses the actual tab-separated output format produced by run_canoes.R /
    CANOES.R (write.table with sep='\t'), with the column order that
    PrintCNVs / CalcCopyNumber generates:
      SAMPLE, CNV, INTERVAL, KB, CHR, MID_BP, TARGETS, NUM_TARG, MLCN, Q_SOME
    INTERVAL uses integer chromosome numbers (no 'chr' prefix) because
    RUN_CANOES strips the prefix with sed before invoking the R script.
    """
    csv = tmp_path / "canoes.csv"
    csv.write_text(
        "SAMPLE\tCNV\tINTERVAL\tKB\tCHR\tMID_BP\tTARGETS\tNUM_TARG\tMLCN\tQ_SOME\n"
        "SAMPLE1\tDEL\t1:1000-2000\t1.0\t1\t1500\t1..2\t2\t1\t90.5\n"
        "SAMPLE2\tDUP\t2:5000-8000\t3.0\t2\t6500\t3..5\t3\t3\t75.2\n"
    )
    return str(csv)


class TestCanoesCSVtoVCF:
    """Tests for bin/canoes_csv_to_vcf.py."""

    def test_individual_vcfs_created(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Individual per-sample VCFs are written to the output directory."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        assert os.path.isfile(os.path.join(out, "SAMPLE2_CANOES_output.vcf"))

    def test_only_one_sample_per_vcf(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Each VCF file contains exactly one sample column."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        for sample, fname in (
            ("SAMPLE1", "SAMPLE1_CANOES_output.vcf"),
            ("SAMPLE2", "SAMPLE2_CANOES_output.vcf"),
        ):
            header, _ = _read_vcf_records(os.path.join(out, fname))
            col_header = [l for l in header if l.startswith("#CHROM")][0]
            columns = col_header.split("\t")
            # FORMAT is column index 8; sample columns start at index 9
            sample_cols = columns[9:]
            assert sample_cols == [sample], (
                f"{fname}: expected sample column [{sample}], got {sample_cols}"
            )

    def test_tool_field_in_info(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Every data record in individual VCFs contains TOOL=CANOES in INFO."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CANOES_output.vcf", "SAMPLE2_CANOES_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                info = record.split("\t")[7]
                assert "TOOL=CANOES" in info, (
                    f"TOOL=CANOES missing from INFO in {fname}: {info}"
                )

    def test_vcf_header_contains_tool_meta(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The VCF header declares the TOOL INFO field."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        tool_meta = [l for l in header if "##INFO=<ID=TOOL" in l]
        assert tool_meta, "##INFO=<ID=TOOL,...> not found in VCF header"

    def test_fileformat_line_present(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The VCF output starts with ##fileformat=VCFv4.x."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        assert header[0].startswith("##fileformat=VCFv4"), (
            f"First header line is not ##fileformat: {header[0]}"
        )

    def test_source_line_is_canoes(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The ##source header line identifies CANOES as the caller."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        source_lines = [l for l in header if l.startswith("##source=")]
        assert source_lines, "No ##source= line in VCF header"
        assert "CANOES" in source_lines[0], (
            f"##source does not mention CANOES: {source_lines[0]}"
        )

    def test_contig_lines_present(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Contig lines derived from the FAI file appear in the header."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert contig_lines, "No ##contig= lines found in VCF header"

    def test_chrom_field_has_chr_prefix(self, canoes_csv, sample_file, fai_file, tmp_path):
        """CHROM column in every data record has the 'chr' prefix.

        CANOES.R outputs integer chromosome numbers (1, 2, …) because
        RUN_CANOES strips the prefix with sed before invoking R.
        The converter must re-add 'chr' so that the VCF is compliant with
        the GRCh3x reference naming convention.
        """
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CANOES_output.vcf", "SAMPLE2_CANOES_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                chrom = record.split("\t")[0]
                assert chrom.startswith("chr"), (
                    f"CHROM field does not start with 'chr' in {fname}: got '{chrom}'"
                )

    def test_malformed_interval_row_is_skipped(self, sample_file, fai_file, tmp_path):
        """Malformed INTERVAL values must be skipped instead of writing POS=0."""
        import canoes_csv_to_vcf as mod

        bad_csv = tmp_path / "canoes_bad_interval.csv"
        bad_csv.write_text(
            "SAMPLE\tCNV\tINTERVAL\tKB\tCHR\tMID_BP\tTARGETS\tNUM_TARG\tMLCN\tQ_SOME\n"
            "SAMPLE1\tDEL\tBAD_INTERVAL\t1.0\t1\t1500\t1..2\t2\t1\t90.5\n"
        )
        out = str(tmp_path / "out")
        mod.process_canoes_data(str(bad_csv), sample_file, fai_file, out)

        vcf = os.path.join(out, "SAMPLE1_CANOES_output.vcf")
        assert os.path.isfile(vcf)
        _, data = _read_vcf_records(vcf)
        assert data == [], "Malformed INTERVAL row should be skipped from VCF output"

    def test_malformed_fai_lines_are_ignored(self, canoes_csv, sample_file, tmp_path):
        """Malformed FAI lines should be ignored without crashing conversion."""
        import canoes_csv_to_vcf as mod

        bad_fai = tmp_path / "bad.fai"
        bad_fai.write_text(
            "chr1\t248956422\t112\t70\t71\n"
            "chr2\n"
            "chr3\tNOT_INT\t0\t70\t71\n"
        )
        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, str(bad_fai), out)

        vcf = os.path.join(out, "SAMPLE1_CANOES_output.vcf")
        assert os.path.isfile(vcf)
        header, _ = _read_vcf_records(vcf)
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert any("ID=chr1" in l for l in contig_lines)
        assert not any("ID=chr2" in l for l in contig_lines)
        assert not any("ID=chr3" in l for l in contig_lines)


# ===========================================================================
# CLAMMS
# ===========================================================================

@pytest.fixture()
def clamms_bed(tmp_path):
    """
    Minimal CLAMMS BED with 18 tab-separated fields per row (one call each
    for SAMPLE1 and SAMPLE2).
    """
    bed = tmp_path / "clamms.bed"
    rows = [
        # CHROM START  END    INTERVAL             SAMPLE   CNV  MLCN  NUM_WIN  Q_SOME  Q_EXACT  QL_EXT  L_EXT  QR_EXT  R_EXT  QL_CONT  L_CONT  QR_CONT  R_CONT
        "chr1\t1000\t2000\tchr1:1000-2000\tSAMPLE1\tDEL\t1\t5\t600\t50\t30\t800\t25\t1200\t15\t1100\t20\t1300",
        "chr2\t5000\t8000\tchr2:5000-8000\tSAMPLE2\tDUP\t3\t8\t450\t40\t20\t4000\t15\t9000\t10\t5200\t12\t7800",
    ]
    bed.write_text("\n".join(rows) + "\n")
    return str(bed)


class TestCLAMMSBEDtoVCF:
    """Tests for bin/clamms_bed_to_vcf.py."""

    def test_individual_vcfs_created(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Individual per-sample VCFs are written to the output directory."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        assert os.path.isfile(os.path.join(out, "SAMPLE2_CLAMMS_output.vcf"))

    def test_only_one_sample_per_vcf(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Each VCF file contains exactly one sample column."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        for sample, fname in (
            ("SAMPLE1", "SAMPLE1_CLAMMS_output.vcf"),
            ("SAMPLE2", "SAMPLE2_CLAMMS_output.vcf"),
        ):
            header, _ = _read_vcf_records(os.path.join(out, fname))
            col_header = [l for l in header if l.startswith("#CHROM")][0]
            columns = col_header.split("\t")
            sample_cols = columns[9:]
            assert sample_cols == [sample], (
                f"{fname}: expected sample column [{sample}], got {sample_cols}"
            )

    def test_tool_field_in_info(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Every data record in individual VCFs contains TOOL=CLAMMS in INFO."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CLAMMS_output.vcf", "SAMPLE2_CLAMMS_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                info = record.split("\t")[7]
                assert "TOOL=CLAMMS" in info, (
                    f"TOOL=CLAMMS missing from INFO in {fname}: {info}"
                )

    def test_vcf_header_contains_tool_meta(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The VCF header declares the TOOL INFO field."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        tool_meta = [l for l in header if "##INFO=<ID=TOOL" in l]
        assert tool_meta, "##INFO=<ID=TOOL,...> not found in VCF header"

    def test_fileformat_line_present(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The VCF output starts with ##fileformat=VCFv4.x."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        assert header[0].startswith("##fileformat=VCFv4"), (
            f"First header line is not ##fileformat: {header[0]}"
        )

    def test_source_line_is_clamms(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The ##source header line identifies CLAMMS as the caller."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        source_lines = [l for l in header if l.startswith("##source=")]
        assert source_lines, "No ##source= line in VCF header"
        assert "CLAMMS" in source_lines[0], (
            f"##source does not mention CLAMMS: {source_lines[0]}"
        )

    def test_contig_lines_present(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Contig lines derived from the FAI file appear in the header."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert contig_lines, "No ##contig= lines found in VCF header"

    def test_chrom_field_has_chr_prefix(self, clamms_bed, sample_file, fai_file, tmp_path):
        """CHROM column in every data record has the 'chr' prefix.

        When the input BED already contains 'chr'-prefixed chromosome names
        (e.g. chr1, chr2) the converter must preserve the prefix without
        doubling it (i.e. must not produce 'chrchr1').
        """
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CLAMMS_output.vcf", "SAMPLE2_CLAMMS_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                chrom = record.split("\t")[0]
                assert chrom.startswith("chr"), (
                    f"CHROM field does not start with 'chr' in {fname}: got '{chrom}'"
                )
                assert not chrom.startswith("chrchr"), (
                    f"CHROM field has double 'chr' prefix in {fname}: got '{chrom}'"
                )

    def test_chrom_field_chr_added_when_missing(self, sample_file, fai_file, tmp_path):
        """When the input BED has no 'chr' prefix (e.g. '1', '2'), the
        converter must add 'chr' so that the VCF is GRCh3x-compliant.

        This matches real CLAMMS workflow output: NORMALIZE_DOC runs
        ``sed 's/^chr//g'`` and TRAIN_MODELS runs ``sed 's/chr//g'``,
        so call_cnv produces chromosome names without the prefix.
        """
        import clamms_bed_to_vcf as mod

        # Build a BED file whose CHROM column has no 'chr' prefix
        bed = tmp_path / "clamms_no_chr.bed"
        rows = [
            # CHROM  START  END   INTERVAL       SAMPLE   CNV  MLCN  NUM_WIN  Q_SOME  Q_EXACT  QL_EXT  L_EXT  QR_EXT  R_EXT  QL_CONT  L_CONT  QR_CONT  R_CONT
            "1\t1000\t2000\t1:1000-2000\tSAMPLE1\tDEL\t1\t5\t600\t50\t30\t800\t25\t1200\t15\t1100\t20\t1300",
            "2\t5000\t8000\t2:5000-8000\tSAMPLE2\tDUP\t3\t8\t450\t40\t20\t4000\t15\t9000\t10\t5200\t12\t7800",
        ]
        bed.write_text("\n".join(rows) + "\n")

        out = str(tmp_path / "out")
        mod.process_clamms_data(str(bed), sample_file, fai_file, out)

        for fname in ("SAMPLE1_CLAMMS_output.vcf", "SAMPLE2_CLAMMS_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                chrom = record.split("\t")[0]
                assert chrom.startswith("chr"), (
                    f"CHROM field does not start with 'chr' in {fname}: "
                    f"got '{chrom}' (converter must add 'chr' when input lacks it)"
                )

    def test_chrX_chrY_records_are_retained(self, sample_file, tmp_path):
        """Sex chromosome CNVs should be retained in output VCFs.

        Uses mixed input formatting (`chrX` and `Y`) to verify the converter
        both retains sex chromosomes and normalizes bare `Y` to `chrY`.
        """
        fai = tmp_path / "reference.fa.fai"
        fai.write_text(
            "chr1\t248956422\t112\t70\t71\n"
            "chrX\t156040895\t0\t70\t71\n"
            "chrY\t57227415\t0\t70\t71\n"
        )

        bed = tmp_path / "clamms_xy.bed"
        rows = [
            "chrX\t1000\t2000\tchrX:1000-2000\tSAMPLE1\tDEL\t1\t5\t600\t50\t30\t800\t25\t1200\t15\t1100\t20\t1300",
            "Y\t5000\t8000\tY:5000-8000\tSAMPLE2\tDUP\t3\t8\t450\t40\t20\t4000\t15\t9000\t10\t5200\t12\t7800",
        ]
        bed.write_text("\n".join(rows) + "\n")

        out = str(tmp_path / "out")
        clamms_converter.process_clamms_data(str(bed), sample_file, str(fai), out)

        sample1_header, sample1_data = _read_vcf_records(
            os.path.join(out, "SAMPLE1_CLAMMS_output.vcf")
        )
        sample2_header, sample2_data = _read_vcf_records(
            os.path.join(out, "SAMPLE2_CLAMMS_output.vcf")
        )

        assert sample1_data, "Expected SAMPLE1 VCF to retain chrX CNV record"
        assert sample2_data, "Expected SAMPLE2 VCF to retain chrY CNV record"
        assert any(line.startswith("#CHROM") for line in sample1_header)
        assert any(line.startswith("#CHROM") for line in sample2_header)
        sample1_fields = sample1_data[0].split("\t")
        sample2_fields = sample2_data[0].split("\t")
        assert len(sample1_fields) > 1, "Expected SAMPLE1 VCF data row to be tab-delimited"
        assert len(sample2_fields) > 1, "Expected SAMPLE2 VCF data row to be tab-delimited"
        assert sample1_fields[0] == "chrX", f"Expected CHROM field to be chrX, got {sample1_fields[0]}"
        assert sample2_fields[0] == "chrY", f"Expected CHROM field to be chrY, got {sample2_fields[0]}"

    def test_bed_sample_not_in_sample_file_is_ignored(self, sample_file, fai_file, tmp_path):
        """BED rows for samples not in sample list are ignored without creating VCFs."""
        import clamms_bed_to_vcf as mod

        bed = tmp_path / "clamms_unknown_sample.bed"
        bed.write_text(
            "chr1\t1000\t2000\tchr1:1000-2000\tUNKNOWN_SAMPLE\tDEL\t1\t5\t600\t50\t30\t800\t25\t1200\t15\t1100\t20\t1300\n"
        )

        out = str(tmp_path / "out")
        mod.process_clamms_data(str(bed), sample_file, fai_file, out)

        assert not os.path.exists(os.path.join(out, "UNKNOWN_SAMPLE_CLAMMS_output.vcf"))
        # No known sample has calls in this input, so no VCFs should be created.
        assert not os.path.exists(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        assert not os.path.exists(os.path.join(out, "SAMPLE2_CLAMMS_output.vcf"))

    def test_clamms_malformed_fai_lines_are_ignored(self, clamms_bed, sample_file, tmp_path):
        """Malformed FAI lines should be ignored while valid lines are retained."""
        import clamms_bed_to_vcf as mod

        bad_fai = tmp_path / "bad.fai"
        bad_fai.write_text(
            "chr1\t248956422\t112\t70\t71\n"
            "chr2\n"
            "chr3\tNOT_INT\t0\t70\t71\n"
        )
        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, str(bad_fai), out)

        vcf = os.path.join(out, "SAMPLE1_CLAMMS_output.vcf")
        assert os.path.isfile(vcf)
        header, _ = _read_vcf_records(vcf)
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert any("ID=chr1" in l for l in contig_lines)
        assert not any("ID=chr2" in l for l in contig_lines)
        assert not any("ID=chr3" in l for l in contig_lines)
