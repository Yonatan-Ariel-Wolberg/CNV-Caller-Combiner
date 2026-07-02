#!/usr/bin/env python3
"""
Tests for the 'indelible' workflow pipeline.

Validates the complete indelible workflow at the Python level:
  1. The FILTER_INDELIBLE step (awk: columns 39 & 40 < 2) correctly filters
     annotated TSV rows.
  2. The bin/indelible_tsv_to_vcf.py conversion script handles the full
     annotated TSV format produced by 'indelible.py annotate'.
  3. The end-to-end flow (annotated TSV → filter → VCF) produces correct
     output for multiple samples.
"""

import os
import sys
import pytest

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))


# ---------------------------------------------------------------------------
# Full annotated TSV column layout (40 columns, matching indelible.py annotate)
# The first 37 columns follow the actual INDELIBLE output format; the last 3
# are added by 'indelible.py database' + 'indelible.py annotate'.
#
# FILTER_INDELIBLE awk command:
#   awk '{ if (($39 < 2) && ($40 < 2)) { print } }' annotated > filtered.tsv
#   $39 (1-based) = 0-based index 38 = AF_freq
#   $40 (1-based) = 0-based index 39 = BP_freq
#
# Column order (1-based):
#   1  chrom                 20 seq_longest           38 n_samples
#   2  position              21 predicted             39 AF_freq
#   3  coverage              22 prob_N                40 BP_freq
#   4  insertion_context     23 prob_Y
#   5  deletion_context      24 ddg2p
#   6  sr_total              25 hgnc
#   7  sr_total_long         26 hgnc_constrained
#   8  sr_total_short        27 exonic
#   9  sr_long_5             28 transcripts
#  10  sr_short_5            29 exon_numbers
#  11  sr_long_3             30 maf
#  12  sr_short_3            31 blast_hit
#  13  sr_entropy            32 blast_strand
#  14  context_entropy       33 blast_identity
#  15  entropy_upstream      34 blast_dist
#  16  entropy_downstream    35 blast_hgnc
#  17  sr_sw_similarity      36 blast_hgnc_constrained
#  18  avg_avg_sr_qual       37 blast_ddg2p
#  19  avg_mapq
# ---------------------------------------------------------------------------
_ANNOTATED_HEADER = (
    "chrom\tposition\tcoverage\tinsertion_context\tdeletion_context\t"
    "sr_total\tsr_total_long\tsr_total_short\tsr_long_5\tsr_short_5\t"
    "sr_long_3\tsr_short_3\tsr_entropy\tcontext_entropy\t"
    "entropy_upstream\tentropy_downstream\tsr_sw_similarity\t"
    "avg_avg_sr_qual\tavg_mapq\tseq_longest\tpredicted\tprob_N\tprob_Y\t"
    "ddg2p\thgnc\thgnc_constrained\texonic\ttranscripts\texon_numbers\t"
    "maf\tblast_hit\tblast_strand\tblast_identity\tblast_dist\t"
    "blast_hgnc\tblast_hgnc_constrained\tblast_ddg2p\t"
    "n_samples\tAF_freq\tBP_freq"
)


def _make_annotated_row(chrom, position, predicted, prob_Y, af_freq, bp_freq,
                        hgnc="GENE1", exonic="False"):
    """Return a 40-column annotated TSV row (tab-separated string).

    Columns follow the actual INDELIBLE output format (37 base columns) plus
    the 3 database-derived columns appended by 'indelible.py annotate':
      38: n_samples  (number of database samples carrying this variant)
      39: AF_freq    (allele frequency in database; awk $39)
      40: BP_freq    (breakpoint frequency in database; awk $40)
    """
    prob_N = round(1.0 - float(prob_Y), 6)
    return (
        f"{chrom}\t{position}\t50\t0\t0\t"
        f"10\t7\t3\t4\t2\t3\t1\t"
        f"1.5\t1.2\t1.3\t1.1\t1.0\t35.0\t29.0\t"
        f"ACGT\t{predicted}\t{prob_N}\t{prob_Y}\t"
        f"NA\t{hgnc}\tNA\t{exonic}\tENST001\t3\t0.01\t"
        f"NA\tNA\tNA\tNA\tNA\tNA\tNA\t"
        f"5\t{af_freq}\t{bp_freq}"
    )


def _apply_filter(input_tsv, output_tsv):
    """
    Python equivalent of the FILTER_INDELIBLE awk command:
        awk '{ if (($39 < 2) && ($40 < 2)) { print } }' annotated > filtered.tsv

    In awk, non-numeric strings evaluate to 0 in numeric context, so the header
    row is always kept (0 < 2 is true for both columns).
    """
    with open(input_tsv) as fin, open(output_tsv, 'w') as fout:
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 40:
                try:
                    col39 = float(parts[38])  # AF_freq (1-based col 39)
                    col40 = float(parts[39])  # BP_freq (1-based col 40)
                    if col39 < 2 and col40 < 2:
                        fout.write(line)
                except ValueError:
                    # Non-numeric (header) → 0 < 2 in awk → kept
                    fout.write(line)
            else:
                fout.write(line)


def _read_vcf_records(vcf_path):
    """Return (header_lines, data_lines) from a VCF file."""
    header, data = [], []
    with open(vcf_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                header.append(line)
            else:
                data.append(line)
    return header, data


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def fai_file(tmp_path):
    """Minimal FASTA index covering chromosomes used in tests."""
    fai = tmp_path / "reference.fa.fai"
    fai.write_text(
        "chr1\t248956422\t112\t70\t71\n"
        "chr2\t242193529\t252513167\t70\t71\n"
    )
    return str(fai)


@pytest.fixture()
def full_annotated_tsv(tmp_path):
    """
    Realistic INDELIBLE annotated TSV (40 columns) with four data rows:
      pos 1000 – predicted=Y, AF_freq=1, BP_freq=1  → PASS filter
      pos 2000 – predicted=Y, AF_freq=3, BP_freq=1  → filtered OUT (AF_freq >= 2)
      pos 5000 – predicted=N, AF_freq=1, BP_freq=5  → filtered OUT (BP_freq >= 2)
      pos 6000 – predicted=N, AF_freq=0, BP_freq=0  → PASS filter
    """
    tsv = tmp_path / "SAMPLE1.counts.scored.annotated"
    rows = [
        _ANNOTATED_HEADER,
        _make_annotated_row("chr1", "1000", "Y", "0.95", "1", "1"),
        _make_annotated_row("chr1", "2000", "Y", "0.90", "3", "1"),
        _make_annotated_row("chr2", "5000", "N", "0.20", "1", "5"),
        _make_annotated_row("chr2", "6000", "N", "0.10", "0", "0"),
    ]
    tsv.write_text("\n".join(rows) + "\n")
    return str(tsv)


# ===========================================================================
# FILTER_INDELIBLE step
# ===========================================================================

class TestFilterIndelible:
    """Tests for the FILTER_INDELIBLE awk logic (columns 39 & 40 < 2)."""

    def test_passing_rows_are_kept(self, full_annotated_tsv, tmp_path):
        """Rows where AF_freq < 2 AND BP_freq < 2 are retained."""
        filtered = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered)

        with open(filtered) as f:
            data_lines = [l for l in f if not l.startswith("chrom\t") and l.strip()]
        positions = [l.split('\t')[1] for l in data_lines]
        assert "1000" in positions, "Row at position 1000 should pass filter"
        assert "6000" in positions, "Row at position 6000 should pass filter"

    def test_high_af_freq_rows_are_removed(self, full_annotated_tsv, tmp_path):
        """Rows with AF_freq >= 2 are removed."""
        filtered = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered)

        with open(filtered) as f:
            data_lines = [l for l in f if not l.startswith("chrom\t") and l.strip()]
        positions = [l.split('\t')[1] for l in data_lines]
        assert "2000" not in positions, "Row with AF_freq=3 should be filtered out"

    def test_high_bp_freq_rows_are_removed(self, full_annotated_tsv, tmp_path):
        """Rows with BP_freq >= 2 are removed."""
        filtered = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered)

        with open(filtered) as f:
            data_lines = [l for l in f if not l.startswith("chrom\t") and l.strip()]
        positions = [l.split('\t')[1] for l in data_lines]
        assert "5000" not in positions, "Row with BP_freq=5 should be filtered out"

    def test_correct_number_of_rows_after_filter(self, full_annotated_tsv, tmp_path):
        """Exactly 2 of the 4 data rows survive the filter."""
        filtered = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered)

        with open(filtered) as f:
            data_lines = [l for l in f if not l.startswith("chrom\t") and l.strip()]
        assert len(data_lines) == 2, f"Expected 2 rows after filter, got {len(data_lines)}"

    def test_boundary_values_at_threshold(self, tmp_path):
        """Rows with AF_freq or BP_freq exactly equal to 2 are filtered out (not < 2)."""
        tsv = tmp_path / "boundary.tsv"
        rows = [
            _ANNOTATED_HEADER,
            _make_annotated_row("chr1", "100", "Y", "0.9", "2", "1"),  # AF_freq=2 → out
            _make_annotated_row("chr1", "200", "Y", "0.9", "1", "2"),  # BP_freq=2 → out
            _make_annotated_row("chr1", "300", "Y", "0.9", "1", "1"),  # both 1 → in
        ]
        tsv.write_text("\n".join(rows) + "\n")

        filtered = str(tmp_path / "filtered.tsv")
        _apply_filter(str(tsv), filtered)

        with open(filtered) as f:
            data_lines = [l for l in f if not l.startswith("chrom\t") and l.strip()]
        positions = [l.split('\t')[1] for l in data_lines]
        assert "300" in positions
        assert "100" not in positions
        assert "200" not in positions


# ===========================================================================
# Full annotated TSV → VCF conversion
# ===========================================================================

class TestFullAnnotatedTsvToVcf:
    """Tests for indelible_tsv_to_vcf.py with the full 40-column annotated format."""

    def test_vcf_created_from_full_annotated_tsv(self, full_annotated_tsv, fai_file, tmp_path):
        """VCF is produced from a full-column annotated TSV without errors."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(full_annotated_tsv, out, "SAMPLE1", fai_file)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))

    def test_correct_record_count(self, full_annotated_tsv, fai_file, tmp_path):
        """All four data rows in the annotated TSV produce four VCF records."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(full_annotated_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert len(data) == 4, f"Expected 4 VCF records, got {len(data)}"

    def test_optional_hgnc_field_in_info(self, full_annotated_tsv, fai_file, tmp_path):
        """The HGNC optional annotation field appears in the VCF INFO column."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(full_annotated_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert data, "No VCF data records"
        first_info = data[0].split('\t')[7]
        assert "HGNC=GENE1" in first_info, f"HGNC=GENE1 not found in INFO: {first_info}"

    def test_optional_exonic_field_in_info(self, full_annotated_tsv, fai_file, tmp_path):
        """The EXONIC optional annotation field appears in the VCF INFO column."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(full_annotated_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert data, "No VCF data records"
        first_info = data[0].split('\t')[7]
        assert "EXONIC=False" in first_info, f"EXONIC=False not found in INFO: {first_info}"


# ===========================================================================
# End-to-end: filter → convert
# ===========================================================================

class TestIndelibleWorkflowEndToEnd:
    """End-to-end tests: annotated TSV → FILTER_INDELIBLE → VCF conversion."""

    def test_filtered_tsv_produces_correct_vcf_record_count(
            self, full_annotated_tsv, fai_file, tmp_path):
        """After filtering (2 of 4 rows pass), the VCF contains exactly 2 records."""
        import indelible_tsv_to_vcf as mod

        filtered_tsv = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered_tsv)

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(filtered_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert len(data) == 2, f"Expected 2 VCF records after filter, got {len(data)}"

    def test_filtered_vcf_contains_correct_positions(
            self, full_annotated_tsv, fai_file, tmp_path):
        """After filtering, the VCF contains only the positions that passed the filter."""
        import indelible_tsv_to_vcf as mod

        filtered_tsv = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered_tsv)

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(filtered_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        positions = {r.split('\t')[1] for r in data}
        assert positions == {"1000", "6000"}, (
            f"Expected positions {{1000, 6000}}, got {positions}"
        )

    def test_multiple_samples_each_get_separate_vcf(self, tmp_path, fai_file):
        """Running the conversion for each sample produces one VCF per sample."""
        import indelible_tsv_to_vcf as mod

        # Create a separate annotated TSV for each of three samples
        samples = ["CHILD_01_1", "MOM_02_2", "DAD_03_3"]
        out = str(tmp_path / "out")

        for sample in samples:
            tsv = tmp_path / f"{sample}.counts.scored.annotated"
            tsv.write_text(
                _ANNOTATED_HEADER + "\n" +
                _make_annotated_row("chr1", "1000", "Y", "0.95", "1", "1") + "\n"
            )
            mod.convert_indelible_tsv_to_vcf(str(tsv), out, sample, fai_file)

        for sample in samples:
            assert os.path.isfile(os.path.join(out, f"{sample}_INDELIBLE_output.vcf")), (
                f"VCF not created for sample {sample}"
            )

    def test_pass_filter_preserved_after_end_to_end(
            self, full_annotated_tsv, fai_file, tmp_path):
        """Rows that pass the filter retain their PASS/LowQuality status in the VCF."""
        import indelible_tsv_to_vcf as mod

        filtered_tsv = str(tmp_path / "filtered.tsv")
        _apply_filter(full_annotated_tsv, filtered_tsv)

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(filtered_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        filters = {r.split('\t')[1]: r.split('\t')[6] for r in data}
        # position 1000 → predicted=Y → PASS
        assert filters.get("1000") == "PASS", (
            f"Expected PASS for position 1000, got {filters.get('1000')}"
        )
        # position 6000 → predicted=N → LowQuality
        assert filters.get("6000") == "LowQuality", (
            f"Expected LowQuality for position 6000, got {filters.get('6000')}"
        )
