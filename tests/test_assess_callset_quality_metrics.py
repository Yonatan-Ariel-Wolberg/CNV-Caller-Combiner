#!/usr/bin/env python3
"""
Tests for bin/assess_callset_quality_metrics.py.
"""

import os
import subprocess
import pytest
import pysam

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
ASSESS_SCRIPT = os.path.join(BIN_DIR, "assess_callset_quality_metrics.py")


def create_test_vcf(vcf_path, records, sample_id="SAMPLE_01"):
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 10000000)
    header.info.add("SVTYPE", 1, "String", "Variant type")
    header.info.add("END", 1, "Integer", "End position")
    header.info.add("SVLEN", 1, "Integer", "Length")
    header.formats.add("GT", 1, "String", "Genotype")
    header.samples.add(sample_id)

    with pysam.VariantFile(vcf_path, "w", header=header) as vf:
        for chrom, start, end, svtype, qual in records:
            rec = header.new_record(contig=chrom, start=start, stop=end, alleles=("N", f"<{svtype}>"))
            rec.qual = qual
            rec.info["SVTYPE"] = svtype
            rec.info["SVLEN"] = end - start
            rec.samples[sample_id]["GT"] = (0, 1)
            vf.write(rec)


def test_assess_callset_quality_metrics(tmp_path):
    vcf1 = str(tmp_path / "sample1.vcf")
    vcf2 = str(tmp_path / "sample2.vcf")
    bed = str(tmp_path / "probes.bed")
    summary = str(tmp_path / "summary.txt")
    csv = str(tmp_path / "metrics.csv")

    create_test_vcf(vcf1, [
        ("chr1", 1000, 3000, "DEL", 80.0),
        ("chr1", 5000, 10000, "DUP", 65.0)
    ], "SAMPLE_01")

    create_test_vcf(vcf2, [
        ("chr1", 20000, 25000, "DEL", 95.0)
    ], "SAMPLE_02")

    with open(bed, "w") as f:
        f.write("chr1\t1200\t1500\n")
        f.write("chr1\t1800\t2100\n")
        f.write("chr1\t6000\t7000\n")

    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        ASSESS_SCRIPT,
        "--vcf-files", vcf1, vcf2,
        "--stage-name", "Test_Stage",
        "--probes-bed", bed,
        "--output-summary", summary,
        "--output-csv", csv
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}\nOutput: {res.stdout}"

    assert os.path.exists(summary)
    assert os.path.exists(csv)

    with open(summary) as f:
        content = f.read()
        assert "Total CNV Calls:" in content
        assert "3" in content
        assert "Deletions (DEL):" in content
        assert "Duplications (DUP):" in content
        assert "Median Length:" in content
        assert "Test_Stage" in content
        assert "3-WAY QUALITY METRIC COMPARISON: TOTAL CNVS vs DELETIONS (DEL) vs DUPLICATIONS (DUP)" in content
        assert "Median Length (DEL):" in content
        assert "Median Length (DUP):" in content
        assert "Total Genomic Footprint (DEL):" in content
        assert "Total Genomic Footprint (DUP):" in content
        assert "Median Probes per CNV (DEL):" in content
        assert "Median Probes per CNV (DUP):" in content
        assert "Mean Calls per Sample (DEL):" in content
        assert "Mean Calls per Sample (DUP):" in content


def test_assess_callset_empty_fallback(tmp_path):
    empty_vcf = str(tmp_path / "empty.vcf")
    create_test_vcf(empty_vcf, [], "SAMPLE_EMPTY")
    summary = str(tmp_path / "empty_summary.txt")
    csv = str(tmp_path / "empty_metrics.csv")

    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        ASSESS_SCRIPT,
        "--vcf-files", empty_vcf,
        "--stage-name", "Empty_Stage",
        "--output-summary", summary,
        "--output-csv", csv
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    with open(summary) as f:
        content = f.read()
        assert "Total CNVs: 0" in content
