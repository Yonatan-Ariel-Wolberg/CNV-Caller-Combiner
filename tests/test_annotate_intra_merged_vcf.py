#!/usr/bin/env python3
import os
import subprocess
import pysam
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
ANNOTATE_SCRIPT = os.path.join(BIN_DIR, "annotate_intra_merged_vcf.py")


def test_annotate_intra_merged_vcf_creates_format_tools(tmp_path):
    """Test that annotate_intra_merged_vcf.py decorates VCF with FORMAT/TOOLS and QS fields."""
    in_vcf_path = str(tmp_path / "survivor_merged.vcf")
    canoes_vcf_path = str(tmp_path / "sample.canoes.normalised.vcf.gz")
    gcnv_vcf_path = str(tmp_path / "sample.gatk_gcnv.normalised.vcf.gz")
    out_vcf_path = str(tmp_path / "sample.annotated.vcf.gz")

    # 1. Create mock canoes normalised VCF
    header_canoes = pysam.VariantHeader()
    header_canoes.contigs.add("chr1", 50000000)
    header_canoes.info.add("SVTYPE", 1, "String", "Type")
    header_canoes.info.add("END", 1, "Integer", "End")
    header_canoes.formats.add("GT", 1, "String", "Genotype")
    header_canoes.formats.add("OQ", 1, "Float", "Original Quality")
    header_canoes.formats.add("OAS", 1, "String", "Original Algorithm")
    header_canoes.samples.add("sample_1")
    with pysam.VariantFile(canoes_vcf_path, "wz", header=header_canoes) as f:
        rec = header_canoes.new_record(contig="chr1", start=100000, stop=105000, alleles=("N", "<DEL>"))
        rec.qual = 92.5
        rec.info["SVTYPE"] = "DEL"
        rec.samples["sample_1"]["GT"] = (0, 1)
        rec.samples["sample_1"]["OQ"] = 74.0
        rec.samples["sample_1"]["OAS"] = "Q_SOME"
        f.write(rec)
    pysam.tabix_index(canoes_vcf_path, preset="vcf", force=True)

    # 2. Create mock GATK-gCNV normalised VCF
    header_gcnv = pysam.VariantHeader()
    header_gcnv.contigs.add("chr1", 50000000)
    header_gcnv.info.add("SVTYPE", 1, "String", "Type")
    header_gcnv.info.add("END", 1, "Integer", "End")
    header_gcnv.formats.add("GT", 1, "String", "Genotype")
    header_gcnv.formats.add("OQ", 1, "Float", "Original Quality")
    header_gcnv.formats.add("OAS", 1, "String", "Original Algorithm")
    header_gcnv.samples.add("sample_1")
    with pysam.VariantFile(gcnv_vcf_path, "wz", header=header_gcnv) as f:
        rec = header_gcnv.new_record(contig="chr1", start=100005, stop=105010, alleles=("N", "<DEL>"))
        rec.qual = 98.0
        rec.info["SVTYPE"] = "DEL"
        rec.samples["sample_1"]["GT"] = (0, 1)
        rec.samples["sample_1"]["OQ"] = 28.0
        rec.samples["sample_1"]["OAS"] = "QS"
        f.write(rec)
    pysam.tabix_index(gcnv_vcf_path, preset="vcf", force=True)

    # 3. Create mock SURVIVOR merged VCF
    header_survivor = pysam.VariantHeader()
    header_survivor.contigs.add("chr1", 50000000)
    header_survivor.add_line("##SAMPLE=<ID=0,File=sample.canoes.normalised.vcf.gz>")
    header_survivor.add_line("##SAMPLE=<ID=1,File=sample.gatk_gcnv.normalised.vcf.gz>")
    header_survivor.info.add("SVTYPE", 1, "String", "Type")
    header_survivor.info.add("END", 1, "Integer", "End")
    header_survivor.info.add("SUPP", 1, "Integer", "Supp")
    header_survivor.info.add("SUPP_VEC", 1, "String", "Supp Vec")
    header_survivor.formats.add("GT", 1, "String", "Genotype")
    header_survivor.samples.add("sample_canoes")
    header_survivor.samples.add("sample_gcnv")
    with pysam.VariantFile(in_vcf_path, "w", header=header_survivor) as f:
        rec = header_survivor.new_record(contig="chr1", start=100000, stop=105000, alleles=("N", "<DEL>"))
        rec.qual = 98.0
        rec.info["SVTYPE"] = "DEL"
        rec.info["SUPP"] = 2
        rec.info["SUPP_VEC"] = "11"
        rec.samples["sample_canoes"]["GT"] = (0, 1)
        rec.samples["sample_gcnv"]["GT"] = (0, 1)
        f.write(rec)

    # 4. Run annotate_intra_merged_vcf.py
    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        ANNOTATE_SCRIPT,
        "--input-vcf", in_vcf_path,
        "--output-vcf", out_vcf_path,
        "--sample-id", "SAMPLE_001",
        "--tool-vcfs", f"CANOES={canoes_vcf_path},GATK_GCNV={gcnv_vcf_path}",
        "--merger-mode", "survivor"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    # 5. Verify output VCF
    with pysam.VariantFile(out_vcf_path) as vcf_out:
        samples = list(vcf_out.header.samples)
        assert samples == ["SAMPLE_001"], f"Expected single sample column 'SAMPLE_001', got {samples}"

        # Check FORMAT definitions
        assert "TOOLS" in vcf_out.header.formats
        assert "CONCORDANCE" in vcf_out.header.formats
        assert "QUAL_NORM" in vcf_out.header.formats
        assert "MEAN_QUAL_NORM" in vcf_out.header.formats
        assert "BEST_TOOL" in vcf_out.header.formats
        assert "BEST_TOOL_RAW_SCORE" in vcf_out.header.formats
        assert "BEST_TOOL_METRIC" in vcf_out.header.formats
        assert "QS_CANOES" in vcf_out.header.formats
        assert "QS_GATK_GCNV" in vcf_out.header.formats
        assert "OQ_CANOES" in vcf_out.header.formats
        assert "OAS_CANOES" in vcf_out.header.formats

        # Check records
        recs = list(vcf_out)
        assert len(recs) == 1
        r = recs[0]
        s = r.samples["SAMPLE_001"]

        assert s["TOOLS"] == ("CANOES", "GATK_GCNV")
        assert s["CONCORDANCE"] == 2
        assert s["QUAL_NORM"] == 98.0
        assert s["MEAN_QUAL_NORM"] == 95.25
        assert s["BEST_TOOL"] == "GATK_GCNV"
        assert s["BEST_TOOL_RAW_SCORE"] == 28.0
        assert s["BEST_TOOL_METRIC"] == "QS"
        assert s["QS_CANOES"] == 92.5
        assert s["QS_GATK_GCNV"] == 98.0
        assert s["OQ_CANOES"] == 74.0
        assert s["OAS_CANOES"] == "Q_SOME"
        assert s["OQ_GATK_GCNV"] == 28.0
        assert s["OAS_GATK_GCNV"] == "QS"
