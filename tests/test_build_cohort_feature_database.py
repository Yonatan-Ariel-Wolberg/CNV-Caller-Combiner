#!/usr/bin/env python3
import os
import sqlite3
import subprocess
import pysam
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_DIR = os.path.join(REPO_ROOT, "bin")
DB_SCRIPT = os.path.join(BIN_DIR, "build_cohort_feature_database.py")


def test_build_cohort_feature_database(tmp_path):
    cohort_vcf = str(tmp_path / "cohort.cnv.vcf.gz")
    features_tsv = str(tmp_path / "SAMPLE_001_features.tsv")
    probes_bed = str(tmp_path / "probes.bed")
    ref_fasta = str(tmp_path / "ref.fa")
    out_dir = str(tmp_path / "db_out")

    # 1. Mock FASTA (long enough: 200,000 bp)
    with open(ref_fasta, "w") as f:
        f.write(">chr1\n" + "GCGCGCGCGC" * 20000 + "\n")
    pysam.faidx(ref_fasta)

    # 2. Mock BED
    with open(probes_bed, "w") as f:
        f.write("chr1\t100000\t100500\tPROBE1\n")
        f.write("chr1\t102000\t102500\tPROBE2\n")

    # 3. Mock Sample Features TSV
    with open(features_tsv, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tis_canoes\tis_gcnv\tqual_norm_canoes\tqual_norm_gatk_gcnv\tconcordance\trd_ratio\tl2r_mean\tmax_qual_norm\tgc_dev\tis_extreme_gc\tis_low_mappability\tprobe_density\toverhang_ratio\n")
        f.write("SAMPLE_001\t1\t100001\t105001\t1\t1\t92.5\t98.0\t2\t0.52\t-0.95\t98.0\t0.02\t0\t0\t0.4\t0.9434\n")

    # 4. Mock Cohort VCF (with decorated FORMAT fields)
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Type")
    header.info.add("END", 1, "Integer", "End")
    header.info.add("SVLEN", 1, "Integer", "Length")
    header.info.add("AF", 1, "Float", "Allele Frequency")
    header.info.add("AC", 1, "Integer", "Allele Count")
    header.info.add("AN", 1, "Integer", "Total Alleles")
    header.info.add("PRECISE", 0, "Flag", "Precise")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("TOOLS", ".", "String", "Tools")
    header.formats.add("MEAN_QUAL_NORM", 1, "Float", "Mean Normalized Quality")
    header.formats.add("BEST_TOOL", 1, "String", "Best Tool")
    header.formats.add("BEST_TOOL_RAW_SCORE", 1, "Float", "Best Tool Raw Score")
    header.formats.add("BEST_TOOL_METRIC", 1, "String", "Best Tool Metric")
    header.formats.add("QS_CANOES", 1, "Float", "Canoes QS")
    header.formats.add("QS_GATK_GCNV", 1, "Float", "GCNV QS")
    header.formats.add("OQ_CANOES", 1, "Float", "Canoes OQ")
    header.formats.add("OQ_GATK_GCNV", 1, "Float", "GCNV OQ")
    header.formats.add("OAS_CANOES", 1, "String", "Canoes OAS")
    header.formats.add("OAS_GATK_GCNV", 1, "String", "GCNV OAS")
    header.samples.add("SAMPLE_001")
    header.samples.add("SAMPLE_002")

    with pysam.VariantFile(cohort_vcf, "wz", header=header) as f:
        rec = header.new_record(contig="chr1", start=100000, stop=105000, alleles=("N", "<DEL>"))
        rec.qual = 98.0
        rec.info["SVTYPE"] = "DEL"
        rec.info["SVLEN"] = -5000
        rec.info["AF"] = 0.5
        rec.info["AC"] = 1
        rec.info["AN"] = 4
        rec.info["PRECISE"] = True
        rec.samples["SAMPLE_001"]["GT"] = (0, 1)
        rec.samples["SAMPLE_001"]["TOOLS"] = ("CANOES", "GATK_GCNV")
        rec.samples["SAMPLE_001"]["MEAN_QUAL_NORM"] = 95.25
        rec.samples["SAMPLE_001"]["BEST_TOOL"] = "GATK_GCNV"
        rec.samples["SAMPLE_001"]["BEST_TOOL_RAW_SCORE"] = 28.0
        rec.samples["SAMPLE_001"]["BEST_TOOL_METRIC"] = "QS"
        rec.samples["SAMPLE_001"]["QS_CANOES"] = 92.5
        rec.samples["SAMPLE_001"]["QS_GATK_GCNV"] = 98.0
        rec.samples["SAMPLE_001"]["OQ_CANOES"] = 74.0
        rec.samples["SAMPLE_001"]["OQ_GATK_GCNV"] = 28.0
        rec.samples["SAMPLE_001"]["OAS_CANOES"] = "Q_SOME"
        rec.samples["SAMPLE_001"]["OAS_GATK_GCNV"] = "QS"
        rec.samples["SAMPLE_002"]["GT"] = (0, 0)
        f.write(rec)
    pysam.tabix_index(cohort_vcf, preset="vcf", force=True)

    # 5. Run build_cohort_feature_database.py
    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        DB_SCRIPT,
        "--cohort-vcf", cohort_vcf,
        "--sample-features-list", features_tsv,
        "--target-bed", probes_bed,
        "--fasta", ref_fasta,
        "--outdir", out_dir,
        "--output-prefix", "test_cohort"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    # 6. Verify SQLite Database
    db_path = os.path.join(out_dir, "test_cohort.db")
    assert os.path.exists(db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Check cohort variants
    cursor.execute("SELECT variant_id, cohort_af, is_precise FROM cohort_variants")
    vars_row = cursor.fetchone()
    assert vars_row is not None
    assert vars_row[0] == "1:100001-105001:DEL"
    assert vars_row[1] == 0.5
    assert vars_row[2] == 1

    # Check carrier calls with tool-specific quality scores, best tool attributes, and caller flags
    cursor.execute("SELECT sample_id, genotype, supporting_callers, qs_canoes, qs_gatk_gcnv, rd_ratio, mean_qual_norm, best_tool, best_tool_raw_score, best_tool_metric, number_of_callers, is_canoes, is_gatk_gcnv, is_gcnv, gc_dev, is_extreme_gc, is_low_mappability, probe_density, overhang_ratio FROM carrier_calls")
    car_row = cursor.fetchone()
    assert car_row is not None
    assert car_row[0] == "SAMPLE_001"
    assert car_row[1] == "0/1"
    assert "CANOES" in car_row[2] and "GATK_GCNV" in car_row[2]
    assert car_row[3] == 92.5
    assert car_row[4] == 98.0
    assert car_row[5] == 0.52
    assert car_row[6] == 95.25
    assert car_row[7] == "GATK_GCNV"
    assert car_row[8] == 28.0
    assert car_row[9] == "QS"
    assert car_row[10] == 2
    assert car_row[11] == 1
    assert car_row[12] == 1
    assert car_row[13] == 1
    assert car_row[14] == 0.02
    assert car_row[15] == 0
    assert car_row[16] == 0
    assert car_row[17] == 0.4
    assert car_row[18] == 0.9434

    # Check denormalized carrier matrix
    cursor.execute("SELECT sample_id, variant_id, cohort_af, gc_content, n_target_probes, qs_canoes, qs_gatk_gcnv, best_tool, best_tool_raw_score, number_of_callers, is_canoes, is_gatk_gcnv FROM cohort_carrier_matrix")
    mat_row = cursor.fetchone()
    assert mat_row is not None
    assert mat_row[0] == "SAMPLE_001"
    assert mat_row[1] == "1:100001-105001:DEL"
    assert mat_row[2] == 0.5
    assert mat_row[3] == 1.0
    assert mat_row[4] == 2
    assert mat_row[5] == 92.5
    assert mat_row[6] == 98.0
    assert mat_row[7] == "GATK_GCNV"
    assert mat_row[8] == 28.0
    assert mat_row[9] == 2
    assert mat_row[10] == 1
    assert mat_row[11] == 1

    conn.close()

    # 7. Verify compressed TSVs
    assert os.path.exists(os.path.join(out_dir, "test_cohort.variants.tsv.gz"))
    assert os.path.exists(os.path.join(out_dir, "test_cohort.carriers.tsv.gz"))


def test_build_cohort_feature_database_single_caller_exomedepth(tmp_path):
    """Test building cohort database when only 1 caller (ExomeDepth) was run."""
    cohort_vcf = str(tmp_path / "cohort_single.cnv.vcf.gz")
    features_tsv = str(tmp_path / "SAMPLE_ED_features.tsv")
    probes_bed = str(tmp_path / "probes.bed")
    ref_fasta = str(tmp_path / "ref.fa")
    out_dir = str(tmp_path / "db_single_out")

    # 1. Mock FASTA
    with open(ref_fasta, "w") as f:
        f.write(">chr1\n" + "GCGCGCGCGC" * 10000 + "\n")
    pysam.faidx(ref_fasta)

    # 2. Mock BED
    with open(probes_bed, "w") as f:
        f.write("chr1\t100000\t100500\tPROBE1\n")

    # 3. Mock Sample Features TSV with only ExomeDepth
    with open(features_tsv, "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tis_exomedepth\tqual_norm_exomedepth\tconcordance\tnumber_of_callers\trd_ratio\tl2r_mean\tmax_qual_norm\n")
        f.write("SAMPLE_ED\t1\t100001\t102001\t1\t94.0\t1\t1\t0.45\t-1.12\t94.0\n")

    # 4. Mock Cohort VCF (from single caller)
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Type")
    header.info.add("END", 1, "Integer", "End")
    header.info.add("SVLEN", 1, "Integer", "Length")
    header.info.add("AF", 1, "Float", "Allele Frequency")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("TOOLS", ".", "String", "Tools")
    header.formats.add("QS_EXOMEDEPTH", 1, "Float", "ED QS")
    header.formats.add("OQ_EXOMEDEPTH", 1, "Float", "ED OQ")
    header.formats.add("OAS_EXOMEDEPTH", 1, "String", "ED OAS")
    header.samples.add("SAMPLE_ED")

    with pysam.VariantFile(cohort_vcf, "wz", header=header) as f:
        rec = header.new_record(contig="chr1", start=100000, stop=102000, alleles=("N", "<DEL>"))
        rec.qual = 94.0
        rec.info["SVTYPE"] = "DEL"
        rec.info["SVLEN"] = -2000
        rec.info["AF"] = 1.0
        rec.samples["SAMPLE_ED"]["GT"] = (0, 1)
        rec.samples["SAMPLE_ED"]["TOOLS"] = "EXOMEDEPTH"
        rec.samples["SAMPLE_ED"]["QS_EXOMEDEPTH"] = 94.0
        rec.samples["SAMPLE_ED"]["OQ_EXOMEDEPTH"] = 45.0
        rec.samples["SAMPLE_ED"]["OAS_EXOMEDEPTH"] = "BF"
        f.write(rec)
    pysam.tabix_index(cohort_vcf, preset="vcf", force=True)

    # 5. Run build_cohort_feature_database.py
    cmd = [
        os.path.join(REPO_ROOT, "test_env/bin/python3"),
        DB_SCRIPT,
        "--cohort-vcf", cohort_vcf,
        "--sample-features-list", features_tsv,
        "--target-bed", probes_bed,
        "--fasta", ref_fasta,
        "--outdir", out_dir,
        "--output-prefix", "single_cohort"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    # 6. Verify SQLite Database
    db_path = os.path.join(out_dir, "single_cohort.db")
    assert os.path.exists(db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT sample_id, number_of_callers, is_exomedepth, is_canoes, qs_exomedepth, oq_exomedepth, oas_exomedepth, best_tool, best_tool_raw_score FROM carrier_calls")
    row = cursor.fetchone()
    assert row is not None
    assert row[0] == "SAMPLE_ED"
    assert row[1] == 1
    assert row[2] == 1
    assert row[3] == 0
    assert row[4] == 94.0
    assert row[5] == 45.0
    assert row[6] == "BF"
    assert row[7] == "EXOMEDEPTH"
    assert row[8] == 45.0

    conn.close()

