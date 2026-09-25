#!/usr/bin/env python3
"""
test_breakpoint_refinement.py

Comprehensive test suite verifying:
1. Size-aware coordinate displacement and reciprocal overlap constraints for focal CNVs (< 10 kb).
2. Distant intronic SV rejection for focal calls despite large global max_delta.
3. Successful intronic bridging for macro-CNVs (>= 10 kb).
4. Extraction of SvABA / DRAGEN-SV microhomology (HOMSEQ, HOMLEN) and insertion (INSERTION, INSLEN) signatures.
5. Biophysical mutational mechanism classification (NHEJ, MMBIR_FOSTES, NAHR_CANDIDATE, COMPLEX_INS).
6. Ingestion into SQLite cohort feature database.
"""

import os
import sqlite3
import subprocess
import pysam
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REFINE_SCRIPT = os.path.join(REPO_ROOT, "bin", "refine_breakpoints.py")
DB_SCRIPT = os.path.join(REPO_ROOT, "bin", "build_cohort_feature_database.py")
PYTHON_EXEC = os.path.join(REPO_ROOT, "test_env", "bin", "python3")


def create_consensus_vcf(filepath, records):
    """Create a test consensus VCF file."""
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Variant type")
    header.info.add("END", 1, "Integer", "End coordinate")
    header.info.add("SVLEN", 1, "Integer", "Variant length")
    header.info.add("AF", "A", "Float", "Allele Frequency")
    header.formats.add("GT", 1, "String", "Genotype")
    header.samples.add("SAMPLE1")

    with pysam.VariantFile(filepath, "wz", header=header) as vcf_out:
        for r in records:
            rec = vcf_out.header.new_record(
                contig=r["chrom"],
                start=r["start"] - 1,  # pysam 0-based start
                stop=r["end"],
                alleles=("N", f"<{r['svtype']}>"),
            )
            rec.info["SVTYPE"] = r["svtype"]
            rec.info["SVLEN"] = r["end"] - r["start"] if r["svtype"] == "DUP" else -(r["end"] - r["start"])
            rec.info["AF"] = (0.5,)
            rec.samples["SAMPLE1"]["GT"] = (0, 1)
            vcf_out.write(rec)

    pysam.tabix_index(filepath, preset="vcf", force=True)


def create_sv_vcf(filepath, records):
    """Create a test SV VCF with SvABA / DRAGEN microhomology and insertion annotations."""
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Variant type")
    header.info.add("END", 1, "Integer", "End coordinate")
    header.info.add("SVLEN", 1, "Integer", "Variant length")
    header.info.add("HOMSEQ", 1, "String", "Microhomology sequence")
    header.info.add("HOMLEN", 1, "Integer", "Microhomology length")
    header.info.add("INSERTION", 1, "String", "Inserted sequence")
    header.info.add("INSLEN", 1, "Integer", "Insertion length")
    header.formats.add("GT", 1, "String", "Genotype")
    header.samples.add("SAMPLE1")

    with pysam.VariantFile(filepath, "wz", header=header) as vcf_out:
        for r in records:
            rec = vcf_out.header.new_record(
                contig=r["chrom"],
                start=r["start"] - 1,
                stop=r["end"],
                alleles=("N", f"<{r['svtype']}>"),
            )
            rec.info["SVTYPE"] = r["svtype"]
            rec.info["SVLEN"] = r["end"] - r["start"] if r["svtype"] == "DUP" else -(r["end"] - r["start"])
            rec.filter.add("PASS")
            if "HOMSEQ" in r:
                rec.info["HOMSEQ"] = r["HOMSEQ"]
            if "HOMLEN" in r:
                rec.info["HOMLEN"] = r["HOMLEN"]
            if "INSERTION" in r:
                rec.info["INSERTION"] = r["INSERTION"]
            if "INSLEN" in r:
                rec.info["INSLEN"] = r["INSLEN"]
            rec.samples["SAMPLE1"]["GT"] = (0, 1)
            vcf_out.write(rec)

    pysam.tabix_index(filepath, preset="vcf", force=True)


def test_focal_cnv_size_aware_rejection_and_acceptance(tmp_path):
    """
    Verify that:
    - A focal CNV (1,000 bp) REJECTS a distant intronic SV (15 kb away) even with --max-delta 25000.
    - A focal CNV (1,000 bp) ACCEPTS a closely matched intronic SV (400 bp delta, encapsulating).
    - A macro CNV (50,000 bp) ACCEPTS an SV with 16,000 bp intronic delta.
    """
    consensus_vcf = str(tmp_path / "consensus.vcf.gz")
    sv_vcf = str(tmp_path / "sv.vcf.gz")
    out_vcf = str(tmp_path / "refined.vcf.gz")
    summary_file = str(tmp_path / "summary.txt")

    consensus_records = [
        # Call 1: Focal DEL (1,000 bp) - targets exon at 100,000 - 101,000
        {"chrom": "chr1", "start": 100000, "end": 101000, "svtype": "DEL"},
        # Call 2: Focal DUP (800 bp) - targets exon at 200,000 - 200,800
        # Will attempt to match distant SV at 215,000 - 216,000 -> MUST BE REJECTED
        {"chrom": "chr1", "start": 200000, "end": 200800, "svtype": "DUP"},
        # Call 3: Macro DEL (50,000 bp) - multi-exon at 300,000 - 350,000
        # Will match SV with deep intronic breakpoints at 292,000 - 358,000 -> MUST BE ACCEPTED
        {"chrom": "chr1", "start": 300000, "end": 350000, "svtype": "DEL"},
    ]
    create_consensus_vcf(consensus_vcf, consensus_records)

    sv_records = [
        # SV 1: Closely spans Call 1 with 6 bp microhomology (MMBIR/FoSTeS)
        {
            "chrom": "chr1",
            "start": 99800,
            "end": 101200,
            "svtype": "DEL",
            "HOMSEQ": "AGCTAG",
            "HOMLEN": 6,
        },
        # SV 2: Distant intronic DUP 14,200 bp away from Call 2
        {
            "chrom": "chr1",
            "start": 215000,
            "end": 215800,
            "svtype": "DUP",
            "HOMLEN": 1,
        },
        # SV 3: Deep intronic breakpoints spanning Call 3 (delta = 8,000 + 8,000 = 16,000 bp)
        # Non-templated complex insertion of 14 bp
        {
            "chrom": "chr1",
            "start": 292000,
            "end": 358000,
            "svtype": "DEL",
            "INSERTION": "ATCGATCGATCGAT",
            "INSLEN": 14,
        },
    ]
    create_sv_vcf(sv_vcf, sv_records)

    cmd = [
        PYTHON_EXEC,
        REFINE_SCRIPT,
        "--consensus-vcf", consensus_vcf,
        "--sv-vcf", sv_vcf,
        "--output-vcf", out_vcf,
        "--summary", summary_file,
        "--max-delta", "25000",
        "--focal-threshold", "10000",
        "--min-focal-delta", "1000",
        "--focal-slack-ratio", "1.0",
        "--min-recip-overlap", "0.20",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    # Verify refined VCF records
    with pysam.VariantFile(out_vcf) as vf:
        recs = list(vf)
        assert len(recs) == 3

        # Call 1: Refined base-pair exact
        rec1 = recs[0]
        assert rec1.pos == 99800
        assert rec1.stop == 101200
        assert "PRECISE" in rec1.info
        assert "IMPRECISE" not in rec1.info
        assert rec1.info["ORIG_POS"] == 100000
        assert rec1.info["ORIG_END"] == 101000
        assert rec1.info["SLACK_DIST"] == 400
        assert rec1.info["HOMSEQ"] == "AGCTAG"
        assert rec1.info["HOMLEN"] == 6
        assert rec1.info["MECH"] == "MMBIR_FOSTES"
        assert float(rec1.info["REFINED_RECIP_OVERLAP"]) > 0.5

        # Call 2: Distant SV was rejected due to size-aware constraint
        rec2 = recs[1]
        assert rec2.pos == 200000
        assert rec2.stop == 200800
        assert "IMPRECISE" in rec2.info
        assert "PRECISE" not in rec2.info
        assert rec2.info["MECH"] == "UNREFINED"
        assert int(rec2.info["SIZE_AWARE_DELTA"]) == 1000  # min_focal_delta for 800 bp

        # Call 3: Macro CNV was refined using intronic SV
        rec3 = recs[2]
        assert rec3.pos == 292000
        assert rec3.stop == 358000
        assert "PRECISE" in rec3.info
        assert rec3.info["ORIG_POS"] == 300000
        assert rec3.info["ORIG_END"] == 350000
        assert rec3.info["SLACK_DIST"] == 16000
        assert rec3.info["INSERTION"] == "ATCGATCGATCGAT"
        assert rec3.info["INSLEN"] == 14
        assert rec3.info["MECH"] == "COMPLEX_INS"

    # Verify summary report contents
    with open(summary_file, "r") as f:
        report = f.read()
        assert "Total Input CNVs:         3" in report
        assert "Focal (< 10000 bp):  2" in report
        assert "Macro (>= 10000 bp): 1" in report
        assert "Focal Refined:        1 (50.00%)" in report
        assert "Macro Refined:        1 (100.00%)" in report
        assert "MMBIR/FoSTeS (4-24 bp homology):     1" in report
        assert "Complex Insertion (> 10 bp ins):     1" in report
        assert "Unrefined (Exon-Bounds):  1 (33.33%)" in report


def test_mutational_mechanism_classification_full_spectrum(tmp_path):
    """Test all 4 primary mutational mechanisms (NHEJ, MMBIR_FOSTES, NAHR_CANDIDATE, COMPLEX_INS)."""
    consensus_vcf = str(tmp_path / "mech_consensus.vcf.gz")
    sv_vcf = str(tmp_path / "mech_sv.vcf.gz")
    out_vcf = str(tmp_path / "mech_refined.vcf.gz")

    consensus_records = [
        {"chrom": "chr1", "start": 10000, "end": 12000, "svtype": "DEL"},
        {"chrom": "chr1", "start": 30000, "end": 32000, "svtype": "DEL"},
        {"chrom": "chr1", "start": 50000, "end": 52000, "svtype": "DEL"},
        {"chrom": "chr1", "start": 70000, "end": 72000, "svtype": "DEL"},
    ]
    create_consensus_vcf(consensus_vcf, consensus_records)

    sv_records = [
        # Call 1: Blunt join (HOMLEN=0) -> NHEJ
        {"chrom": "chr1", "start": 9950, "end": 12050, "svtype": "DEL", "HOMLEN": 0},
        # Call 2: Microhomology 12 bp -> MMBIR_FOSTES
        {"chrom": "chr1", "start": 29900, "end": 32100, "svtype": "DEL", "HOMSEQ": "ATCGGCTAATCG", "HOMLEN": 12},
        # Call 3: Extended homology 35 bp -> NAHR_CANDIDATE
        {"chrom": "chr1", "start": 49950, "end": 52050, "svtype": "DEL", "HOMSEQ": "A" * 35, "HOMLEN": 35},
        # Call 4: Non-templated insertion 18 bp -> COMPLEX_INS
        {"chrom": "chr1", "start": 69900, "end": 72100, "svtype": "DEL", "INSERTION": "GC" * 9, "INSLEN": 18},
    ]
    create_sv_vcf(sv_vcf, sv_records)

    cmd = [
        PYTHON_EXEC,
        REFINE_SCRIPT,
        "--consensus-vcf", consensus_vcf,
        "--sv-vcf", sv_vcf,
        "--output-vcf", out_vcf,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0

    with pysam.VariantFile(out_vcf) as vf:
        recs = list(vf)
        assert recs[0].info["MECH"] == "NHEJ"
        assert recs[1].info["MECH"] == "MMBIR_FOSTES"
        assert recs[2].info["MECH"] == "NAHR_CANDIDATE"
        assert recs[3].info["MECH"] == "COMPLEX_INS"


def test_cohort_database_ingestion_of_refined_features(tmp_path):
    """Verify that build_cohort_feature_database.py stores refined signatures into SQLite."""
    refined_vcf = str(tmp_path / "cohort_refined.vcf.gz")
    out_dir = str(tmp_path / "cohort_db")

    # Mock refined cohort VCF
    header = pysam.VariantHeader()
    header.contigs.add("chr1", 50000000)
    header.info.add("SVTYPE", 1, "String", "Type")
    header.info.add("END", 1, "Integer", "End")
    header.info.add("SVLEN", 1, "Integer", "Length")
    header.info.add("AF", "A", "Float", "AF")
    header.info.add("PRECISE", 0, "Flag", "Precise")
    header.info.add("SV_REFINED", 1, "String", "Refined Caller")
    header.info.add("SLACK_DIST", 1, "Integer", "Slack")
    header.info.add("SIZE_AWARE_DELTA", 1, "Integer", "Delta")
    header.info.add("REFINED_RECIP_OVERLAP", 1, "Float", "Overlap")
    header.info.add("REFINED_SIZE_SIM", 1, "Float", "Sim")
    header.info.add("HOMSEQ", 1, "String", "Homseq")
    header.info.add("HOMLEN", 1, "Integer", "Homlen")
    header.info.add("INSERTION", 1, "String", "Insertion")
    header.info.add("INSLEN", 1, "Integer", "Inslen")
    header.info.add("MECH", 1, "String", "Mechanism")
    header.formats.add("GT", 1, "String", "Genotype")
    header.samples.add("SAMPLE1")

    with pysam.VariantFile(refined_vcf, "wz", header=header) as vf:
        rec = header.new_record(contig="chr1", start=99800 - 1, stop=101200, alleles=("N", "<DEL>"))
        rec.info["SVTYPE"] = "DEL"
        rec.info["SVLEN"] = -1400
        rec.info["AF"] = (0.5,)
        rec.info["PRECISE"] = True
        rec.info["SV_REFINED"] = "SVABA"
        rec.info["SLACK_DIST"] = 400
        rec.info["SIZE_AWARE_DELTA"] = 1000
        rec.info["REFINED_RECIP_OVERLAP"] = 0.7138
        rec.info["REFINED_SIZE_SIM"] = 0.7143
        rec.info["HOMSEQ"] = "AGCTAG"
        rec.info["HOMLEN"] = 6
        rec.info["MECH"] = "MMBIR_FOSTES"
        rec.samples["SAMPLE1"]["GT"] = (0, 1)
        vf.write(rec)
    pysam.tabix_index(refined_vcf, preset="vcf", force=True)

    cmd = [
        PYTHON_EXEC,
        DB_SCRIPT,
        "--cohort-vcf", refined_vcf,
        "--outdir", out_dir,
        "--output-prefix", "test_refinement_db",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}"

    db_file = os.path.join(out_dir, "test_refinement_db.db")
    assert os.path.exists(db_file)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT variant_id, is_precise, sv_refined_caller, homlen, homseq, mutational_mechanism, refined_recip_overlap
        FROM cohort_variants
    """)
    row = cursor.fetchone()
    assert row is not None
    assert row[0] == "1:99800-101200:DEL"
    assert row[1] == 1
    assert row[2] == "SVABA"
    assert row[3] == 6
    assert row[4] == "AGCTAG"
    assert row[5] == "MMBIR_FOSTES"
    assert abs(row[6] - 0.7138) < 1e-4

    conn.close()


def test_dragen_exome_maxdepth_and_cnv_type_resolution(tmp_path):
    """
    Verify that:
    1. DRAGEN CNVs with INFO/SVTYPE=CNV and ALT=<DEL> or ALT=<DUP> are resolved correctly.
    2. DRAGEN SVs flagged with MaxDepth (common in targeted exome sequencing) are accepted.
    3. Microhomology annotations and NAHR_CANDIDATE / NHEJ mechanisms are properly populated.
    """
    consensus_vcf = str(tmp_path / "dragen_cnv.vcf.gz")
    sv_vcf = str(tmp_path / "dragen_sv.vcf.gz")
    out_vcf = str(tmp_path / "dragen_refined.vcf.gz")
    summary_file = str(tmp_path / "dragen_summary.txt")

    # Manually write raw DRAGEN-style VCF lines to test exact parsing
    raw_consensus = str(tmp_path / "dragen_cnv.vcf")
    with open(raw_consensus, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr19,length=58617616>\n")
        f.write("##contig=<ID=chr20,length=64444167>\n")
        f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
        f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")
        f.write("chr19\t7766015\tDRAGEN:LOSS:chr19:7766015-7766146\tN\t<DEL>\t42\tPASS\tSVLEN=-132;SVTYPE=CNV;END=7766146\tGT\t0/1\n")
        f.write("chr20\t10000000\tDRAGEN:GAIN:chr20:10000000-10005000\tN\t<DUP>\t50\tPASS\tSVLEN=5000;SVTYPE=CNV;END=10005000\tGT\t0/1\n")

    pysam.tabix_compress(raw_consensus, consensus_vcf, force=True)
    pysam.tabix_index(consensus_vcf, preset="vcf", force=True)

    raw_sv = str(tmp_path / "dragen_sv.vcf")
    with open(raw_sv, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr19,length=58617616>\n")
        f.write("##FILTER=<ID=MaxDepth,Description=\"Exome high coverage\">\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
        f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n")
        f.write("##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Homology length\">\n")
        f.write("##INFO=<ID=HOMSEQ,Number=1,Type=String,Description=\"Homology sequence\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")
        f.write("chr19\t7765914\tDRAGEN:DEL:313\tN\t<DEL>\t999\tMaxDepth\tEND=7766052;SVTYPE=DEL;HOMLEN=68;HOMSEQ=GCTGAAGGCTGCAGTGGGTGAGTTGCCAGAGAAATCCAAGCTGCAGGAGATCTACCAGGAGCTGACCC\tGT\t0/1\n")

    pysam.tabix_compress(raw_sv, sv_vcf, force=True)
    pysam.tabix_index(sv_vcf, preset="vcf", force=True)

    cmd = [
        PYTHON_EXEC,
        REFINE_SCRIPT,
        "--consensus-vcf", consensus_vcf,
        "--sv-vcf", sv_vcf,
        "--output-vcf", out_vcf,
        "--summary", summary_file,
        "--engine", "DRAGEN_SV",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Error: {res.stderr}\n{res.stdout}"

    with pysam.VariantFile(out_vcf) as vf:
        recs = list(vf.fetch())
        assert len(recs) == 2, "100% variant retention required"
        rec1, rec2 = recs[0], recs[1]

        # First record refined to base-pair exact DRAGEN SV coordinates
        assert rec1.chrom == "chr19"
        assert rec1.pos == 7765914
        assert rec1.stop == 7766052
        assert rec1.info["PRECISE"] is True
        assert rec1.info["SV_REFINED"] == "DRAGEN_SV"
        assert rec1.info["HOMLEN"] == 68
        assert rec1.info["MECH"] == "NAHR_CANDIDATE"

        # Second record retained as IMPRECISE (no SV match)
        assert rec2.chrom == "chr20"
        assert rec2.pos == 10000000
        assert rec2.stop == 10005000
        assert rec2.info["IMPRECISE"] is True
        assert rec2.info["MECH"] == "UNREFINED"
