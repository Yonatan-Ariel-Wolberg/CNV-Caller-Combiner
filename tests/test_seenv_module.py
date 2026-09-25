#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

class TestSeeNVModule(unittest.TestCase):
    def test_seenv_cli_help(self):
        script_path = os.path.abspath("bin/render_cnv_infographics.py")
        res = subprocess.run(["python3", script_path, "--help"], capture_output=True, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("Render SeeNV-style CNV infographics", res.stdout)

    def test_seenv_infographic_generation(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            # Create mock VCF with ACMG 5 and ACMG 4 calls
            vcf_content = (
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">\n"
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
                "##INFO=<ID=CALLERS,Number=1,Type=String,Description=\"Callers\">\n"
                "##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description=\"Number of callers\">\n"
                "##INFO=<ID=ACMG_CLASS,Number=1,Type=String,Description=\"ACMG classification\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t1000000\tCNV_1\tN\t<DEL>\t99\tPASS\tEND=1050000;SVTYPE=DEL;CALLERS=CANOES,CLAMMS,XHMM;NUM_CALLERS=3;ACMG_CLASS=5\n"
                "chr2\t2000000\tCNV_2\tN\t<DUP>\t99\tPASS\tEND=2080000;SVTYPE=DUP;CALLERS=GCNV,CNVKIT;NUM_CALLERS=2;ACMG_CLASS=4\n"
            )
            vcf_file = tmp_path / "test.vcf"
            vcf_file.write_text(vcf_content)

            # Create mock BAM file using samtools
            sam_content = (
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:50000000\n"
                "@SQ\tSN:chr2\tLN:50000000\n"
                "@RG\tID:TEST\tSM:TEST_SAMPLE\n"
                "read1\t99\tchr1\t1000100\t60\t75M\t=\t1000300\t275\t" + "A"*75 + "\t" + "F"*75 + "\tRG:Z:TEST\n"
                "read1\t147\tchr1\t1000300\t60\t75M\t=\t1000100\t-275\t" + "A"*75 + "\t" + "F"*75 + "\tRG:Z:TEST\n"
            )
            sam_file = tmp_path / "test.sam"
            sam_file.write_text(sam_content)
            bam_file = tmp_path / "test.bam"
            if shutil.which("samtools"):
                subprocess.run(["samtools", "view", "-bS", str(sam_file), "-o", str(bam_file)], check=True)
                subprocess.run(["samtools", "index", str(bam_file)], check=True)
            else:
                import pysam
                pysam.view("-bS", str(sam_file), "-o", str(bam_file), catch_stdout=False)
                pysam.index(str(bam_file))

            out_dir = tmp_path / "infographics_out"
            script_path = os.path.abspath("bin/render_cnv_infographics.py")
            
            res = subprocess.run([
                sys.executable, script_path,
                "--vcf", str(vcf_file),
                "--bam", str(bam_file),
                "--sample_id", "TEST_SAMPLE",
                "--outdir", str(out_dir),
                "--max_calls", "5"
            ], capture_output=True, text=True)

            self.assertEqual(res.returncode, 0)
            self.assertTrue((out_dir / "infographics_summary.tsv").exists())
            
            summary_lines = (out_dir / "infographics_summary.tsv").read_text().strip().split("\n")
            self.assertEqual(len(summary_lines), 3) # Header + 2 calls
            
            pngs = list(out_dir.glob("*.png"))
            self.assertEqual(len(pngs), 2)

if __name__ == "__main__":
    unittest.main()
