#!/usr/bin/env python3
"""
validate_run.py

Pre-flight diagnostic validator for CNV-Caller-Combiner.
Ensures that all input files, reference assets, caller resources,
and cluster software requirements are in place before launching a run.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"
    BOLD = "\033[1m"
    RESET = "\033[0m"


def check_file(path_str: str, desc: str, required: bool = True) -> tuple[bool, str]:
    if not path_str:
        if required:
            return False, f"Missing configuration for {desc}"
        return True, "Optional (not configured)"
    p = Path(path_str).resolve()
    if not p.exists():
        return False, f"File does not exist: {p}"
    if not p.is_file():
        return False, f"Path is not a regular file: {p}"
    if not os.access(p, os.R_OK):
        return False, f"File is not readable (permission denied): {p}"
    size_mb = p.stat().st_size / (1024 * 1024)
    return True, f"Found ({size_mb:.1f} MB)"


def check_dir(path_str: str, desc: str) -> tuple[bool, str]:
    if not path_str:
        return False, f"Missing configuration for {desc}"
    p = Path(path_str).resolve()
    if not p.exists():
        return False, f"Directory does not exist: {p}"
    if not p.is_dir():
        return False, f"Path is not a directory: {p}"
    return True, f"Found ({p})"


def check_command(cmd: str) -> tuple[bool, str]:
    loc = shutil.which(cmd)
    if not loc:
        return False, f"'{cmd}' not found in PATH"
    try:
        if cmd == "java":
            res = subprocess.run(["java", "-version"], capture_output=True, text=True, timeout=5)
            version = res.stderr.splitlines()[0] if res.stderr else "installed"
            return True, f"{version}"
        elif cmd == "nextflow":
            res = subprocess.run(["nextflow", "-v"], capture_output=True, text=True, timeout=5)
            version = res.stdout.strip()
            return True, f"{version}"
        elif cmd in ["apptainer", "singularity"]:
            res = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=5)
            version = res.stdout.strip()
            return True, f"{version}"
    except Exception:
        pass
    return True, f"Available ({loc})"


def main():
    parser = argparse.ArgumentParser(
        description="Pre-flight diagnostic validator for CNV-Caller-Combiner configuration."
    )
    parser.add_argument(
        "--params", "-p",
        required=True,
        help="Path to the JSON parameters file (e.g. params-local-2-test-balanced.json)."
    )
    args = parser.parse_args()

    params_path = Path(args.params).resolve()
    if not params_path.is_file():
        sys.stderr.write(f"❌ Error: Parameters file does not exist: {params_path}\n")
        sys.exit(1)

    try:
        with open(params_path) as f:
            params = json.load(f)
    except Exception as e:
        sys.stderr.write(f"❌ Error: Failed to parse JSON configuration: {e}\n")
        sys.exit(1)

    print("=" * 72)
    print(f" 🛡️  CNV-CALLER-COMBINER PRE-FLIGHT DIAGNOSTIC VALIDATOR")
    print(f" Config: {params_path}")
    print("=" * 72)

    checks = []

    # 1. System Software
    for cmd in ["java", "nextflow", "apptainer"]:
        alt = "singularity" if cmd == "apptainer" else None
        ok, msg = check_command(cmd)
        if not ok and alt:
            ok, msg = check_command(alt)
            cmd = alt
        checks.append((f"Software: {cmd}", ok, msg, True))

    # 2. Samplesheet
    ss_bams = params.get("samplesheet_bams", "")
    ss_crams = params.get("samplesheet_crams", "")

    if ss_bams and Path(ss_bams).is_file():
        ss_path = ss_bams
        ss_type = "BAM"
    elif ss_crams and Path(ss_crams).is_file():
        ss_path = ss_crams
        ss_type = "CRAM"
    elif ss_bams and ss_crams and not Path(ss_bams).exists() and Path(ss_crams).is_file():
        ss_path = ss_crams
        ss_type = "CRAM"
    else:
        ss_path = ss_bams or ss_crams
        ss_type = "Unknown"

    ok, msg = check_file(ss_path, f"Samplesheet")
    if ok:
        try:
            with open(ss_path) as f:
                header = f.readline().strip().split("\t")
                rows = [line.strip().split("\t") for line in f if line.strip()]
            sample_count = len(rows)
            msg = f"{sample_count} samples declared in {Path(ss_path).name} ({ss_type} input)"
            missing_align = 0
            for r in rows:
                if len(r) > 1 and not Path(r[1]).is_file():
                    missing_align += 1
            if missing_align > 0:
                ok = False
                msg += f" (⚠️ {missing_align} alignment files missing on disk!)"
        except Exception as e:
            ok = False
            msg = f"Failed to read samplesheet: {e}"
    elif ss_bams and ss_crams and not Path(ss_bams).exists() and Path(ss_crams).is_file():
        ok = True
        msg = f"CRAM samplesheet exists; downstream BAM samplesheet will be generated."
    checks.append((f"Cohort Samplesheet", ok, msg, True))

    # 3. Reference Genome & Dictionaries
    ref_fa = params.get("ref", params.get("reference_fasta", params.get("fasta", "")))
    ok, msg = check_file(ref_fa, "Reference FASTA")
    checks.append(("Reference FASTA", ok, msg, True))

    fai = params.get("fai", f"{ref_fa}.fai" if ref_fa else "")
    ok, msg = check_file(fai, "Reference FASTA Index (.fai)")
    checks.append(("Reference Index (.fai)", ok, msg, True))

    dict_file = params.get("dict", f"{ref_fa.rsplit('.', 1)[0]}.dict" if ref_fa else "")
    ok, msg = check_file(dict_file, "Reference Sequence Dict (.dict)", required=False)
    checks.append(("Sequence Dict (.dict)", ok, msg, False))

    # 4. Exome Targets & Intervals
    targets = params.get("targets", params.get("probes", ""))
    ok, msg = check_file(targets, "Exome Targets / Probes BED")
    checks.append(("Exome Targets BED", ok, msg, True))

    intervals = params.get("interval_list", params.get("exome_targets", ""))
    ok, msg = check_file(intervals, "Exome Interval List", required=False)
    checks.append(("Interval List", ok, msg, False))

    # 5. Caller-specific resources
    mappability = params.get("mappability", "")
    ok, msg = check_file(mappability, "CLAMMS Mappability Track", required=False)
    checks.append(("CLAMMS Mappability", ok, msg, False))

    special_reg = params.get("special_reg", "")
    ok, msg = check_file(special_reg, "CLAMMS Special Regions", required=False)
    checks.append(("CLAMMS Special Regions", ok, msg, False))

    xhmm_conf = params.get("xhmm_conf", "")
    ok, msg = check_file(xhmm_conf, "XHMM Configuration", required=False)
    checks.append(("XHMM Config", ok, msg, False))

    # 6. AnnotSV Database
    annotsv_db = params.get("annotsv_db", "")
    ok, msg = check_dir(annotsv_db, "AnnotSV Knowledgebase Directory")
    checks.append(("AnnotSV Database", ok, msg, True))

    # 7. Output directory & storage
    outdir = params.get("outdir", "")
    if outdir:
        out_p = Path(outdir).resolve()
        try:
            out_p.mkdir(parents=True, exist_ok=True)
            test_file = out_p / ".write_test"
            test_file.touch()
            test_file.unlink()
            stat = shutil.disk_usage(out_p)
            free_gb = stat.free / (1024 * 1024 * 1024)
            ok = free_gb > 10
            msg = f"Writable, {free_gb:.1f} GB free disk space"
        except Exception as e:
            ok = False
            msg = f"Cannot write to output dir: {e}"
    else:
        ok = False
        msg = "Missing 'outdir' parameter"
    checks.append(("Output Storage", ok, msg, True))

    # Print Report
    all_passed = True
    critical_failed = False

    for name, passed, info, critical in checks:
        if passed:
            icon = f"{Colors.GREEN}[✓]{Colors.RESET}"
        elif critical:
            icon = f"{Colors.RED}[✗]{Colors.RESET}"
            all_passed = False
            critical_failed = True
        else:
            icon = f"{Colors.YELLOW}[-]{Colors.RESET}"
            all_passed = False

        print(f" {icon} {name:<26} : {info}")

    print("-" * 72)
    if critical_failed:
        print(f"❌ {Colors.RED}{Colors.BOLD}PRE-FLIGHT AUDIT FAILED!{Colors.RESET}")
        print("Please resolve the red [✗] critical requirements before starting the run.\n")
        sys.exit(1)
    elif not all_passed:
        print(f"⚠️  {Colors.YELLOW}{Colors.BOLD}PRE-FLIGHT PASSED WITH NOTICES.{Colors.RESET}")
        print("Optional assets are missing, but the core pipeline is ready to execute.\n")
    else:
        print(f"🎉 {Colors.GREEN}{Colors.BOLD}ALL PRE-FLIGHT CHECKS PASSED! READY TO EXECUTE.{Colors.RESET}\n")

if __name__ == "__main__":
    main()
