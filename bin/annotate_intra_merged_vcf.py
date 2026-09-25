#!/usr/bin/env python3
"""
annotate_intra_merged_vcf.py

Harmonizes and decorates an intra-sample merged CNV VCF (from SURVIVOR or Truvari).
Collapses any multi-caller columns into a single unified sample column named after
the patient/sample ID, and annotates:
  1. FORMAT/TOOLS: Names of callers supporting the CNV event in this sample.
  2. FORMAT/CONCORDANCE: Integer count of agreeing callers.
  3. FORMAT/QUAL_NORM: Maximum normalized quality score (0-1000 scale).
  4. FORMAT/QS_<CALLER>: Per-caller normalized quality scores.
  5. FORMAT/OQ_<CALLER>: Per-caller original raw quality scores (from FORMAT/OQ).
  6. FORMAT/OAS_<CALLER>: Per-caller original algorithm metric names (from FORMAT/OAS).
  7. INFO/TOOLS & INFO/NUM_TOOLS: Summary caller names and count.
"""

import argparse
from collections import defaultdict
import logging
import os
import re
import sys
import numpy as np
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

_CALLER_FILENAME_PATTERNS = [
    ('CANOES',     re.compile(r'CANOES', re.IGNORECASE)),
    ('CLAMMS',     re.compile(r'CLAMMS', re.IGNORECASE)),
    ('XHMM',       re.compile(r'XHMM', re.IGNORECASE)),
    ('GATK_GCNV',  re.compile(r'GCN[V]|GCNV|GATK', re.IGNORECASE)),
    ('CNVKIT',     re.compile(r'CNVKit|CNVKIT', re.IGNORECASE)),
    ('DRAGEN',     re.compile(r'DRAGEN|\.cnv\.vcf', re.IGNORECASE)),
    ('EXOMEDEPTH', re.compile(r'EXOMEDEPTH', re.IGNORECASE)),
]

SUPPORTED_CALLERS = [name for name, _ in _CALLER_FILENAME_PATTERNS]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Decorate intra-merged VCF with FORMAT/TOOLS and per-caller quality subfields."
    )
    parser.add_argument("--input-vcf", required=True, help="Path to input intra-merged VCF (SURVIVOR or Truvari).")
    parser.add_argument("--output-vcf", required=True, help="Path to output decorated VCF (.vcf or .vcf.gz).")
    parser.add_argument("--sample-id", required=True, help="Target sample ID for the unified sample column.")
    parser.add_argument("--tool-vcfs", default="", help="Comma-separated list of caller=path or VCF paths.")
    parser.add_argument("--merger-mode", default="auto", choices=["auto", "survivor", "truvari"], help="Merger used.")
    parser.add_argument("--collapsed-vcf", default=None, help="Optional Truvari collapsed VCF (truvari collapse -c).")
    return parser.parse_args()


def parse_tool_vcfs(tool_vcfs_arg):
    """Parse tool_vcfs argument into a dict of CALLER_NAME -> path."""
    if not tool_vcfs_arg:
        return {}
    res = {}
    tokens = [t.strip() for t in tool_vcfs_arg.split(",") if t.strip()]
    for token in tokens:
        if "=" in token:
            k, v = token.split("=", 1)
            k_upper = k.strip().upper().replace("-", "_")
            res[k_upper] = v.strip()
        else:
            # Infer caller from path
            v = token.strip()
            for caller, pat in _CALLER_FILENAME_PATTERNS:
                if pat.search(v):
                    res[caller] = v
                    break
    return res


def caller_order_from_survivor_header(vcf_header):
    sample_entries = {}
    for rec in vcf_header.records:
        if str(rec).startswith('##SAMPLE='):
            raw = str(rec).strip()
            id_match = re.search(r'ID=(\d+)', raw)
            file_match = re.search(r'File=([^,>]+)', raw)
            if id_match and file_match:
                sample_entries[int(id_match.group(1))] = file_match.group(1)
    if not sample_entries:
        return []
    order = []
    for idx in sorted(sample_entries):
        fpath = sample_entries[idx]
        name = f"TOOL_{idx}"
        for caller, pat in _CALLER_FILENAME_PATTERNS:
            if pat.search(fpath):
                name = caller
                break
        order.append(name)
    return order


def callers_from_supp_vec(supp_vec, caller_order):
    callers = []
    for i, bit in enumerate(str(supp_vec)):
        if bit == '1':
            c_name = caller_order[i] if i < len(caller_order) else f"TOOL_{i}"
            callers.append(c_name)
    return callers


def safe_info_get(info, key, default=None):
    """Safely get an INFO field from a pysam VariantRecord without throwing ValueError on undeclared headers."""
    try:
        if key in info:
            val = info[key]
            return val if val is not None else default
        return default
    except (KeyError, ValueError):
        return default


def build_matchid_caller_map(collapsed_vcf_path):
    if not collapsed_vcf_path or not os.path.exists(collapsed_vcf_path):
        return {}
    matchid_map = defaultdict(set)
    vcf = pysam.VariantFile(collapsed_vcf_path)
    if 'MatchId' not in vcf.header.info:
        vcf.close()
        return {}
    for rec in vcf:
        m_id = safe_info_get(rec.info, 'MatchId')
        if m_id is None:
            continue
        tool_val = safe_info_get(rec.info, 'TOOL', safe_info_get(rec.info, 'TOOLS', ''))
        matched_caller = None
        for caller, pat in _CALLER_FILENAME_PATTERNS:
            if pat.search(str(tool_val)):
                matched_caller = caller
                break
        if matched_caller:
            if isinstance(m_id, (tuple, list)):
                for single_id in m_id:
                    matchid_map[int(float(single_id))].add(matched_caller)
            else:
                matchid_map[int(float(m_id))].add(matched_caller)
    vcf.close()
    return matchid_map


def sanitize_and_clean_vcf(filepath):
    """Remove empty lines and sanitize malformed FORMAT/sample columns (e.g. SURVIVOR haploid genotype bugs)."""
    if not filepath or not os.path.exists(filepath):
        return
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = []
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                lines.append(raw_line)
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) >= 10:
                fmt_keys = fields[8].split(":")
                n_keys = len(fmt_keys)
                for s_idx in range(9, len(fields)):
                    s_vals = fields[s_idx].split(":")
                    if len(s_vals) != n_keys:
                        # Handle SURVIVOR bug where haploid GT ('1:99') is copied as 3 chars ('1:9') splitting GT into two fields
                        if len(s_vals) == n_keys + 1 and s_vals[0] in ("0", "1", ".") and s_vals[1].isdigit():
                            gt = "./." if s_vals[0] == "." else f"{s_vals[0]}/{s_vals[0]}"
                            s_vals = [gt] + s_vals[2:]
                        elif len(s_vals) > n_keys:
                            s_vals = s_vals[:n_keys]
                        elif len(s_vals) < n_keys:
                            s_vals = s_vals + ["."] * (n_keys - len(s_vals))
                        fields[s_idx] = ":".join(s_vals)
                lines.append("\t".join(fields) + "\n")
            else:
                lines.append(raw_line)
    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(lines)


def main():
    args = parse_args()
    tool_vcfs = parse_tool_vcfs(args.tool_vcfs)

    # Sanitize input VCFs against empty lines and malformed FORMAT/sample fields
    sanitize_and_clean_vcf(args.input_vcf)
    if args.collapsed_vcf:
        sanitize_and_clean_vcf(args.collapsed_vcf)

    in_vcf = pysam.VariantFile(args.input_vcf)

    # Determine merger mode
    merger_mode = args.merger_mode
    if merger_mode == "auto":
        if "SUPP_VEC" in in_vcf.header.info:
            merger_mode = "survivor"
        elif "MatchId" in in_vcf.header.info:
            merger_mode = "truvari"
        else:
            merger_mode = "survivor"

    caller_order = caller_order_from_survivor_header(in_vcf.header)
    if not caller_order:
        caller_order = list(tool_vcfs.keys())

    matchid_caller_map = build_matchid_caller_map(args.collapsed_vcf) if merger_mode == "truvari" else {}

    # Open tool VCFs for quality score lookup
    open_tools = {}
    for c_name, c_path in tool_vcfs.items():
        if os.path.exists(c_path):
            try:
                open_tools[c_name] = pysam.VariantFile(c_path)
            except Exception as e:
                logging.warning(f"Could not open tool VCF {c_path}: {e}")

    # Build new header
    new_header = pysam.VariantHeader()
    # Copy contigs
    for c in in_vcf.header.contigs.values():
        new_header.contigs.add(c.name, c.length)

    # Copy existing INFO headers except tools if already present
    for info in in_vcf.header.info.values():
        if info.name not in ["TOOLS", "NUM_TOOLS"]:
            new_header.info.add(info.name, info.number, info.type, info.description)

    # Add INFO headers
    new_header.info.add("TOOLS", ".", "String", "Names of callers supporting this CNV event in this sample")
    new_header.info.add("NUM_TOOLS", 1, "Integer", "Number of callers agreeing on this event in this sample")

    # Add FORMAT headers
    new_header.formats.add("GT", 1, "String", "Genotype")
    new_header.formats.add("TOOLS", ".", "String", "Names of callers supporting this CNV event in this sample")
    new_header.formats.add("CONCORDANCE", 1, "Integer", "Number of callers agreeing on this event in this sample")
    new_header.formats.add("QUAL_NORM", 1, "Float", "Maximum normalized quality score among supporting callers (0-1000 scale)")
    new_header.formats.add("MEAN_QUAL_NORM", 1, "Float", "Mean normalized quality score across supporting callers (0-1000 scale)")
    new_header.formats.add("BEST_TOOL", 1, "String", "Name of the strongest caller by normalized quality score")
    new_header.formats.add("BEST_TOOL_RAW_SCORE", 1, "Float", "Raw pre-normalized quality score of the strongest caller")
    new_header.formats.add("BEST_TOOL_METRIC", 1, "String", "Metric name of the strongest caller's raw score")

    for c_name in SUPPORTED_CALLERS:
        new_header.formats.add(f"QS_{c_name}", 1, "Float", f"{c_name} normalized quality score (0-1000)")
        new_header.formats.add(f"OQ_{c_name}", 1, "Float", f"{c_name} original raw quality score")
        new_header.formats.add(f"OAS_{c_name}", 1, "String", f"{c_name} original algorithm metric name")

    # Single unified sample column
    new_header.samples.add(args.sample_id)

    # Output VCF
    out_mode = "wz" if args.output_vcf.endswith(".gz") else "w"
    tmp_out = args.output_vcf if not args.output_vcf.endswith(".gz") else args.output_vcf.replace(".gz", ".tmp.vcf")
    out_vcf = pysam.VariantFile(tmp_out, "w", header=new_header)

    records_to_write = []
    for rec in in_vcf:
        chrom = rec.chrom
        pos = rec.pos
        end = rec.stop
        info = rec.info

        # 1. Resolve callers
        supporting = []
        if merger_mode == "survivor":
            supp_vec = str(safe_info_get(info, "SUPP_VEC", "1"))
            supporting = callers_from_supp_vec(supp_vec, caller_order)
        elif merger_mode == "truvari":
            # Check representative TOOL
            tool_val = safe_info_get(info, "TOOL", safe_info_get(info, "TOOLS", ""))
            for caller, pat in _CALLER_FILENAME_PATTERNS:
                if pat.search(str(tool_val)):
                    supporting.append(caller)
                    break
            # Check collapsed MatchId
            m_id = safe_info_get(info, "MatchId")
            if m_id is not None:
                if isinstance(m_id, (tuple, list)):
                    for single_id in m_id:
                        supporting.extend(matchid_caller_map.get(int(float(single_id)), []))
                else:
                    supporting.extend(matchid_caller_map.get(int(float(m_id)), []))

        # Deduplicate while preserving order
        seen = set()
        dedup_callers = []
        for c in supporting:
            if c not in seen:
                seen.add(c)
                dedup_callers.append(c)

        if not dedup_callers:
            dedup_callers = ["UNKNOWN"]

        # 2. Extract per-caller quality scores
        qual_norm_vals = {}
        oq_vals = {}
        oas_vals = {}

        for c_name in dedup_callers:
            if c_name in open_tools:
                t_vcf = open_tools[c_name]
                if chrom in t_vcf.header.contigs:
                    try:
                        matches = list(t_vcf.fetch(chrom, max(0, pos - 20), pos + 20))
                        if matches:
                            m_rec = matches[0]
                            if m_rec.qual is not None:
                                qual_norm_vals[c_name] = float(m_rec.qual)
                            # Check FORMAT in matched record
                            for s_sample in m_rec.samples.values():
                                if "OQ" in s_sample and s_sample["OQ"] is not None:
                                    oq_vals[c_name] = float(s_sample["OQ"])
                                if "OAS" in s_sample and s_sample["OAS"] is not None:
                                    oas_vals[c_name] = str(s_sample["OAS"])
                    except Exception:
                        pass

        if qual_norm_vals:
            best_tool = max(qual_norm_vals.keys(), key=lambda k: qual_norm_vals[k])
            max_qual_norm = qual_norm_vals[best_tool]
            mean_qual_norm = sum(qual_norm_vals.values()) / len(qual_norm_vals)
            best_tool_raw = oq_vals.get(best_tool)
            best_tool_metric = oas_vals.get(best_tool, "")
        else:
            best_tool = dedup_callers[0] if dedup_callers else "UNKNOWN"
            max_qual_norm = float(rec.qual) if rec.qual is not None else 100.0
            mean_qual_norm = max_qual_norm
            best_tool_raw = None
            best_tool_metric = ""

        # 3. Create new record
        new_rec = new_header.new_record()
        new_rec.chrom = chrom
        new_rec.pos = pos
        new_rec.id = rec.id
        new_rec.ref = rec.ref
        new_rec.alts = rec.alts
        new_rec.qual = max_qual_norm
        new_rec.filter.add("PASS")

        # Copy existing info fields
        for k, v in info.items():
            if k in new_header.info and k not in ["TOOLS", "NUM_TOOLS"]:
                try:
                    new_rec.info[k] = v
                except Exception:
                    pass

        # Set INFO fields
        new_rec.info["TOOLS"] = ",".join(dedup_callers)
        new_rec.info["NUM_TOOLS"] = len(dedup_callers)

        # Set single SAMPLE FORMAT fields
        sample = new_rec.samples[args.sample_id]
        sample["GT"] = (0, 1)
        sample["TOOLS"] = ",".join(dedup_callers)
        sample["CONCORDANCE"] = len(dedup_callers)
        sample["QUAL_NORM"] = round(max_qual_norm, 2)
        sample["MEAN_QUAL_NORM"] = round(mean_qual_norm, 2)
        sample["BEST_TOOL"] = best_tool
        if best_tool_raw is not None:
            sample["BEST_TOOL_RAW_SCORE"] = round(float(best_tool_raw), 2)
        if best_tool_metric:
            sample["BEST_TOOL_METRIC"] = best_tool_metric

        for c_name in SUPPORTED_CALLERS:
            if c_name in qual_norm_vals:
                sample[f"QS_{c_name}"] = round(qual_norm_vals[c_name], 2)
            if c_name in oq_vals:
                sample[f"OQ_{c_name}"] = round(oq_vals[c_name], 2)
            if c_name in oas_vals:
                sample[f"OAS_{c_name}"] = oas_vals[c_name]

        records_to_write.append(new_rec)

    in_vcf.close()
    for t_vcf in open_tools.values():
        t_vcf.close()

    # Sort records by contig order defined in header, then start position
    contig_order = {name: i for i, name in enumerate(new_header.contigs)}
    records_to_write.sort(key=lambda r: (contig_order.get(r.chrom, 999999), r.pos, r.stop))

    for r in records_to_write:
        out_vcf.write(r)
    out_vcf.close()

    # Compress and index if target is .gz
    if args.output_vcf.endswith(".gz"):
        pysam.tabix_compress(tmp_out, args.output_vcf, force=True)
        pysam.tabix_index(args.output_vcf, preset="vcf", force=True)
        if os.path.exists(tmp_out):
            os.remove(tmp_out)

    logging.info(f"Wrote decorated intra-merged VCF with FORMAT/TOOLS to {args.output_vcf}")


if __name__ == "__main__":
    main()
