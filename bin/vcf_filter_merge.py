#!/usr/bin/env python3

# =============================================================================
# Detailed documentation: vcf_filter_merge.py
# Purpose: This pipeline helper script supports CNV calling and/or post-processing tasks in
#   this project, with parameters and I/O conventions aligned to repository workflows.
# Scope: General utility execution.
# Inputs: CLI arguments, repository config/params files, and sample-level CNV data.
# Outputs: Generated intermediate artifacts, final result files, and/or run metadata.
# Operational notes: Keep parameter names synchronized with related workflow/process
#   definitions and update this file alongside upstream interface changes.
# =============================================================================
import sys
import argparse
import os

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            info_dict[k] = v
        else:
            info_dict[item] = True
    return info_dict

def build_info(info_dict):
    parts = []
    for k, v in info_dict.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f"{k}={v}")
    return ';'.join(parts)

def get_quality_score(var):
    format_keys = var['FORMAT'].split(':')
    format_vals = var['SAMPLE_VALS'].split(':')
    for k, v in zip(format_keys, format_vals):
        if k in ['Q_SOME', 'GQ']:
            try:
                return float(v)
            except ValueError:
                return 0.0
    try:
        return float(var['QUAL']) if var['QUAL'] != '.' else 0.0
    except ValueError:
        return 0.0

def filter_record(var, caller):
    format_keys = var['FORMAT'].split(':')
    format_vals = var['SAMPLE_VALS'].split(':')
    format_dict = dict(zip(format_keys, format_vals))
    
    gt = format_dict.get('GT', './.')
    if gt in ['0', '0/0', '0|0']:
        return False
        
    if caller.upper() == 'CLAMMS':
        try:
            q_some = float(format_dict.get('Q_SOME', 0))
            q_exact = float(format_dict.get('Q_EXACT', 0))
            return q_some > 500 and q_exact > 0
        except ValueError:
            return False
        
    elif caller.upper() == 'CANOES':
        try:
            q_some = float(format_dict.get('Q_SOME', 0))
            return q_some > 80
        except ValueError:
            return False
        
    return True

def main():
    parser = argparse.ArgumentParser(description="Filter and merge adjacent VCF variants within a sample.")
    parser.add_argument('--input_vcf', required=True, help="Input VCF file path")
    parser.add_argument('--output_vcf', required=True, help="Output VCF file path")
    parser.add_argument('--caller', required=True, choices=['CLAMMS', 'CANOES'], help="Caller name for filtering logic")
    parser.add_argument('--distance', type=int, default=1500, help="Max distance to merge adjacent variants")
    parser.add_argument('--new_sample', help="New sample ID to write in the VCF header")
    args = parser.parse_args()
    
    header_lines = []
    records = []
    
    with open(args.input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    if args.new_sample and len(parts) >= 10:
                        parts[9] = args.new_sample
                        line = '\t'.join(parts) + '\n'
                header_lines.append(line)
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
                
            chrom = parts[0]
            pos = int(parts[1])
            vid = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filt = parts[6]
            info_str = parts[7]
            fmt = parts[8]
            sample_vals = parts[9]
            
            info_dict = parse_info(info_str)
            end_val = info_dict.get('END')
            if end_val is not None:
                end = int(end_val)
            else:
                end = pos + len(ref) - 1
                
            format_keys = fmt.split(':')
            format_vals = sample_vals.split(':')
            format_dict = dict(zip(format_keys, format_vals))
            gt_val = format_dict.get('GT', '0')
            
            svtype = info_dict.get('SVTYPE', 'CNV')
            if svtype == 'CNV':
                if gt_val == '1':
                    svtype = 'DEL'
                elif gt_val == '2':
                    svtype = 'DUP'
                elif alt == '<DEL>':
                    svtype = 'DEL'
                elif alt == '<DUP>':
                    svtype = 'DUP'
                
            info_dict['SVTYPE'] = svtype
            info_dict['END'] = end
            
            var = {
                'CHROM': chrom,
                'POS': pos,
                'END': end,
                'ID': vid,
                'REF': ref,
                'ALT': alt,
                'QUAL': qual,
                'FILTER': filt,
                'INFO_DICT': info_dict,
                'FORMAT': fmt,
                'SAMPLE_VALS': sample_vals,
                'SVTYPE': svtype
            }
            
            if filter_record(var, args.caller):
                records.append(var)
                
    groups = {}
    for var in records:
        key = (var['CHROM'], var['SVTYPE'])
        if key not in groups:
            groups[key] = []
        groups[key].append(var)
        
    merged_records = []
    
    for key, var_list in groups.items():
        var_list.sort(key=lambda x: x['POS'])
        
        current = None
        for var in var_list:
            if current is None:
                current = var
            else:
                gap = var['POS'] - current['END']
                if gap <= args.distance:
                    current['END'] = max(current['END'], var['END'])
                    current['INFO_DICT']['END'] = current['END']
                    
                    q_curr = get_quality_score(current)
                    q_var = get_quality_score(var)
                    if q_var > q_curr:
                        current['FORMAT'] = var['FORMAT']
                        current['SAMPLE_VALS'] = var['SAMPLE_VALS']
                        current['QUAL'] = var['QUAL']
                else:
                    merged_records.append(current)
                    current = var
                    
        if current is not None:
            merged_records.append(current)
            
    def get_sort_key(var):
        chr_name = var['CHROM'].replace('chr', '')
        try:
            return (0, int(chr_name), var['POS'])
        except ValueError:
            return (1, chr_name, var['POS'])
            
    merged_records.sort(key=get_sort_key)
    
    with open(args.output_vcf, 'w') as out:
        out.writelines(header_lines)
        for var in merged_records:
            svlen = var['END'] - var['POS'] + 1
            if var['SVTYPE'] == 'DEL':
                var['INFO_DICT']['SVLEN'] = -svlen
            else:
                var['INFO_DICT']['SVLEN'] = svlen
                
            info_str = build_info(var['INFO_DICT'])
            
            row = [
                var['CHROM'],
                str(var['POS']),
                var['ID'],
                var['REF'],
                var['ALT'],
                var['QUAL'],
                var['FILTER'],
                info_str,
                var['FORMAT'],
                var['SAMPLE_VALS']
            ]
            out.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    main()
