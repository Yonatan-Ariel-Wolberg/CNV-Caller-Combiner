# =============================================================================
# Detailed documentation: evaluate_random_forest_performance.py
# Purpose: Evaluate the classification performance of a Random Forest model
#   or caller predictions at the probe level using a capture-target probes BED,
#   a truth-set BED, and a call-set BED.
# Inputs:
#   - truth_bed (str): Path to truth set BED file. Must contain at least 4 tab-separated
#                      columns: CHR, START, STOP, SAMPLE
#   - callset_bed (str): Path to call set BED file. Must contain at least 4 tab-separated
#                      columns: CHR, START, STOP, SAMPLE
#   - probe_bed (str): Path to capture target probes BED file. Must contain at least 3
#                      tab-separated columns: CHR, START, STOP
# Outputs:
#   - Prints the Confusion Matrix, F1 Score (F_beta with beta=1), MCC, Sensitivity,
#     and Specificity to stdout.
# Operational notes:
#   - Leverages scikit-learn's confusion_matrix to generate confusion metrics.
#   - Categorizes probe coordinates using overlapping criteria.
# =============================================================================
import pandas as pd
from sklearn.metrics import confusion_matrix

def load_bed_file(bed_file):
    """
    Load a capture-target BED file containing chromosome coordinates.
    
    Args:
        bed_file (str): Path to the BED file.
        
    Returns:
        pd.DataFrame: DataFrame with columns 'chr', 'start', and 'end'.
    """
    return pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end'])

def load_cnv_file(cnv_file):
    """
    Load a 4-column CNV BED file with sample information.
    
    Args:
        cnv_file (str): Path to the CNV BED file.
        
    Returns:
        pd.DataFrame: DataFrame with columns 'chr', 'start', 'end', and 'sample'.
    """
    return pd.read_csv(cnv_file, sep='\t', header=None, names=['chr', 'start', 'end', 'sample'])

def find_overlaps(probes, cnvs):
    """
    Find genomic overlaps between capture targets (probes) and CNVs.
    
    Args:
        probes (pd.DataFrame): Dataframe of probe coordinates.
        cnvs (pd.DataFrame): Dataframe of CNV coordinates.
        
    Returns:
        list: List of tuples of overlapping (probe, cnv) rows.
    """
    overlaps = []
    for index, probe in probes.iterrows():
        for _, cnv in cnvs.iterrows():
            if (probe['chr'] == cnv['chr'] and 
                probe['start'] < cnv['end'] and 
                probe['end'] > cnv['start']):
                overlaps.append((probe, cnv))
    return overlaps

def categorize_probes(probes, truth_cnv, callset_cnv):
    """
    Categorizes each probe into TP, FN, FP, TN based on truth-set and callset overlaps.
    
    Args:
        probes (pd.DataFrame): Probes coordinates dataframe.
        truth_cnv (pd.DataFrame): Truth-set CNV dataframe.
        callset_cnv (pd.DataFrame): Callset CNV dataframe.
        
    Returns:
        dict: A dictionary mapping 'TP', 'FN', 'FP', 'TN' to lists of probe tuples.
    """
    categorized = {
        'TP': [],
        'FN': [],
        'FP': [],
        'TN': []
    }
    
    # Find overlaps with truth set
    truth_overlaps = find_overlaps(probes, truth_cnv)
    for probe in probes.itertuples(index=False):
        is_in_truth = any((probe.chr == overlap[0]['chr'] and 
                           probe.start < overlap[0]['end'] and 
                           probe.end > overlap[0]['start']) for overlap in truth_overlaps)
        
        is_called = any((probe.chr == cnv['chr'] and 
                         probe.start < cnv['end'] and 
                         probe.end > cnv['start']) for cnv in callset_cnv.itertuples(index=False))
        
        if is_in_truth and is_called:
            categorized['TP'].append(probe)
        elif is_in_truth and not is_called:
            categorized['FN'].append(probe)
        elif not is_in_truth and is_called:
            categorized['FP'].append(probe)
        else:
            categorized['TN'].append(probe)
    
    return categorized

def compute_metrics(categorized):
    """
    Compute standard performance metrics (Confusion Matrix, F_beta, MCC, Sensitivity, Specificity)
    from categorized probe counts.
    
    Args:
        categorized (dict): Categorized probes dictionary containing TP, FN, FP, TN lists.
        
    Returns:
        tuple: (confusion_matrix_ndarray, F_beta, MCC, sensitivity, specificity)
    """
    TP = len(categorized['TP'])
    FN = len(categorized['FN'])
    FP = len(categorized['FP'])
    TN = len(categorized['TN'])

    # Confusion Matrix
    confusion = confusion_matrix(
        [1] * (TP + FN) + [0] * (FP + TN),
        [1] * TP + [0] * FN + [1] * FP + [0] * TN
    )

    # Sensitivity (Recall)
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0

    # Precision
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0

    # Specificity
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0

    # F_beta Calculation
    beta = 1  # F1 Score, can be adjusted for other beta values
    if precision + sensitivity > 0:
        F_beta = (1 + beta**2) * (precision * sensitivity) / ((beta**2 * precision) + sensitivity)
    else:
        F_beta = 0

    # Matthews Correlation Coefficient (MCC)
    MCC = (TP * TN - FP * FN) / (((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))**0.5) if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0 else 1

    return confusion, F_beta, MCC, sensitivity, specificity

def main(truth_bed, callset_bed, probe_bed):
    """
    Run evaluation comparing predictions against the truth set.
    
    Args:
        truth_bed (str): Path to truth set BED file.
        callset_bed (str): Path to call set BED file.
        probe_bed (str): Path to probes BED file.
    """
    probes = load_bed_file(probe_bed)
    truth_cnv = load_cnv_file(truth_bed)
    callset_cnv = load_cnv_file(callset_bed)

    categorized = categorize_probes(probes, truth_cnv, callset_cnv)
    confusion, F_beta, MCC, sensitivity, specificity = compute_metrics(categorized)

    print(f"Confusion Matrix:\n{confusion}")
    print(f"F_Beta: {F_beta}")
    print(f"MCC: {MCC}")
    print(f"Sensitivity: {sensitivity}")
    print(f"Specificity: {specificity}")

# Example usage
# main('truthset_cnv.bed', 'callset_cnv.bed', 'probes_sanger.bed')
