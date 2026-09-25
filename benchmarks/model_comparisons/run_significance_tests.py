import numpy as np
from scipy import stats

def run_fisher_exact(tp1, fp1, tp2, fp2, label1, label2):
    # Contingency table:
    #             Success (TP)   Failure (FP)
    # Workflow 1      tp1           fp1
    # Workflow 2      tp2           fp2
    table = np.array([[tp1, fp1], [tp2, fp2]])
    oddsratio, p_value = stats.fisher_exact(table)
    print(f"Fisher's Exact Test: {label1} vs {label2}")
    print(f"  Contingency Table: {table.tolist()}")
    print(f"  Odds Ratio: {oddsratio:.4f}")
    print(f"  p-value: {p_value:e}")
    print(f"  Significant (alpha=0.05): {p_value < 0.05}\n")
    return p_value

def run_mcnemar_binomial(b, c, label1, label2):
    # When b+c is small (like with n=59 truthset variants), the binomial test is exact.
    # We test if the discordant calls (b vs c) deviate from a 50/50 distribution.
    n = b + c
    k = min(b, c)
    # Two-sided binomial test p-value:
    p_value = stats.binom.cdf(k, n, 0.5) * 2
    # Cap at 1.0
    p_value = min(1.0, p_value)
    print(f"McNemar's Exact Binomial Test: {label1} vs {label2}")
    # b: detected only by label1, c: detected only by label2
    print(f"  Discordant: {label1}_only={b}, {label2}_only={c}")
    print(f"  p-value: {p_value:e}")
    print(f"  Significant (alpha=0.05): {p_value < 0.05}\n")
    return p_value

print("=== 1. MERGER COMPARISON: TRUVARI VS SURVIVOR (Stage 3 Consensus) ===")
# Truvari: TP = 4404, FP = 8959 -> Precision = 32.96%
# SURVIVOR: TP = 4404, FP = 8959 -> (At stage 3 Truvari has 21.5% precision, SURVIVOR has 4.2% precision)
# Let's use the actual TP/FP counts at the consensus stage (Stage 3):
# Truvari: TP=52, FP=190
# SURVIVOR: TP=55, FP=1254
run_fisher_exact(52, 190, 55, 1254, "Truvari Consensus", "SURVIVOR Consensus")

# McNemar's on Sensitivity (Truvari missed 7 variants that SURVIVOR detected; SURVIVOR missed 4 that Truvari detected)
run_mcnemar_binomial(4, 7, "Truvari Sensitivity", "SURVIVOR Sensitivity")


print("=== 2. TARGET COORDINATE BED COMPARISON: MOSDEPTH VS CONTROL ===")
# Local Consensus on Mosdepth: TP = 52, FP = 68 (Precision = 43.3%)
# Local Consensus on Control: TP = 52, FP = 190 (Precision = 21.5%)
run_fisher_exact(52, 68, 52, 190, "Mosdepth Target BED", "Control Target BED")

# McNemar's on Sensitivity (Mosdepth missed 2 variants that Control detected, Control missed 2 that Mosdepth detected)
run_mcnemar_binomial(2, 2, "Mosdepth Sensitivity", "Control Sensitivity")

