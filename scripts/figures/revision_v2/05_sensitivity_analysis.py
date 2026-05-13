#!/usr/bin/env python3
"""
Script 05: Sensitivity Analysis
Stratify by context, leave-one-subtype-out, outlier detection
"""

import pandas as pd
import numpy as np
import json
import statsmodels.api as sm
from pathlib import Path
from scipy import stats

BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"

def run_context_stratification(df):
    """Stratify analysis by context tier"""
    results = {}
    for context in df['context_tier'].dropna().unique():
        sub = df[df['context_tier'] == context].copy()
        if len(sub) >= 3:
            results[context] = {
                'n': len(sub),
                'mean_log10fc': float(sub['log10_FC'].mean()),
                'std_log10fc': float(sub['log10_FC'].std()),
                'mutations': sub['Mutation'].tolist()
            }
    return results

def run_leave_one_subtype_out(df):
    """Leave-one-subtype-out analysis"""
    results = {}
    subtypes = df['subtype'].dropna().unique()

    for subtype in subtypes:
        train = df[df['subtype'] != subtype].copy()
        if len(train) >= 5:
            train_mean = float(train['log10_FC'].mean())
            full_mean = float(df['log10_FC'].mean())
            results[subtype] = {
                'n_removed': int((df['subtype'] == subtype).sum()),
                'mean_with': full_mean,
                'mean_without': train_mean,
                'delta': float(abs(train_mean - full_mean))
            }
    return results

def detect_outliers(df):
    """Outlier detection using 3-sigma rule"""
    mean = df['log10_FC'].mean()
    std = df['log10_FC'].std()
    threshold = 3
    outliers = df[np.abs(df['log10_FC'] - mean) > threshold * std]
    return {
        'threshold_sigma': threshold,
        'n_outliers': len(outliers),
        'outlier_mutations': outliers['Mutation'].tolist() if len(outliers) > 0 else []
    }

def main():
    print("="*60)
    print("Script 05: Sensitivity Analysis")
    print("="*60)

    df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
    df = df[df['log10_FC'].notna()].copy()

    results = {
        'context_stratification': run_context_stratification(df),
        'leave_one_subtype_out': run_leave_one_subtype_out(df),
        'outlier_detection': detect_outliers(df),
        'summary': {
            'n_total': len(df),
            'mean_variation_leave_one_subtype': None
        }
    }

    # Mean variation from LOSO subtype
    if results['leave_one_subtype_out']:
        deltas = [v['delta'] for v in results['leave_one_subtype_out'].values()]
        results['summary']['mean_variation_leave_one_subtype'] = float(np.mean(deltas))
        print(f"Mean variation from leave-one-subtype-out: {np.mean(deltas):.4f} log10FC units")

    out_path = RESULTS_DIR / "sensitivity_results.json"
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n[OK] Saved: {out_path}")
    print("Script 05 complete!")

if __name__ == "__main__":
    main()
