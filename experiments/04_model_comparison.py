#!/usr/bin/env python3
"""
Model Comparison Framework - Replace variance=0 interpretation
Implements M0-M3 models with AIC/BIC comparison and bootstrap stability
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
import statsmodels.api as sm
from statsmodels.regression.mixed_linear_model import MixedLM
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load harmonized phenotype data"""
    df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")

    # Filter to complete observations
    df_complete = df[df['log10_FC'].notna()].copy()

    print(f"Loaded {len(df_complete)} complete observations")
    print(f"  Mutations: {df_complete['Mutation'].nunique()}")
    print(f"  Subtypes: {df_complete['Subtype'].nunique()}")
    print(f"  Studies: {df_complete['study_source'].nunique()}")

    return df_complete

def fit_model_M0(df):
    """M0: Intercept only (null model)"""
    y = df['log10_FC'].values
    X = np.ones(len(y))

    model = sm.OLS(y, X).fit()

    return {
        'model_name': 'M0_intercept_only',
        'aic': model.aic,
        'bic': model.bic,
        'log_likelihood': model.llf,
        'n_params': 1,
        'residual_std': np.std(model.resid)
    }

def fit_model_M1(df):
    """M1: Mutation fixed effects only"""
    # Create dummy variables for mutations
    df_model = df.copy()
    mutation_dummies = pd.get_dummies(df_model['Mutation'], prefix='mut', drop_first=True)

    X = sm.add_constant(mutation_dummies).astype(float)
    y = df_model['log10_FC'].values.astype(float)

    model = sm.OLS(y, X).fit()

    return {
        'model_name': 'M1_mutation_only',
        'aic': model.aic,
        'bic': model.bic,
        'log_likelihood': model.llf,
        'n_params': len(model.params),
        'residual_std': np.std(model.resid),
        'r_squared': model.rsquared
    }

def fit_model_M2(df):
    """M2: Mutation fixed + Subtype random effects"""
    try:
        # Prepare data
        df_model = df.copy()
        df_model['intercept'] = 1.0

        # Create mutation dummies
        mutation_dummies = pd.get_dummies(df_model['Mutation'], prefix='mut', drop_first=True).astype(float)

        # Combine
        X = pd.concat([df_model[['intercept']], mutation_dummies], axis=1).astype(float)

        # Fit mixed model with subtype as random effect
        model = MixedLM(df_model['log10_FC'], X, groups=df_model['Subtype']).fit(method='powell')

        return {
            'model_name': 'M2_mutation_subtype',
            'aic': model.aic,
            'bic': model.bic,
            'log_likelihood': model.llf,
            'n_params': len(model.params) + 1,  # +1 for random effect variance
            'residual_std': np.sqrt(model.scale),
            'random_effect_var': float(model.cov_re.iloc[0, 0]) if hasattr(model, 'cov_re') else 0.0,
            'converged': model.converged
        }
    except Exception as e:
        print(f"Warning: M2 fitting failed: {e}")
        return {
            'model_name': 'M2_mutation_subtype',
            'aic': np.nan,
            'bic': np.nan,
            'error': str(e)
        }

def fit_model_M3(df):
    """M3: Mutation fixed + Study context random effects"""
    try:
        df_model = df.copy()
        df_model['intercept'] = 1.0

        mutation_dummies = pd.get_dummies(df_model['Mutation'], prefix='mut', drop_first=True).astype(float)
        X = pd.concat([df_model[['intercept']], mutation_dummies], axis=1).astype(float)

        # Use study_source as grouping variable
        model = MixedLM(df_model['log10_FC'], X, groups=df_model['study_source']).fit(method='powell')

        return {
            'model_name': 'M3_mutation_study',
            'aic': model.aic,
            'bic': model.bic,
            'log_likelihood': model.llf,
            'n_params': len(model.params) + 1,
            'residual_std': np.sqrt(model.scale),
            'random_effect_var': float(model.cov_re.iloc[0, 0]) if hasattr(model, 'cov_re') else 0.0,
            'converged': model.converged
        }
    except Exception as e:
        print(f"Warning: M3 fitting failed: {e}")
        return {
            'model_name': 'M3_mutation_study',
            'aic': np.nan,
            'bic': np.nan,
            'error': str(e)
        }

def leave_one_study_out_cv(df):
    """Leave-one-study-out cross-validation"""

    studies = df['study_source'].unique()
    cv_results = []

    print(f"\nPerforming leave-one-study-out CV ({len(studies)} studies)...")

    for study in studies:
        train = df[df['study_source'] != study].copy()
        test = df[df['study_source'] == study].copy()

        if len(train) < 5 or len(test) < 1:
            continue

        # Fit M1 on training data
        mutation_dummies_train = pd.get_dummies(train['Mutation'], prefix='mut', drop_first=True).astype(float)
        X_train = sm.add_constant(mutation_dummies_train).astype(float)
        y_train = train['log10_FC'].values.astype(float)

        model = sm.OLS(y_train, X_train).fit()

        # Predict on test data
        mutation_dummies_test = pd.get_dummies(test['Mutation'], prefix='mut', drop_first=True).astype(float)

        # Align columns
        for col in X_train.columns:
            if col not in mutation_dummies_test.columns and col != 'const':
                mutation_dummies_test[col] = 0.0

        X_test = sm.add_constant(mutation_dummies_test[X_train.columns.drop('const')]).astype(float)

        try:
            y_pred = model.predict(X_test)
            y_true = test['log10_FC'].values

            mse = np.mean((y_true - y_pred) ** 2)
            mae = np.mean(np.abs(y_true - y_pred))

            cv_results.append({
                'held_out_study': study,
                'n_train': len(train),
                'n_test': len(test),
                'mse': mse,
                'mae': mae,
                'rmse': np.sqrt(mse)
            })
        except:
            continue

    cv_df = pd.DataFrame(cv_results)

    print(f"[OK] CV completed: {len(cv_df)} folds")
    print(f"  Mean RMSE: {cv_df['rmse'].mean():.3f}")
    print(f"  Mean MAE: {cv_df['mae'].mean():.3f}")

    return cv_df

def bootstrap_ranking_stability(df, n_bootstrap=1000):
    """Bootstrap analysis of mutation ranking stability"""

    print(f"\nBootstrap ranking stability (n={n_bootstrap})...")

    mutations = df['Mutation'].unique()
    rankings = {mut: [] for mut in mutations}

    for i in range(n_bootstrap):
        if (i + 1) % 100 == 0:
            print(f"  Bootstrap {i+1}/{n_bootstrap}")

        # Resample with replacement
        boot_sample = df.sample(n=len(df), replace=True)

        # Calculate mean log10_FC for each mutation
        mut_means = boot_sample.groupby('Mutation')['log10_FC'].mean()

        # Rank mutations (highest to lowest)
        mut_ranks = mut_means.rank(ascending=False)

        for mut in mutations:
            if mut in mut_ranks.index:
                rankings[mut].append(mut_ranks[mut])

    # Summarize rankings
    ranking_summary = []
    for mut in mutations:
        ranks = rankings[mut]
        if len(ranks) > 0:
            ranking_summary.append({
                'mutation': mut,
                'mean_rank': np.mean(ranks),
                'median_rank': np.median(ranks),
                'std_rank': np.std(ranks),
                'rank_95ci_lower': np.percentile(ranks, 2.5),
                'rank_95ci_upper': np.percentile(ranks, 97.5),
                'rank_stability': 1 - (np.std(ranks) / len(mutations))  # Higher = more stable
            })

    ranking_df = pd.DataFrame(ranking_summary)
    ranking_df = ranking_df.sort_values('mean_rank')

    print(f"[OK] Bootstrap complete")
    print(f"  Most stable rankings (low std):")
    for _, row in ranking_df.head(5).iterrows():
        print(f"    {row['mutation']}: rank {row['mean_rank']:.1f} ± {row['std_rank']:.1f}")

    return ranking_df

def main():
    """Main execution"""
    print("="*60)
    print("Model Comparison Framework")
    print("="*60)

    # Load data
    df = load_data()

    # Fit models
    print("\nFitting models...")
    results = []

    print("  M0: Intercept only...")
    results.append(fit_model_M0(df))

    print("  M1: Mutation fixed effects...")
    results.append(fit_model_M1(df))

    print("  M2: Mutation + Subtype random...")
    results.append(fit_model_M2(df))

    print("  M3: Mutation + Study random...")
    results.append(fit_model_M3(df))

    # Save model comparison
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / "model_comparison.csv", index=False)

    # Also save as JSON for easier reading
    with open(OUTPUT_DIR / "model_comparison.json", 'w') as f:
        json.dump(results, f, indent=2)

    print("\n[OK] Model comparison saved")

    # Print comparison
    print("\n" + "="*60)
    print("Model Comparison Results")
    print("="*60)
    for result in results:
        print(f"\n{result['model_name']}:")
        print(f"  AIC: {result.get('aic', 'N/A')}")
        print(f"  BIC: {result.get('bic', 'N/A')}")
        if 'r_squared' in result:
            print(f"  R-squared: {result['r_squared']:.3f}")
        if 'random_effect_var' in result:
            print(f"  Random effect variance: {result['random_effect_var']:.6f}")

    # Determine best model
    valid_results = [r for r in results if not np.isnan(r.get('aic', np.nan))]
    if valid_results:
        best_aic = min(valid_results, key=lambda x: x['aic'])
        best_bic = min(valid_results, key=lambda x: x['bic'])
        print(f"\nBest model by AIC: {best_aic['model_name']}")
        print(f"Best model by BIC: {best_bic['model_name']}")

    # Cross-validation
    cv_results = leave_one_study_out_cv(df)
    cv_results.to_csv(OUTPUT_DIR / "leave_one_study_out_cv.csv", index=False)

    # Bootstrap ranking stability
    bootstrap_ranks = bootstrap_ranking_stability(df, n_bootstrap=1000)
    bootstrap_ranks.to_csv(OUTPUT_DIR / "bootstrap_ranks.csv", index=False)

    print("\n" + "="*60)
    print("Model comparison complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
