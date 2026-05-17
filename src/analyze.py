"""
Core analysis: model comparison, bootstrap ranking, sensitivity, epistasis,
and compensatory mutation analysis.
"""

import pandas as pd
import numpy as np
import json
import statsmodels.api as sm
from statsmodels.regression.mixed_linear_model import MixedLM

from .config import REV2_DATA, REV2_RESULTS


# ── Data loading ──────────────────────────────────────────────────────────

def load_harmonized(clean=True):
    """Load the harmonized phenotype dataset."""
    path = REV2_DATA / ("harmonized_phenotype_data_clean.csv" if clean
                        else "harmonized_phenotype_data.csv")
    df = pd.read_csv(path)
    df = df[df["log10_FC"].notna()].copy()
    return df


# ── Model comparison (M0–M3) ─────────────────────────────────────────────

def _fit_M0(df):
    y = df["log10_FC"].values
    m = sm.OLS(y, np.ones(len(y))).fit()
    return {"aic": m.aic, "bic": m.bic, "log_likelihood": m.llf,
            "n_params": 1, "residual_std": np.std(m.resid)}


def _fit_M1(df):
    dummies = pd.get_dummies(df["Mutation"], prefix="mut", drop_first=True).astype(float)
    X = sm.add_constant(dummies)
    m = sm.OLS(df["log10_FC"].values.astype(float), X).fit()
    return {"aic": m.aic, "bic": m.bic, "log_likelihood": m.llf,
            "n_params": len(m.params), "residual_std": np.std(m.resid),
            "r_squared": m.rsquared}


def _fit_M2(df):
    try:
        dummies = pd.get_dummies(df["Mutation"], prefix="mut", drop_first=True).astype(float)
        X = pd.concat([pd.DataFrame({"intercept": 1.0}, index=df.index), dummies], axis=1)
        m = MixedLM(df["log10_FC"].values, X, groups=df["Subtype"]).fit(method="powell")
        return {"aic": m.aic, "bic": m.bic, "log_likelihood": m.llf,
                "n_params": len(m.params) + 1, "residual_std": np.sqrt(m.scale),
                "random_effect_var": float(m.cov_re.iloc[0, 0]) if m.cov_re is not None else 0.0,
                "converged": m.converged}
    except Exception as e:
        return {"aic": np.nan, "bic": np.nan, "error": str(e)}


def _fit_M3(df):
    try:
        dummies = pd.get_dummies(df["Mutation"], prefix="mut", drop_first=True).astype(float)
        X = pd.concat([pd.DataFrame({"intercept": 1.0}, index=df.index), dummies], axis=1)
        m = MixedLM(df["log10_FC"].values, X, groups=df["study_source"]).fit(method="powell")
        return {"aic": m.aic, "bic": m.bic, "log_likelihood": m.llf,
                "n_params": len(m.params) + 1, "residual_std": np.sqrt(m.scale),
                "random_effect_var": float(m.cov_re.iloc[0, 0]) if m.cov_re is not None else 0.0,
                "converged": m.converged}
    except Exception as e:
        return {"aic": np.nan, "bic": np.nan, "error": str(e)}


def compare_models(df):
    """Fit M0–M3 and return comparison table, save results."""
    results = [
        {"model_name": "M0_intercept_only", **_fit_M0(df)},
        {"model_name": "M1_mutation_only", **_fit_M1(df)},
        {"model_name": "M2_mutation_subtype", **_fit_M2(df)},
        {"model_name": "M3_mutation_study", **_fit_M3(df)},
    ]
    pd.DataFrame(results).to_csv(REV2_RESULTS / "model_comparison.csv", index=False)
    with open(REV2_RESULTS / "model_comparison.json", "w") as f:
        json.dump(results, f, indent=2)
    return results


# ── Leave-one-study-out CV ───────────────────────────────────────────────

def leave_one_study_out_cv(df):
    """LOSO cross-validation for M1 (mutation-only)."""
    studies = df["study_source"].unique()
    rows = []
    for study in studies:
        train = df[df["study_source"] != study]
        test = df[df["study_source"] == study]
        if len(train) < 5 or len(test) < 1:
            continue
        dummies = pd.get_dummies(train["Mutation"], prefix="mut", drop_first=True).astype(float)
        X_train = sm.add_constant(dummies)
        m = sm.OLS(train["log10_FC"].values.astype(float), X_train).fit()

        test_dummies = pd.get_dummies(test["Mutation"], prefix="mut", drop_first=True).astype(float)
        for col in X_train.columns:
            if col not in test_dummies.columns and col != "const":
                test_dummies[col] = 0.0
        X_test = sm.add_constant(test_dummies[X_train.columns.drop("const")])
        try:
            pred = m.predict(X_test)
            err = test["log10_FC"].values - pred
            rows.append({"held_out_study": study, "n_train": len(train),
                         "n_test": len(test), "rmse": np.sqrt(np.mean(err**2)),
                         "mae": np.mean(np.abs(err))})
        except Exception:
            continue
    cv_df = pd.DataFrame(rows)
    cv_df.to_csv(REV2_RESULTS / "leave_one_study_out_cv.csv", index=False)
    return cv_df


# ── Bootstrap ranking stability ──────────────────────────────────────────

def bootstrap_ranks(df, n=1000):
    """Resample with replacement and compute rank mean/SD per mutation."""
    mutations = df["Mutation"].unique()
    rankings = {m: [] for m in mutations}
    for _ in range(n):
        boot = df.sample(n=len(df), replace=True)
        means = boot.groupby("Mutation")["log10_FC"].mean()
        ranks = means.rank(ascending=False)
        for m in mutations:
            if m in ranks.index:
                rankings[m].append(ranks[m])

    rows = []
    for m, r in rankings.items():
        if not r:
            continue
        rows.append({"mutation": m, "mean_rank": np.mean(r), "median_rank": np.median(r),
                      "std_rank": np.std(r),
                      "rank_95ci_lower": np.percentile(r, 2.5),
                      "rank_95ci_upper": np.percentile(r, 97.5)})
    rdf = pd.DataFrame(rows).sort_values("mean_rank")
    rdf.to_csv(REV2_RESULTS / "bootstrap_ranks.csv", index=False)
    return rdf


# ── Sensitivity analysis ─────────────────────────────────────────────────

def context_stratification(df):
    """Stratify by context tier."""
    return {
        ctx: {"n": len(sub), "mean_log10fc": float(sub["log10_FC"].mean()),
              "std_log10fc": float(sub["log10_FC"].std())}
        for ctx in df["context_tier"].dropna().unique()
        if len(sub := df[df["context_tier"] == ctx]) >= 3
    }


def leave_one_subtype_out(df):
    """Leave-one-subtype-out means."""
    full_mean = float(df["log10_FC"].mean())
    results = {}
    for st in df["subtype"].dropna().unique():
        train = df[df["subtype"] != st]
        if len(train) >= 5:
            results[st] = {"n_removed": int((df["subtype"] == st).sum()),
                           "mean_without": float(train["log10_FC"].mean()),
                           "delta": float(abs(train["log10_FC"].mean() - full_mean))}
    return results


def detect_outliers(df, sigma=3):
    """3-sigma outlier detection."""
    mean, std = df["log10_FC"].mean(), df["log10_FC"].std()
    outliers = df[np.abs(df["log10_FC"] - mean) > sigma * std]
    return {"threshold_sigma": sigma, "n_outliers": len(outliers),
            "outlier_mutations": outliers["Mutation"].tolist()}


def run_sensitivity(df):
    """All sensitivity analyses."""
    results = {
        "context_stratification": context_stratification(df),
        "leave_one_subtype_out": leave_one_subtype_out(df),
        "outlier_detection": detect_outliers(df),
    }
    with open(REV2_RESULTS / "sensitivity_results.json", "w") as f:
        json.dump(results, f, indent=2)
    return results


# ── Epistasis analysis ───────────────────────────────────────────────────

def calc_epistasis(df):
    """Additive vs synergistic vs compensatory for double mutants."""
    singles = df[~df["Mutation"].str.contains(r"\+", na=False)]
    singles_mean = singles.groupby("Mutation")["log10_FC"].mean()
    dm = df[df["Mutation"].str.contains(r"\+", na=False)]
    rows = []
    for _, r in dm.iterrows():
        muts = str(r["Mutation"]).split("+")
        if len(muts) != 2:
            continue
        m1, m2 = muts
        m1_fc = singles_mean.get(m1, np.nan)
        m2_fc = singles_mean.get(m2, np.nan)
        if np.isnan(m1_fc) or np.isnan(m2_fc):
            continue
        expected = m1_fc + m2_fc
        observed = r["log10_FC"]
        residual = observed - expected
        itype = ("additive" if abs(residual) < 0.3
                 else "positive_synergy" if residual > 0.3
                 else "negative_synergy")
        rows.append({"combination": r["Mutation"], "observed_log10fc": observed,
                      "expected_log10fc": expected, "interaction_residual": residual,
                      "interaction_type": itype,
                      "fold_amplification": (10**observed) / (10**expected)})
    epi = pd.DataFrame(rows)
    epi.to_csv(REV2_RESULTS / "epistasis_matrix.csv", index=False)
    return epi


def compensatory_patterns(df):
    """M66I-centered compensatory analysis."""
    m66i = df[df["Mutation"] == "M66I"]["log10_FC"].mean()
    m66i_fc = 10 ** m66i
    known = [
        ("M66I+A105T", 111, m66i_fc, "putative_compensatory"),
        ("M66I+T107A", 234, m66i_fc, "putative_compensatory"),
        ("M66I+N74D+A105T", 1337, m66i_fc, "putative_compensatory"),
    ]
    rows = []
    for combo, obs_fc, ref, pat in known:
        rows.append({"combination": combo, "observed_fc": obs_fc,
                      "m66i_alone_fc": ref, "ratio_vs_m66i": obs_fc / ref,
                      "pattern": pat})
    comp = pd.DataFrame(rows)
    comp.to_csv(REV2_RESULTS / "compensatory_patterns.csv", index=False)
    return comp


def context_dependent_combinations(df):
    """List combinations observed in multiple studies/contexts."""
    dm = df[df["Mutation"].str.contains(r"\+", na=False)]
    rows = []
    for combo, grp in dm.groupby("Mutation"):
        if len(grp) > 1:
            for _, r in grp.iterrows():
                rows.append({"combination": combo, "log10_FC": r["log10_FC"],
                              "FC": r.get("FC_numeric", np.nan),
                              "context": r["context_tier"],
                              "study": r["study_source"], "subtype": r["Subtype"]})
    cdf = pd.DataFrame(rows)
    cdf.to_csv(REV2_RESULTS / "context_specific_combinations.csv", index=False)
    return cdf


# ── Main ──────────────────────────────────────────────────────────────────

def run_analysis():
    """Run all downstream analyses."""
    df = load_harmonized()
    print(f"Loaded {len(df)} observations, {df['Mutation'].nunique()} mutations")

    print("\n[1/6] Model comparison...")
    compare_models(df)

    print("[2/6] LOSO cross-validation...")
    leave_one_study_out_cv(df)

    print("[3/6] Bootstrap ranking...")
    bootstrap_ranks(df)

    print("[4/6] Sensitivity analysis...")
    run_sensitivity(df)

    print("[5/6] Epistasis analysis...")
    calc_epistasis(df)
    compensatory_patterns(df)
    context_dependent_combinations(df)

    print(f"\nDone. Results in {REV2_RESULTS}")


if __name__ == "__main__":
    run_analysis()
