"""
Sensitivity Analysis Framework for HIV Capsid Surveillance Analysis

This module provides:
- Leave-one-country-out analysis
- Temporal stability assessment
- Subtype stratification analysis
- Near-full-length vs partial sequence comparison

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Callable
from dataclasses import dataclass, field
from collections import defaultdict
import json

import numpy as np
import pandas as pd
from scipy import stats

from . import config
from .sequence_qc import SequenceRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class SensitivityResult:
    """Container for sensitivity analysis results."""
    analysis_type: str
    parameter: str
    parameter_value: Any

    # Results
    baseline_score: float = 0.0
    test_score: float = 0.0
    score_change: float = 0.0
    score_change_pct: float = 0.0

    # Statistical significance
    p_value: float = 1.0
    is_significant: bool = False

    # Details
    affected_countries: List[str] = field(default_factory=list)
    details: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return {
            'analysis_type': self.analysis_type,
            'parameter': self.parameter,
            'parameter_value': str(self.parameter_value),
            'baseline_score': self.baseline_score,
            'test_score': self.test_score,
            'score_change': self.score_change,
            'score_change_pct': self.score_change_pct,
            'p_value': self.p_value,
            'is_significant': self.is_significant,
            'affected_countries': ','.join(self.affected_countries),
        }


@dataclass
class StabilityResult:
    """Container for stability assessment results."""
    metric_name: str
    stability_score: float = 0.0  # 0-1, higher is more stable
    coefficient_of_variation: float = 0.0
    trend_direction: str = 'stable'
    trend_p_value: float = 1.0

    def to_dict(self) -> Dict:
        return {
            'metric': self.metric_name,
            'stability_score': self.stability_score,
            'coefficient_of_variation': self.coefficient_of_variation,
            'trend_direction': self.trend_direction,
            'trend_p_value': self.trend_p_value,
        }


class LeaveOneCountryOutAnalyzer:
    """
    Analyze sensitivity to individual country data.

    Tests how removing each country affects overall and regional statistics.
    """

    def __init__(self, country_scores: Dict[str, float],
                 regional_mapping: Dict[str, str]):
        """
        Initialize analyzer.

        Args:
            country_scores: Dictionary of country -> CIR-SRI score
            regional_mapping: Dictionary of country -> region
        """
        self.country_scores = country_scores
        self.regional_mapping = regional_mapping
        self.all_countries = list(country_scores.keys())

    def _calculate_regional_means(self, exclude_country: Optional[str] = None) -> Dict[str, float]:
        """Calculate regional mean CIR-SRI scores."""
        regions = defaultdict(list)

        for country, score in self.country_scores.items():
            if country == exclude_country:
                continue
            region = self.regional_mapping.get(country, 'UNKNOWN')
            if region != 'UNKNOWN':
                regions[region].append(score)

        return {region: np.mean(scores) for region, scores in regions.items()}

    def analyze_all(self) -> List[SensitivityResult]:
        """
        Run leave-one-country-out analysis for all countries.

        Returns:
            List of SensitivityResult objects
        """
        # Calculate baseline
        baseline_global = np.mean(list(self.country_scores.values()))
        baseline_regional = self._calculate_regional_means()

        results = []

        for country in self.all_countries:
            # Calculate without this country
            remaining_scores = {c: s for c, s in self.country_scores.items() if c != country}
            test_global = np.mean(list(remaining_scores.values()))

            test_regional = self._calculate_regional_means(exclude_country=country)

            # Calculate changes
            global_change = test_global - baseline_global
            global_change_pct = (global_change / baseline_global * 100) if baseline_global else 0

            # Regional changes
            regional_changes = {}
            for region in baseline_regional:
                baseline_r = baseline_regional.get(region, 0)
                test_r = test_regional.get(region, 0)
                change = test_r - baseline_r
                change_pct = (change / baseline_r * 100) if baseline_r else 0
                regional_changes[region] = {'change': change, 'change_pct': change_pct}

            result = SensitivityResult(
                analysis_type='leave_one_country_out',
                parameter='country',
                parameter_value=country,
                baseline_score=baseline_global,
                test_score=test_global,
                score_change=global_change,
                score_change_pct=global_change_pct,
                affected_countries=[country],
                details={
                    'regional_changes': regional_changes,
                    'countries_remaining': len(remaining_scores),
                    'country_share': self.country_scores[country] / baseline_global
                }
            )

            # Calculate significance (simplified)
            if len(remaining_scores) >= 2:
                _, p_value = stats.ttest_1samp(
                    list(remaining_scores.values()),
                    baseline_global
                )
                result.p_value = p_value
                result.is_significant = p_value < 0.05

            results.append(result)

        # Sort by absolute impact
        results.sort(key=lambda r: abs(r.score_change), reverse=True)

        logger.info(f"Completed LOO analysis for {len(results)} countries")

        return results


class TemporalStabilityAnalyzer:
    """
    Analyze temporal stability of surveillance readiness metrics.
    """

    def __init__(self, sequence_records: List[SequenceRecord]):
        """
        Initialize analyzer.

        Args:
            sequence_records: List of SequenceRecord objects with year info
        """
        self.records = sequence_records
        self.yearly_data = self._aggregate_by_year()

    def _aggregate_by_year(self) -> Dict[int, Dict]:
        """Aggregate data by year."""
        yearly_data = defaultdict(lambda: {
            'n_sequences': 0,
            'n_countries': 0,
            'subtypes': set(),
        })

        for record in self.records:
            if record.year:
                yearly_data[record.year]['n_sequences'] += 1
                yearly_data[record.year]['n_countries'].add(record.country)
                if record.subtype:
                    yearly_data[record.year]['subtypes'].add(record.subtype)

        return dict(yearly_data)

    def analyze_sequence_volume(self, min_year: int = 2000) -> StabilityResult:
        """
        Analyze stability of sequence volume over time.

        Args:
            min_year: Minimum year to include in analysis

        Returns:
            StabilityResult object
        """
        years = sorted(self.yearly_data.keys())
        years = [y for y in years if y >= min_year]

        if len(years) < 2:
            return StabilityResult(metric_name='sequence_volume')

        volumes = [self.yearly_data[y]['n_sequences'] for y in years]

        mean_vol = np.mean(volumes)
        std_vol = np.std(volumes)
        cv = std_vol / mean_vol if mean_vol else 0

        # Trend analysis
        slope, intercept, r_value, p_value, std_err = stats.linregress(years, volumes)

        trend_direction = 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'

        # Stability score (inverse of CV, normalized)
        stability_score = 1 / (1 + cv) if cv is not None else 0

        return StabilityResult(
            metric_name='sequence_volume',
            stability_score=stability_score,
            coefficient_of_variation=cv,
            trend_direction=trend_direction,
            trend_p_value=p_value
        )

    def analyze_country_coverage(self, min_year: int = 2000) -> StabilityResult:
        """Analyze stability of country coverage over time."""
        years = sorted(self.yearly_data.keys())
        years = [y for y in years if y >= min_year]

        if len(years) < 2:
            return StabilityResult(metric_name='country_coverage')

        coverages = [self.yearly_data[y]['n_countries'] for y in years]

        mean_cov = np.mean(coverages)
        std_cov = np.std(coverages)
        cv = std_cov / mean_cov if mean_cov else 0

        slope, _, _, p_value, _ = stats.linregress(years, coverages)
        trend_direction = 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'

        stability_score = 1 / (1 + cv) if cv is not None else 0

        return StabilityResult(
            metric_name='country_coverage',
            stability_score=stability_score,
            coefficient_of_variation=cv,
            trend_direction=trend_direction,
            trend_p_value=p_value
        )

    def analyze_subtype_diversity(self, min_year: int = 2000) -> StabilityResult:
        """Analyze stability of subtype diversity over time."""
        years = sorted(self.yearly_data.keys())
        years = [y for y in years if y >= min_year]

        if len(years) < 2:
            return StabilityResult(metric_name='subtype_diversity')

        diversities = [len(self.yearly_data[y]['subtypes']) for y in years]

        mean_div = np.mean(diversities)
        std_div = np.std(diversities)
        cv = std_div / mean_div if mean_div else 0

        slope, _, _, p_value, _ = stats.linregress(years, diversities)
        trend_direction = 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'

        stability_score = 1 / (1 + cv) if cv is not None else 0

        return StabilityResult(
            metric_name='subtype_diversity',
            stability_score=stability_score,
            coefficient_of_variation=cv,
            trend_direction=trend_direction,
            trend_p_value=p_value
        )

    def analyze_all(self) -> List[StabilityResult]:
        """Run all stability analyses."""
        results = [
            self.analyze_sequence_volume(),
            self.analyze_country_coverage(),
            self.analyze_subtype_diversity(),
        ]

        logger.info(f"Completed temporal stability analysis")

        return results


class SubtypeStratificationAnalyzer:
    """
    Analyze how subtype composition affects surveillance readiness.
    """

    def __init__(self, records: List[SequenceRecord]):
        """
        Initialize analyzer.

        Args:
            records: List of SequenceRecord objects
        """
        self.records = records
        self.subtype_data = self._aggregate_by_subtype()

    def _aggregate_by_subtype(self) -> Dict[str, Dict]:
        """Aggregate data by subtype."""
        subtype_data = defaultdict(lambda: {
            'n_sequences': 0,
            'countries': set(),
            'years': set(),
        })

        for record in self.records:
            subtype = record.subtype or 'UNKNOWN'
            subtype_data[subtype]['n_sequences'] += 1
            if record.country:
                subtype_data[subtype]['countries'].add(record.country)
            if record.year:
                subtype_data[subtype]['years'].add(record.year)

        return dict(subtype_data)

    def analyze_coverage_by_subtype(self) -> pd.DataFrame:
        """
        Analyze RAM site coverage stratified by subtype.

        Returns:
            DataFrame with subtype coverage metrics
        """
        rows = []

        for subtype, data in self.subtype_data.items():
            rows.append({
                'subtype': subtype,
                'n_sequences': data['n_sequences'],
                'n_countries': len(data['countries']),
                'year_range': f"{min(data['years'])}-{max(data['years'])}" if data['years'] else 'N/A',
                'years_covered': len(data['years']),
            })

        df = pd.DataFrame(rows)
        df = df.sort_values('n_sequences', ascending=False)

        return df

    def test_subtype_balance(self) -> Tuple[float, float]:
        """
        Test if subtype distribution is balanced across countries.

        Returns:
            Tuple of (gini_coefficient, chi_square_p_value)
        """
        # Create contingency table
        countries = sorted(set(r.country for r in self.records))
        subtypes = sorted(set(r.subtype for r in self.records if r.subtype))

        contingency = np.zeros((len(countries), len(subtypes)))

        for record in self.records:
            if record.country and record.subtype:
                ci = countries.index(record.country)
                si = subtypes.index(record.subtype)
                contingency[ci, si] += 1

        # Chi-square test
        chi2, p_value, dof, expected = stats.chi2_contingency(contingency)

        # Gini coefficient for subtype distribution
        subtype_totals = contingency.sum(axis=0)
        n = len(subtype_totals)
        if n > 0:
            sorted_totals = np.sort(subtype_totals)
            cumsum = np.cumsum(sorted_totals)
            gini = (2 * np.sum((np.arange(1, n + 1) * sorted_totals))) / (n * cumsum[-1]) - (n + 1) / n
        else:
            gini = 0

        return gini, p_value


class SequenceLengthAnalyzer:
    """
    Compare near-full-length vs partial sequences.
    """

    def __init__(self, records: List[SequenceRecord]):
        """
        Initialize analyzer.

        Args:
            records: List of SequenceRecord objects
        """
        self.records = records
        self._classify_sequences()

    def _classify_sequences(self):
        """Classify sequences by length."""
        self.nfl_sequences = []  # >= 8000 bp
        self.partial_sequences = []  # 7000-8000 bp
        self.short_sequences = []  # < 7000 bp

        for record in self.records:
            if record.length >= 8000:
                self.nfl_sequences.append(record)
            elif record.length >= 7000:
                self.partial_sequences.append(record)
            else:
                self.short_sequences.append(record)

    def compare_ram_coverage(self) -> Dict:
        """
        Compare RAM site coverage between sequence groups.

        Returns:
            Dictionary with comparison statistics
        """
        results = {}

        for group_name, group_records in [
            ('near_full_length', self.nfl_sequences),
            ('partial', self.partial_sequences),
            ('short', self.short_sequences),
        ]:
            if not group_records:
                continue

            # Count sequences with coverage
            site_coverage = defaultdict(int)
            total = len(group_records)

            for record in group_records:
                for site_name, hxb2_pos in config.RAM_SITES.items():
                    nucleotide_pos = hxb2_pos - 1
                    if nucleotide_pos >= 0 and nucleotide_pos + 2 < len(record.sequence):
                        codon = record.sequence[nucleotide_pos:nucleotide_pos + 3]
                        if len(codon) == 3 and all(c in 'ATGC' for c in codon):
                            site_coverage[site_name] += 1

            results[group_name] = {
                'n_sequences': total,
                'site_coverage': {s: c / total for s, c in site_coverage.items()}
            }

        return results

    def analyze_length_distribution(self) -> Dict:
        """Analyze length distribution and impact on coverage."""
        lengths = [r.length for r in self.records]

        return {
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'std_length': np.std(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'n_nfl_sequences': len(self.nfl_sequences),
            'n_partial_sequences': len(self.partial_sequences),
            'n_short_sequences': len(self.short_sequences),
            'pct_nfl': len(self.nfl_sequences) / len(lengths) * 100,
            'pct_partial': len(self.partial_sequences) / len(lengths) * 100,
        }


class ParameterSensitivityAnalyzer:
    """
    Analyze sensitivity to pipeline parameters.
    """

    def __init__(self, analyzer_func: Callable):
        """
        Initialize analyzer.

        Args:
            analyzer_func: Function that takes parameters and returns scores
        """
        self.analyzer_func = analyzer_func

    def run_min_year_sensitivity(self) -> List[SensitivityResult]:
        """
        Test sensitivity to minimum year parameter.

        Returns:
            List of SensitivityResult objects
        """
        results = []

        # Baseline
        baseline_scores = self.analyzer_func(min_year=2000)
        baseline_mean = np.mean(list(baseline_scores.values()))

        for min_year in config.SENSITIVITY_PARAMS['min_year']:
            test_scores = self.analyzer_func(min_year=min_year)
            test_mean = np.mean(list(test_scores.values()))

            change = test_mean - baseline_mean
            change_pct = (change / baseline_mean * 100) if baseline_mean else 0

            result = SensitivityResult(
                analysis_type='min_year_sensitivity',
                parameter='min_year',
                parameter_value=min_year,
                baseline_score=baseline_mean,
                test_score=test_mean,
                score_change=change,
                score_change_pct=change_pct,
                details={
                    'countries_affected': sum(1 for r in self.records if r.year and r.year < min_year),
                }
            )
            results.append(result)

        return results

    def run_length_threshold_sensitivity(self) -> List[SensitivityResult]:
        """Test sensitivity to minimum sequence length."""
        results = []

        baseline_scores = self.analyzer_func(min_sequence_length=7000)
        baseline_mean = np.mean(list(baseline_scores.values()))

        for min_length in config.SENSITIVITY_PARAMS['min_sequence_length']:
            test_scores = self.analyzer_func(min_sequence_length=min_length)
            test_mean = np.mean(list(test_scores.values()))

            change = test_mean - baseline_mean
            change_pct = (change / baseline_mean * 100) if baseline_mean else 0

            result = SensitivityResult(
                analysis_type='length_threshold_sensitivity',
                parameter='min_sequence_length',
                parameter_value=min_length,
                baseline_score=baseline_mean,
                test_score=test_mean,
                score_change=change,
                score_change_pct=change_pct,
            )
            results.append(result)

        return results


def run_complete_sensitivity_analysis(
    country_scores: Dict[str, float],
    regional_mapping: Dict[str, str],
    sequence_records: List[SequenceRecord],
    output_dir: str
) -> Dict[str, Any]:
    """
    Run complete sensitivity analysis suite.

    Args:
        country_scores: Dictionary of country -> CIR-SRI score
        regional_mapping: Dictionary of country -> region
        sequence_records: List of SequenceRecord objects
        output_dir: Output directory

    Returns:
        Dictionary with all results
    """
    logger.info("Starting complete sensitivity analysis...")

    results = {}

    # 1. Leave-one-country-out analysis
    loo_analyzer = LeaveOneCountryOutAnalyzer(country_scores, regional_mapping)
    results['leave_one_country_out'] = loo_analyzer.analyze_all()

    # 2. Temporal stability analysis
    temporal_analyzer = TemporalStabilityAnalyzer(sequence_records)
    results['temporal_stability'] = temporal_analyzer.analyze_all()

    # 3. Subtype stratification
    subtype_analyzer = SubtypeStratificationAnalyzer(sequence_records)
    results['subtype_coverage'] = subtype_analyzer.analyze_coverage_by_subtype()
    results['subtype_balance'] = list(subtype_analyzer.test_subtype_balance())

    # 4. Sequence length analysis
    length_analyzer = SequenceLengthAnalyzer(sequence_records)
    results['length_comparison'] = length_analyzer.compare_ram_coverage()
    results['length_distribution'] = length_analyzer.analyze_length_distribution()

    # Save results
    os.makedirs(output_dir, exist_ok=True)

    # Save LOO results
    loo_df = pd.DataFrame([r.to_dict() for r in results['leave_one_country_out']])
    loo_df.to_csv(os.path.join(output_dir, 'sensitivity_loo.csv'), index=False)

    # Save stability results
    stability_df = pd.DataFrame([r.to_dict() for r in results['temporal_stability']])
    stability_df.to_csv(os.path.join(output_dir, 'sensitivity_temporal_stability.csv'), index=False)

    # Save subtype coverage
    results['subtype_coverage'].to_csv(os.path.join(output_dir, 'sensitivity_subtype_coverage.csv'), index=False)

    # Save length analysis
    length_results = [results['length_comparison'], results['length_distribution']]
    with open(os.path.join(output_dir, 'sensitivity_length_analysis.json'), 'w') as f:
        json.dump({'comparison': results['length_comparison'],
                   'distribution': results['length_distribution']}, f, indent=2)

    logger.info(f"Sensitivity analysis complete. Results saved to {output_dir}")

    return results


def generate_sensitivity_report(results: Dict) -> str:
    """
    Generate text summary of sensitivity analysis.

    Args:
        results: Dictionary from run_complete_sensitivity_analysis

    Returns:
        Formatted report string
    """
    lines = [
        "=" * 80,
        "Sensitivity Analysis Report",
        "=" * 80,
        "",
    ]

    # LOO Analysis
    lines.append("Leave-One-Country-Out Analysis:")
    lines.append("-" * 40)

    loo_results = results.get('leave_one_country_out', [])
    if loo_results:
        most_impactful = sorted(
            loo_results,
            key=lambda r: abs(r.score_change),
            reverse=True
        )[:5]

        lines.append("Most impactful countries (by effect on global mean):")
        for result in most_impactful:
            lines.append(
                f"  {result.parameter_value}: {result.score_change:+.4f} "
                f"({result.score_change_pct:+.2f}%)"
            )

        # Check for instability
        changes = [r.score_change for r in loo_results]
        max_change = max(abs(c) for c in changes)
        lines.append(f"\nMaximum impact: {max_change:.4f}")

    lines.append("")

    # Temporal Stability
    lines.append("Temporal Stability Analysis:")
    lines.append("-" * 40)

    stability_results = results.get('temporal_stability', [])
    for result in stability_results:
        status = "STABLE" if result.stability_score > 0.7 else "MODERATE" if result.stability_score > 0.4 else "UNSTABLE"
        lines.append(f"  {result.metric_name}:")
        lines.append(f"    Stability score: {result.stability_score:.3f} ({status})")
        lines.append(f"    CV: {result.coefficient_of_variation:.3f}")
        lines.append(f"    Trend: {result.trend_direction} (p={result.trend_p_value:.4f})")

    lines.append("")

    # Subtype Analysis
    lines.append("Subtype Stratification:")
    lines.append("-" * 40)

    if 'subtype_balance' in results:
        gini, p_value = results['subtype_balance']
        lines.append(f"  Gini coefficient: {gini:.3f}")
        lines.append(f"  Chi-square p-value: {p_value:.4f}")
        if p_value < 0.05:
            lines.append("  Subtype distribution is significantly unbalanced (p < 0.05)")
        else:
            lines.append("  Subtype distribution is not significantly unbalanced")

    lines.append("")

    # Length Analysis
    lines.append("Sequence Length Analysis:")
    lines.append("-" * 40)

    if 'length_distribution' in results:
        dist = results['length_distribution']
        lines.append(f"  NFL sequences (>=8000): {dist['n_nfl_sequences']} ({dist['pct_nfl']:.1f}%)")
        lines.append(f"  Partial (7000-8000): {dist['n_partial_sequences']} ({dist['pct_partial']:.1f}%)")
        lines.append(f"  Mean length: {dist['mean_length']:.0f}")
        lines.append(f"  Median length: {dist['median_length']:.0f}")

    lines.append("")
    lines.append("=" * 80)

    return "\n".join(lines)


def main():
    """Example usage."""
    import argparse

    parser = argparse.ArgumentParser(description='Run sensitivity analysis')
    parser.add_argument('cir_sri_csv', help='CIR-SRI scores CSV')
    parser.add_argument('sequences_csv', help='Processed sequences CSV')
    parser.add_argument('-o', '--output', default='data/sensitivity', help='Output directory')

    args = parser.parse_args()

    # Load CIR-SRI scores
    cir_sri_df = pd.read_csv(args.cir_sri_csv)
    country_scores = cir_sri_df.set_index('country')['cir_sri_score'].to_dict()
    regional_mapping = cir_sri_df.set_index('country')['region'].to_dict()

    # Load sequence records (simplified)
    seq_df = pd.read_csv(args.sequences_csv)
    records = []
    for _, row in seq_df.iterrows():
        record = SequenceRecord(
            accession=row.get('accession', ''),
            sequence='',
            sequence_aa='',
            country=row.get('country'),
            year=row.get('year'),
            subtype=row.get('subtype'),
            length=row.get('length', 0),
        )
        records.append(record)

    # Run analysis
    results = run_complete_sensitivity_analysis(
        country_scores, regional_mapping, records, args.output
    )

    # Generate report
    report = generate_sensitivity_report(results)
    print(report)

    # Save report
    with open(os.path.join(args.output, 'sensitivity_report.txt'), 'w') as f:
        f.write(report)


if __name__ == '__main__':
    main()