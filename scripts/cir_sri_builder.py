"""
CIR-SRI (Capsid Inhibitor Resistance Surveillance Readiness Index) Builder

This module constructs the composite index with:
- 5-dimension scoring (0-2 each, total 0-10)
- Composite index calculation
- Priority ranking
- Sensitivity analysis framework

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

from . import config
from .cso_calculator import CSOCalculator, CSOResult
from .ros_calculator import ROSCalculator, ROSResult

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class DimensionScore:
    """Score for one dimension of CIR-SRI."""
    dimension_name: str
    raw_score: float = 0.0
    normalized_score: float = 0.0
    percentile: float = 0.0
    category: str = 'low'  # low, medium, high

    def to_dict(self) -> Dict:
        return {
            'dimension': self.dimension_name,
            'raw_score': self.raw_score,
            'normalized_score': self.normalized_score,
            'percentile': self.percentile,
            'category': self.category,
        }


@dataclass
class CIRSRICountryResult:
    """Complete CIR-SRI result for a country."""
    country: str
    region: str

    # Five dimensions (0-2 each)
    sequence_availability: DimensionScore = field(
        default_factory=lambda: DimensionScore('sequence_availability'))
    metadata_completeness: DimensionScore = field(
        default_factory=lambda: DimensionScore('metadata_completeness'))
    ram_site_observability: DimensionScore = field(
        default_factory=lambda: DimensionScore('ram_site_observability'))
    temporal_coverage: DimensionScore = field(
        default_factory=lambda: DimensionScore('temporal_coverage'))
    subtype_diversity: DimensionScore = field(
        default_factory=lambda: DimensionScore('subtype_diversity'))

    # Composite scores
    cso_score: float = 0.0
    ros_score: float = 0.0
    cir_sri_score: float = 0.0  # Total 0-10
    cir_sri_percentile: float = 0.0
    regional_rank: int = 0
    global_rank: int = 0

    # Additional context
    n_sequences: int = 0
    priority_tier: str = 'unknown'  # urgent, high, medium, low
    critical_gaps: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            'country': self.country,
            'region': self.region,
            'sequence_availability': self.sequence_availability.raw_score,
            'metadata_completeness': self.metadata_completeness.raw_score,
            'ram_site_observability': self.ram_site_observability.raw_score,
            'temporal_coverage': self.temporal_coverage.raw_score,
            'subtype_diversity': self.subtype_diversity.raw_score,
            'cso_score': self.cso_score,
            'ros_score': self.ros_score,
            'cir_sri_score': self.cir_sri_score,
            'cir_sri_percentile': self.cir_sri_percentile,
            'regional_rank': self.regional_rank,
            'global_rank': self.global_rank,
            'priority_tier': self.priority_tier,
            'n_sequences': self.n_sequences,
            'critical_gaps': ','.join(self.critical_gaps),
        }

    def get_scores_dict(self) -> Dict[str, float]:
        """Get dimension scores as dictionary."""
        return {
            'sequence_availability': self.sequence_availability.raw_score,
            'metadata_completeness': self.metadata_completeness.raw_score,
            'ram_site_observability': self.ram_site_observability.raw_score,
            'temporal_coverage': self.temporal_coverage.raw_score,
            'subtype_diversity': self.subtype_diversity.raw_score,
        }


class CIRSRIBuilder:
    """
    Build the CIR-SRI (Capsid Inhibitor Resistance Surveillance Readiness Index).

    The index combines 5 dimensions, each scored 0-2:
    1. Sequence Availability (0-2)
    2. Metadata Completeness (0-2)
    3. RAM Site Observability (0-2)
    4. Temporal Coverage (0-2)
    5. Subtype Diversity (0-2)

    Total score: 0-10
    """

    # Dimension score thresholds (raw metric -> score 0, 1, 2)
    # These will be applied based on thresholds in config
    DIMENSION_THRESHOLDS = {
        'sequence_availability': {
            0: ('low', 100),      # < 100 sequences = score 0
            1: ('medium', 1000),  # 100-1000 = score 1
            2: ('high', None),    # > 1000 = score 2
        },
        'metadata_completeness': {
            0: ('low', 0.3),
            1: ('medium', 0.7),
            2: ('high', None),
        },
        'ram_site_observability': {
            0: ('low', 0.5),
            1: ('medium', 0.8),
            2: ('high', None),
        },
        'temporal_coverage': {
            0: ('low', 2),        # < 2 years
            1: ('medium', 5),     # 2-5 years
            2: ('high', None),    # > 5 years
        },
        'subtype_diversity': {
            0: ('low', 2),        # < 2 subtypes
            1: ('medium', 5),     # 2-5 subtypes
            2: ('high', None),    # > 5 subtypes
        },
    }

    # Priority tier thresholds based on CIR-SRI score
    PRIORITY_TIERS = {
        'urgent': (0, 3),
        'high': (3, 5),
        'medium': (5, 7),
        'low': (7, 11),  # 11 to include score of 10
    }

    def __init__(self, cso_results: List[CSOResult], ros_results: List[ROSResult]):
        """
        Initialize CIR-SRI builder with CSO and ROS results.

        Args:
            cso_results: List of CSOResult objects
            ros_results: List of ROSResult objects
        """
        self.cso_results = {r.country: r for r in cso_results}
        self.ros_results = {r.country: r for r in ros_results}

    def _score_sequence_availability(self, n_sequences: int) -> Tuple[float, str]:
        """
        Score sequence availability dimension.

        Args:
            n_sequences: Number of sequences

        Returns:
            Tuple of (score, category)
        """
        if n_sequences < 100:
            return 0.0, 'low'
        elif n_sequences < 1000:
            return 1.0, 'medium'
        else:
            return 2.0, 'high'

    def _score_metadata_completeness(self, completeness: float) -> Tuple[float, str]:
        """
        Score metadata completeness dimension.

        Args:
            completeness: Fraction of metadata fields complete (0-1)

        Returns:
            Tuple of (score, category)
        """
        if completeness < 0.3:
            return 0.0, 'low'
        elif completeness < 0.7:
            return 1.0, 'medium'
        else:
            return 2.0, 'high'

    def _score_ram_site_observability(self, ros_score: float) -> Tuple[float, str]:
        """
        Score RAM site observability dimension.

        Args:
            ros_score: ROS score (0-100)

        Returns:
            Tuple of (score, category)
        """
        # ROS is 0-100, map to 0-2 scale
        if ros_score < 50:
            return 0.0, 'low'
        elif ros_score < 80:
            return 1.0, 'medium'
        else:
            return 2.0, 'high'

    def _score_temporal_coverage(self, n_years: int, recent_year: int) -> Tuple[float, str]:
        """
        Score temporal coverage dimension.

        Args:
            n_years: Number of years with data
            recent_year: Most recent year with data

        Returns:
            Tuple of (score, category)
        """
        current_year = 2024  # Should be dynamically determined

        # Primary: number of years
        if n_years < 2:
            return 0.0, 'low'
        elif n_years < 5:
            return 1.0, 'medium'
        else:
            # Bonus for recency
            years_since_recent = current_year - recent_year if recent_year else 999
            if years_since_recent <= 3:
                return 2.0, 'high'
            else:
                return 1.5, 'medium'

    def _score_subtype_diversity(self, n_subtypes: int) -> Tuple[float, str]:
        """
        Score subtype diversity dimension.

        Args:
            n_subtypes: Number of distinct subtypes observed

        Returns:
            Tuple of (score, category)
        """
        if n_subtypes < 2:
            return 0.0, 'low'
        elif n_subtypes < 5:
            return 1.0, 'medium'
        else:
            return 2.0, 'high'

    def calculate_country_cir_sri(self, country: str) -> Optional[CIRSRICountryResult]:
        """
        Calculate complete CIR-SRI for a country.

        Args:
            country: Country code

        Returns:
            CIRSRICountryResult or None if no data
        """
        cso_result = self.cso_results.get(country)
        ros_result = self.ros_results.get(country)

        if not cso_result and not ros_result:
            return None

        result = CIRSRICountryResult(
            country=country,
            region=config.COUNTRY_TO_REGION.get(country, 'UNKNOWN')
        )

        # Get base metrics
        n_sequences = cso_result.n_sequences if cso_result else 0
        metadata_completeness = (
            cso_result.metadata_completeness_score if cso_result else 0.0
        )
        ros_score = ros_result.ros_score if ros_result else 0.0
        years_covered = cso_result.years_covered if cso_result else []
        n_years = len(years_covered)
        recent_year = max(years_covered) if years_covered else None
        subtypes_observed = cso_result.subtypes_observed if cso_result else []
        n_subtypes = len(subtypes_observed)

        # Calculate dimension scores
        seq_score, seq_cat = self._score_sequence_availability(n_sequences)
        result.sequence_availability = DimensionScore(
            'sequence_availability', seq_score, seq_score, 0, seq_cat
        )

        meta_score, meta_cat = self._score_metadata_completeness(metadata_completeness)
        result.metadata_completeness = DimensionScore(
            'metadata_completeness', meta_score, meta_score, 0, meta_cat
        )

        ram_score, ram_cat = self._score_ram_site_observability(ros_score)
        result.ram_site_observability = DimensionScore(
            'ram_site_observability', ram_score, ram_score, 0, ram_cat
        )

        temp_score, temp_cat = self._score_temporal_coverage(n_years, recent_year)
        result.temporal_coverage = DimensionScore(
            'temporal_coverage', temp_score, temp_score, 0, temp_cat
        )

        subtype_score, subtype_cat = self._score_subtype_diversity(n_subtypes)
        result.subtype_diversity = DimensionScore(
            'subtype_diversity', subtype_score, subtype_score, 0, subtype_cat
        )

        # Calculate composite scores
        result.cso_score = cso_result.cso_score if cso_result else 0.0
        result.ros_score = ros_score
        result.cir_sri_score = sum([
            result.sequence_availability.raw_score,
            result.metadata_completeness.raw_score,
            result.ram_site_observability.raw_score,
            result.temporal_coverage.raw_score,
            result.subtype_diversity.raw_score,
        ])

        result.n_sequences = n_sequences

        # Identify critical gaps
        if ros_result:
            result.critical_gaps = ros_result.missing_sites.copy()
            if result.cir_sri_score < 5:
                result.critical_gaps.append('overall_low_readiness')

        # Assign priority tier
        for tier, (min_score, max_score) in self.PRIORITY_TIERS.items():
            if min_score <= result.cir_sri_score < max_score:
                result.priority_tier = tier
                break

        return result

    def calculate_all_cir_sri(self) -> List[CIRSRICountryResult]:
        """
        Calculate CIR-SRI for all countries with data.

        Returns:
            List of CIRSRICountryResult objects
        """
        # Get all countries with any data
        all_countries = set(self.cso_results.keys()) | set(self.ros_results.keys())

        results = []
        for country in all_countries:
            result = self.calculate_country_cir_sri(country)
            if result:
                results.append(result)

        # Calculate percentiles
        all_scores = [r.cir_sri_score for r in results]
        for result in results:
            result.cir_sri_percentile = stats.percentileofscore(all_scores, result.cir_sri_score)

        # Calculate regional ranks
        by_region = defaultdict(list)
        for result in results:
            by_region[result.region].append(result)

        for region, region_results in by_region.items():
            sorted_results = sorted(region_results, key=lambda r: r.cir_sri_score, reverse=True)
            for rank, result in enumerate(sorted_results, 1):
                result.regional_rank = rank

        # Calculate global ranks
        sorted_results = sorted(results, key=lambda r: r.cir_sri_score, reverse=True)
        for rank, result in enumerate(sorted_results, 1):
            result.global_rank = rank

        logger.info(f"Calculated CIR-SRI for {len(results)} countries")

        return results

    def calculate_regional_summary(self, results: List[CIRSRICountryResult]) -> Dict[str, Dict]:
        """
        Calculate regional CIR-SRI summaries.

        Args:
            results: List of CIRSRICountryResult objects

        Returns:
            Dictionary of regional summaries
        """
        regional_summaries = {}

        for region in config.WHO_REGIONS.keys():
            region_results = [r for r in results if r.region == region]

            if not region_results:
                continue

            scores = [r.cir_sri_score for r in region_results]

            # Count by priority tier
            tier_counts = defaultdict(int)
            for result in region_results:
                tier_counts[result.priority_tier] += 1

            # Dimension averages
            dim_averages = {}
            for dim_name in ['sequence_availability', 'metadata_completeness',
                           'ram_site_observability', 'temporal_coverage', 'subtype_diversity']:
                dim_scores = [getattr(r, dim_name).raw_score for r in region_results]
                dim_averages[dim_name] = np.mean(dim_scores)

            regional_summaries[region] = {
                'n_countries': len(region_results),
                'n_total_countries': len(config.WHO_REGIONS[region]),
                'mean_cir_sri': np.mean(scores),
                'median_cir_sri': np.median(scores),
                'std_cir_sri': np.std(scores),
                'min_cir_sri': min(scores),
                'max_cir_sri': max(scores),
                'tier_counts': dict(tier_counts),
                'dimension_averages': dim_averages,
            }

        return regional_summaries

    def generate_priority_list(self, results: List[CIRSRICountryResult]) -> pd.DataFrame:
        """
        Generate prioritized country list for surveillance improvement.

        Args:
            results: List of CIRSRICountryResult objects

        Returns:
            DataFrame sorted by priority
        """
        rows = []
        for result in results:
            row = result.to_dict()
            row['priority_score'] = (
                (10 - result.cir_sri_score) * 10 +
                len(result.critical_gaps) * 5 +
                (result.regional_rank / max(1, result.n_sequences)) * 0.1
            )
            rows.append(row)

        df = pd.DataFrame(rows)
        df = df.sort_values('priority_score', ascending=False)

        return df

    def save_results(self, results: List[CIRSRICountryResult], output_path: str):
        """
        Save CIR-SRI results to CSV.

        Args:
            results: List of CIRSRICountryResult objects
            output_path: Output file path
        """
        df = pd.DataFrame([r.to_dict() for r in results])
        df.to_csv(output_path, index=False)
        logger.info(f"Saved CIR-SRI results to {output_path}")

    def save_dimension_details(self, results: List[CIRSRICountryResult], output_path: str):
        """
        Save detailed dimension scores to CSV.

        Args:
            results: List of CIRSRICountryResult objects
            output_path: Output file path
        """
        rows = []
        for result in results:
            for dim_name in ['sequence_availability', 'metadata_completeness',
                           'ram_site_observability', 'temporal_coverage', 'subtype_diversity']:
                dim = getattr(result, dim_name)
                rows.append({
                    'country': result.country,
                    'region': result.region,
                    'dimension': dim.dimension_name,
                    'raw_score': dim.raw_score,
                    'category': dim.category,
                })

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved dimension details to {output_path}")

    def save_priority_list(self, df: pd.DataFrame, output_path: str):
        """
        Save priority list to CSV.

        Args:
            df: Priority DataFrame
            output_path: Output file path
        """
        df.to_csv(output_path, index=False)
        logger.info(f"Saved priority list to {output_path}")

    def save_regional_summary(self, summaries: Dict, output_path: str):
        """
        Save regional summaries to CSV.

        Args:
            summaries: Regional summaries dictionary
            output_path: Output file path
        """
        rows = []
        for region, stats in summaries.items():
            row = {
                'region': region,
                'n_countries': stats['n_countries'],
                'n_total_countries': stats['n_total_countries'],
                'coverage_pct': stats['n_countries'] / stats['n_total_countries'] * 100,
                'mean_cir_sri': stats['mean_cir_sri'],
                'median_cir_sri': stats['median_cir_sri'],
                'std_cir_sri': stats['std_cir_sri'],
                'min_cir_sri': stats['min_cir_sri'],
                'max_cir_sri': stats['max_cir_sri'],
                'urgent_tier_n': stats['tier_counts'].get('urgent', 0),
                'high_tier_n': stats['tier_counts'].get('high', 0),
                'medium_tier_n': stats['tier_counts'].get('medium', 0),
                'low_tier_n': stats['tier_counts'].get('low', 0),
            }

            # Add dimension averages
            for dim_name, dim_avg in stats['dimension_averages'].items():
                row[f'dim_{dim_name}_avg'] = dim_avg

            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved regional summary to {output_path}")

    def generate_summary_report(self, results: List[CIRSRICountryResult],
                               regional_summaries: Dict) -> str:
        """
        Generate a text summary report.

        Args:
            results: List of CIRSRICountryResult objects
            regional_summaries: Regional summaries dictionary

        Returns:
            Formatted report string
        """
        report_lines = [
            "=" * 80,
            "CIR-SRI (Capsid Inhibitor Resistance Surveillance Readiness Index) Summary Report",
            "=" * 80,
            "",
            f"Total countries analyzed: {len(results)}",
            "",
        ]

        # Global statistics
        scores = [r.cir_sri_score for r in results]
        report_lines.extend([
            "Global Statistics:",
            f"  Mean CIR-SRI: {np.mean(scores):.2f}",
            f"  Median CIR-SRI: {np.median(scores):.2f}",
            f"  Std Dev: {np.std(scores):.2f}",
            f"  Range: {min(scores):.1f} - {max(scores):.1f}",
            "",
        ])

        # Priority tier distribution
        tier_counts = defaultdict(int)
        for result in results:
            tier_counts[result.priority_tier] += 1

        report_lines.extend([
            "Priority Tier Distribution:",
            f"  Urgent (0-3): {tier_counts.get('urgent', 0)} countries",
            f"  High (3-5): {tier_counts.get('high', 0)} countries",
            f"  Medium (5-7): {tier_counts.get('medium', 0)} countries",
            f"  Low (7-10): {tier_counts.get('low', 0)} countries",
            "",
        ])

        # Regional summaries
        report_lines.append("Regional Summaries:")
        report_lines.append("-" * 40)
        for region, stats in regional_summaries.items():
            report_lines.append(f"\n{region}:")
            report_lines.append(f"  Countries: {stats['n_countries']}/{stats['n_total_countries']}")
            report_lines.append(f"  Mean CIR-SRI: {stats['mean_cir_sri']:.2f}")
            report_lines.append(f"  Tier distribution: Urgent={stats['tier_counts'].get('urgent', 0)}, "
                              f"High={stats['tier_counts'].get('high', 0)}, "
                              f"Medium={stats['tier_counts'].get('medium', 0)}, "
                              f"Low={stats['tier_counts'].get('low', 0)}")

        # Top and bottom countries
        sorted_results = sorted(results, key=lambda r: r.cir_sri_score, reverse=True)

        report_lines.extend([
            "",
            "Top 10 Countries by CIR-SRI:",
            "-" * 40,
        ])
        for result in sorted_results[:10]:
            report_lines.append(f"  {result.country:5s} ({result.region:4s}): {result.cir_sri_score:.1f} "
                              f"[{result.priority_tier}]")

        report_lines.extend([
            "",
            "Bottom 10 Countries by CIR-SRI:",
            "-" * 40,
        ])
        for result in sorted_results[-10:]:
            report_lines.append(f"  {result.country:5s} ({result.region:4s}): {result.cir_sri_score:.1f} "
                              f"[{result.priority_tier}]")

        report_lines.append("")
        report_lines.append("=" * 80)

        return "\n".join(report_lines)


def main():
    """Example usage of CIR-SRI builder."""
    import argparse

    parser = argparse.ArgumentParser(description='Calculate CIR-SRI')
    parser.add_argument('cso_csv', help='CSO scores CSV file')
    parser.add_argument('ros_csv', help='ROS scores CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    args = parser.parse_args()

    # Load CSO results (simplified conversion)
    cso_df = pd.read_csv(args.cso_csv)
    cso_results = []
    for _, row in cso_df.iterrows():
        cso_result = CSOResult(
            country=row['country'],
            region=config.COUNTRY_TO_REGION.get(row['country'], 'UNKNOWN'),
            n_sequences=row.get('n_sequences', 0),
            cso_score=row.get('cso_score', 0),
            metadata_completeness_score=row.get('metadata_completeness_score', 0),
            years_covered=[int(y) for y in str(row.get('years_covered', '')).split(',') if y.strip()],
            subtypes_observed=str(row.get('subtypes_observed', '')).split(','),
        )
        cso_results.append(cso_result)

    # Load ROS results (simplified conversion)
    ros_df = pd.read_csv(args.ros_csv)
    ros_results = []
    for _, row in ros_df.iterrows():
        ros_result = ROSResult(
            country=row['country'],
            region=config.COUNTRY_TO_REGION.get(row['country'], 'UNKNOWN'),
            ros_score=row.get('ros_score', 0),
            n_sequences=row.get('n_sequences', 0),
        )
        ros_results.append(ros_result)

    # Build CIR-SRI
    builder = CIRSRIBuilder(cso_results, ros_results)
    results = builder.calculate_all_cir_sri()
    regional_summaries = builder.calculate_regional_summary(results)
    priority_df = builder.generate_priority_list(results)

    # Save results
    os.makedirs(args.output, exist_ok=True)
    builder.save_results(results, os.path.join(args.output, 'cir_sri_scores.csv'))
    builder.save_dimension_details(results, os.path.join(args.output, 'dimension_details.csv'))
    builder.save_priority_list(priority_df, os.path.join(args.output, 'priority_list.csv'))
    builder.save_regional_summary(regional_summaries, os.path.join(args.output, 'regional_summary.csv'))

    # Generate and save report
    report = builder.generate_summary_report(results, regional_summaries)
    with open(os.path.join(args.output, 'cir_sri_report.txt'), 'w') as f:
        f.write(report)

    print(report)


if __name__ == '__main__':
    main()