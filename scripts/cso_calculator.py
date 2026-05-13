"""
Capsid Sequence Observability (CSO) Calculator

This module calculates the CSO score for countries/regions, measuring:
- Sequence availability and coverage
- Metadata completeness
- Temporal trend analysis
- Subtype coverage assessment

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
import json

import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d

from . import config
from .sequence_qc import SequenceRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class CSOResult:
    """Container for CSO calculation results."""
    country: str
    region: str
    n_sequences: int = 0
    n_countries_in_region: int = 0

    # Component scores (0-1 scale)
    sequence_availability_score: float = 0.0
    metadata_completeness_score: float = 0.0
    temporal_coverage_score: float = 0.0
    subtype_diversity_score: float = 0.0
    regional_proportion_score: float = 0.0

    # Overall CSO (weighted average)
    cso_score: float = 0.0

    # Detailed metrics
    years_covered: List[int] = field(default_factory=list)
    subtypes_observed: List[str] = field(default_factory=list)
    metadata_fields_complete: Dict[str, float] = field(default_factory=dict)
    sequence_length_mean: float = 0.0
    sequence_length_std: float = 0.0

    # Percentile ranks
    cso_percentile: float = 0.0
    regional_rank: int = 0

    def to_dict(self) -> Dict:
        """Convert to dictionary for export."""
        return {
            'country': self.country,
            'region': self.region,
            'n_sequences': self.n_sequences,
            'sequence_availability_score': self.sequence_availability_score,
            'metadata_completeness_score': self.metadata_completeness_score,
            'temporal_coverage_score': self.temporal_coverage_score,
            'subtype_diversity_score': self.subtype_diversity_score,
            'regional_proportion_score': self.regional_proportion_score,
            'cso_score': self.cso_score,
            'cso_percentile': self.cso_percentile,
            'regional_rank': self.regional_rank,
            'years_covered': ','.join(map(str, self.years_covered)),
            'subtypes_observed': ','.join(self.subtypes_observed),
        }


class CSOCalculator:
    """
    Calculate Capsid Sequence Observability scores for countries.

    The CSO score quantifies how well a country's sequence data supports
    capsid inhibitor surveillance, based on multiple dimensions.
    """

    # Component weights for final CSO score
    COMPONENT_WEIGHTS = {
        'sequence_availability': 0.30,
        'metadata_completeness': 0.25,
        'temporal_coverage': 0.20,
        'subtype_diversity': 0.15,
        'regional_proportion': 0.10,
    }

    # Reference values for normalization
    REFERENCE_VALUES = {
        'max_sequences_per_country': 10000,
        'min_years_for_full_temporal': 10,
        'min_subtypes_for_full_diversity': 10,
        'max_regional_proportion': 0.5,  # Max 50% of region's sequences
    }

    def __init__(self, records: List[SequenceRecord]):
        """
        Initialize CSO calculator with sequence data.

        Args:
            records: List of SequenceRecord objects from QC pipeline
        """
        self.records = records
        self.country_data = self._aggregate_by_country()
        self.regional_totals = self._calculate_regional_totals()

    def _aggregate_by_country(self) -> Dict[str, Dict]:
        """
        Aggregate sequence data by country.

        Returns:
            Dictionary with country-level statistics
        """
        country_data = defaultdict(lambda: {
            'sequences': [],
            'years': [],
            'subtypes': set(),
            'metadata_fields': {
                'country': 0,
                'year': 0,
                'subtype': 0,
                'date': 0,
            },
            'sequence_lengths': [],
        })

        for record in self.records:
            country = record.country or 'UNKNOWN'

            country_data[country]['sequences'].append(record)
            country_data[country]['sequence_lengths'].append(record.length)

            if record.year:
                country_data[country]['years'].append(record.year)
            if record.subtype:
                country_data[country]['subtypes'].add(record.subtype)

            # Track metadata completeness
            if record.country:
                country_data[country]['metadata_fields']['country'] += 1
            if record.year:
                country_data[country]['metadata_fields']['year'] += 1
            if record.subtype:
                country_data[country]['metadata_fields']['subtype'] += 1
            if record.collection_date:
                country_data[country]['metadata_fields']['date'] += 1

        return dict(country_data)

    def _calculate_regional_totals(self) -> Dict[str, Dict]:
        """
        Calculate regional sequence totals.

        Returns:
            Dictionary with region-level statistics
        """
        regional_totals = defaultdict(lambda: {
            'n_countries': 0,
            'n_sequences': 0,
            'countries': set(),
        })

        for country, data in self.country_data.items():
            region = config.COUNTRY_TO_REGION.get(country, 'UNKNOWN')
            if region != 'UNKNOWN':
                regional_totals[region]['countries'].add(country)
                regional_totals[region]['n_countries'] += 1
                regional_totals[region]['n_sequences'] += len(data['sequences'])

        return dict(regional_totals)

    def _calculate_sequence_availability_score(self, country: str) -> float:
        """
        Calculate sequence availability score (0-1).

        Scores based on number of sequences relative to maximum observed.

        Args:
            country: Country code

        Returns:
            Score from 0 to 1
        """
        if country not in self.country_data:
            return 0.0

        n_sequences = len(self.country_data[country]['sequences'])
        n_total = sum(len(d['sequences']) for d in self.country_data.values())

        if n_total == 0:
            return 0.0

        # Proportion of global sequences
        proportion = n_sequences / n_total

        # Scale by maximum observed proportion
        max_proportion = max(
            len(d['sequences']) / n_total
            for d in self.country_data.values()
        ) if self.country_data else 1.0

        if max_proportion > 0:
            score = min(1.0, proportion / max_proportion)
        else:
            score = 0.0

        return score

    def _calculate_metadata_completeness_score(self, country: str) -> float:
        """
        Calculate metadata completeness score (0-1).

        Args:
            country: Country code

        Returns:
            Score from 0 to 1
        """
        if country not in self.country_data:
            return 0.0

        data = self.country_data[country]
        n_sequences = len(data['sequences'])

        if n_sequences == 0:
            return 0.0

        # Calculate completeness for each field
        completeness = {}
        for field_name, count in data['metadata_fields'].items():
            completeness[field_name] = count / n_sequences

        # Weighted average of fields
        weights = {'country': 0.3, 'year': 0.4, 'subtype': 0.2, 'date': 0.1}
        score = sum(completeness.get(f, 0) * w for f, w in weights.items())

        return min(1.0, score)

    def _calculate_temporal_coverage_score(self, country: str) -> float:
        """
        Calculate temporal coverage score (0-1).

        Based on:
        - Number of years with data
        - Consistency of sampling
        - Recency of data

        Args:
            country: Country code

        Returns:
            Score from 0 to 1
        """
        if country not in self.country_data:
            return 0.0

        years = self.country_data[country]['years']
        if not years:
            return 0.0

        unique_years = sorted(set(years))
        n_years = len(unique_years)

        # Score for number of years
        max_years = self.REFERENCE_VALUES['min_years_for_full_temporal']
        year_score = min(1.0, n_years / max_years)

        # Score for recency
        current_year = 2024  # Should be dynamically determined
        max_year_gap = max(current_year - max(years), 0) if unique_years else current_year
        recency_score = max(0.0, 1.0 - (max_year_gap / 10))

        # Score for consistency (coefficient of variation)
        if len(years) > 1:
            year_counts = pd.Series(years).value_counts()
            # Ideal: uniform distribution
            expected_per_year = len(years) / n_years
            cv = year_counts.std() / expected_per_year if expected_per_year > 0 else 0
            consistency_score = max(0.0, 1.0 - cv)
        else:
            consistency_score = 0.5

        # Combine scores
        score = (year_score * 0.4 + recency_score * 0.4 + consistency_score * 0.2)

        return min(1.0, score)

    def _calculate_subtype_diversity_score(self, country: str) -> float:
        """
        Calculate subtype diversity score (0-1).

        Args:
            country: Country code

        Returns:
            Score from 0 to 1
        """
        if country not in self.country_data:
            return 0.0

        subtypes = self.country_data[country]['subtypes']
        n_subtypes = len(subtypes)

        # Reference: globally observed subtypes
        all_subtypes = set()
        for data in self.country_data.values():
            all_subtypes.update(data['subtypes'])

        global_n_subtypes = len(all_subtypes) if all_subtypes else 1

        # Score based on proportion of global diversity
        diversity_score = min(1.0, n_subtypes / global_n_subtypes)

        # Bonus for presence of common subtypes
        common_subtypes = {'B', 'C', 'AE', 'AG', 'G', 'D', 'F', 'H', 'J', 'K'}
        common_observed = subtypes & common_subtypes
        common_bonus = len(common_observed) / len(common_subtypes) * 0.2

        score = diversity_score * 0.8 + common_bonus

        return min(1.0, score)

    def _calculate_regional_proportion_score(self, country: str) -> float:
        """
        Calculate regional proportion score (0-1).

        Measures country's contribution to regional sequence pool.

        Args:
            country: Country code

        Returns:
            Score from 0 to 1
        """
        if country not in self.country_data:
            return 0.0

        region = config.COUNTRY_TO_REGION.get(country)
        if not region or region not in self.regional_totals:
            return 0.0

        country_n = len(self.country_data[country]['sequences'])
        region_n = self.regional_totals[region]['n_sequences']

        if region_n == 0:
            return 0.0

        proportion = country_n / region_n
        max_prop = self.REFERENCE_VALUES['max_regional_proportion']

        score = min(1.0, proportion / max_prop)

        return score

    def calculate_country_cso(self, country: str) -> CSOResult:
        """
        Calculate complete CSO for a single country.

        Args:
            country: Country code

        Returns:
            CSOResult object with all scores
        """
        result = CSOResult(
            country=country,
            region=config.COUNTRY_TO_REGION.get(country, 'UNKNOWN')
        )

        if country not in self.country_data:
            return result

        data = self.country_data[country]
        result.n_sequences = len(data['sequences'])
        result.n_countries_in_region = len(
            self.regional_totals.get(result.region, {}).get('countries', set())
        )

        # Calculate component scores
        result.sequence_availability_score = self._calculate_sequence_availability_score(country)
        result.metadata_completeness_score = self._calculate_metadata_completeness_score(country)
        result.temporal_coverage_score = self._calculate_temporal_coverage_score(country)
        result.subtype_diversity_score = self._calculate_subtype_diversity_score(country)
        result.regional_proportion_score = self._calculate_regional_proportion_score(country)

        # Calculate overall CSO
        components = [
            ('sequence_availability', result.sequence_availability_score),
            ('metadata_completeness', result.metadata_completeness_score),
            ('temporal_coverage', result.temporal_coverage_score),
            ('subtype_diversity', result.subtype_diversity_score),
            ('regional_proportion', result.regional_proportion_score),
        ]

        result.cso_score = sum(
            score * self.COMPONENT_WEIGHTS[component]
            for component, score in components
        )

        # Store detailed metrics
        result.years_covered = sorted(set(data['years']))
        result.subtypes_observed = sorted(list(data['subtypes']))
        result.sequence_length_mean = np.mean(data['sequence_lengths'])
        result.sequence_length_std = np.std(data['sequence_lengths'])

        # Metadata completeness details
        n = result.n_sequences
        for field_name, count in data['metadata_fields'].items():
            result.metadata_fields_complete[field_name] = count / n if n > 0 else 0

        return result

    def calculate_all_cso(self) -> List[CSOResult]:
        """
        Calculate CSO for all countries with sequences.

        Returns:
            List of CSOResult objects
        """
        results = []

        for country in self.country_data.keys():
            result = self.calculate_country_cso(country)
            results.append(result)

        # Calculate percentile ranks
        all_scores = [r.cso_score for r in results]
        for result in results:
            result.cso_percentile = stats.percentileofscore(all_scores, result.cso_score)

        # Calculate regional ranks
        by_region = defaultdict(list)
        for result in results:
            by_region[result.region].append(result)

        for region, region_results in by_region.items():
            sorted_results = sorted(region_results, key=lambda r: r.cso_score, reverse=True)
            for rank, result in enumerate(sorted_results, 1):
                result.regional_rank = rank

        # Sort by overall CSO score
        results.sort(key=lambda r: r.cso_score, reverse=True)

        logger.info(f"Calculated CSO for {len(results)} countries")

        return results

    def calculate_regional_cso(self) -> Dict[str, Dict]:
        """
        Calculate aggregated CSO by region.

        Returns:
            Dictionary of regional statistics
        """
        regional_stats = {}

        for region, countries in config.WHO_REGIONS.items():
            region_records = [
                r for r in self.records
                if config.COUNTRY_TO_REGION.get(r.country) == region
            ]

            if not region_records:
                continue

            # Aggregate metrics
            n_sequences = len(region_records)
            years = sorted(set(r.year for r in region_records if r.year))
            subtypes = set(r.subtype for r in region_records if r.subtype)

            # Calculate mean CSO for countries in region
            country_cso = [
                self.calculate_country_cso(c).cso_score
                for c in countries
                if c in self.country_data
            ]

            regional_stats[region] = {
                'n_countries': len([c for c in countries if c in self.country_data]),
                'n_total_countries': len(countries),
                'n_sequences': n_sequences,
                'n_countries_with_data': len(country_cso),
                'mean_cso': np.mean(country_cso) if country_cso else 0,
                'median_cso': np.median(country_cso) if country_cso else 0,
                'std_cso': np.std(country_cso) if country_cso else 0,
                'years_covered': years,
                'year_range': (min(years), max(years)) if years else (None, None),
                'n_subtypes': len(subtypes),
                'subtypes': sorted(subtypes),
            }

        return regional_stats

    def save_results(self, results: List[CSOResult], output_path: str):
        """
        Save CSO results to CSV file.

        Args:
            results: List of CSOResult objects
            output_path: Output file path
        """
        df = pd.DataFrame([r.to_dict() for r in results])
        df.to_csv(output_path, index=False)
        logger.info(f"Saved CSO results to {output_path}")

    def save_regional_summary(self, regional_stats: Dict, output_path: str):
        """
        Save regional summary to CSV file.

        Args:
            regional_stats: Regional statistics dictionary
            output_path: Output file path
        """
        rows = []
        for region, stats in regional_stats.items():
            row = {
                'region': region,
                'n_countries': stats['n_countries'],
                'n_total_countries': stats['n_total_countries'],
                'coverage_pct': stats['n_countries'] / stats['n_total_countries'] * 100,
                'n_sequences': stats['n_sequences'],
                'mean_cso': stats['mean_cso'],
                'median_cso': stats['median_cso'],
                'std_cso': stats['std_cso'],
                'year_min': stats['year_range'][0],
                'year_max': stats['year_range'][1],
                'n_subtypes': stats['n_subtypes'],
                'subtypes': ','.join(stats['subtypes']),
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved regional summary to {output_path}")

    def generate_temporal_trends(self, country: str) -> Dict[str, Any]:
        """
        Generate temporal trend analysis for a country.

        Args:
            country: Country code

        Returns:
            Dictionary with trend statistics
        """
        if country not in self.country_data:
            return {}

        data = self.country_data[country]
        years = data['years']

        if not years:
            return {}

        # Count sequences per year
        year_counts = pd.Series(years).value_counts().sort_index()

        # Calculate trend (linear regression)
        if len(year_counts) >= 2:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                year_counts.index, year_counts.values
            )
            trend_direction = 'increasing' if slope > 0 else 'decreasing' if slope < 0 else 'stable'
        else:
            slope = intercept = r_value = p_value = std_err = 0
            trend_direction = 'insufficient_data'

        return {
            'country': country,
            'year_counts': year_counts.to_dict(),
            'trend_slope': slope,
            'trend_direction': trend_direction,
            'r_squared': r_value ** 2,
            'p_value': p_value,
            'mean_annual_sequences': np.mean(year_counts.values),
            'std_annual_sequences': np.std(year_counts.values),
        }


def main():
    """Example usage of CSO calculator."""
    import argparse

    parser = argparse.ArgumentParser(description='Calculate CSO scores')
    parser.add_argument('input_csv', help='Processed sequences CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    args = parser.parse_args()

    # Load processed sequences
    df = pd.read_csv(args.input_csv)

    # Convert to SequenceRecord objects (simplified)
    records = []
    for _, row in df.iterrows():
        record = SequenceRecord(
            accession=row.get('accession', ''),
            sequence=row.get('sequence', ''),
            sequence_aa='',
            country=row.get('country'),
            year=row.get('year'),
            subtype=row.get('subtype'),
            length=row.get('length', 0),
        )
        records.append(record)

    # Calculate CSO
    calculator = CSOCalculator(records)
    results = calculator.calculate_all_cso()
    regional_stats = calculator.calculate_regional_cso()

    # Save results
    os.makedirs(args.output, exist_ok=True)
    calculator.save_results(results, os.path.join(args.output, 'cso_scores.csv'))
    calculator.save_regional_summary(regional_stats, os.path.join(args.output, 'regional_summary.csv'))

    print(f"\nCSO Summary:")
    print(f"  Countries analyzed: {len(results)}")
    print(f"  Top 5 countries by CSO:")
    for result in results[:5]:
        print(f"    {result.country}: {result.cso_score:.3f}")


if __name__ == '__main__':
    main()