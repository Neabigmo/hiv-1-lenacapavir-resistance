"""
RAM-site Observability Score (ROS) Calculator

This module calculates ROS scores measuring:
- RAM site coverage per country/region
- Per-sequence RAM observability
- Heatmap data generation for visualization
- Missing site identification

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Set
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np
import pandas as pd

from . import config
from .sequence_qc import SequenceRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class RAMSiteResult:
    """Container for RAM site observability result."""
    site_name: str
    hxb2_position: int
    ca_position: int

    # Observability metrics
    sequences_with_site: int = 0
    total_sequences: int = 0
    observability_score: float = 0.0

    # Quality metrics
    ambiguous_observations: int = 0
    stop_codon_at_site: int = 0
    amino_acid_counts: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'site_name': self.site_name,
            'hxb2_position': self.hxb2_position,
            'ca_position': self.ca_position,
            'sequences_with_site': self.sequences_with_site,
            'total_sequences': self.total_sequences,
            'observability_score': self.observability_score,
            'ambiguous_observations': self.ambiguous_observations,
            'stop_codon_at_site': self.stop_codon_at_site,
        }


@dataclass
class ROSResult:
    """Container for ROS calculation results."""
    country: str
    region: str

    # Overall ROS
    ros_score: float = 0.0

    # Site-level results
    site_results: List[RAMSiteResult] = field(default_factory=list)

    # Missing sites
    missing_sites: List[str] = field(default_factory=list)
    partially_observed_sites: List[str] = field(default_factory=list)
    fully_observed_sites: List[str] = field(default_factory=list)

    # Statistics
    n_sequences: int = 0
    n_sequences_full_ram_coverage: int = 0

    # Composite metrics
    coverage_completeness: float = 0.0
    site_quality_mean: float = 0.0

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'country': self.country,
            'region': self.region,
            'ros_score': self.ros_score,
            'n_sequences': self.n_sequences,
            'n_sequences_full_ram_coverage': self.n_sequences_full_ram_coverage,
            'full_coverage_pct': self.n_sequences_full_ram_coverage / max(1, self.n_sequences) * 100,
            'missing_sites': ','.join(self.missing_sites),
            'partially_observed_sites': ','.join(self.partially_observed_sites),
            'fully_observed_sites': ','.join(self.fully_observed_sites),
            'coverage_completeness': self.coverage_completeness,
            'site_quality_mean': self.site_quality_mean,
        }


class ROSCalculator:
    """
    Calculate RAM-site Observability Scores.

    ROS measures how well sequence data covers the specific resistance-
    associated mutation sites in the capsid protein.
    """

    # Codon position within RAM site (0-indexed, for 3-base codon)
    CODON_OFFSET = 0

    def __init__(self, records: List[SequenceRecord]):
        """
        Initialize ROS calculator.

        Args:
            records: List of SequenceRecord objects
        """
        self.records = records
        self.country_data = self._aggregate_by_country()

    def _aggregate_by_country(self) -> Dict[str, List[SequenceRecord]]:
        """Aggregate records by country."""
        country_data = defaultdict(list)
        for record in self.records:
            country = record.country or 'UNKNOWN'
            country_data[country].append(record)
        return dict(country_data)

    def _get_codon_at_position(self, sequence: str, position: int) -> Tuple[str, bool]:
        """
        Extract codon at a specific nucleotide position.

        Args:
            sequence: DNA sequence
            position: 0-indexed nucleotide position

        Returns:
            Tuple of (codon_string, is_ambiguous)
        """
        if position < 0 or position + 2 >= len(sequence):
            return '', True

        codon = sequence[position:position + 3].upper()

        # Check for ambiguous bases
        ambiguous_bases = {'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}
        is_ambiguous = any(base in codon for base in ambiguous_bases)

        # Validate it's a valid codon (only A, T, G, C)
        if not all(base in 'ATGC' for base in codon):
            is_ambiguous = True

        return codon, is_ambiguous

    def _codon_to_amino_acid(self, codon: str) -> str:
        """
        Convert DNA codon to amino acid.

        Args:
            codon: DNA codon string (3 bases)

        Returns:
            Single letter amino acid code, or 'X' if invalid
        """
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        }

        return codon_table.get(codon.upper(), 'X')

    def _analyze_sequence_for_sites(self, record: SequenceRecord) -> Dict[str, Dict]:
        """
        Analyze a single sequence for RAM site coverage.

        Args:
            record: SequenceRecord to analyze

        Returns:
            Dictionary of site_name -> site_analysis dict
        """
        site_analyses = {}

        for site_name, hxb2_pos in config.RAM_SITES.items():
            # Convert to 0-indexed
            nucleotide_pos = hxb2_pos - 1

            # Get codon at this position
            codon, is_ambiguous = self._get_codon_at_position(
                record.sequence, nucleotide_pos
            )

            if codon:
                aa = self._codon_to_amino_acid(codon)
                is_stop = aa == '*'
            else:
                aa = 'X'
                is_stop = False

            site_analyses[site_name] = {
                'observed': len(codon) == 3,
                'ambiguous': is_ambiguous,
                'amino_acid': aa,
                'is_stop': is_stop,
                'is_missing': len(codon) != 3,
            }

        return site_analyses

    def calculate_site_observability(self, country: str) -> List[RAMSiteResult]:
        """
        Calculate observability for each RAM site within a country.

        Args:
            country: Country code

        Returns:
            List of RAMSiteResult objects
        """
        if country not in self.country_data:
            return []

        records = self.country_data[country]
        site_results = []

        for site_name, hxb2_pos in config.RAM_SITES.items():
            ca_position = hxb2_pos - config.CA_START + 1

            result = RAMSiteResult(
                site_name=site_name,
                hxb2_position=hxb2_pos,
                ca_position=ca_position,
                total_sequences=len(records),
            )

            aa_counts = defaultdict(int)

            for record in records:
                site_analyses = self._analyze_sequence_for_sites(record)
                analysis = site_analyses.get(site_name, {})

                if analysis.get('observed') and not analysis.get('is_missing'):
                    result.sequences_with_site += 1

                    if analysis.get('ambiguous'):
                        result.ambiguous_observations += 1

                    if analysis.get('is_stop'):
                        result.stop_codon_at_site += 1

                    aa = analysis.get('amino_acid', 'X')
                    aa_counts[aa] += 1

            # Calculate observability score
            if result.total_sequences > 0:
                result.observability_score = result.sequences_with_site / result.total_sequences

            result.amino_acid_counts = dict(aa_counts)
            site_results.append(result)

        return site_results

    def calculate_country_ros(self, country: str) -> ROSResult:
        """
        Calculate complete ROS for a country.

        Args:
            country: Country code

        Returns:
            ROSResult object
        """
        result = ROSResult(
            country=country,
            region=config.COUNTRY_TO_REGION.get(country, 'UNKNOWN')
        )

        if country not in self.country_data:
            return result

        records = self.country_data[country]
        result.n_sequences = len(records)

        # Calculate site-level observability
        site_results = self.calculate_site_observability(country)
        result.site_results = site_results

        # Identify missing, partially observed, and fully observed sites
        for site_result in site_results:
            if site_result.sequences_with_site == 0:
                result.missing_sites.append(site_result.site_name)
            elif site_result.observability_score < 0.8:
                result.partially_observed_sites.append(site_result.site_name)
            else:
                result.fully_observed_sites.append(site_result.site_name)

        # Calculate overall ROS score
        # Weight: mean observability * coverage completeness
        if site_results:
            mean_observability = np.mean([s.observability_score for s in site_results])
            site_quality = np.mean([
                s.observability_score * (1 - s.ambiguous_observations / max(1, s.sequences_with_site))
                for s in site_results
            ])

            result.site_quality_mean = site_quality
            result.coverage_completeness = len(result.fully_observed_sites) / len(site_results)

            # ROS = coverage * quality * 100
            result.ros_score = result.coverage_completeness * site_quality * 100
        else:
            result.ros_score = 0.0

        # Count sequences with full RAM coverage
        for record in records:
            site_analyses = self._analyze_sequence_for_sites(record)
            all_sites_observed = all(
                analysis.get('observed') and not analysis.get('is_stop')
                for analysis in site_analyses.values()
            )
            if all_sites_observed:
                result.n_sequences_full_ram_coverage += 1

        return result

    def calculate_all_ros(self) -> List[ROSResult]:
        """
        Calculate ROS for all countries.

        Returns:
            List of ROSResult objects
        """
        results = []

        for country in self.country_data.keys():
            result = self.calculate_country_ros(country)
            results.append(result)

        # Sort by ROS score
        results.sort(key=lambda r: r.ros_score, reverse=True)

        logger.info(f"Calculated ROS for {len(results)} countries")

        return results

    def calculate_regional_ros(self) -> Dict[str, Dict]:
        """
        Calculate aggregated ROS by region.

        Returns:
            Dictionary of regional ROS statistics
        """
        regional_stats = {}

        all_results = self.calculate_all_ros()

        for region in config.WHO_REGIONS.keys():
            region_results = [r for r in all_results if r.region == region]

            if not region_results:
                continue

            # Aggregate statistics
            ros_scores = [r.ros_score for r in region_results]
            ros_scores_valid = [s for s in ros_scores if s > 0]

            # Aggregate site coverage
            all_missing = set()
            all_partially_observed = set()
            all_fully_observed = set()

            for result in region_results:
                all_missing.update(result.missing_sites)
                all_partially_observed.update(result.partially_observed_sites)
                all_fully_observed.update(result.fully_observed_sites)

            regional_stats[region] = {
                'n_countries': len(region_results),
                'n_countries_with_data': len(ros_scores_valid),
                'mean_ros': np.mean(ros_scores_valid) if ros_scores_valid else 0,
                'median_ros': np.median(ros_scores_valid) if ros_scores_valid else 0,
                'std_ros': np.std(ros_scores_valid) if ros_scores_valid else 0,
                'min_ros': min(ros_scores_valid) if ros_scores_valid else 0,
                'max_ros': max(ros_scores_valid) if ros_scores_valid else 0,
                'missing_sites': sorted(all_missing),
                'partially_observed_sites': sorted(all_partially_observed),
                'fully_observed_sites': sorted(all_fully_observed),
            }

        return regional_stats

    def generate_heatmap_data(self, results: List[ROSResult]) -> pd.DataFrame:
        """
        Generate heatmap data for visualization.

        Args:
            results: List of ROSResult objects

        Returns:
            DataFrame with site x country matrix
        """
        site_names = list(config.RAM_SITES.keys())
        countries = sorted(set(r.country for r in results))

        # Create matrix
        matrix = np.zeros((len(countries), len(site_names)))

        for i, country in enumerate(countries):
            country_result = next((r for r in results if r.country == country), None)
            if country_result:
                for j, site_name in enumerate(site_names):
                    site_result = next(
                        (s for s in country_result.site_results if s.site_name == site_name),
                        None
                    )
                    if site_result:
                        matrix[i, j] = site_result.observability_score * 100

        df = pd.DataFrame(matrix, index=countries, columns=site_names)
        return df

    def generate_regional_heatmap_data(self) -> pd.DataFrame:
        """
        Generate regional-level heatmap data.

        Returns:
            DataFrame with site x region matrix
        """
        site_names = list(config.RAM_SITES.keys())
        regions = list(config.WHO_REGIONS.keys())

        matrix = np.zeros((len(regions), len(site_names)))

        for i, region in enumerate(regions):
            region_records = [
                r for r in self.records
                if config.COUNTRY_TO_REGION.get(r.country) == region
            ]

            if region_records:
                # Calculate aggregated observability for this region
                for j, site_name in enumerate(site_names):
                    hxb2_pos = config.RAM_SITES[site_name]
                    nucleotide_pos = hxb2_pos - 1

                    observed_count = 0
                    for record in region_records:
                        codon, is_ambiguous = self._get_codon_at_position(
                            record.sequence, nucleotide_pos
                        )
                        if codon and len(codon) == 3 and not is_ambiguous:
                            observed_count += 1

                    matrix[i, j] = observed_count / len(region_records) * 100

        df = pd.DataFrame(matrix, index=regions, columns=site_names)
        return df

    def identify_critical_gaps(self, results: List[ROSResult],
                               threshold: float = 50.0) -> List[Dict]:
        """
        Identify critical gaps in RAM site coverage.

        Args:
            results: List of ROSResult objects
            threshold: Observability threshold for critical gaps

        Returns:
            List of gap analysis dictionaries
        """
        gaps = []

        for site_name in config.RAM_SITES.keys():
            # Count countries with low coverage
            low_coverage_countries = []
            no_coverage_countries = []

            for result in results:
                site_result = next(
                    (s for s in result.site_results if s.site_name == site_name),
                    None
                )
                if site_result:
                    if site_result.observability_score * 100 < threshold:
                        if site_result.observability_score == 0:
                            no_coverage_countries.append(result.country)
                        else:
                            low_coverage_countries.append(result.country)

            if no_coverage_countries or low_coverage_countries:
                gaps.append({
                    'site': site_name,
                    'hxb2_position': config.RAM_SITES[site_name],
                    'no_coverage_countries': no_coverage_countries,
                    'low_coverage_countries': low_coverage_countries,
                    'no_coverage_n': len(no_coverage_countries),
                    'low_coverage_n': len(low_coverage_countries),
                    'total_countries': len(results),
                })

        # Sort by number of countries with no coverage
        gaps.sort(key=lambda g: g['no_coverage_n'], reverse=True)

        return gaps

    def save_results(self, results: List[ROSResult], output_path: str):
        """
        Save ROS results to CSV.

        Args:
            results: List of ROSResult objects
            output_path: Output file path
        """
        df = pd.DataFrame([r.to_dict() for r in results])
        df.to_csv(output_path, index=False)
        logger.info(f"Saved ROS results to {output_path}")

    def save_site_details(self, results: List[ROSResult], output_path: str):
        """
        Save detailed site-level results to CSV.

        Args:
            results: List of ROSResult objects
            output_path: Output file path
        """
        rows = []
        for result in results:
            for site_result in result.site_results:
                row = site_result.to_dict()
                row['country'] = result.country
                row['region'] = result.region
                rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved site details to {output_path}")

    def save_heatmap(self, df: pd.DataFrame, output_path: str):
        """
        Save heatmap data to CSV.

        Args:
            df: Heatmap DataFrame
            output_path: Output file path
        """
        df.to_csv(output_path)
        logger.info(f"Saved heatmap data to {output_path}")


def main():
    """Example usage of ROS calculator."""
    import argparse

    parser = argparse.ArgumentParser(description='Calculate ROS scores')
    parser.add_argument('input_csv', help='Processed sequences CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    args = parser.parse_args()

    # Load processed sequences (simplified)
    df = pd.read_csv(args.input_csv)

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

    # Calculate ROS
    calculator = ROSCalculator(records)
    results = calculator.calculate_all_ros()
    regional_stats = calculator.calculate_regional_ros()
    heatmap_data = calculator.generate_heatmap_data(results)
    regional_heatmap = calculator.generate_regional_heatmap_data()

    # Save results
    os.makedirs(args.output, exist_ok=True)
    calculator.save_results(results, os.path.join(args.output, 'ros_scores.csv'))
    calculator.save_site_details(results, os.path.join(args.output, 'site_details.csv'))
    calculator.save_heatmap(heatmap_data, os.path.join(args.output, 'country_heatmap.csv'))
    calculator.save_heatmap(regional_heatmap, os.path.join(args.output, 'regional_heatmap.csv'))

    # Save critical gaps
    gaps = calculator.identify_critical_gaps(results)
    gaps_df = pd.DataFrame(gaps)
    gaps_df.to_csv(os.path.join(args.output, 'critical_gaps.csv'), index=False)

    print(f"\nROS Summary:")
    print(f"  Countries analyzed: {len(results)}")
    print(f"  Top 5 countries by ROS:")
    for result in results[:5]:
        print(f"    {result.country}: {result.ros_score:.1f}")
    print(f"\n  Critical gaps identified: {len(gaps)}")


if __name__ == '__main__':
    main()