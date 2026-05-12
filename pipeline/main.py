"""
HIV Capsid Surveillance Readiness Analysis Pipeline

Main orchestrator that runs the complete analysis pipeline:
1. Sequence quality control
2. CSO (Capsid Sequence Observability) calculation
3. ROS (RAM-site Observability Score) calculation
4. CIR-SRI (Composite Index) construction
5. Visualization generation
6. Sensitivity analysis
7. Report generation

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

import numpy as np
import pandas as pd

from . import config
from .sequence_qc import SequenceQC, QCReport
from .cso_calculator import CSOCalculator, CSOResult
from .ros_calculator import ROSCalculator, ROSResult
from .cir_sri_builder import CIRSRIBuilder, CIRSRICountryResult
from .figures import generate_all_figures
from .sensitivity_analysis import run_complete_sensitivity_analysis, generate_sensitivity_report

# Configure logging
def setup_logging(log_dir: str, log_level: int = logging.INFO) -> logging.Logger:
    """Set up logging configuration."""
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'pipeline_{timestamp}.log')

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    logger = logging.getLogger('pipeline')
    logger.info(f"Pipeline started at {datetime.now()}")
    logger.info(f"Log file: {log_file}")

    return logger


class Pipeline:
    """
    Main pipeline orchestrator for HIV capsid surveillance analysis.
    """

    def __init__(self, output_dir: str, data_dir: str, logger: Optional[logging.Logger] = None):
        """
        Initialize pipeline.

        Args:
            output_dir: Base output directory
            data_dir: Input data directory
            logger: Logger instance
        """
        self.output_dir = Path(output_dir)
        self.data_dir = Path(data_dir)
        self.logger = logger or logging.getLogger(__name__)

        # Create subdirectories
        self.subdirs = {
            'processed': self.output_dir / 'processed_sequences',
            'qc_reports': self.output_dir / 'qc_reports',
            'scores': self.output_dir / 'scores',
            'figures': self.output_dir / 'figures',
            'tables': self.output_dir / 'tables',
            'sensitivity': self.output_dir / 'sensitivity',
            'logs': self.output_dir / 'logs',
        }

        for subdir in self.subdirs.values():
            subdir.mkdir(parents=True, exist_ok=True)

        # Results storage
        self.sequences: List = []
        self.qc_report: Optional[QCReport] = None
        self.cso_results: List[CSOResult] = []
        self.ros_results: List[ROSResult] = []
        self.cir_sri_results: List[CIRSRICountryResult] = []
        self.regional_summary: Dict = {}

    def run_step(self, step_name: str, func, *args, **kwargs):
        """Run a pipeline step with error handling."""
        self.logger.info(f"Running step: {step_name}")
        try:
            result = func(*args, **kwargs)
            self.logger.info(f"Step {step_name} completed successfully")
            return result
        except Exception as e:
            self.logger.error(f"Step {step_name} failed: {e}")
            raise

    def step1_qc(self, input_files: List[str]) -> List:
        """
        Step 1: Quality control of input sequences.

        Args:
            input_files: List of input sequence files

        Returns:
            List of processed SequenceRecord objects
        """
        self.logger.info(f"Step 1: Quality Control - Processing {len(input_files)} files")

        qc = SequenceQC()
        all_records = []

        for file_path in input_files:
            if not os.path.exists(file_path):
                self.logger.warning(f"File not found: {file_path}")
                continue

            records, report = qc.process_file(file_path)
            all_records.extend(records)

            # Save per-file QC report
            filename = os.path.basename(file_path)
            qc.save_qc_report(
                report,
                self.subdirs['qc_reports'] / f"qc_{filename}.csv"
            )

        self.sequences = all_records
        self.qc_report = report

        # Save processed sequences
        qc.save_processed_sequences(
            all_records,
            self.subdirs['processed'] / 'all_processed.fasta'
        )

        self.logger.info(f"QC complete: {len(all_records)} sequences passed")

        return all_records

    def step2_cso(self) -> List[CSOResult]:
        """
        Step 2: Calculate CSO (Capsid Sequence Observability) scores.

        Returns:
            List of CSOResult objects
        """
        self.logger.info("Step 2: Calculating CSO scores")

        calculator = CSOCalculator(self.sequences)
        results = calculator.calculate_all_cso()
        regional_stats = calculator.calculate_regional_cso()

        # Save results
        calculator.save_results(results, self.subdirs['scores'] / 'cso_scores.csv')
        calculator.save_regional_summary(regional_stats, self.subdirs['scores'] / 'cso_regional_summary.csv')

        self.cso_results = results

        self.logger.info(f"CSO calculation complete: {len(results)} countries")

        return results

    def step3_ros(self) -> List[ROSResult]:
        """
        Step 3: Calculate ROS (RAM-site Observability Score) scores.

        Returns:
            List of ROSResult objects
        """
        self.logger.info("Step 3: Calculating ROS scores")

        calculator = ROSCalculator(self.sequences)
        results = calculator.calculate_all_ros()
        regional_stats = calculator.calculate_regional_ros()

        # Save results
        calculator.save_results(results, self.subdirs['scores'] / 'ros_scores.csv')
        calculator.save_site_details(results, self.subdirs['scores'] / 'ros_site_details.csv')

        # Generate heatmaps
        heatmap = calculator.generate_heatmap_data(results)
        regional_heatmap = calculator.generate_regional_heatmap_data()
        calculator.save_heatmap(heatmap, self.subdirs['scores'] / 'ros_country_heatmap.csv')
        calculator.save_heatmap(regional_heatmap, self.subdirs['scores'] / 'ros_regional_heatmap.csv')

        # Identify critical gaps
        gaps = calculator.identify_critical_gaps(results)
        gaps_df = pd.DataFrame(gaps)
        gaps_df.to_csv(self.subdirs['scores'] / 'ros_critical_gaps.csv', index=False)

        self.ros_results = results

        self.logger.info(f"ROS calculation complete: {len(results)} countries")

        return results

    def step4_cir_sri(self) -> List[CIRSRICountryResult]:
        """
        Step 4: Build CIR-SRI (Composite Index).

        Returns:
            List of CIRSRICountryResult objects
        """
        self.logger.info("Step 4: Building CIR-SRI")

        builder = CIRSRIBuilder(self.cso_results, self.ros_results)
        results = builder.calculate_all_cir_sri()
        regional_summaries = builder.calculate_regional_summary(results)

        # Save results
        builder.save_results(results, self.subdirs['scores'] / 'cir_sri_scores.csv')
        builder.save_dimension_details(results, self.subdirs['scores'] / 'cir_sri_dimensions.csv')
        builder.save_priority_list(
            builder.generate_priority_list(results),
            self.subdirs['scores'] / 'priority_list.csv'
        )
        builder.save_regional_summary(regional_summaries, self.subdirs['scores'] / 'cir_sri_regional_summary.csv')

        # Generate report
        report = builder.generate_summary_report(results, regional_summaries)
        with open(self.subdirs['scores'] / 'cir_sri_report.txt', 'w') as f:
            f.write(report)

        self.cir_sri_results = results
        self.regional_summary = regional_summaries

        self.logger.info(f"CIR-SRI calculation complete: {len(results)} countries")

        return results

    def step5_figures(self):
        """Step 5: Generate all figures."""
        self.logger.info("Step 5: Generating figures")

        # Load data for figures
        cir_sri_df = pd.DataFrame([r.to_dict() for r in self.cir_sri_results])
        regional_data = self.regional_summary

        # Load heatmap data
        heatmap_path = self.subdirs['scores'] / 'ros_country_heatmap.csv'
        site_heatmap = pd.read_csv(heatmap_path, index_col=0) if heatmap_path.exists() else pd.DataFrame()

        # Generate temporal data
        temporal_data = self._generate_temporal_data()

        # Generate all figures
        generate_all_figures(
            cir_sri_df=cir_sri_df,
            regional_summary=regional_data,
            site_heatmap=site_heatmap,
            temporal_data=temporal_data,
            output_dir=str(self.subdirs['figures'])
        )

        self.logger.info("Figure generation complete")

    def _generate_temporal_data(self) -> pd.DataFrame:
        """Generate temporal coverage data for heatmap."""
        yearly_country_counts = {}

        for record in self.sequences:
            if record.year and record.country:
                key = (record.year, record.country)
                yearly_country_counts[key] = yearly_country_counts.get(key, 0) + 1

        # Convert to DataFrame
        data_dict = defaultdict(dict)
        for (year, country), count in yearly_country_counts.items():
            data_dict[country][year] = count

        df = pd.DataFrame(data_dict).fillna(0)
        df.index.name = 'country'

        return df

    def step6_sensitivity(self):
        """Step 6: Run sensitivity analysis."""
        self.logger.info("Step 6: Running sensitivity analysis")

        # Prepare data
        country_scores = {r.country: r.cir_sri_score for r in self.cir_sri_results}
        regional_mapping = {r.country: r.region for r in self.cir_sri_results}

        # Run analysis
        results = run_complete_sensitivity_analysis(
            country_scores=country_scores,
            regional_mapping=regional_mapping,
            sequence_records=self.sequences,
            output_dir=str(self.subdirs['sensitivity'])
        )

        # Generate report
        report = generate_sensitivity_report(results)
        with open(self.subdirs['sensitivity'] / 'sensitivity_report.txt', 'w') as f:
            f.write(report)

        self.logger.info("Sensitivity analysis complete")

    def step7_generate_summary(self) -> str:
        """Step 7: Generate final summary report."""
        self.logger.info("Step 7: Generating summary report")

        lines = [
            "=" * 80,
            "HIV CAPSID SURVEILLANCE READINESS ANALYSIS",
            "=" * 80,
            "",
            f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Output Directory: {self.output_dir}",
            "",
            "Pipeline Summary:",
            "-" * 40,
            f"  Step 1 - QC: {len(self.sequences)} sequences processed",
            f"  Step 2 - CSO: {len(self.cso_results)} countries scored",
            f"  Step 3 - ROS: {len(self.ros_results)} countries scored",
            f"  Step 4 - CIR-SRI: {len(self.cir_sri_results)} countries indexed",
            "",
        ]

        # Overall statistics
        if self.cir_sri_results:
            scores = [r.cir_sri_score for r in self.cir_sri_results]
            lines.extend([
                "Global Statistics:",
                f"  Mean CIR-SRI: {np.mean(scores):.2f}",
                f"  Median CIR-SRI: {np.median(scores):.2f}",
                f"  Std Dev: {np.std(scores):.2f}",
                f"  Min: {min(scores):.1f}, Max: {max(scores):.1f}",
                "",
            ])

            # Priority tier distribution
            tier_counts = {}
            for result in self.cir_sri_results:
                tier = result.priority_tier
                tier_counts[tier] = tier_counts.get(tier, 0) + 1

            lines.extend([
                "Priority Tier Distribution:",
                f"  Urgent (0-3): {tier_counts.get('urgent', 0)} countries",
                f"  High (3-5): {tier_counts.get('high', 0)} countries",
                f"  Medium (5-7): {tier_counts.get('medium', 0)} countries",
                f"  Low (7-10): {tier_counts.get('low', 0)} countries",
                "",
            ])

        # Regional summary
        if self.regional_summary:
            lines.extend([
                "Regional Summary:",
                "-" * 40,
            ])
            for region, stats in self.regional_summary.items():
                lines.append(
                    f"  {region}: Mean={stats.get('mean_cir_sri', 0):.2f}, "
                    f"Countries={stats.get('n_countries', 0)}"
                )
            lines.append("")

        # Top countries
        if self.cir_sri_results:
            sorted_results = sorted(self.cir_sri_results, key=lambda r: r.cir_sri_score, reverse=True)
            lines.extend([
                "Top 10 Countries:",
                "-" * 40,
            ])
            for i, result in enumerate(sorted_results[:10], 1):
                lines.append(
                    f"  {i}. {result.country} ({result.region}): {result.cir_sri_score:.1f} "
                    f"[{result.priority_tier}]"
                )
            lines.append("")

        # Output files
        lines.extend([
            "Output Files:",
            "-" * 40,
        ])
        for name, subdir in self.subdirs.items():
            files = list(subdir.glob('*'))
            if files:
                lines.append(f"  {name}/: {len(files)} files")

        lines.extend([
            "",
            "=" * 80,
            "Pipeline completed successfully",
            "=" * 80,
        ])

        report = "\n".join(lines)

        # Save report
        with open(self.output_dir / 'pipeline_summary.txt', 'w') as f:
            f.write(report)

        self.logger.info("Summary report generated")

        return report

    def run(self, input_files: List[str], skip_steps: List[str] = None) -> Dict:
        """
        Run the complete pipeline.

        Args:
            input_files: List of input sequence files
            skip_steps: Steps to skip (for re-running partial analyses)

        Returns:
            Dictionary with all results
        """
        skip_steps = skip_steps or []

        self.logger.info("Starting complete pipeline run")

        start_time = datetime.now()

        try:
            # Step 1: QC
            if 'qc' not in skip_steps:
                self.step1_qc(input_files)

            # Step 2: CSO
            if 'cso' not in skip_steps:
                self.step2_cso()

            # Step 3: ROS
            if 'ros' not in skip_steps:
                self.step3_ros()

            # Step 4: CIR-SRI
            if 'cir_sri' not in skip_steps:
                self.step4_cir_sri()

            # Step 5: Figures
            if 'figures' not in skip_steps:
                self.step5_figures()

            # Step 6: Sensitivity
            if 'sensitivity' not in skip_steps:
                self.step6_sensitivity()

            # Step 7: Summary
            summary = self.step7_generate_summary()

            end_time = datetime.now()
            elapsed = end_time - start_time

            self.logger.info(f"Pipeline completed in {elapsed}")

            return {
                'success': True,
                'elapsed_time': str(elapsed),
                'sequences_processed': len(self.sequences),
                'countries_analyzed': len(self.cir_sri_results),
                'output_dir': str(self.output_dir),
            }

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            return {
                'success': False,
                'error': str(e),
                'elapsed_time': str(datetime.now() - start_time),
            }

    def load_previous_results(self, results_dir: str):
        """
        Load results from a previous run for re-analysis.

        Args:
            results_dir: Directory with previous results
        """
        results_dir = Path(results_dir)

        # Load CIR-SRI results
        cir_sri_path = results_dir / 'scores' / 'cir_sri_scores.csv'
        if cir_sri_path.exists():
            cir_sri_df = pd.read_csv(cir_sri_path)
            # Convert to result objects (simplified)
            for _, row in cir_sri_df.iterrows():
                result = CIRSRICountryResult(
                    country=row['country'],
                    region=row['region'],
                    cir_sri_score=row['cir_sri_score'],
                    priority_tier=row['priority_tier'],
                )
                self.cir_sri_results.append(result)

        # Load regional summary
        regional_path = results_dir / 'scores' / 'cir_sri_regional_summary.csv'
        if regional_path.exists():
            regional_df = pd.read_csv(regional_path)
            self.regional_summary = regional_df.set_index('region').to_dict('index')


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description='HIV Capsid Surveillance Readiness Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python -m scripts.main --input data/sequences.fasta --output results/

  # Re-run from ROS step (assuming CSO results exist)
  python -m scripts.main --input data/sequences.fasta --output results/ --skip qc cso

  # Run with multiple input files
  python -m scripts.main --input data/*.fasta --output results/
        """
    )

    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input sequence files (FASTA or GenBank format)'
    )
    parser.add_argument(
        '-o', '--output',
        default='results',
        help='Output directory (default: results)'
    )
    parser.add_argument(
        '--skip',
        nargs='+',
        choices=['qc', 'cso', 'ros', 'cir_sri', 'figures', 'sensitivity'],
        help='Steps to skip'
    )
    parser.add_argument(
        '--load-results',
        help='Load previous results directory for re-analysis'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(os.path.join(args.output, 'logs'), log_level)

    # Initialize pipeline
    pipeline = Pipeline(
        output_dir=args.output,
        data_dir=os.path.dirname(args.input[0]) if args.input else '.',
        logger=logger
    )

    # Load previous results if specified
    if args.load_results:
        logger.info(f"Loading previous results from {args.load_results}")
        pipeline.load_previous_results(args.load_results)

    # Run pipeline
    results = pipeline.run(
        input_files=args.input,
        skip_steps=args.skip
    )

    # Print summary
    if results['success']:
        print("\n" + "=" * 60)
        print("Pipeline completed successfully!")
        print("=" * 60)
        print(f"Sequences processed: {results['sequences_processed']}")
        print(f"Countries analyzed: {results['countries_analyzed']}")
        print(f"Output directory: {results['output_dir']}")
        print(f"Total time: {results['elapsed_time']}")
        print("=" * 60)

        # Print summary report
        summary_path = os.path.join(args.output, 'pipeline_summary.txt')
        if os.path.exists(summary_path):
            print("\n--- Pipeline Summary ---\n")
            with open(summary_path, 'r') as f:
                print(f.read())
    else:
        print("\n" + "=" * 60)
        print("Pipeline FAILED!")
        print("=" * 60)
        print(f"Error: {results.get('error', 'Unknown error')}")
        print(f"Time elapsed: {results.get('elapsed_time', 'N/A')}")
        print("=" * 60)
        sys.exit(1)


if __name__ == '__main__':
    main()