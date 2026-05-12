"""
Sequence Quality Control Pipeline for HIV Capsid Surveillance Analysis

This module handles:
- FASTA/GenBank file parsing
- Metadata extraction from sequence headers
- HXB2 coordinate mapping
- Hypermutation detection
- Stop codon and frameshift identification
- Sequence deduplication

Author: HIV Capsid Surveillance Analysis Pipeline
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict
import hashlib

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from . import config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class SequenceRecord:
    """Container for processed sequence data."""
    accession: str
    sequence: str
    sequence_aa: str
    country: Optional[str] = None
    year: Optional[int] = None
    month: Optional[int] = None
    subtype: Optional[str] = None
    collection_date: Optional[str] = None
    lineage: Optional[str] = None
    length: int = 0
    ca_start_pos: Optional[int] = None
    ca_end_pos: Optional[int] = None
    has_stop_codons: bool = False
    is_hypermutated: bool = False
    qc_passed: bool = True
    qc_failures: List[str] = field(default_factory=list)


@dataclass
class QCReport:
    """Quality control summary report."""
    total_sequences: int = 0
    passed_sequences: int = 0
    failed_sequences: int = 0
    failure_reasons: Dict[str, int] = field(default_factory=dict)
    sequences_by_country: Dict[str, int] = field(default_factory=dict)
    sequences_by_year: Dict[str, int] = field(default_factory=dict)
    sequences_by_subtype: Dict[str, int] = field(default_factory=dict)
    duplicate_count: int = 0
    hypermutated_count: int = 0
    stop_codon_count: int = 0


class MetadataParser:
    """Parse metadata from various sequence header formats."""

    # Common patterns for HIV sequence headers
    HEADER_PATTERNS = {
        # LANL-style headers
        'lanl': r'(?P<accession>[A-Z]{2,3}\d+)\s+(?P<country>[A-Z]{2})\s+(?P<year>\d{4})\s+(?P< subtype>[A-Z0-9]+)',
        # GenBank-style
        'genbank': r'(?P<accession>\w+)\s+(?P<description>.+)',
        # Los Alamos HIV Sequence Database
        'los_alamos': r'(?P<accession>[A-Z]{2,3}\d+)\|(?P<name>[^|]+)\|(?P<country>[A-Z]{2})\|(?P<year>\d{4})',
        # Simplified: accession country year subtype
        'simple': r'(?P<accession>\S+)\s+(?P<country>[A-Z]{2})\s+(?P<year>\d{4})',
    }

    # Countries that may appear in headers
    VALID_COUNTRIES = set(config.COUNTRY_TO_REGION.keys())

    # Subtype pattern
    SUBTYPE_PATTERN = re.compile(
        r'(?:subtype|clade|CRF\d+_[a-z]+|CRF\d+)\s*[:\s]*([A-Z0-9_]+)',
        re.IGNORECASE
    )

    def __init__(self):
        self.compiled_patterns = {
            name: re.compile(pattern)
            for name, pattern in self.HEADER_PATTERNS.items()
        }

    def parse_header(self, header: str) -> Dict:
        """
        Extract metadata from sequence header.

        Args:
            header: Sequence header string

        Returns:
            Dictionary of extracted metadata
        """
        metadata = {
            'accession': header.split()[0] if header else None,
            'country': None,
            'year': None,
            'subtype': None,
            'raw_header': header,
        }

        # Try each pattern
        for pattern_name, pattern in self.compiled_patterns.items():
            match = pattern.search(header)
            if match:
                matched = match.groupdict()
                metadata.update({k: v for k, v in matched.items() if v is not None})

                # Validate country code
                if metadata.get('country') and metadata['country'] not in self.VALID_COUNTRIES:
                    # Try to find a valid country code nearby
                    possible_codes = re.findall(r'[A-Z]{2}', header)
                    for code in possible_codes:
                        if code in self.VALID_COUNTRIES:
                            metadata['country'] = code
                            break

                # Parse year
                if metadata.get('year'):
                    try:
                        metadata['year'] = int(metadata['year'])
                    except (ValueError, TypeError):
                        metadata['year'] = None

                break

        # Extract subtype from description if not found
        if not metadata.get('subtype'):
            subtype_match = self.SUBTYPE_PATTERN.search(header)
            if subtype_match:
                metadata['subtype'] = subtype_match.group(1)

        return metadata

    def parse_date(self, date_str: str) -> Tuple[Optional[int], Optional[int]]:
        """
        Parse date string into year and month.

        Args:
            date_str: Date string in various formats

        Returns:
            Tuple of (year, month)
        """
        if not date_str:
            return None, None

        year, month = None, None

        # Try YYYY-MM-DD format
        match = re.match(r'(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?', date_str)
        if match:
            year = int(match.group(1))
            month = int(match.group(2)) if match.group(2) else None

        return year, month


class HypermutationDetector:
    """Detect APOBEC-mediated hypermutation in HIV sequences."""

    def __init__(self, threshold: float = config.QC_THRESHOLDS['hypermutation_threshold']):
        self.threshold = threshold
        # G→A mutation patterns characteristic of APOBEC
        self.g_to_a_pattern = re.compile(r'G[AG]A|G[AG]G')

    def detect(self, sequence: str) -> Tuple[bool, float]:
        """
        Detect hypermutation in a sequence.

        Args:
            sequence: DNA sequence string

        Returns:
            Tuple of (is_hypermutated, g_to_a_ratio)
        """
        if not sequence:
            return False, 0.0

        # Count G nucleotides and G→A mutations
        g_count = sequence.count('G')
        total_nucs = len(sequence)

        if total_nucs == 0:
            return False, 0.0

        # Count G→A patterns in codons
        g_to_a_count = len(self.g_to_a_pattern.findall(sequence))

        # Calculate ratio
        g_to_a_ratio = g_to_a_count / total_nucs if total_nucs > 0 else 0.0

        is_hypermutated = g_to_a_ratio > self.threshold

        return is_hypermutated, g_to_a_ratio

    def find_mutation_clusters(self, sequence: str, window_size: int = 100) -> List[int]:
        """
        Find positions with clusters of potential hypermutation.

        Args:
            sequence: DNA sequence
            window_size: Size of sliding window

        Returns:
            List of start positions of mutation clusters
        """
        clusters = []
        sequence_upper = sequence.upper()

        for i in range(0, len(sequence_upper) - window_size, 10):
            window = sequence_upper[i:i + window_size]
            g_to_a_in_window = len(self.g_to_a_pattern.findall(window))
            if g_to_a_in_window / window_size > self.threshold:
                clusters.append(i)

        return clusters


class SequenceDeduplicator:
    """Remove duplicate sequences while preserving best representative."""

    def __init__(self):
        self.seen_sequences = {}

    def add_sequence(self, record: SequenceRecord) -> bool:
        """
        Check if sequence is duplicate and track if not.

        Args:
            record: SequenceRecord to check/add

        Returns:
            True if sequence is new, False if duplicate
        """
        seq_hash = hashlib.md5(record.sequence.encode()).hexdigest()

        if seq_hash in self.seen_sequences:
            # Check if this record has better metadata
            existing = self.seen_sequences[seq_hash]
            if self._is_better_metadata(record, existing):
                self.seen_sequences[seq_hash] = record
            return False
        else:
            self.seen_sequences[seq_hash] = record
            return True

    def _is_better_metadata(self, new: SequenceRecord, existing: SequenceRecord) -> bool:
        """Compare two records, return True if new has better metadata."""
        # Prefer records with more metadata fields filled
        new_fields = sum([
            new.country is not None,
            new.year is not None,
            new.subtype is not None
        ])
        existing_fields = sum([
            existing.country is not None,
            existing.year is not None,
            existing.subtype is not None
        ])

        if new_fields != existing_fields:
            return new_fields > existing_fields

        # Prefer longer sequences
        return len(new.sequence) > len(existing.sequence)

    def get_deduplicated(self) -> List[SequenceRecord]:
        """Return list of unique sequences."""
        return list(self.seen_sequences.values())


class SequenceQC:
    """Main quality control pipeline for HIV sequences."""

    STOP_CODONS = {'TAA', 'TAG', 'TGA'}

    def __init__(self, config_dict: Optional[Dict] = None):
        self.config = config_dict or config.QC_THRESHOLDS
        self.metadata_parser = MetadataParser()
        self.hypermutation_detector = HypermutationDetector()
        self.deduplicator = SequenceDeduplicator()

    def parse_fasta(self, file_path: str) -> List[SequenceRecord]:
        """
        Parse FASTA file and extract sequences with metadata.

        Args:
            file_path: Path to FASTA file

        Returns:
            List of SequenceRecord objects
        """
        logger.info(f"Parsing FASTA file: {file_path}")

        records = []
        for seq_record in SeqIO.parse(file_path, 'fasta'):
            metadata = self.metadata_parser.parse_header(seq_record.description)

            seq_record_obj = SequenceRecord(
                accession=metadata.get('accession', seq_record.id),
                sequence=str(seq_record.seq),
                sequence_aa=str(seq_record.seq.translate()),
                country=metadata.get('country'),
                year=metadata.get('year'),
                subtype=metadata.get('subtype'),
                collection_date=metadata.get('collection_date'),
                length=len(seq_record.seq),
            )

            records.append(seq_record_obj)

        logger.info(f"Parsed {len(records)} sequences from FASTA")
        return records

    def parse_genbank(self, file_path: str) -> List[SequenceRecord]:
        """
        Parse GenBank file and extract sequences with metadata.

        Args:
            file_path: Path to GenBank file

        Returns:
            List of SequenceRecord objects
        """
        logger.info(f"Parsing GenBank file: {file_path}")

        records = []
        for seq_record in SeqIO.parse(file_path, 'genbank'):
            # Extract metadata from qualifiers
            qualifiers = seq_record.annotations.get('taxonomy', [])

            # Get country from features
            country = None
            year = None
            subtype = None

            for feature in seq_record.features:
                if feature.type == 'source':
                    qual = feature.qualifiers
                    country = qual.get('country', [None])[0]
                    if country:
                        # Extract 2-letter code if present
                        country_match = re.search(r'([A-Z]{2}):', country)
                        if country_match:
                            country = country_match.group(1)

                elif feature.type == 'region' and feature.qualifiers.get('note'):
                    note = str(feature.qualifiers['note'])
                    subtype_match = re.search(r'subtype[:\s]*([A-Z0-9]+)', note, re.I)
                    if subtype_match:
                        subtype = subtype_match.group(1)

            # Check for date annotation
            dates = seq_record.annotations.get('dates', '')
            year_match = re.search(r'(\d{4})', dates)
            if year_match:
                year = int(year_match.group(1))

            seq_record_obj = SequenceRecord(
                accession=seq_record.id,
                sequence=str(seq_record.seq),
                sequence_aa=str(seq_record.seq.translate()),
                country=country,
                year=year,
                subtype=subtype,
                collection_date=seq_record.annotations.get('date'),
                length=len(seq_record.seq),
            )

            records.append(seq_record_obj)

        logger.info(f"Parsed {len(records)} sequences from GenBank")
        return records

    def map_to_hxb2(self, record: SequenceRecord) -> Tuple[Optional[int], Optional[int]]:
        """
        Determine CA region coordinates relative to HXB2.

        Args:
            record: SequenceRecord with sequence

        Returns:
            Tuple of (ca_start_pos, ca_end_pos) in 0-indexed coordinates
        """
        # This requires alignment to HXB2 reference
        # For now, assume input sequences are already aligned or
        # we use a simple heuristic based on sequence length

        seq_length = len(record.sequence)

        # HXB2 reference is 9718 bp
        # CA region is 1338-1872 (535 aa = 1605 bp)
        # If sequence is near full-length, map positions

        if seq_length >= 8000:
            # Near full-length: assume good mapping
            # Position 1 corresponds to first base
            ca_start = config.CA_START - 1  # Convert to 0-indexed
            ca_end = config.CA_END - 1

            if ca_start >= 0 and ca_end < seq_length:
                record.ca_start_pos = ca_start
                record.ca_end_pos = ca_end
                return ca_start, ca_end

        return None, None

    def extract_ca_region(self, record: SequenceRecord) -> str:
        """
        Extract CA (capsid) region from sequence.

        Args:
            record: SequenceRecord with mapped coordinates

        Returns:
            CA region sequence or empty string if not mappable
        """
        if record.ca_start_pos is not None and record.ca_end_pos is not None:
            return record.sequence[record.ca_start_pos:record.ca_end_pos + 1]
        return ''

    def check_stop_codons(self, record: SequenceRecord) -> bool:
        """
        Check for stop codons in CA region.

        Args:
            record: SequenceRecord

        Returns:
            True if stop codons found in CA region
        """
        ca_region = self.extract_ca_region(record)
        if not ca_region:
            return False

        # Check each codon
        for i in range(0, len(ca_region) - 2, 3):
            codon = ca_region[i:i + 3].upper()
            if codon in self.STOP_CODONS:
                record.has_stop_codons = True
                return True

        return False

    def check_hypermutation(self, record: SequenceRecord) -> bool:
        """
        Check for hypermutation in sequence.

        Args:
            record: SequenceRecord

        Returns:
            True if hypermutation detected
        """
        is_hypermutated, _ = self.hypermutation_detector.detect(record.sequence)
        record.is_hypermutated = is_hypermutated
        return is_hypermutated

    def run_qc(self, record: SequenceRecord) -> bool:
        """
        Run full QC checks on a sequence record.

        Args:
            record: SequenceRecord to QC

        Returns:
            True if sequence passes all QC checks
        """
        record.qc_failures = []

        # Check sequence length
        if record.length < self.config['min_sequence_length']:
            record.qc_failures.append(f"short_length:{record.length}")
            record.qc_passed = False

        # Check collection year
        if record.year and record.year < self.config['min_year']:
            record.qc_failures.append(f"old_sequence:{record.year}")
            record.qc_passed = False

        # Check ambiguous characters
        ambiguous = sum(1 for c in record.sequence if c.upper() not in 'ACGTN')
        ambiguous_ratio = ambiguous / len(record.sequence) if record.sequence else 1
        if ambiguous_ratio > self.config['max_ambiguous_chars']:
            record.qc_failures.append(f"high_ambiguity:{ambiguous_ratio:.3f}")
            record.qc_passed = False

        # Map to HXB2 coordinates
        self.map_to_hxb2(record)

        # Check for stop codons
        if self.config['stop_codon_tolerance'] == 0:
            if self.check_stop_codons(record):
                record.qc_failures.append("stop_codon_in_CA")
                record.qc_passed = False

        # Check for hypermutation
        if self.config['hypermutation_detection']:
            if self.check_hypermutation(record):
                record.qc_failures.append("hypermutated")
                record.qc_passed = False

        return record.qc_passed

    def process_file(self, file_path: str) -> Tuple[List[SequenceRecord], QCReport]:
        """
        Process a sequence file through the full QC pipeline.

        Args:
            file_path: Path to sequence file (FASTA or GenBank)

        Returns:
            Tuple of (list of SequenceRecords, QCReport)
        """
        file_ext = Path(file_path).suffix.lower()

        if file_ext in ['.fasta', '.fa', '.faa']:
            records = self.parse_fasta(file_path)
        elif file_ext in ['.gb', '.gbk', '.genbank']:
            records = self.parse_genbank(file_path)
        else:
            logger.error(f"Unsupported file format: {file_ext}")
            return [], QCReport()

        # Run QC on each record
        passed_records = []
        report = QCReport(total_sequences=len(records))

        for record in records:
            if self.run_qc(record):
                passed_records.append(record)

                # Track for deduplication
                self.deduplicator.add_sequence(record)
            else:
                # Count failures
                for failure in record.qc_failures:
                    reason = failure.split(':')[0]
                    report.failure_reasons[reason] = report.failure_reasons.get(reason, 0) + 1

        # Get deduplicated sequences
        unique_records = self.deduplicator.get_deduplicated()
        report.duplicate_count = len(passed_records) - len(unique_records)

        # Update report
        report.passed_sequences = len(unique_records)
        report.failed_sequences = report.total_sequences - report.passed_sequences
        report.hypermutated_count = sum(1 for r in records if r.is_hypermutated)
        report.stop_codon_count = sum(1 for r in records if r.has_stop_codons)

        # Count by categories
        for record in unique_records:
            if record.country:
                report.sequences_by_country[record.country] = \
                    report.sequences_by_country.get(record.country, 0) + 1
            if record.year:
                report.sequences_by_year[str(record.year)] = \
                    report.sequences_by_year.get(str(record.year), 0) + 1
            if record.subtype:
                report.sequences_by_subtype[record.subtype] = \
                    report.sequences_by_subtype.get(record.subtype, 0) + 1

        logger.info(f"QC complete: {report.passed_sequences}/{report.total_sequences} passed")
        logger.info(f"Removed {report.duplicate_count} duplicates")

        return unique_records, report

    def save_processed_sequences(self, records: List[SequenceRecord],
                                  output_path: str):
        """
        Save processed sequences to FASTA file.

        Args:
            records: List of SequenceRecords
            output_path: Output file path
        """
        seq_records = []
        for rec in records:
            description = f"{rec.accession}|{rec.country or 'UNK'}|{rec.year or 'UNK'}|{rec.subtype or 'UNK'}"
            seq_records.append(
                SeqRecord(
                    Seq(rec.sequence),
                    id=rec.accession,
                    description=description
                )
            )

        SeqIO.write(seq_records, output_path, 'fasta')
        logger.info(f"Saved {len(records)} sequences to {output_path}")

    def save_qc_report(self, report: QCReport, output_path: str):
        """
        Save QC report to file.

        Args:
            report: QCReport object
            output_path: Output file path
        """
        report_dict = {
            'total_sequences': report.total_sequences,
            'passed_sequences': report.passed_sequences,
            'failed_sequences': report.failed_sequences,
            'duplicate_count': report.duplicate_count,
            'hypermutated_count': report.hypermutated_count,
            'stop_codon_count': report.stop_codon_count,
        }

        # Add failure reasons
        for reason, count in report.failure_reasons.items():
            report_dict[f'failure_{reason}'] = count

        # Add country counts
        for country, count in report.sequences_by_country.items():
            report_dict[f'country_{country}'] = count

        df = pd.DataFrame([report_dict])
        df.to_csv(output_path, index=False)
        logger.info(f"Saved QC report to {output_path}")


def main():
    """Example usage of the QC pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description='HIV Capsid Sequence QC Pipeline')
    parser.add_argument('input_file', help='Input FASTA or GenBank file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--prefix', default='processed', help='Output file prefix')

    args = parser.parse_args()

    qc = SequenceQC()
    records, report = qc.process_file(args.input_file)

    os.makedirs(args.output, exist_ok=True)

    # Save outputs
    base_prefix = os.path.join(args.output, args.prefix)
    qc.save_processed_sequences(records, f"{base_prefix}.fasta")
    qc.save_qc_report(report, f"{base_prefix}_qc_report.csv")

    print(f"\nQC Summary:")
    print(f"  Total: {report.total_sequences}")
    print(f"  Passed: {report.passed_sequences}")
    print(f"  Failed: {report.failed_sequences}")
    print(f"  Duplicates removed: {report.duplicate_count}")


if __name__ == '__main__':
    main()