# Configuration for HIV capsid surveillance readiness analysis
# Part of: submission_hiv_medicine pipeline

# =============================================================================
# HXB2 Coordinates for CA Region
# =============================================================================
# Capsid (CA) region in HXB2 genome: positions 1338-1872 (1-indexed)
CA_START = 1338
CA_END = 1872
CA_LENGTH = CA_END - CA_START + 1  # 535 amino acids

# =============================================================================
# RAM Target Sites (Resistance-Associated Mutations in Capsid)
# =============================================================================
# HXB2 coordinates for RAM sites (1-indexed)
# Sites identified in literature as important for capsid inhibitor resistance

RAM_SITES = {
    'L56': 1338 + 55,   # HXB2 ~1393 (CA position 56)
    'N57': 1338 + 56,   # HXB2 ~1394 (CA position 57)
    'M66': 1338 + 65,   # HXB2 ~1403 (CA position 66)
    'Q67': 1338 + 66,   # HXB2 ~1404 (CA position 67)
    'K70': 1338 + 69,   # HXB2 ~1407 (CA position 70)
    'N74': 1338 + 73,   # HXB2 ~1411 (CA position 74)
    'A105': 1338 + 104, # HXB2 ~1442 (CA position 105)
    'T107': 1338 + 106, # HXB2 ~1444 (CA position 107)
}

# RAM amino acid positions in CA protein (1-indexed)
RAM_AA_POSITIONS = ['56', '57', '66', '67', '70', '74', '105', '107']

# RAM amino acid alternatives (for resistance detection)
RAM_ALTERNATIVES = {
    'L56': ['L56V', 'L56I', 'L56M'],
    'N57': ['N57K', 'N57S', 'N57H', 'N57D'],
    'M66': ['M66I', 'M66V', 'M66K'],
    'Q67': ['Q67H', 'Q67K', 'Q67R', 'Q67N'],
    'K70': ['K70N', 'K70E', 'K70R', 'K70T'],
    'N74': ['N74D', 'N74Y', 'N74H'],
    'A105': ['A105T', 'A105V', 'A105S'],
    'T107': ['T107I', 'T107N', 'T107S'],
}

# =============================================================================
# WHO Regions
# =============================================================================
WHO_REGIONS = {
    'AFRO': ['ZA', 'KE', 'BW', 'MW', 'ZW', 'UG', 'ZM', 'TZ', 'NG', 'ET', 'GH', 'MW', 'MZ', 'NA', 'RW', 'SN', 'SZ', 'TZ', 'ZA', 'ZM', 'ZW'],
    'AMRO': ['US', 'BR', 'PE', 'CA', 'AR', 'MX', 'CO', 'CL', 'VE', 'EC', 'BO', 'PY', 'UY'],
    'SEARO': ['TH', 'IN', 'ID', 'VN', 'MM', 'BD', 'BT', 'LK', 'MV', 'NP'],
    'WPRO': ['CN', 'AU', 'JP', 'KR', 'PH', 'TH', 'VN', 'MY', 'SG', 'NZ', 'TW', 'HK'],
    'EURO': ['GB', 'DE', 'FR', 'BE', 'ES', 'IT', 'NL', 'RU', 'PL', 'SE', 'CH', 'AT', 'DK', 'FI', 'NO', 'IE', 'PT'],
    'EMRO': ['EG', 'IR', 'PK', 'SD', 'SA', 'YE', 'JO', 'LB', 'SY', 'IQ', 'AF'],
}

# Reverse lookup: country to region
COUNTRY_TO_REGION = {}
for region, countries in WHO_REGIONS.items():
    for country in countries:
        COUNTRY_TO_REGION[country] = region

# =============================================================================
# CIR-SRI Scoring Thresholds
# =============================================================================
CIR_SRI_THRESHOLDS = {
    'sequence_availability': {
        'low': 100,      # < 100 sequences = low
        'medium': 1000,  # 100-1000 = medium
        # > 1000 = high
    },
    'metadata_completeness': {
        'low': 0.3,      # < 30% complete = low
        'medium': 0.7,   # 30-70% = medium
        # > 70% = high
    },
    'ram_site_observability': {
        'low': 0.5,      # < 50% sites observed = low
        'medium': 0.8,   # 50-80% = medium
        # > 80% = high
    },
    'temporal_coverage': {
        'low': 2,        # < 2 years of data = low
        'medium': 5,     # 2-5 years = medium
        # > 5 years = high
    },
    'subtype_diversity': {
        'low': 2,        # < 2 subtypes = low
        'medium': 5,     # 2-5 subtypes = medium
        # > 5 subtypes = high
    },
}

# =============================================================================
# Sensitivity Analysis Parameters
# =============================================================================
SENSITIVITY_PARAMS = {
    'min_year': [2000, 2010, 2015, 2018],
    'min_sequence_length': [7000, 8000, 9000],
    'exclude_subtype_b': [True, False],
    'min_ram_coverage': [0.5, 0.7, 0.9],
    'temporal_window_years': [3, 5, 10],
}

# =============================================================================
# Data Quality Thresholds
# =============================================================================
QC_THRESHOLDS = {
    'min_sequence_length': 7000,       # Minimum length for inclusion
    'min_year': 2000,                  # Earliest year for inclusion
    'max_ambiguous_chars': 0.05,       # Max 5% ambiguous characters
    'stop_codon_tolerance': 0,         # No stop codons allowed in CA region
    'hypermutation_detection': True,   # Enable APOBEC detection
    'hypermutation_threshold': 0.15,   # 15% G→A ratio threshold
}

# =============================================================================
# Visualization Settings
# =============================================================================
PLOT_SETTINGS = {
    'dpi': 300,
    'fig_size': (10, 8),
    'font_size': 12,
    'title_font_size': 14,
    'color_palette': {
        'AFRO': '#E64A19',   # Deep Orange
        'AMRO': '#1976D2',   # Blue
        'SEARO': '#388E3C',  # Green
        'WPRO': '#7B1FA2',   # Purple
        'EURO': '#00796B',   # Teal
        'EMRO': '#C2185B',   # Pink
    },
    'cmap_cso': 'YlOrRd',
    'cmap_ros': 'RdYlGn',
}

# =============================================================================
# Output Paths
# =============================================================================
OUTPUT_DIRS = {
    'processed_sequences': 'data/processed',
    'qc_reports': 'data/qc_reports',
    'scores': 'data/scores',
    'figures': 'figures',
    'tables': 'tables',
    'logs': 'logs',
}

# =============================================================================
# Analysis Parameters
# =============================================================================
ANALYSIS_PARAMS = {
    'bootstrap_iterations': 1000,
    'confidence_level': 0.95,
    'correlation_threshold': 0.7,
    'outlier_iqr_multiplier': 1.5,
}

# =============================================================================
# File Patterns
# =============================================================================
FILE_PATTERNS = {
    'fasta': '*.fasta',
    'genbank': '*.gb',
    'metadata': 'metadata*.csv',
    'output_prefix': 'hiv_capsid_analysis',
}