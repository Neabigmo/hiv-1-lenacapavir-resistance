"""Shared paths and configuration for the lenacapavir resistance analysis."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
REV2_DATA = PROCESSED_DIR / "revision_v2"

RESULTS_DIR = ROOT / "results"
REV2_RESULTS = RESULTS_DIR / "revision_v2"

STRUCT_DIR = DATA_DIR / "structures"
FOLDX_DIR = ROOT / "tools" / "foldx"

# Ensure output dirs exist
for d in [PROCESSED_DIR, REV2_DATA, RESULTS_DIR, REV2_RESULTS, STRUCT_DIR]:
    d.mkdir(parents=True, exist_ok=True)
