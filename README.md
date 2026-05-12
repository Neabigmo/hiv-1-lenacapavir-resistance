# Lenacapavir HIV-1 Resistance: Three-Tier Evidence-Based Interpretive Framework

A systematic evidence synthesis and surveillance framework for clinical interpretation of lenacapavir resistance-associated mutations (RAMs) in HIV-1.

## Core Scientific Insights

1. **Mutation Identity Drives First-Pass Triage**: The mutation-only model (M1) explains 94.6% of variance in resistance levels (AIC = 40.9, R² = 0.946). Adding subtype as a random effect (M2) does not improve fit, supporting universal first-pass surveillance rules.

2. **N57H and M66I Dominate**: These two mutations show 4,890-fold and 3,200-fold resistance respectively — undisputed Tier 1 candidates. Bootstrap ranking (1,000 resamples) confirms their top positions with high stability (SD < 0.5 ranks).

3. **Subtype Contributes Negligibly**: At current data depth (16 mutations, 26 observations, 4 subtypes), subtype explains no detectable additional variation. However, this conclusion is bounded by data sparsity — not a universal negative.

4. **Q67H+K70R Is the Key Edge Case**: This double mutant exhibits 76-fold context-dependent variability (0.59×–46.3×) across three independent studies, making it the highest-priority candidate for future subtype-stratified phenotyping.

5. **Three-Tier Surveillance Framework**:
   - **Tier 1** (Strong): N57H, M66I — universal first-pass surveillance regardless of subtype
   - **Tier 2** (Moderate): Q67H+K70R — enhanced monitoring for context-sensitive combinations
   - **Tier 3** (Hypothesis): Other mutations — research priority, need more data

6. **Data Sparsity Is the Principal Limitation**: No mutation has ≥3 independent subtype-stratified observations. Evidence supports major mutation-level ranking but is underpowered for formal subtype interaction tests.

## Code Structure

### `pipeline/` — Core Analysis Pipeline

| Script | Description |
|--------|-------------|
| `config.py` | Global configuration: RAM site definitions (L56/N57/M66/Q67/K70/N74/A105/T107), WHO region mappings, CIR-SRI scoring thresholds, QC parameters, color palettes |
| `sequence_qc.py` | HIV-1 sequence quality control: FASTA/GenBank parsing, HXB2 coordinate mapping, hypermutation detection, stop codon/frameshift filtering |
| `cso_calculator.py` | Capsid Sequence Observability (CSO) scoring: evaluates sequence availability per country/region, metadata completeness, temporal trends, subtype coverage |
| `ros_calculator.py` | RAM-site Observability Score (ROS): per-site and per-sequence coverage of lenacapavir RAM positions across global HIV-1 sequences |
| `cir_sri_builder.py` | CIR-SRI (Capsid Inhibitor Resistance Surveillance Readiness Index): 5-dimension composite index (sequence availability, metadata quality, RAM coverage, temporal depth, subtype diversity), each 0-2, total 0-10 |
| `sensitivity_analysis.py` | Sensitivity and robustness testing: leave-one-country-out cross-validation, temporal stability analysis, subtype stratification, near-full-length vs partial sequence comparison |
| `figures.py` | Publication-quality visualization module: world maps, bar charts, heatmaps, radar charts (300+ DPI, journal-ready output) |
| `main.py` | Pipeline orchestrator: runs complete analysis from sequence QC through CIR-SRI calculation, sensitivity analysis, and visualization generation |

### `experiments/` — Figure & Analysis Scripts

| # | Script | Corresponding Analysis |
|---|--------|----------------------|
| 01 | `evidence_synthesis.py` | PRISMA 2020 evidence flow: 87 records → 11 studies included |
| 02 | `data_harmonization.py` | Data standardization: all FC values → log₁₀ scale, context tier assignment |
| 03 | `sequence_conservation.py` | Capsid amino acid conservation scores at each RAM position |
| 04 | `model_comparison.py` | **Core statistical model**: M0–M3 nested model comparison (AIC/BIC), leave-one-study-out CV, bootstrap ranking stability (1,000 resamples) |
| 05 | `sensitivity_analysis.py` | Context stratification, leave-one-subtype-out validation |
| 06 | `epistasis_analysis.py` | Epistatic interaction analysis of double mutants vs single mutant additivity |
| 07 | `compensatory_analysis.py` | M66I-centered compensatory mutation identification |
| 08 | `prepare_structures.py` | Automated PDB download: 6VKV (WT) and 7RAO (M66I) |
| 09 | `foldx_analysis.py` | FoldX computational mutagenesis: ΔΔG calculation for key mutations |
| 10 | `structural_metrics.py` | Structural perturbation metrics from published literature |
| 11 | `generate_figure1_v2.py` | **Fig 1**: PRISMA 2020 flow diagram + Mutation×Subtype×Context heatmap |
| 12 | `generate_figure2.py` | **Fig 2**: Core phenotypic evidence (raincloud, forest plot, ΔAIC comparison, bootstrap ranking) |
| 13 | `generate_figure3.py` | **Fig 3**: Epistasis network (Q67H+K70R context sensitivity, M66I combination network) |
| 14 | `generate_figure4.py` | **Fig 4**: Structural context (hexamer localization, mechanism panel, structure-phenotype correlation r=0.886) |
| 15 | `generate_figure5.py` | **Fig 5**: Evolutionary constraints (conservation, subtype frequencies, fitness cost) |
| 16 | `generate_figure6_v2.py` | **Fig 6**: Three-tier framework (evidence matrix, surveillance algorithm, validation roadmap) |
| 17 | `generate_graphical_abstract.py` | Graphical abstract: complete study workflow overview |
| 18 | `genetic_barrier.py` | Genetic barrier analysis: nucleotide substitution requirements per RAM |
| 19 | `phylogenetic_mapping.py` | Phylogenetic mapping: RAM distribution across HIV-1 subtype tree |
| 20 | `generate_supplementary.py` | Supplementary figures S1–S6 |
| 21 | `generate_figure_s7.py` | Supplementary Figure S7: ranking stability sensitivity analysis |

## Environment & Dependencies

### Conda Environment

```bash
conda create -n sw_mgli python=3.10
conda activate sw_mgli
pip install -r requirements.txt
```

### Key Dependencies (`requirements.txt`)

- **Data science**: `pandas`, `numpy`, `scipy`, `statsmodels`
- **Visualization**: `matplotlib` (≥3.6), `seaborn`, `networkx`
- **Bioinformatics**: `biopython`
- **Structure**: `pymol` (PyMOL open-source, for 3D structural figures)
- **Computational mutagenesis**: `FoldX` (optional, for ΔΔG calculations)
- **PDF/text extraction**: `PyPDF2` (optional)

### LaTeX Compilation (for manuscript)

TeX Live 2025 required. Compile with:

```bash
cd manuscript/
xelatex -interaction=nonstopmode -synctex=1 lenacapavir_revised_v3.tex
# Run twice to resolve cross-references
```

## Data Overview

### `data/curated/`

| File | Description |
|------|-------------|
| `unified_database.csv` | Core unified dataset: 10 columns including mutation, fold-change (numeric), subtype, context tier, source PMID, fitness (WT%), quality score, and provenance notes |

### `data/processed/revision_v2/`

| File | Description |
|------|-------------|
| `harmonized_phenotype_data.csv` | **Primary analysis dataset**: 26 harmonized quantitative observations with context tier, assay type, study source, quality score, observation ID, mutation type, and data provenance tracking |
| `availability_matrix.csv` | Observation counts per mutation × (subtype × context) combination |
| `observation_summary.csv` | Per-mutation summary statistics (n observations, n studies, mean/median log₁₀FC, context tier breakdown) |
| `study_metadata.csv` | Study-level metadata for all 11 included sources |
| `evidence_flow.csv` / `.json` | PRISMA 2020 flow diagram counts (87→78→22→11) |
| `column_documentation.json` | Complete column definitions and data provenance for harmonized dataset |
| `inclusion_criteria.json` | Study inclusion/exclusion criteria applied during evidence synthesis |
| `evidence_metadata.json` | Metadata about the evidence synthesis process and data sources |
| `harmonization_validation.json` | Cross-validation results of the data harmonization procedure |

## How to Reproduce

### Step 1: Data Preparation

```bash
# Literature-extracted data is provided in data/curated/unified_database.csv
# To re-extract from raw literature CSVs:
python experiments/01_evidence_synthesis.py
python experiments/02_data_harmonization.py
```

### Step 2: Core Analysis

```bash
# Run model comparison (M0-M3, bootstrap):
python experiments/04_model_comparison.py

# Run epistasis analysis:
python experiments/06_epistasis_analysis.py

# Run structural analysis (requires FoldX):
python experiments/09_foldx_analysis.py
python experiments/10_structural_metrics.py
```

### Step 3: Generate Figures

```bash
# Generate all main figures:
python experiments/11_generate_figure1_v2.py
python experiments/12_generate_figure2.py
# ... through 16_generate_figure6_v2.py
```

### Step 4: Full Pipeline (optional, requires HIV sequence data)

```bash
# Requires hiv-db.fasta (590 MB, not included — see Data Access below)
python pipeline/main.py
```

## Data Access

### Literature Data (Included)
All manually extracted quantitative data from 11 published studies is provided in `data/curated/unified_database.csv`. Sources include:

- CAPELLA 2-year resistance data (PMID 39873394)
- CALIBRATE week 28 resistance data
- mBio 2022 structural mechanism (PMID 36000858)
- PMC12077089 (L56V/N57H high-fold resistance)
- PMC11995365 (RevLun phenotypic analysis)
- NATAP 2022 lenacapavir resistance
- PMC9039614 (capsid diversity/conservation)
- Additional PMC and clinical trial sources

### PDB Structures (Auto-downloaded)

```bash
python experiments/08_prepare_structures.py
```

Downloads: 6VKV (WT lenacapavir-CA hexamer) and 7RAO (M66I mutant capsid) from rcsb.org.

### HIV Sequence Data (Requires Registration)

The full HIV sequence dataset (hiv-db.fasta, 590 MB, 130,704 sequences) is **not included** due to size and licensing. To obtain:

1. Visit [LANL HIV Sequence Database](https://www.hiv.lanl.gov/)
2. Request access to PR/RT fragment sequences
3. Download to `data/raw/hiv-db.fasta`

**Note**: The current hiv-db.fasta contains only 693nt PR/RT fragments and does **not** include capsid region sequences. The CIR-SRI pipeline can process any HIV-1 sequence data with capsid coverage.

## Citation

If you use this code or data in your research, please cite:

> Yang Y. Clinical Interpretation Framework for Lenacapavir Resistance in HIV-1: A Three-Tier Evidence-Based Interpretive Framework. 2025.

## License

This project is intended for academic and research use. Data extracted from published literature remains subject to the original publishers' terms.
