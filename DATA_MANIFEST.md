# Data Manifest — Lenacapavir Resistance Evidence Synthesis

This document catalogs all data sources used in this project, their provenance, contents, and which analyses consume them.

---

## 1. Literature-Extracted Quantitative Data

### 1.1 Unified Database

**File**: `data/curated/unified_database.csv`

| Column | Description |
|--------|-------------|
| `Mutation` | Capsid mutation name (e.g., N57H, M66I, Q67H+K70R) |
| `FC_numeric` | Fold-change (EC50_mutant / EC50_WT), numeric |
| `Subtype` | HIV-1 subtype (B, C, CRF02_AG, D, or untyped) |
| `Context` | Evidence context tier (clinical, in_vitro, natural_polymorphism) |
| `Source_PMID` | PubMed ID or study source identifier |
| `Fitness_WT_pct` | Relative fitness as percentage of wild-type |
| `Notes` | Free-text notes about data quality or special conditions |
| `Quality_score` | 1–3 quality rating (3 = highest rigor) |

**Used by**: 01_evidence_synthesis, 02_data_harmonization, 04_model_comparison, 06_epistasis_analysis

### 1.2 Harmonized Phenotype Dataset

**File**: `data/processed/revision_v2/harmonized_phenotype_data.csv`

The primary analysis-ready dataset, created by `02_data_harmonization.py`. All fold-change values are transformed to log₁₀ scale.

| Column | Description |
|--------|-------------|
| `mutation` | Mutation name |
| `log10_fc` | log₁₀(fold-change) — harmonized primary endpoint |
| `fc_numeric` | Raw fold-change value |
| `subtype` | HIV-1 subtype |
| `context_tier` | Context tier: clinical, in_vitro, natural_polymorphism |
| `assay_type` | Assay platform used (e.g., PhenoSense, MT4) |
| `study_source` | Source study name |
| `quality_score` | 1–3 quality rating |
| `observation_id` | Unique observation identifier (format: OBS_####) |
| `mutation_type` | single / combination |
| `data_provenance` | Full provenance chain |
| `tier_label` | Claim tier: Strong / Moderate / Hypothesis |
| `nucleotide_changes` | Nucleotide substitutions required |

**Used by**: ALL experiment scripts (04–21)

### 1.3 Observation Summary

**File**: `data/processed/revision_v2/observation_summary.csv`

Per-mutation aggregated statistics:

| Column | Description |
|--------|-------------|
| `mutation` | Mutation name |
| `n_observations` | Number of quantitative observations |
| `n_studies` | Number of independent studies |
| `mean_log10fc` | Mean log₁₀ fold-change |
| `median_log10fc` | Median log₁₀ fold-change |
| `min_log10fc` | Minimum log₁₀ fold-change |
| `max_log10fc` | Maximum log₁₀ fold-change |
| `clinical_tier_n` | Number of clinical-tier observations |
| `in_vitro_tier_n` | Number of in vitro observations |

**Used by**: figure2, figure_s2

### 1.4 Availability Matrix

**File**: `data/processed/revision_v2/availability_matrix.csv`

Observation counts per mutation × (subtype × context) combination. Used to visualize data sparsity patterns across subtypes and experimental contexts.

| Column | Description |
|--------|-------------|
| `mutation` | Mutation name |
| `subtype` | HIV-1 subtype |
| `context` | Context tier |
| `n_observations` | Number of observations |

**Used by**: figure1 (Panel B), figure_s2

### 1.5 Study Metadata

**File**: `data/processed/revision_v2/study_metadata.csv`

| Column | Description |
|--------|-------------|
| `study_name` | Study identifier |
| `pmid` | PubMed ID |
| `year` | Publication year |
| `assay_type` | Phenotypic assay type |
| `subtypes_tested` | Subtypes included |
| `mutations_reported` | Mutations with quantitative data |
| `context_tier` | Evidence context |

**Used by**: evidence_flow, figure1

### 1.6 Inclusion Criteria

**File**: `data/processed/revision_v2/inclusion_criteria.json`

Study inclusion/exclusion criteria used during evidence synthesis.

| Field | Description |
|-------|-------------|
| `inclusion_criteria` | Criteria for study inclusion (e.g., quantitative FC data, clinical isolates) |
| `exclusion_criteria` | Criteria for study exclusion (e.g., review articles, no original data) |
| `applied_at` | Stage of PRISMA flow where criteria were applied |

**Used by**: 01_evidence_synthesis

### 1.7 Evidence Metadata

**File**: `data/processed/revision_v2/evidence_metadata.json`

Metadata about the evidence synthesis process and data sources.

| Field | Description |
|-------|-------------|
| `total_records` | Number of records identified (87) |
| `studies_included` | Number of studies included (11) |
| `total_observations` | Total quantitative observations (26) |
| `synthesis_date` | Date of evidence synthesis |

**Used by**: 01_evidence_synthesis

### 1.8 Harmonization Validation

**File**: `data/processed/revision_v2/harmonization_validation.json`

Cross-validation results of the data harmonization procedure.

| Field | Description |
|-------|-------------|
| `validation_method` | Approach used for validation |
| `consistency_metrics` | Metrics assessing cross-study consistency |
| `flagged_observations` | Observations requiring special attention |

**Used by**: 02_data_harmonization, 05_sensitivity_analysis

---

## 2. Source Data by Study (11 Included Studies)

| # | Source | Year | PMID / DOI | Mutations | Subtypes | Context |
|---|--------|------|------------|-----------|----------|---------|
| 1 | CAPELLA 2-year | 2024 | 39873394 | Q67H+K70R | B, Mixed | clinical |
| 2 | CALIBRATE wk28 | 2024 | (trial) | M66I, Q67H, Q67H+K70R, N57H | B | clinical |
| 3 | mBio 2022 | 2022 | 36000858 | Q67H, N74D, Q67H+N74D | B | in_vitro |
| 4 | PMC12077089 | 2025 | 12077089 | L56V, N57H | D, CRF02_AG | in_vitro |
| 5 | PMC11995365 | 2025 | 11995365 | M66I, A105T, T107A, M66I+A105T, M66I+T107A | B | in_vitro |
| 6 | NATAP 2022 | 2022 | — | L56V, N57H, M66I, Q67H, N74D, K70R, K70N | B | in_vitro |
| 7 | PMC9039614 | 2022 | 9039614 | Capsid diversity & conservation | Multiple | natural |
| 8 | Uganda study | 2025 | — | Q67H+K70R, K70N+N74K, Q67K+K70H | A1, D | clinical |
| 9 | Clinical isolates | 2024 | — | M66I+N74D+A105T, M66I+A105T, M66I+T107A | Clinical | clinical |
| 10 | Subtype B panel | 2024 | — | Q67H+T107N, Q67H+N74S, K436E+I437T | Subtype_B | in_vitro |
| 11 | Primary T cells | 2024 | — | L56V+N57H | Primary | in_vitro |

---

## 3. Structural Data

### 3.1 PDB Structures

| File | Source | Description | Used By |
|------|--------|-------------|---------|
| `6VKV.pdb` | rcsb.org (auto-downloaded) | WT lenacapavir-CA hexamer complex | 08_prepare_structures, 09_foldx_analysis, figure4 |
| `7RAO.pdb` | rcsb.org (auto-downloaded) | M66I mutant capsid | 08_prepare_structures, 09_foldx_analysis, figure4 |

### 3.2 FoldX Analysis Output

| File | Description |
|------|-------------|
| `6VKV_Repair.pdb` | FoldX-repaired WT structure |
| `6VKV_Repair.fxout` | FoldX repair energy output |
| `rotabase.txt` | FoldX rotamer library (2.6 MB) |

---

## 4. HIV Sequence Data (Not Included in Repository)

**File**: `data/raw/hiv-db.fasta` (590 MB, 130,704 sequences, 693nt PR/RT fragments)
**Source**: LANL HIV Sequence Database
**Access**: Requires registration at https://www.hiv.lanl.gov/
**Note**: Contains only PR/RT fragments — does NOT include capsid region sequences.
**Pipeline usage**: sequence_qc.py, cso_calculator.py, ros_calculator.py, cir_sri_builder.py

---

## 5. Derived Data Files

### 5.1 Genetic Barrier Analysis

Generated by `experiments/18_genetic_barrier.py`:

| File | Description |
|------|-------------|
| `genetic_barrier_scores.csv` | Per-subtype genetic barrier: nucleotide changes, barrier class (low/medium/high), fitness cost, relative escape probability |
| `genetic_barrier_summary.csv` | Aggregated genetic barrier by mutation across subtypes |

### 5.2 Phylogenetic Analysis

Generated by `experiments/19_phylogenetic_mapping.py`:

| File | Description |
|------|-------------|
| `phylogenetic_distances.csv` | Pairwise genetic distances between subtypes |
| `phylogenetic_mutation_matrix.csv` | Mutation × subtype matrix with codon, fold-change, tier, geographic distribution |
| `phylogenetic_tree.csv` | Distance matrix for tree construction |

### 5.3 Evidence Flow

Generated by `01_evidence_synthesis.py`:

| File | Description |
|------|-------------|
| `evidence_flow.csv` | PRISMA flow counts: 87 identified, 78 screened, 22 full-text, 11 included |
| `evidence_flow.json` | Same data in JSON format |
| `study_metadata.csv` | Study-level metadata for 11 included studies |

---

## 6. Data Flow Diagram

```
raw/papers/ (41 CSV files)
    │
    ▼
01_evidence_synthesis.py ──► evidence_flow.csv, study_metadata.csv
    │
    ▼
02_data_harmonization.py ──► harmonized_phenotype_data.csv
    │                              │
    ▼                              ▼
03_sequence_conservation.py    availability_matrix.csv
04_model_comparison.py         observation_summary.csv
05_sensitivity_analysis.py
06_epistasis_analysis.py            │
07_compensatory_analysis.py         ▼
08_prepare_structures.py     Figure generation scripts (11-21)
09_foldx_analysis.py               │
10_structural_metrics.py           ▼
                             manuscript/figures/*.pdf
```
