# Data Acquisition

This study is a secondary analysis of published data. No primary experiments were conducted.

## Phenotypic Resistance Data

| Source | Data Type | Obtained From | Access |
|--------|-----------|---------------|--------|
| Andreatta et al. 2024 | In vitro selection FC | Nature Communications, PMID: PMC12077089 | Open access |
| Yant et al. 2019 / mBio 2022 | SDM FC values | Nature Medicine / mBio, PMID: 31501601 / PMC9600929 | Open access |
| CAPELLA trial investigators 2025 | Clinical isolate FC | J Infect Dis, PMID: 38153326 | Subscription |
| Margot et al. 2022 | Clinical isolate FC | J Antimicrob Chemother, DOI: 10.1093/jac/dkab470 | Subscription |
| Uganda A1/D study (2025) | Natural polymorphism FC | J Antimicrob Chemother | Subscription |
| Link et al. 2020 | SDM FC values | Nature, PMID: 32612233 | Subscription |

**Extraction method**: Fold-change (FC) values were extracted from tables, figures (via WebPlotDigitizer where needed), and supplementary materials. All values recorded as EC50(mutant)/EC50(wild-type).

## Sequence Conservation Data

| Source | Content | Access |
|--------|---------|--------|
| mBio 2025 capsid polymorphism study | >10,000 sequences conservation scores | Open access |
| FDA SUNLENCA / YEZTUGO prescribing labels | Natural polymorphism frequencies | FDA website |
| Durand Nka et al. 2023 | Subtype-specific frequencies | J Antimicrob Chemother |

## Structural Data

| PDB ID | Description | Source |
|--------|-------------|--------|
| 6VKV | LEN–CA hexamer complex (WT) | RCSB PDB (https://www.rcsb.org/structure/6VKV) |

Structural perturbation data (H-bond loss, steric clash, ΔΔG binding) were extracted from mBio 2022 supplementary materials.

## Literature Search

- **Database**: PubMed
- **Date**: 2024-11-15
- **Query**: `(lenacapavir OR GS-6207) AND (resistance OR capsid OR mutation OR fold-change)`
- **Coverage**: 2020–2026
- **Additional sources**: WHO HIVDR reports, FDA prescribing labels, conference abstracts (CAPELLA, CALIBRATE)

## Data Processing

All extracted values were harmonized to log₁₀(fold-change) scale. Context tiers assigned as:
- **Clinical**: patient-derived isolates
- **In vitro**: site-directed mutagenesis or in vitro selection
- **Natural polymorphism**: baseline prevalence surveys

## File Format

Place raw data CSVs in `data/raw/` with columns: `Mutation, FC_numeric, Subtype, Context, Source, Quality`.
The harmonization pipeline (`src/harmonize.py`) reads from `data/processed/revision/` by default.
