# HIV-1 Lenacapavir Resistance — Evidence-Based Surveillance Prioritization

A systematic evidence synthesis and three-tier framework for prioritizing HIV-1 capsid resistance mutations for surveillance in the context of lenacapavir (LEN) rollout for treatment and PrEP.

## Key Finding

**Mutation identity alone drives resistance ranking** — subtype and study effects are negligible in available data (24 observations, 13 genotypes). N57H (4,890-fold) and M66I (3,200-fold) are the most stable single-mutation resistance signals.

## Three-Tier Framework

| Tier | Mutations | Recommendation |
|------|-----------|---------------|
| 1A | N57H, M66I | Subtype-agnostic confirmatory testing |
| 1B | K70N+N74K, Q67H+N74D | Enhanced monitoring, document genetic context |
| 2 | Q67H+K70R | Context-aware phenotyping (78.5× range) |
| 3 | A105T, T107A, compensatory | Research priority |

## Structure

```
├── src/                # Core analysis (4 modules)
│   ├── harmonize.py    # Evidence synthesis → harmonized dataset
│   ├── analyze.py      # Models, bootstrap, CV, sensitivity, epistasis
│   ├── conservation.py # Conservation & polymorphism tables
│   └── structure.py    # Structural perturbation + FoldX wrapper
├── data/               # Data directory (see DATA_ACQUISITION.md)
├── results/            # Generated output tables (CSV, JSON)
├── DATA_ACQUISITION.md # Full data provenance
├── requirements.txt    # Python dependencies
└── .gitignore
```

## Quick Start

```bash
pip install -r requirements.txt
python -m src.harmonize       # → data/processed/revision_v2/
python -m src.analyze          # → results/revision_v2/
python -m src.conservation     # (optional)
python -m src.structure        # (optional)
```

## Methods

- **Evidence**: PRISMA-informed synthesis; 87 records → 11 sources → 24 quantitative observations
- **Models**: Nested linear models M0–M3 (intercept → mutation → mutation+subtype → mutation+study)
- **Bootstrap**: 1,000-resample ranking stability
- **Epistasis**: Additive vs synergistic vs compensatory classification
- **Structural**: Perturbation scoring from published mBio 2022 data

All data from published literature. See [DATA_ACQUISITION.md](DATA_ACQUISITION.md).

## Citation

Yang Y. Evidence-Based Prioritization of HIV-1 Lenacapavir Resistance Mutations for Surveillance Triage. *Journal of Medical Virology* (submitted).
