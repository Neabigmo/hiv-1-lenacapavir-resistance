# Key Scientific Insights

## 1. Mutation Identity Alone Explains ~95% of Resistance Variance

**Finding**: The mutation-only linear model (M1: log₁₀FC = β₀ + Σβₖ·mutₖ) achieves AIC = 40.9, R² = 0.946. Adding subtype as a random effect (M2) does not improve fit (ΔAIC = +2.1). Adding study-level random effects (M3) does not meaningfully change results (ΔAIC = +1.8).

**Implication**: For first-pass surveillance triage, mutation identity is sufficient. Clinical laboratories do not need to wait for subtype information before flagging a high-resistance mutation.

**Caveat**: This conclusion is explicitly bounded by the small sample (16 mutations, 26 observations). It means "subtype does not contribute detectably at this data depth," not "subtype never contributes."

## 2. N57H and M66I: Undisputed Tier 1 Mutations

| Mutation | Median FC | Log₁₀FC | Bootstrap Rank (mean ± SD) | N Studies |
|----------|-----------|---------|---------------------------|-----------|
| N57H | 4,890× | 3.69 | 1.0 ± 0.0 | 3 |
| M66I | 3,200× | 3.51 | 2.0 ± 0.0 | 3 |
| L56V | 7× | 0.85 | 4.2 ± 1.1 | 2 |

**Key**: N57H and M66I hold rank 1 and 2 with zero standard deviation across 1,000 bootstrap resamples. They are the only mutations with replicated observations across multiple independent studies.

## 3. Q67H+K70R: The Critical Edge Case

The double mutant Q67H+K70R displays an extraordinary 76-fold range across three independent studies:

| Study | Subtype | Fold-Change | Context |
|-------|---------|-------------|---------|
| CAPELLA 2yr (PMID 39873394) | B subtype | 0.59× | clinical |
| CALIBRATE wk28 | Mixed | 10.4× | clinical |
| Uganda clinical isolates | A1/D | 46.3× | clinical |

This range cannot be explained by mutation identity alone. It represents the highest-priority target for future subtype-stratified phenotyping.

## 4. Three-Tier Interpretive Framework

### Tier 1 — Strong Evidence (Universal First-Pass)

| Mutation | Criteria |
|----------|----------|
| N57H | ≥3 independent replications, bootstrap rank SD = 0, structural coherence (direct binding pocket), highly conserved (99.5%) |
| M66I | ≥3 independent replications, bootstrap rank SD = 0, structural mechanism validated (steric hindrance, PMID 36000858), M66 is the most conserved position (99.95%) |

**Action**: Flag immediately in any subtype. No confirmatory testing needed before clinical alert.

### Tier 2 — Moderate Evidence (Enhanced Monitoring)

| Mutation | Criteria |
|----------|----------|
| Q67H+K70R | 76-fold context-dependent range; requires subtype-stratified interpretation |
| Q67H+N74D | Structural synergy (ΔΔG > single mutants additive); limited replication |

**Action**: Enhanced monitoring. Consider genetic background when interpreting.

### Tier 3 — Hypothesis-Generating (Research Priority)

| Mutation | Reason |
|----------|--------|
| L56V | Only 2 observations, high variability |
| N74D | Only 1 observation in isolation |
| K70R | Only combination data |
| A105T, T107A | Single observations, low fold-change |
| All other combinations | n ≤ 1 per combination |

**Action**: Targeted validation needed. Do not use for clinical decisions.

## 5. Data Sparsity Quantified

| Metric | Value |
|--------|-------|
| Total observations | 26 |
| Unique mutations | 16 |
| Subtypes with data | 4 (B, C, CRF02_AG, D) out of >10 major subtypes |
| Complete cases (mutation + subtype + FC) | 23 / 26 |
| Mutations with ≥3 observations | 2 (N57H, M66I) |
| Subtype-stratified observations per mutation | 0–3 |

**Bottom line**: Current evidence is sufficient for ranking major mutation-level effects but severely underpowered for formal subtype interaction tests. The near-zero subtype variance should be interpreted as "unresolved" rather than "absent."

## 6. Structural & Evolutionary Corroboration

The mutation-first hierarchy is supported by independent structural and evolutionary evidence:

- **Binding pocket location**: N57H and M66I directly contact the LEN binding site; Q67H and K70R are at the pocket rim; N74D and A105T are distal
- **Conservation**: M66 (99.95%) > N57 (99.5%) > Q67 (98%) > L56 (96%) > N74 (97%) > K70 (92%)
- **Structural mechanism**: M66I causes steric hindrance (validated by 7RAO structure); N57H alters pocket electrostatics
- **Structure-phenotype correlation**: Residue-level structural perturbation correlates with resistance phenotype (r = 0.886, p < 0.01)

## 7. Recommendations for Future Research

1. **Expand subtype coverage**: Prioritize Q67H+K70R testing in non-B subtypes (particularly C, CRF01_AE, A1)
2. **Standardize assay platforms**: Current data spans PhenoSense, MT4, and primary cell assays — harmonization adds uncertainty
3. **Publish negative results**: Studies testing LEN RAMs and finding <2-fold changes are underrepresented in the literature
4. **Dual-report fitness and FC**: Only 3 observations currently report both — insufficient for fitness-adjusted resistance thresholds
