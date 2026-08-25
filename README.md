# PD-L1 TPS & CPS Immunoscore Calculator

> **PD-L1 Scoring for Immunotherapy Eligibility**
> Reference: Herbst RS et al. Nature. 2014;515(7528):563-567

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10+-3776AB.svg)

## Overview

Real implementation of PD-L1 scoring for immunotherapy eligibility assessment:

- **TPS (Tumor Proportion Score)**: PD-L1 positive tumor cells / total tumor cells × 100
- **CPS (Combined Positive Score)**: (tumor + lymphocytes + macrophages) / total tumor cells × 100
- **IC Score**: Immune cell area-based scoring for atezolizumab (SP142 assay)
- **Cancer-type interpretation** for NSCLC, gastric, TNBC, urothelial, cervical, HNSCC, esophageal
- **Assay comparison**: 22C3, 28-8, SP142, SP263

## Quick Start

```bash
# Calculate TPS
python pdl1_scorer.py tps --positive 50 --total 100

# Calculate CPS
python pdl1_scorer.py cps --tumor-pos 10 --lymph-pos 5 --macro-pos 3 --total-tumor 100

# Interpret for cancer type
python pdl1_scorer.py interpret --score 60 --cancer-type nsclc

# Get assay info
python pdl1_scorer.py assay --name 22C3

# Batch processing
python pdl1_scorer.py batch -i cases.csv -o results.csv
```

## TPS Categories

| TPS | Category | Clinical Significance |
|-----|----------|----------------------|
| <1% | Negative | No PD-L1 expression |
| 1-49% | Low expression | May benefit from combo therapy |
| ≥50% | High expression | Pembrolizumab monotherapy eligible (NSCLC) |

## CPS Thresholds by Cancer Type

| Cancer Type | Threshold | Therapy |
|-------------|-----------|---------|
| NSCLC | TPS ≥50% | Pembrolizumab monotherapy (1st line) |
| Gastric/GEJ | CPS ≥5 | Pembrolizumab (3rd line) |
| TNBC | CPS ≥10 | Pembrolizumab + chemo (1st line) |
| Urothelial | CPS ≥10 | Pembrolizumab (2nd line) |
| Cervical | CPS ≥1 | Pembrolizumab + chemo (1st line) |
| HNSCC | CPS ≥20 | Pembrolizumab monotherapy (1st line) |

## Assay Comparison

| Assay | Clone | Platform | Score Types | Primary Use |
|-------|-------|----------|-------------|-------------|
| 22C3 | Dako | Dako Autostainer | TPS, CPS | Pembrolizumab companion Dx |
| 28-8 | Dako | Dako Autostainer | TPS | Nivolumab complementary Dx |
| SP142 | Ventana | Ventana BenchMark | TPS, IC | Atezolizumab companion Dx |
| SP263 | Ventana | Ventana BenchMark | TPS | Durvalumab companion Dx |

## Python API

```python
from pdl1_scorer import calculate_tps, calculate_cps, interpret_for_cancer_type

# TPS
result = calculate_tps(pd_l1_positive_tumor_cells=50, total_viable_tumor_cells=100)
print(result["tps"])  # 50.0
print(result["tps_category"])  # "High expression"

# CPS
result = calculate_cps(tumor_pos=10, lymph_pos=5, macro_pos=3, total_tumor=100)
print(result["cps"])  # 18.0

# Cancer-type interpretation
result = interpret_for_cancer_type(60.0, "nsclc")
print(result["therapy"])  # "Pembrolizumab monotherapy (1st line)"
```

## Running Tests

```bash
python -m pytest test_pdl1_scorer.py -v
```

## License

MIT License. See [LICENSE](LICENSE) for details.
