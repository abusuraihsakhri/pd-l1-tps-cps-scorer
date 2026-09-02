# Pd L1 TPS CPS Scorer

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

PD-L1 TPS & CPS Immunoscore Calculator

Implements PD-L1 scoring for immunotherapy eligibility assessment:
  - Tumor Proportion Score (TPS)
  - Combined Positive Score (CPS)
  - Immune Cell (IC) score

References:
  - Herbst RS et al. Nature. 2014;515(7528):563-567
  - Kulangara K et al. Appl Immunohistochem Mol Morphol. 2019;27(2):99-106
  - FDA pembrolizumab prescribing information

Zero-dependency Python implementation (stdlib only).
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`calculate_tps()`**: Calculate Tumor Proportion Score (TPS).

TPS = (PD-L1 positive tumor cells / Total viable tumor cells) x 100

Parameters
----------
pd_l1_positive_tumor_cells : int
    Number of viable tumor cells showing partial or complete membrane staining.
total_viable_tumor_cells : int
    Total number of viable tumor cells evaluated.

Returns
-------
dict with tps (float), tps_category, description, interpretation
- **`calculate_cps()`**: Calculate Combined Positive Score (CPS).

CPS = (PD-L1 positive cells [tumor + lymphocytes + macrophages] /
       Total viable tumor cells) x 100

The denominator is ONLY tumor cells (not total cells).

Parameters
----------
pd_l1_positive_tumor_cells : int
    Number of PD-L1 positive tumor cells.
pd_l1_positive_lymphocytes : int
    Number of PD-L1 positive lymphocytes.
pd_l1_positive_macrophages : int
    Number of PD-L1 positive macrophages.
total_viable_tumor_cells : int
    Total number of viable tumor cells in the denominator.

Returns
-------
dict with cps (float), cps_category, component counts, interpretation
- **`calculate_ic_score()`**: Calculate Immune Cell (IC) score for atezolizumab (SP142 assay).

IC score = % of tumor area covered by PD-L1 positive immune cells.

Parameters
----------
pd_l1_positive_immune_cells_area : float
    Area occupied by PD-L1 positive immune cells (IC+).
total_tumor_area : float
    Total tumor area evaluated.

Returns
-------
dict with ic_score, ic_category, description
- **`interpret_for_cancer_type()`**: Interpret a PD-L1 score for a specific cancer type.

Parameters
----------
score_value : float
    The TPS or CPS value.
cancer_type : str
    Cancer type key (e.g., 'nsclc', 'gastric', 'tnbc', 'urothelial').

Returns
-------
dict with cancer_type, score_type, score_value, threshold_met, therapy
- **`get_assay_info()`**: Get information about a specific PD-L1 assay.

---

## 📐 Mathematical Formulation & Logic

```text
  Calculate Tumor Proportion Score (TPS).
  Calculate Combined Positive Score (CPS).
  Calculate Immune Cell (IC) score for atezolizumab (SP142 assay).
  IC score = % of tumor area covered by PD-L1 positive immune cells.
  tps_result = calculate_tps(tumor_pos, total)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t pd-l1-tps-cps-scorer .
docker run -p 8000:8000 pd-l1-tps-cps-scorer
```
