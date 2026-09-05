# PD-L1 TPS & CPS Immunoscore Calculator

> **Domain:** Clinical Decision Support & Biomedical Computing
> **Reference Guidelines & Standards:** Standard Clinical Formulations & ISO/IEC Quality Frameworks

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## What It Does

PD-L1 TPS & CPS Immunoscore Calculator implements PD-L1 scoring for immunotherapy eligibility assessment:

- **Tumor Proportion Score (TPS)**
- **Combined Positive Score (CPS)**
- **Immune Cell (IC) score**

### References

- Herbst RS et al. Nature. 2014;515(7528):563-567
- Kulangara K et al. Appl Immunohistochem Mol Morphol. 2019;27(2):99-106
- FDA pembrolizumab prescribing information

Zero-dependency Python implementation (stdlib only for core scoring).
License: MIT

---

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/pd-l1-tps-cps-scorer.git
cd pd-l1-tps-cps-scorer

# Install dependencies (for enterprise features and testing)
pip install fastapi uvicorn pydantic pytest
```

---

## Key Capabilities & Algorithmic Modules

### Analytical Functions

- **`calculate_tps()`**: Calculate Tumor Proportion Score (TPS).
  - TPS = (PD-L1 positive tumor cells / Total viable tumor cells) x 100
  - Returns: dict with tps, tps_category, description, interpretation

- **`calculate_cps()`**: Calculate Combined Positive Score (CPS).
  - CPS = (PD-L1 positive cells [tumor + lymphocytes + macrophages] / Total viable tumor cells) x 100
  - Returns: dict with cps, cps_category, component counts, interpretation

- **`calculate_ic_score()`**: Calculate Immune Cell (IC) score for atezolizumab (SP142 assay).
  - IC score = % of tumor area covered by PD-L1 positive immune cells
  - Returns: dict with ic_score, ic_category, description

- **`interpret_for_cancer_type()`**: Interpret a PD-L1 score for a specific cancer type.
  - Supports: nsclc, gastric, tnbc, urothelial, cervical, esophageal, hnscc
  - Returns: dict with cancer_type, score_type, score_value, threshold_met, therapy

- **`get_assay_info()`**: Get information about a specific PD-L1 assay.
  - Supports: 22C3, 28-8, SP142, SP263

---

## CLI Quickstart & Usage

### 1. Calculate TPS
```bash
python cli.py tps --positive 50 --total 100
```

### 2. Calculate CPS
```bash
python cli.py cps --tumor-pos 10 --lymph-pos 5 --macro-pos 3 --total-tumor 100
```

### 3. Calculate IC Score
```bash
python cli.py ic --ic-area 5.0 --tumor-area 100.0
```

### 4. Interpret for Cancer Type
```bash
python cli.py interpret --score 60 --cancer-type nsclc
```

### 5. Get Assay Information
```bash
python cli.py assay --name 22C3
```

### 6. Batch Process CSV
```bash
python cli.py batch -i sample.csv -o results.csv
```

### 7. Enterprise Supervisor Audit
```bash
python cli.py audit --task-id TASK-001 --primary 25.0 --secondary 8.0
```

### 8. Supervisory Chat
```bash
python cli.py chat "Explain PD-L1 scoring"
```

### 9. Verify Audit Trail
```bash
python cli.py verify-audit
```

### 10. Start REST API Server
```bash
python cli.py serve --host 0.0.0.0 --port 8000
```

---

## REST API Endpoints

When running the server (`python cli.py serve`), the following endpoints are available:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health and metadata check |
| `/metrics` | GET | Prometheus operational metrics |
| `/api/audit` | POST | Dispatch task payload across specialized workers |
| `/api/chat` | POST | Air-gapped supervisory conversational assistant |
| `/api/audit/logs` | GET | Retrieve and verify HMAC audit trail |

---

## Security & Enterprise Architecture

- **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
- **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation.
- **Secure Key Management:** Audit keys generated via `secrets.token_hex()` or configured via `AUDIT_SECRET_KEY` environment variable.
- **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics.

---

## Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py 1000
```

---

## Container Deployment

```bash
docker build -t pd-l1-tps-cps-scorer .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key pd-l1-tps-cps-scorer
```

Or using Docker Compose:

```bash
docker-compose up -d
```

---

## Project Structure

```
pd-l1-tps-cps-scorer/
├── pdl1_scorer.py      # Core scoring algorithms & CLI
├── cli.py              # CLI entry point
├── enrichment.py       # Enrichment feature engines
├── simulator.py        # High-throughput simulation
├── agents/             # Enterprise agent framework
│   ├── base.py         # Security, PHI guard, audit trail
│   ├── models.py       # Pydantic schemas
│   ├── supervisor.py   # Supervisor orchestrator
│   ├── workers.py      # Specialized worker agents
│   ├── api.py          # FastAPI REST server
│   ├── learning.py     # Bayesian calibration engine
│   ├── metrics.py      # Prometheus metrics
│   └── streamer.py     # WebSocket telemetry
├── web/                # Operations console (HTML)
├── tests/              # Test suite
├── Dockerfile          # Container definition
├── docker-compose.yml  # Compose configuration
└── openapi_spec.json   # OpenAPI 3.1 specification
```
