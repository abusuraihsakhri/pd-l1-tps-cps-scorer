# PD-L1 TPS & CPS Scorer

A numerical PD-L1 immunohistochemistry scoring calculator for research, teaching and manual pathology review. Computes Tumor Proportion Score (TPS), Combined Positive Score (CPS) and area-based immune-cell (IC) percentage.

## Browser application

Run `python -m http.server 8001` at the repository root and visit `http://localhost:8001/`; the root page redirects to `web/`. The static site is also structured for GitHub Pages deployment from the `master` branch. Confirm live publishing in the repository's Pages settings; a working public endpoint is not assumed.

The browser calculator:
- Uses TPS = 100 × PD-L1-positive viable tumor cells / total viable tumor cells.
- Uses CPS = 100 × (positive tumor cells + positive lymphocytes + positive macrophages) / total viable tumor cells, **capped at 100**.
- Uses IC = 100 × PD-L1-positive immune-cell area / total evaluated tumor area.
- Validates inputs, displays descriptive expression strata and flags counts below 100 viable tumor cells. Categories use **unrounded** scores.
- Supports local CSV import and scored CSV export (1,000 cases / 2 MB maximum).

**Privacy:** All browser calculations and CSV imports stay within the browser. No backend submission, tracking API or persistent browser storage is used. Do not include patient identifiers.

## Command line

Python 3.10 or later. Core calculations use the standard library.

```bash
python cli.py tps --positive 50 --total 100
python cli.py cps --tumor-pos 10 --lymph-pos 5 --macro-pos 3 --total-tumor 100
python cli.py ic --ic-area 5 --tumor-area 100
python cli.py interpret --score 50 --cancer-type nsclc
python cli.py assay --name 28-8
python cli.py batch --input sample.csv --output results.csv
```

Batch CSV requires columns `pd_l1_positive_tumor_cells` and `total_viable_tumor_cells`. Optional columns: `pd_l1_positive_lymphocytes`, `pd_l1_positive_macrophages` and `cancer_type`. For NSCLC, the descriptive score strata use **TPS**, not CPS. No treatment recommendations are generated.

## Optional REST API and supervisor

```bash
python -m pip install "fastapi>=0.110,<1" "uvicorn>=0.29,<1" "pydantic>=2,<3"
python cli.py serve --host 127.0.0.1 --port 8000
```

The endpoints are `/health`, `/metrics`, `/api/audit`, `/api/chat` and `/api/audit/logs`. The supervisor is a generic demonstration system, **not** a validated PD-L1 scoring engine or clinical decision-support system. The chat response is deterministic mock output. Do not use this API to process protected health information.

Audit log retrieval requires an environment-configured `API_TOKEN` and an `Authorization: Bearer <API_TOKEN>` header. For Docker Compose, copy `.env.example` to a private `.env`, configure strong random `AUDIT_SECRET_KEY` and `API_TOKEN` values, and run `docker compose up --build`. Do not commit .env.

## Tests

```bash
python -m pip install "fastapi>=0.110,<1" "uvicorn>=0.29,<1" "pydantic>=2,<3" "pytest>=8,<10"
pytest -v
node --check web/app.js
node --test web/app.test.cjs
```

GitHub Actions runs Python tests on 3.10, 3.11 and 3.12 and JavaScript tests under Node.js 22. The website uses standard HTML/CSS/JavaScript and is intended for modern mobile and desktop browsers; it does not require server-side Python or Pyodide.

## Clinical limitations

The numerical results rely on qualified pathological assessment, validated staining, correct cell identification and assay-specific specimen adequacy requirements. This tool does **not** review images, validate assays, determine treatment eligibility or issue therapy recommendations. Cutoffs and prescribing indications depend on tumor type, assay, line of treatment and current approved labeling.

Scoring references: [Agilent PD-L1 IHC 22C3 pharmDx overview](https://www.agilent.com/en/product/pharmdx/pd-l1-ihc-22c3-pharmdx-overview) and [official intended-use details](https://www.agilent.com/en-us/product/pharmdx/pd-l1-ihc-22c3-pharmdx/pd-l1-ihc-22c3-pharmdx-product-page).

## License

[MIT License](LICENSE).
