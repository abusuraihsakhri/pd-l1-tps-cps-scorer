#!/usr/bin/env python3
"""
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
"""

import argparse
import csv
import json
import sys
from typing import Dict, Any, Optional


# ---------------------------------------------------------------------------
# TPS Calculation
# ---------------------------------------------------------------------------

def calculate_tps(pd_l1_positive_tumor_cells: int,
                  total_viable_tumor_cells: int) -> Dict[str, Any]:
    """
    Calculate Tumor Proportion Score (TPS).

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
    """
    if pd_l1_positive_tumor_cells < 0:
        raise ValueError(f"pd_l1_positive_tumor_cells must be >= 0; got {pd_l1_positive_tumor_cells}")
    if total_viable_tumor_cells <= 0:
        raise ValueError(f"total_viable_tumor_cells must be > 0; got {total_viable_tumor_cells}")
    if pd_l1_positive_tumor_cells > total_viable_tumor_cells:
        raise ValueError("Positive cells cannot exceed total cells")

    tps = round((pd_l1_positive_tumor_cells / total_viable_tumor_cells) * 100, 1)

    if tps < 1.0:
        category = "Negative"
        description = "TPS <1%: No PD-L1 expression detected."
    elif tps < 50.0:
        category = "Low expression"
        description = "TPS 1-49%: Low PD-L1 expression."
    else:
        category = "High expression"
        description = "TPS >=50%: High PD-L1 expression (strong positive)."

    return {
        "score_type": "TPS",
        "tps": tps,
        "pd_l1_positive_tumor_cells": pd_l1_positive_tumor_cells,
        "total_viable_tumor_cells": total_viable_tumor_cells,
        "tps_category": category,
        "description": description,
    }


# ---------------------------------------------------------------------------
# CPS Calculation
# ---------------------------------------------------------------------------

def calculate_cps(pd_l1_positive_tumor_cells: int,
                  pd_l1_positive_lymphocytes: int,
                  pd_l1_positive_macrophages: int,
                  total_viable_tumor_cells: int) -> Dict[str, Any]:
    """
    Calculate Combined Positive Score (CPS).

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
    """
    if pd_l1_positive_tumor_cells < 0:
        raise ValueError("pd_l1_positive_tumor_cells must be >= 0")
    if pd_l1_positive_lymphocytes < 0:
        raise ValueError("pd_l1_positive_lymphocytes must be >= 0")
    if pd_l1_positive_macrophages < 0:
        raise ValueError("pd_l1_positive_macrophages must be >= 0")
    if total_viable_tumor_cells <= 0:
        raise ValueError(f"total_viable_tumor_cells must be > 0; got {total_viable_tumor_cells}")

    total_positive = (pd_l1_positive_tumor_cells +
                      pd_l1_positive_lymphocytes +
                      pd_l1_positive_macrophages)

    cps = round((total_positive / total_viable_tumor_cells) * 100, 1)

    if cps < 1.0:
        category = "Negative"
        description = "CPS <1: No PD-L1 expression detected."
    elif cps < 10.0:
        category = "Low expression"
        description = "CPS 1-9: Low PD-L1 expression."
    elif cps < 50.0:
        category = "Moderate expression"
        description = "CPS 10-49: Moderate PD-L1 expression."
    else:
        category = "High expression"
        description = "CPS >=50: High PD-L1 expression."

    return {
        "score_type": "CPS",
        "cps": cps,
        "pd_l1_positive_tumor_cells": pd_l1_positive_tumor_cells,
        "pd_l1_positive_lymphocytes": pd_l1_positive_lymphocytes,
        "pd_l1_positive_macrophages": pd_l1_positive_macrophages,
        "total_positive_cells": total_positive,
        "total_viable_tumor_cells": total_viable_tumor_cells,
        "cps_category": category,
        "description": description,
    }


# ---------------------------------------------------------------------------
# IC Score Calculation
# ---------------------------------------------------------------------------

def calculate_ic_score(pd_l1_positive_immune_cells_area: float,
                       total_tumor_area: float) -> Dict[str, Any]:
    """
    Calculate Immune Cell (IC) score for atezolizumab (SP142 assay).

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
    """
    if pd_l1_positive_immune_cells_area < 0:
        raise ValueError("pd_l1_positive_immune_cells_area must be >= 0")
    if total_tumor_area <= 0:
        raise ValueError(f"total_tumor_area must be > 0; got {total_tumor_area}")
    if pd_l1_positive_immune_cells_area > total_tumor_area:
        raise ValueError("IC area cannot exceed total tumor area")

    ic = round((pd_l1_positive_immune_cells_area / total_tumor_area) * 100, 1)

    if ic < 1.0:
        category = "IC0"
        description = "IC <1%: No PD-L1 positive immune cells."
    elif ic < 5.0:
        category = "IC1"
        description = "IC 1-4%: Low immune cell PD-L1 expression."
    else:
        category = "IC2/3"
        description = "IC >=5%: High immune cell PD-L1 expression."

    return {
        "score_type": "IC",
        "ic_score": ic,
        "ic_category": category,
        "description": description,
    }


# ---------------------------------------------------------------------------
# Interpretation by Cancer Type
# ---------------------------------------------------------------------------

CANCER_TYPE_THRESHOLDS = {
    "nsclc": {
        "name": "Non-Small Cell Lung Cancer (NSCLC)",
        "score_type": "TPS",
        "thresholds": [
            {"min": 50.0, "label": "TPS >=50%", "therapy": "Pembrolizumab monotherapy (1st line)"},
            {"min": 1.0, "label": "TPS 1-49%", "therapy": "Pembrolizumab + chemotherapy (1st line)"},
            {"min": 0.0, "label": "TPS <1%", "therapy": "PD-L1 negative; consider chemotherapy"},
        ],
    },
    "gastric": {
        "name": "Gastric/GEJ Adenocarcinoma",
        "score_type": "CPS",
        "thresholds": [
            {"min": 10.0, "label": "CPS >=10", "therapy": "Pembrolizumab (3rd line or 1st line if CPS >=10)"},
            {"min": 5.0, "label": "CPS >=5", "therapy": "Pembrolizumab (3rd line)"},
            {"min": 1.0, "label": "CPS 1-4", "therapy": "Limited data; consider clinical trial"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated"},
        ],
    },
    "tnbc": {
        "name": "Triple-Negative Breast Cancer (TNBC)",
        "score_type": "CPS",
        "thresholds": [
            {"min": 10.0, "label": "CPS >=10", "therapy": "Pembrolizumab + chemotherapy (1st line)"},
            {"min": 1.0, "label": "CPS 1-9", "therapy": "Consider pembrolizumab + chemo (1st line)"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated"},
        ],
    },
    "urothelial": {
        "name": "Urothelial Carcinoma",
        "score_type": "CPS",
        "thresholds": [
            {"min": 10.0, "label": "CPS >=10", "therapy": "Pembrolizumab (2nd line or 1st line cisplatin-ineligible)"},
            {"min": 1.0, "label": "CPS 1-9", "therapy": "Pembrolizumab (2nd line if cisplatin-ineligible)"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated for 1st line"},
        ],
    },
    "cervical": {
        "name": "Cervical Cancer",
        "score_type": "CPS",
        "thresholds": [
            {"min": 1.0, "label": "CPS >=1", "therapy": "Pembrolizumab + chemotherapy (1st line)"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated"},
        ],
    },
    "esophageal": {
        "name": "Esophageal Squamous Cell Carcinoma",
        "score_type": "CPS",
        "thresholds": [
            {"min": 10.0, "label": "CPS >=10", "therapy": "Pembrolizumab + chemotherapy (1st line)"},
            {"min": 1.0, "label": "CPS 1-9", "therapy": "Limited data; consider clinical trial"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated"},
        ],
    },
    "hnscc": {
        "name": "Head and Neck Squamous Cell Carcinoma",
        "score_type": "CPS",
        "thresholds": [
            {"min": 20.0, "label": "CPS >=20", "therapy": "Pembrolizumab monotherapy (1st line recurrent/metastatic)"},
            {"min": 1.0, "label": "CPS 1-19", "therapy": "Pembrolizumab + chemotherapy (1st line)"},
            {"min": 0.0, "label": "CPS <1", "therapy": "PD-L1 negative; pembrolizumab not indicated"},
        ],
    },
}


def interpret_for_cancer_type(score_value: float, cancer_type: str) -> Dict[str, Any]:
    """
    Interpret a PD-L1 score for a specific cancer type.

    Parameters
    ----------
    score_value : float
        The TPS or CPS value.
    cancer_type : str
        Cancer type key (e.g., 'nsclc', 'gastric', 'tnbc', 'urothelial').

    Returns
    -------
    dict with cancer_type, score_type, score_value, threshold_met, therapy
    """
    cancer_type = cancer_type.lower().strip()
    if cancer_type not in CANCER_TYPE_THRESHOLDS:
        valid = ", ".join(sorted(CANCER_TYPE_THRESHOLDS.keys()))
        raise ValueError(f"Unknown cancer type '{cancer_type}'. Valid: {valid}")

    info = CANCER_TYPE_THRESHOLDS[cancer_type]

    for threshold in info["thresholds"]:
        if score_value >= threshold["min"]:
            return {
                "cancer_type": info["name"],
                "score_type": info["score_type"],
                "score_value": score_value,
                "threshold_met": threshold["label"],
                "therapy": threshold["therapy"],
            }

    return {
        "cancer_type": info["name"],
        "score_type": info["score_type"],
        "score_value": score_value,
        "threshold_met": "Below all thresholds",
        "therapy": "No PD-L1-directed therapy indicated",
    }


# ---------------------------------------------------------------------------
# Assay Information
# ---------------------------------------------------------------------------

ASSAY_INFO = {
    "22C3": {
        "clone": "22C3 (Dako)",
        "platform": "Dako Autostainer",
        "primary_use": "Pembrolizumab companion diagnostic",
        "score_types": ["TPS", "CPS"],
        "approved_cancers": ["NSCLC", "Gastric/GEJ", "TNBC", "Urothelial", "Cervical", "HNSCC", "Esophageal"],
    },
    "28-8": {
        "clone": "28-8 (Dako)",
        "platform": "Dako Autostainer",
        "primary_use": "Nivolumab complementary diagnostic",
        "score_types": ["TPS"],
        "approved_cancers": ["NSCLC", "Melanoma", "Urothelial"],
    },
    "SP142": {
        "clone": "SP142 (Ventana)",
        "platform": "Ventana BenchMark",
        "primary_use": "Atezolizumab companion diagnostic",
        "score_types": ["TPS", "IC"],
        "approved_cancers": ["TNBC", "NSCLC", "Urothelial"],
    },
    "SP263": {
        "clone": "SP263 (Ventana)",
        "platform": "Ventana BenchMark",
        "primary_use": "Durvalumab companion diagnostic",
        "score_types": ["TPS"],
        "approved_cancers": ["NSCLC", "Urothelial"],
    },
}


def get_assay_info(assay_name: str) -> Dict[str, Any]:
    """Get information about a specific PD-L1 assay."""
    key = assay_name.upper().replace(" ", "").replace("-", "")
    if key not in ASSAY_INFO:
        valid = ", ".join(sorted(ASSAY_INFO.keys()))
        raise ValueError(f"Unknown assay '{assay_name}'. Valid: {valid}")
    return ASSAY_INFO[key]


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def process_batch(input_csv: str, output_csv: str) -> int:
    """
    Process a CSV of PD-L1 cases.

    Expected columns for TPS: pd_l1_positive_tumor_cells, total_viable_tumor_cells
    Optional for CPS: pd_l1_positive_lymphocytes, pd_l1_positive_macrophages
    Optional: cancer_type
    """
    with open(input_csv, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + ["tps", "tps_category", "cps", "cps_category",
                               "therapy_recommendation"]
    out_rows = []

    for row in rows:
        tumor_pos = int(row.get("pd_l1_positive_tumor_cells", 0))
        total = int(row.get("total_viable_tumor_cells", 1))
        lymph = int(row.get("pd_l1_positive_lymphocytes", 0))
        macro = int(row.get("pd_l1_positive_macrophages", 0))
        cancer = row.get("cancer_type", "nsclc")

        tps_result = calculate_tps(tumor_pos, total)
        cps_result = calculate_cps(tumor_pos, lymph, macro, total)

        therapy = ""
        try:
            interp = interpret_for_cancer_type(cps_result["cps"], cancer)
            therapy = interp["therapy"]
        except ValueError:
            therapy = "Unknown cancer type"

        row_dict = dict(row)
        row_dict["tps"] = tps_result["tps"]
        row_dict["tps_category"] = tps_result["tps_category"]
        row_dict["cps"] = cps_result["cps"]
        row_dict["cps_category"] = cps_result["cps_category"]
        row_dict["therapy_recommendation"] = therapy
        out_rows.append(row_dict)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return len(out_rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="PD-L1 TPS & CPS Immunoscore Calculator"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # TPS
    p_tps = subparsers.add_parser("tps", help="Calculate TPS")
    p_tps.add_argument("--positive", type=int, required=True,
                        help="PD-L1 positive tumor cells")
    p_tps.add_argument("--total", type=int, required=True,
                        help="Total viable tumor cells")

    # CPS
    p_cps = subparsers.add_parser("cps", help="Calculate CPS")
    p_cps.add_argument("--tumor-pos", type=int, required=True,
                        help="PD-L1 positive tumor cells")
    p_cps.add_argument("--lymph-pos", type=int, default=0,
                        help="PD-L1 positive lymphocytes")
    p_cps.add_argument("--macro-pos", type=int, default=0,
                        help="PD-L1 positive macrophages")
    p_cps.add_argument("--total-tumor", type=int, required=True,
                        help="Total viable tumor cells")

    # IC
    p_ic = subparsers.add_parser("ic", help="Calculate IC score")
    p_ic.add_argument("--ic-area", type=float, required=True,
                       help="PD-L1 positive immune cell area")
    p_ic.add_argument("--tumor-area", type=float, required=True,
                       help="Total tumor area")

    # Interpret
    p_interp = subparsers.add_parser("interpret", help="Interpret for cancer type")
    p_interp.add_argument("--score", type=float, required=True,
                           help="TPS or CPS value")
    p_interp.add_argument("--cancer-type", required=True,
                           help="Cancer type (nsclc, gastric, tnbc, urothelial, etc.)")

    # Assay info
    p_assay = subparsers.add_parser("assay", help="Get assay information")
    p_assay.add_argument("--name", required=True,
                          help="Assay name (22C3, 28-8, SP142, SP263)")

    # Batch
    p_batch = subparsers.add_parser("batch", help="Batch process CSV")
    p_batch.add_argument("-i", "--input", required=True)
    p_batch.add_argument("-o", "--output", default="results.csv")

    # Audit (enterprise supervisor)
    p_audit = subparsers.add_parser("audit", help="Run supervisor audit on a task")
    p_audit.add_argument("--task-id", required=True, help="Task identifier")
    p_audit.add_argument("--target", default="KEY-001", help="Target identifier")
    p_audit.add_argument("--primary", type=float, default=10.0, help="Primary metric")
    p_audit.add_argument("--secondary", type=float, default=3.0, help="Secondary metric")
    p_audit.add_argument("--descriptor", default="NOMINAL", help="Status descriptor")
    p_audit.add_argument("--critical", action="store_true", help="Critical flag")

    # Chat (enterprise supervisory assistant)
    p_chat = subparsers.add_parser("chat", help="Supervisory chat interface")
    p_chat.add_argument("query", nargs="+", help="Query text")

    # Verify audit trail
    p_verify = subparsers.add_parser("verify-audit", help="Verify HMAC audit trail integrity")

    # Serve (FastAPI REST server)
    p_serve = subparsers.add_parser("serve", help="Start FastAPI REST API server")
    p_serve.add_argument("--host", default="0.0.0.0", help="Host to bind")
    p_serve.add_argument("--port", type=int, default=8000, help="Port to bind")

    args = parser.parse_args(argv)

    if args.command == "tps":
        result = calculate_tps(args.positive, args.total)
        print(json.dumps(result, indent=2))

    elif args.command == "cps":
        result = calculate_cps(args.tumor_pos, args.lymph_pos,
                               args.macro_pos, args.total_tumor)
        print(json.dumps(result, indent=2))

    elif args.command == "ic":
        result = calculate_ic_score(args.ic_area, args.tumor_area)
        print(json.dumps(result, indent=2))

    elif args.command == "interpret":
        result = interpret_for_cancer_type(args.score, args.cancer_type)
        print(json.dumps(result, indent=2))

    elif args.command == "assay":
        result = get_assay_info(args.name)
        print(json.dumps(result, indent=2))

    elif args.command == "batch":
        count = process_batch(args.input, args.output)
        print(f"Processed {count} records -> {args.output}")

    elif args.command == "audit":
        from agents.supervisor import SystemSupervisor
        from agents.models import SystemTaskPayload
        supervisor = SystemSupervisor(model_provider="mock")
        payload = SystemTaskPayload(
            task_id=args.task_id,
            target_identifier=args.target,
            primary_metric=args.primary,
            secondary_metric=args.secondary,
            status_descriptor=args.descriptor,
            is_critical_flag=args.critical,
        )
        dossier = supervisor.process_task(payload)
        print(json.dumps(dossier.to_dict(), indent=2, default=str))

    elif args.command == "chat":
        from agents.supervisor import SystemSupervisor
        supervisor = SystemSupervisor(model_provider="mock")
        query = " ".join(args.query)
        response = supervisor.query_supervisory_chat(query)
        print(json.dumps({"query": query, "response": response}, indent=2))

    elif args.command == "verify-audit":
        from agents.base import AuditLogger
        valid = AuditLogger.verify_integrity()
        trail_len = len(AuditLogger.get_trail())
        print(json.dumps({"audit_valid": valid, "trail_length": trail_len}, indent=2))

    elif args.command == "serve":
        try:
            import uvicorn
            from agents.api import app
            uvicorn.run(app, host=args.host, port=args.port)
        except ImportError:
            print("ERROR: uvicorn and fastapi are required. Install with: pip install fastapi uvicorn")
            sys.exit(1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
