#!/usr/bin/env python3
"""Tests for PD-L1 TPS & CPS Scorer - 20 real test cases."""
import csv
import pytest
from pdl1_scorer import (
    calculate_tps, calculate_cps, calculate_ic_score,
    interpret_for_cancer_type, get_assay_info, process_batch
)


# ---------------------------------------------------------------------------
# TPS Tests
# ---------------------------------------------------------------------------

class TestTPS:

    def test_tps_zero(self):
        result = calculate_tps(0, 100)
        assert result["tps"] == 0.0
        assert result["tps_category"] == "Negative"

    def test_tps_below_1_negative(self):
        result = calculate_tps(5, 1000)
        assert result["tps"] == 0.5
        assert result["tps_category"] == "Negative"

    def test_tps_exactly_1(self):
        result = calculate_tps(1, 100)
        assert result["tps"] == 1.0
        assert result["tps_category"] == "Low expression"

    def test_tps_low_expression(self):
        result = calculate_tps(30, 100)
        assert result["tps"] == 30.0
        assert result["tps_category"] == "Low expression"

    def test_tps_exactly_50(self):
        result = calculate_tps(50, 100)
        assert result["tps"] == 50.0
        assert result["tps_category"] == "High expression"

    def test_tps_high_expression(self):
        result = calculate_tps(90, 100)
        assert result["tps"] == 90.0
        assert result["tps_category"] == "High expression"

    def test_tps_invalid_zero_total(self):
        with pytest.raises(ValueError, match="must be > 0"):
            calculate_tps(10, 0)

    def test_tps_positive_exceeds_total(self):
        with pytest.raises(ValueError, match="cannot exceed total"):
            calculate_tps(101, 100)


# ---------------------------------------------------------------------------
# CPS Tests
# ---------------------------------------------------------------------------

class TestCPS:

    def test_cps_basic(self):
        result = calculate_cps(10, 5, 3, 100)
        assert result["cps"] == 18.0
        assert result["total_positive_cells"] == 18

    def test_cps_zero(self):
        result = calculate_cps(0, 0, 0, 100)
        assert result["cps"] == 0.0
        assert result["cps_category"] == "Negative"

    def test_cps_low(self):
        result = calculate_cps(2, 1, 0, 500)
        assert result["cps"] == 0.6
        assert result["cps_category"] == "Negative"

    def test_cps_moderate(self):
        result = calculate_cps(5, 3, 2, 100)
        assert result["cps"] == 10.0
        assert result["cps_category"] == "Moderate expression"

    def test_cps_high(self):
        result = calculate_cps(30, 15, 5, 100)
        assert result["cps"] == 50.0
        assert result["cps_category"] == "High expression"

    def test_cps_invalid_negative_lymph(self):
        with pytest.raises(ValueError, match="must be >= 0"):
            calculate_cps(10, -1, 0, 100)


# ---------------------------------------------------------------------------
# IC Score Tests
# ---------------------------------------------------------------------------

class TestICScore:

    def test_ic_zero(self):
        result = calculate_ic_score(0.0, 100.0)
        assert result["ic_score"] == 0.0
        assert result["ic_category"] == "IC0"

    def test_ic_low(self):
        result = calculate_ic_score(2.0, 100.0)
        assert result["ic_score"] == 2.0
        assert result["ic_category"] == "IC1"

    def test_ic_high(self):
        result = calculate_ic_score(10.0, 100.0)
        assert result["ic_score"] == 10.0
        assert result["ic_category"] == "IC2/3"


# ---------------------------------------------------------------------------
# Cancer Type Interpretation Tests
# ---------------------------------------------------------------------------

class TestCancerTypeInterpretation:

    def test_nsclc_tps_high(self):
        result = interpret_for_cancer_type(60.0, "nsclc")
        assert result["cancer_type"] == "Non-Small Cell Lung Cancer (NSCLC)"
        assert result["threshold_met"] == "TPS >=50%"
        assert "No treatment recommendation" in result["therapy"]

    def test_nsclc_tps_low(self):
        result = interpret_for_cancer_type(20.0, "nsclc")
        assert result["threshold_met"] == "TPS 1-<50%"

    def test_gastric_cps_high(self):
        result = interpret_for_cancer_type(15.0, "gastric")
        assert "No treatment recommendation" in result["therapy"]

    def test_tnbc_cps_threshold(self):
        result = interpret_for_cancer_type(10.0, "tnbc")
        assert "No treatment recommendation" in result["therapy"]

    def test_urothelial_cps_low(self):
        result = interpret_for_cancer_type(5.0, "urothelial")
        assert "No treatment recommendation" in result["therapy"]

    def test_invalid_cancer_type(self):
        with pytest.raises(ValueError, match="Unknown cancer type"):
            interpret_for_cancer_type(10.0, "unknown_cancer")


# ---------------------------------------------------------------------------
# Assay Info Tests
# ---------------------------------------------------------------------------

class TestAssayInfo:

    def test_22c3_assay(self):
        result = get_assay_info("22C3")
        assert result["clone"] == "22C3 (Dako)"
        assert "TPS" in result["score_types"]

    def test_sp142_assay(self):
        result = get_assay_info("SP142")
        assert "IC" in result["score_types"]

    def test_invalid_assay(self):
        with pytest.raises(ValueError, match="Unknown assay"):
            get_assay_info("INVALID")


# ---------------------------------------------------------------------------
# Batch Processing Tests
# ---------------------------------------------------------------------------

class TestBatchProcessing:

    def test_batch_basic(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "pd_l1_positive_tumor_cells,total_viable_tumor_cells,"
            "pd_l1_positive_lymphocytes,pd_l1_positive_macrophages,cancer_type\n"
            "50,100,5,3,nsclc\n"
            "0,100,0,0,gastric\n",
            encoding="utf-8"
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 2
        content = csv_out.read_text(encoding="utf-8")
        assert "tps" in content
        assert "cps" in content


def test_cps_caps_at_100():
    assert calculate_cps(100, 40, 0, 100)["cps"] == 100.0


def test_cps_rejects_impossible_tumor_count():
    with pytest.raises(ValueError, match="cannot exceed"):
        calculate_cps(101, 0, 0, 100)


@pytest.mark.parametrize("fn,args", [
    (calculate_tps, (True, 100)),
    (calculate_tps, (1.5, 100)),
    (calculate_cps, (1, 0, 0, 0)),
    (calculate_ic_score, (float("nan"), 100.0)),
    (calculate_ic_score, (5.0, float("inf"))),
])
def test_invalid_inputs_rejected(fn, args):
    with pytest.raises(ValueError):
        fn(*args)


def test_classification_not_based_on_display_rounding():
    # 0.96% rounds to 1.0%, but remains below the 1% threshold.
    assert calculate_tps(12, 1250)["tps_category"] == "Negative"
    assert calculate_cps(12, 0, 0, 1250)["cps_category"] == "Negative"
    assert calculate_ic_score(0.96, 100)["ic_category"] == "IC0"


def test_assay_name_28_8():
    assert get_assay_info("28-8")["clone"] == "28-8 (Dako)"


def test_batch_uses_tps_for_nsclc(tmp_path):
    input_path = tmp_path / "cases.csv"
    output_path = tmp_path / "results.csv"
    input_path.write_text(
        "pd_l1_positive_tumor_cells,total_viable_tumor_cells,"
        "pd_l1_positive_lymphocytes,pd_l1_positive_macrophages,cancer_type\n"
        "0,100,50,0,nsclc\n", encoding="utf8"
    )
    process_batch(str(input_path), str(output_path))
    with output_path.open(newline="", encoding="utf8") as fh:
        row = next(csv.DictReader(fh))
    assert row["tps"] == "0.0"
    assert row["cps"] == "50.0"
    assert "No treatment recommendation" in row["therapy_recommendation"]


def test_batch_requires_columns(tmp_path):
    input_path = tmp_path / "wrong.csv"
    input_path.write_text("Patient_ID,v1\nP001,3\n", encoding="utf8")
    with pytest.raises(ValueError, match="Missing required CSV columns"):
        process_batch(str(input_path), str(tmp_path / "out.csv"))
