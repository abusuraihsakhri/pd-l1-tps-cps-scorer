"""
Automated Pytest for new CLI commands (audit, chat, verify-audit).
"""
import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from cli import main


class TestAuditCLI:
    def test_audit_basic(self):
        result = main(["audit", "--task-id", "TEST-001"])
        assert result == 0

    def test_audit_with_critical_flag(self):
        result = main(["audit", "--task-id", "TEST-002", "--critical"])
        assert result == 0

    def test_audit_with_custom_metrics(self):
        result = main([
            "audit", "--task-id", "TEST-003",
            "--primary", "35.5",
            "--secondary", "15.0",
            "--descriptor", "DISCORDANT_ANOMALY"
        ])
        assert result == 0


class TestChatCLI:
    def test_chat_basic(self):
        result = main(["chat", "Explain", "PD-L1", "scoring"])
        assert result == 0

    def test_chat_single_word(self):
        result = main(["chat", "status"])
        assert result == 0


class TestVerifyAuditCLI:
    def test_verify_audit(self):
        result = main(["verify-audit"])
        assert result == 0


class TestScoringCLI:
    def test_tps_command(self):
        result = main(["tps", "--positive", "50", "--total", "100"])
        assert result == 0

    def test_cps_command(self):
        result = main([
            "cps", "--tumor-pos", "10",
            "--lymph-pos", "5", "--macro-pos", "3",
            "--total-tumor", "100"
        ])
        assert result == 0

    def test_ic_command(self):
        result = main(["ic", "--ic-area", "5.0", "--tumor-area", "100.0"])
        assert result == 0

    def test_interpret_command(self):
        result = main(["interpret", "--score", "60", "--cancer-type", "nsclc"])
        assert result == 0

    def test_assay_command(self):
        result = main(["assay", "--name", "22C3"])
        assert result == 0
