"""Shared test configuration.

Passing tests mean the software checks passed. They are not mathematical proofs.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True  # keep the frozen legacy tree free of caches
for p in (ROOT / "src", ROOT / "legacy" / "pair_correlation_sensitivity_v0" / "src", ROOT / "tools"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def pytest_terminal_summary(terminalreporter):
    terminalreporter.write_line(
        "NOTE: passing tests show that software checks passed; they are not mathematical proofs. "
        "See docs/status-policy.md.")


@pytest.fixture(scope="session")
def root():
    return ROOT


def has_cvxpy():
    try:
        import cvxpy  # noqa: F401
        return True
    except ImportError:
        return False


requires_cvxpy = pytest.mark.skipif(not has_cvxpy(), reason="cvxpy not installed (optional extra 'sdp')")
