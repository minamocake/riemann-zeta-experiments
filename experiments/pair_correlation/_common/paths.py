"""Path helpers so experiment scripts run from a plain checkout (no installation needed)."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
LEGACY_SRC = ROOT / "legacy" / "pair_correlation_sensitivity_v0" / "src"
CERTIFICATES = ROOT / "certificates"
OUTPUTS = ROOT / "outputs"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
