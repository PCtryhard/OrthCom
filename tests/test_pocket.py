"""Basic tests for OrthCom pocket geometry."""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.pocket_geometry import Pocket


def test_framework_runs():
    """Verify that pytest executes and basic imports succeed."""
    pocket = Pocket()
    assert pocket.volume is None
    assert pocket.residues == []
