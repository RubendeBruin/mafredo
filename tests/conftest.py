from __future__ import annotations
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parent.parent  # repo root
TEST_DATA = ROOT / "tests" / "files"


@pytest.fixture(scope="session")
def data_path() -> Path:
    return TEST_DATA
