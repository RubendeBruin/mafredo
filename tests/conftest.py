from pathlib import Path
import pytest


here = Path(__file__).parent

@pytest.fixture
def files():
    """Fixture to provide the path to the files directory."""
    return here / 'files'
