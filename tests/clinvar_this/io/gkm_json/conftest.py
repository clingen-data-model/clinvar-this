"""Shared paths and input formats for GKM tests."""

import pathlib

import pytest

DATA_DIR = pathlib.Path(__file__).parents[2] / "data" / "io_gkm_json"

NATIVE_DIR = DATA_DIR / "native"
BUNDLE_DIR = DATA_DIR / "bundles"


@pytest.fixture(params=["native", "bundle"])
def gkm_input_format(request: pytest.FixtureRequest) -> str:
    """Run shared transformation expectations against each input format."""
    return request.param
