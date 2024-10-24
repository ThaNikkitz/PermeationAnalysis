"""
Unit and regression test for the PermeationAnalysis package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import PermeationAnalysis


def test_PermeationAnalysis_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "PermeationAnalysis" in sys.modules
