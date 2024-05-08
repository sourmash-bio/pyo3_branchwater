import os
import pytest

from .sourmash_tst_utils import TempDirectory, RunnerContext

@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)

# Set environment variable PYTEST_RUNNING
def pytest_configure(config):
    os.environ["PYTEST_RUNNING"] = "1"
