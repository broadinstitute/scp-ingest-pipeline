import pytest


@pytest.fixture(autouse=True)
def set_testing(monkeypatch):
    """Sets the TESTING environment variable."""
    monkeypatch.setenv("TESTING", "")
