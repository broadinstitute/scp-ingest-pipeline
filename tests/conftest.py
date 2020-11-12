import pytest


@pytest.fixture(autouse=True)
def set_testing(monkeypatch):
    """Sets the TESTING environment variable."""
    monkeypatch.setenv("TESTING", "")


@pytest.fixture
def delete_testing(monkeypatch):
    """Deletes the TESTING environment variable."""
    monkeypatch.delenv("TESTING", raising=False)
