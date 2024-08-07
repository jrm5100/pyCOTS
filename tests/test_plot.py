from collections import Counter
from pathlib import Path

import pytest

from pycots.plot import plot_counts


@pytest.fixture(scope="function")
def plot_file(tmp_path_factory):
    plot_file = tmp_path_factory.mktemp("output") / "test.png"
    return Path(plot_file)


def test_plot_none(plot_file):
    counts = Counter()
    plot_counts(counts, plot_file)


def test_plot_one(plot_file):
    counts = Counter({"ACGT": 1})
    plot_counts(counts, plot_file)


def test_plot_multiple(plot_file):
    counts = Counter({"ACGT": 13, "ACGG": 19})
    plot_counts(counts, plot_file)
