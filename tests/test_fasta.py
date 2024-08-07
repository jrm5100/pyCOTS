from pathlib import Path

import pytest

from pycots.fasta import iterate_fasta


@pytest.fixture(scope="function")
def single_line_fasta(tmp_path_factory):
    content = """>seq1
ACGTACGT
>seq2
TGACTGAC
"""
    fasta = tmp_path_factory.mktemp("data") / "test.fasta"
    with open(fasta, "w") as o:
        o.write(content)
    return Path(fasta)


@pytest.fixture(scope="function")
def multiline_fasta(tmp_path_factory):
    content = """>seq1
ACGTACGT
CGTAC
>seq2
ACGTacgt
ACGT
>seq3
ACGT
"""
    fasta = tmp_path_factory.mktemp("data") / "test.fasta"
    with open(fasta, "w") as o:
        o.write(content)
    return Path(fasta)


def test_single_line(single_line_fasta):
    actual = list(iterate_fasta(single_line_fasta))
    expected = [("seq1", "ACGTACGT"), ("seq2", "TGACTGAC")]
    assert actual == expected


def test_multi_line(multiline_fasta):
    actual = list(iterate_fasta(multiline_fasta))
    expected = [("seq1", "ACGTACGTCGTAC"), ("seq2", "ACGTACGTACGT"), ("seq3", "ACGT")]
    assert actual == expected
