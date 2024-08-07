from cyvcf2 import VCF
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from pycots.matching import seq2matrix, add_variant, find_matches, find_compatible_seqs


@pytest.fixture(scope="function")
def gnomad_test_vcf(shared_datadir):
    return VCF((shared_datadir / "gnomad_test_chr21.vcf.gz"))


class TestSeq2matrix:

    def test_simple(self):
        seq = "ACGT"
        actual = seq2matrix(seq=seq, ignore_n=False)
        expected = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        assert_array_equal(actual, expected)

    def test_softmasked(self):
        seq = "AcgT"
        actual = seq2matrix(seq=seq, ignore_n=False)
        expected = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        assert_array_equal(actual, expected)

    def test_not_ignore_n(self):
        seq = "ANGT"
        actual = seq2matrix(seq=seq, ignore_n=False)
        expected = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]])
        assert_array_equal(actual, expected)

    def test_ignore_n(self):
        seq = "ANGT"
        actual = seq2matrix(seq=seq, ignore_n=True)
        expected = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        assert_array_equal(actual, expected)

    def test_ambiguous(self):
        seq = "RYSWKMBDHV"
        actual = seq2matrix(seq=seq, ignore_n=True)
        expected = np.array(
            [
                [1, 0, 0, 1, 0, 1, 0, 1, 1, 1],
                [0, 1, 1, 0, 0, 1, 1, 0, 1, 1],
                [1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
                [0, 1, 0, 1, 1, 0, 1, 1, 1, 0],
            ]
        )
        assert_array_equal(actual, expected)

    def test_invalid_nt(self):
        seq = "!@#"
        actual = seq2matrix(seq)
        expected = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        assert_array_equal(actual, expected)


class TestAddVariant:

    def test_add_non_base(self):
        matrix = seq2matrix("ACGTA")
        with pytest.raises(
            ValueError,
            match="Alternate allele must be one of A,C,T, or G. Tried to add '.'",
        ):
            add_variant(matrix, 5, ".")

    def test_add_ref(self):
        matrix = seq2matrix("ACGTA")
        expected = matrix.copy()
        add_variant(matrix, 1, "A")
        assert_array_equal(matrix, expected)

    def test_add_pos_1(self):
        matrix = seq2matrix("ACGTA")
        add_variant(matrix, 1, "T")
        expected = np.array(
            [[1, 0, 0, 0, 1], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 0, 1, 0]]
        )
        assert_array_equal(matrix, expected)

    def test_add_pos_last(self):
        matrix = seq2matrix("ACGTA")
        add_variant(matrix, 5, "T")
        expected = np.array(
            [[1, 0, 0, 0, 1], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 1]]
        )
        assert_array_equal(matrix, expected)

    def test_add_multiple(self):
        matrix = seq2matrix("ACGTA")
        add_variant(matrix, 5, "T")
        add_variant(matrix, 5, "C")
        expected = np.array(
            [[1, 0, 0, 0, 1], [0, 1, 0, 0, 1], [0, 0, 1, 0, 0], [0, 0, 0, 1, 1]]
        )
        assert_array_equal(matrix, expected)


class TestFindMatches:

    def test_match_1bp(self):
        ref = seq2matrix("A")
        pattern = seq2matrix("A")
        assert find_matches(ref, pattern) == [0]

    def test_mismatch_1bp(self):
        ref = seq2matrix("A")
        pattern = seq2matrix("C")
        assert_array_equal(find_matches(ref, pattern), [])

    def test_1bp_pattern(self):
        ref = seq2matrix("ACGTACGTA")
        pattern = seq2matrix("A")
        assert_array_equal(find_matches(ref, pattern), [0, 4, 8])

    def test_long_pattern(self):
        ref = seq2matrix("ACGT")
        pattern = seq2matrix("ACGTACGT")
        assert_array_equal(find_matches(ref, pattern), [])

    def test_with_variants(self):
        ref = seq2matrix("ACGTACGTACGT")
        add_variant(ref, 1, "C")
        add_variant(ref, 3, "C")
        add_variant(ref, 7, "C")
        pattern = seq2matrix("CCC")
        assert_array_equal(find_matches(ref, pattern), [0])
        assert_array_equal(find_matches(ref, pattern, max_mismatches=1), [0, 1, 4, 5])
        assert_array_equal(
            find_matches(ref, pattern, max_mismatches=2), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        )

    def test_with_ambiguous_pattern(self):
        ref = seq2matrix("ACGTACGTACGT")
        pattern = seq2matrix("CRD")
        assert_array_equal(find_matches(ref, pattern), [1, 5, 9])


class TestCompatibleSeqs:

    def test_1bp(self):
        ref = seq2matrix("A")
        pattern = seq2matrix("A")
        assert find_compatible_seqs(
            ref,
            pattern,
        ) == ["A"]

    def test_1bp_ambiguous(self):
        ref = seq2matrix("A")
        pattern = seq2matrix("N")
        assert find_compatible_seqs(
            ref,
            pattern,
        ) == ["A"]
        pattern = seq2matrix("M")
        assert find_compatible_seqs(
            ref,
            pattern,
        ) == ["A"]

    def test_1bp_with_variant(self):
        ref = seq2matrix("A")
        add_variant(ref, 1, "T")
        assert find_compatible_seqs(ref, seq2matrix("A")) == ["A"]
        assert find_compatible_seqs(ref, seq2matrix("T")) == ["T"]
        assert find_compatible_seqs(ref, seq2matrix("N")) == ["A", "T"]

    def test_multibp_with_variant(self):
        ref = seq2matrix("CT")
        add_variant(ref, 2, "C")
        assert find_compatible_seqs(ref, seq2matrix("CT")) == ["CT"]
        assert find_compatible_seqs(ref, seq2matrix("CC")) == ["CC"]
        assert find_compatible_seqs(ref, seq2matrix("CN")) == ["CC", "CT"]
        assert find_compatible_seqs(ref, seq2matrix("NN")) == ["CC", "CT"]
