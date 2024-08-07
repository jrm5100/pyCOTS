import itertools

import numpy as np
from numpy.typing import NDArray
from scipy.signal import correlate

NT_MAP = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "R": [1, 0, 1, 0],
    "Y": [0, 1, 0, 1],
    "S": [0, 1, 1, 0],
    "W": [1, 0, 0, 1],
    "K": [0, 0, 1, 1],
    "M": [1, 1, 0, 0],
    "B": [0, 1, 1, 1],
    "D": [1, 0, 1, 1],
    "H": [1, 1, 0, 1],
    "V": [1, 1, 1, 0],
    "N": [1, 1, 1, 1],
}

IDX_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

BASE_TO_IDX = {v: k for k, v in IDX_TO_BASE.items()}


def seq2matrix(seq: str, ignore_n=False) -> NDArray[np.bool_]:
    """Convert an n-length nucleotide sequence to a 4xn bit array

    :param seq: The sequence, optionally using IUPAC ambiguity codes
    :param ignore_n: If True, 'N' is represented by all zeros (no base), otherwise all ones (any base)
    :yields: A binary array with 4 rows and n columns where n is the sequence length
    """
    matrix = np.zeros((4, len(seq)), dtype=bool)
    for idx, nt in enumerate(seq):
        if nt == "N" and ignore_n:
            continue  # Leave N as no base in the reference sequence rather than any base
        matrix[:, idx] = NT_MAP.get(nt.upper(), [0, 0, 0, 0])
    return matrix


def add_variant(seq_matrix: NDArray[np.bool_], position: int, base: str):
    """Modify the sequence matrix in-place using SNV information

    :param seq_matrix: The matrix to be modified
    :param position: The 1-based (VCF style) position
    :param base: The alt base from the variant
    """
    if base not in {"A", "C", "G", "T"}:
        raise ValueError(
            f"Alternate allele must be one of A,C,T, or G. Tried to add '{base}'"
        )
    seq_matrix[BASE_TO_IDX[base], position - 1] = True


def find_matches(
    ref_matrix: NDArray[np.bool_],
    pattern_matrix: NDArray[np.bool_],
    max_mismatches: int = 0,
):
    """Return 0-based positions in the ref matrix where the pattern matches

    :param ref_matrix: The matrix of the sequence to be searched
    :param pattern_matrix: The matrix of the pattern used to search
    :param max_mismatches: Only report positions where the pattern matches with less than this number of mismatches
    :returns: A 1 dimensional numpy array listing 0-based positions where the pattern matches according to the given criteria
    """
    minimum_matched = pattern_matrix.shape[1] - max_mismatches
    assert minimum_matched > 0, ValueError(
        "Too many mismatches allowed- anything would match"
    )
    result = correlate(ref_matrix.astype(int), pattern_matrix.astype(int), mode="valid")
    return np.where(result >= (minimum_matched))[1]


def find_compatible_seqs(
    ref_matrix: NDArray[np.bool_], pattern_matrix: NDArray[np.bool_]
):
    """Return sequences generated from the ref_matrix that are compatible with the pattern matrix

    Positions are 0-based. Both matricies must have the same shape.
    It is possible that no matching sequences exist.

    :param ref_matrix: The matrix searched for a match
    :param pattern_matrix: The matrix used to define a match
    """
    assert ref_matrix.shape == pattern_matrix.shape
    matched_nt_matrix = np.logical_and(ref_matrix, pattern_matrix)
    if np.any(np.sum(matched_nt_matrix, axis=0) == 0):
        return []  # No matched sequences

    valid_bases_indices = [np.nonzero(col)[0] for col in matched_nt_matrix.T]

    return [
        "".join([IDX_TO_BASE[int(idx)] for idx in combination])
        for combination in itertools.product(*valid_bases_indices)
    ]
