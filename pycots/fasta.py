from pathlib import Path
from typing import Iterator, Tuple
import itertools
import gzip


def iterate_fasta(fasta: Path) -> Iterator[Tuple[str, str]]:
    """Yield FASTA names and sequences from an optionally compressed file.

    Ignores soft-masking.

    :param fasta: The path to the FASTA file (optionally compressed).
    :type fasta: Path
    :yields: A tuple containing the FASTA header (without the '>') and the corresponding sequence.
    :rtype: Iterator[Tuple[str, str]]
    """
    open_fn = gzip.open if fasta.suffix == ".gz" else open
    with open_fn(fasta, "rt") as f:
        for is_header, group in itertools.groupby(
            f, key=lambda line: line.startswith(">")
        ):
            if is_header:
                contig_name = next(group).strip()[1:]
            else:
                seq = "".join(line.strip() for line in group).upper()
                yield (contig_name, seq)
