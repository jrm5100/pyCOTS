from pathlib import Path
from collections import Counter

import click
from cyvcf2 import VCF

from pycots.fasta import iterate_fasta
from pycots.matching import add_variant, find_matches, find_compatible_seqs, seq2matrix
from pycots.plot import plot_counts


@click.command()
@click.option("--spacer", "-s", required=True, help="Spacer sequence")
@click.option("--pam", "-p", required=True, help="Pam sequence")
@click.option(
    "--ref",
    "-r",
    required=True,
    type=click.Path(exists=True),
    help="Optionally compressed reference FASTA file",
)
@click.option(
    "--vcf",
    "-v",
    required=True,
    type=click.Path(exists=True),
    help="Compressed and indexed VCF file with SNVs to consider when matching",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(exists=False, writable=True),
    help="Output file location for the file listing matched locations and PAM sequences",
)
@click.option(
    "--plot",
    "-l",
    required=True,
    type=click.Path(exists=False, writable=True),
    help="Histogram of PAM sequence matches",
)
def cli(spacer, pam, ref, vcf, output, plot) -> None:
    assert 16 <= len(spacer) <= 25, ValueError(f"Invalid Spacer Length {len(spacer)}")
    assert 2 <= len(pam) <= 6, ValueError(f"Invalid PAM Length {len(pam)}")
    spacer_matrix = seq2matrix(spacer)
    pam_matrix = seq2matrix(pam)
    vcf = VCF(vcf)
    with Path(output).open("w") as o:
        o.write("locus\tpam\n")
        counter = Counter()
        for seq_name, seq in iterate_fasta(Path(ref)):
            ref_matrix = seq2matrix(seq, ignore_n=True)
            # Add variants to the reference sequence
            for var in vcf(seq_name):
                if len(var.ALT) == 1:
                    add_variant(ref_matrix, var.POS, var.ALT[0])
            match_starts = find_matches(ref_matrix, spacer_matrix, max_mismatches=1)
            for m in match_starts:
                for pam_seq in find_compatible_seqs(
                    ref_matrix[:, m + len(spacer) : m + len(spacer) + len(pam)],
                    pam_matrix,
                ):
                    o.write(f"{seq_name}:{m+1}\t{pam_seq}\n")
                    counter[pam_seq] += 1
        if len(counter) == 0:
            click.echo("No Spacer-PAM matches found", color="green")
        else:
            click.echo("Spacer-PAM match(es) found.", color="yellow")
            plot_counts(counter, plot)


if __name__ == "__main__":
    cli()
