from click.testing import CliRunner
import pytest

from pycots.cli import cli


@pytest.fixture()
def runner():
    return CliRunner()


def test_cli_missing_args(runner):
    result = runner.invoke(cli, [])
    assert result.exit_code == 2
    assert "Error: Missing option" in result.output


def test_cli_with_sample_data(shared_datadir, tmp_path_factory):
    runner = CliRunner()
    ref_file = shared_datadir / "chrM.fa"
    vcf_file = shared_datadir / "chrM.vcf.gz"
    output_file = tmp_path_factory.mktemp("output") / "test.txt"
    output_plot = tmp_path_factory.mktemp("output") / "test.png"

    result = runner.invoke(
        cli,
        [
            "--spacer",
            "TAGCTTTTATTCCAGT",
            "--pam",
            "ACN",
            "--ref",
            ref_file,
            "--vcf",
            vcf_file,
            "--output",
            str(output_file),
            "--plot",
            str(output_plot),
        ],
    )
    assert result.exit_code == 0

    # Should find one match, which utilized a variant
    with open(output_file, "r") as f:
        content = f.read()
        assert content == "locus\tpam\nchrM:4582\tACT\n"

    # Check plot
    # TODO
