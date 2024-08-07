from collections import Counter
from pathlib import Path

import click
import matplotlib.pyplot as plt
import matplotlib.ticker
import pandas as pd


def plot_counts(counts: Counter, plot_file: Path):
    """Plot a histogram to a file

    :param counts: A Counter object holding the number of times each PAM sequence was found
    :param plot_file: Name of a png file to write the plot to
    """
    if len(counts) == 0:
        click.echo("Warning: No PAM matches found, nothing plotted", color="yellow")
        return
    counts = pd.Series(counts).sort_values(ascending=False)
    fig, ax = plt.subplots(1, figsize=(10, 6))
    counts.plot(kind="bar", ax=ax)
    # Styling
    ax.set_xlabel("Sequence", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_title("Number of matches for each unique PAM", fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    plt.tight_layout()
    # Save plot
    fig.savefig(plot_file)
