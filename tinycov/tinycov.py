# Compute and plot coverage in rolling windows along the genome of input bam file.
# cmdoret, 20190920

import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
import click


def parse_bam(bam, chromlist, res):
    """
    Parse input bam file and yield chromosomes one by one along with a
    rolling window mean of coverage.

    Parameters
    ----------
    bam : pysam.Samfile
        A pysam file handle to an alignment file.
    chromlist : list of str
        A list of chromosome names to include in the analysis.

    Returns
    -------
    generator of str, int, pandas.Series
        For each element in the generator, there are 3 values: The chromosome
        name, its length and an array of rolling coverage values.
    """
    for chromo, length in zip(bam.references, bam.lengths):
        if chromo in chromlist:
            depths = np.zeros(length + 1)
            for base in bam.pileup(chromo, setpper="all"):
                try:
                    depths[base.reference_pos] += base.nsegments
                except AttributeError:
                    depths[base.pos] += len(base.pileups)
            df = pd.DataFrame(depths, index=np.arange(length + 1), columns=["depth"])
            yield chromo, length, df.rolling(window=res, center=True).mean()


def aneuploidy_thresh(depths, ploidy=2):
    """Compute coverage thresholds for aneuploidies based on default ploidy."""

    med = np.nanmedian(depths)
    cn_values = np.array(range(1, (2 * ploidy) + 1), dtype=np.float)
    cov_mult = cn_values / ploidy
    ltypes = ["-", "--", "-.", ":"]
    cn_cov = {
        f"{int(cn_values[i])}N": [med * mult, ltypes[i % len(ltypes)]]
        for i, mult in enumerate(cov_mult)
    }
    return cn_cov


@click.command()
@click.option(
    "--res",
    "-r",
    default=10000,
    help="Size of windows in which to compute coverage, in basepairs.",
    show_default=True,
)
@click.option(
    "--skip",
    "-s",
    default=1000,
    help="Stride between windows, in basepairs.",
    show_default=True,
)
@click.option(
    "--name",
    "-n",
    default="",
    help="Name of the sample (plot title). Base name of input file by default",
)
@click.option(
    "--blacklist",
    "-b",
    default="",
    help="Exclude those chromosomes from the plot. List of comma-separated chromosome names.",
)
@click.option(
    "--whitelist",
    "-w",
    default="",
    help="Only include those chromosomes in the plot. List of comma-separated chromosome names.",
)
@click.option(
    "--out",
    "-o",
    default="",
    help="Output file where to write the plot. If not provided, the plot is shown interactively",
    type=click.Path(exists=False),
)
@click.option('-t',
    '--text',
    default='',
    help="Output file where to write the raw data table.",
    type=click.Path(exists=False),
)
@click.option(
    "--ploidy",
    "-p",
    default=2,
    help="Ploidy of input sample, used to estimate coverage threshold for aneuploidies",
)
@click.argument("bam", type=click.Path(exists=True))
def covplot(bam, out, res, skip, name, blacklist, whitelist, ploidy, text):
    click.echo("Visualise read coverage in rolling windows from a bam file.")
    sns.set_style("white")
    # Factor for representing basepairs (1000 for kb)
    scale = 1000
    bam_handle = ps.Samfile(bam)
    blacklist = blacklist.split(",")
    whitelist = whitelist.split(",")
    chromlist = []
    if len(whitelist[0]):
        chromlist = whitelist
    else:
        chromlist = list(bam_handle.references)
        if len(blacklist[0]):
            for chrom in blacklist:
                chromlist.remove(chrom)
    all_depths = []
    if text:
        text_out = open(text, 'w')
    with sns.color_palette("husl", bam_handle.nreferences):
        min_count, max_count = 0, 0
        offset, chrom_id = np.zeros(len(chromlist) + 1), 0
        for chrom, length, counts in parse_bam(bam_handle, chromlist, res):
            coverage = counts[counts.columns[0]].values[::skip] 
            centers = counts.index.values[::skip]
            plt.plot(
                ( centers + offset[chrom_id]) / scale, coverage,
                ".",
            )
            # Write data as text, if requested
            if text:
                for bp, cov in zip(centers, coverage):
                    if not np.isnan(cov):
                        text_out.write(f'{chrom}\t{bp-res//2}\t{bp+res//2}\t{cov}\n')
            highest = np.max(counts.iloc[::skip, 0])
            lowest = np.min(counts.iloc[::skip, 0])
            if lowest < min_count:
                min_count = lowest
            if highest > max_count:
                max_count = highest

            plt.axvline(offset[chrom_id] / scale)
            offset[chrom_id + 1] = offset[chrom_id] + length
            chrom_id += 1
            all_depths.append(counts)
    # Close text file, if opened
    if text:
        text_out.close()
    all_depths = np.concatenate(all_depths, axis=0).flatten()
    for n, chrom in enumerate(chromlist):
        plt.text(
            ((offset[n + 1] - offset[n]) / 2 + offset[n]) / scale,
            1.05 * max_count,
            chrom,
            size=10,
        )
    aneuploidies = aneuploidy_thresh(all_depths, ploidy)
    for aneup, cov in aneuploidies.items():
        if aneup == str(ploidy) + "N":
            lw = 2
            color = "green"
        else:
            lw = 1
            color = "grey"
        plt.axhline(y=cov[0], label=aneup, ls=cov[1], lw=lw, color=color, alpha=0.5)
    plt.xlabel("kb")
    # plt.legend()
    plt.gca().set_ylim([min_count, 1.1 * max_count])
    plt.ylabel(f"coverage ({res}-bp averaged)")
    plt.gca().set_xlim([0, offset[-1] / scale])
    if len(name) == 0:
        name = os.path.splitext(os.path.basename(bam))[0]
    plt.title(name)
    plt.legend()
    if len(out):
        plt.savefig(out)
    else:
        plt.show()

if __name__ == "__main__":
    covplot()
