from typing import List
import click
from tinycov.version import __version__
import tinycov as tc
import functools


def common_options(f):
    """Define options shared between subcommands to avoid code duplication"""
    options = [
        click.option(
            "--res",
            "-r",
            default=10000,
            help="Size of windows in which to compute coverage, in basepairs.",
            show_default=True,
        ),
        click.option(
            "--skip",
            "-s",
            default=1000,
            help="Stride between windows, in basepairs.",
            show_default=True,
        ),
        click.option(
            "--bins",
            "-B",
            default=None,
            help="Tab-separated file of three columns (chromosome, start, end) without header containing a custom binning to use. Overrides --res and --skip, optional.",
            show_default=False,
        ),
        click.option(
            "--name",
            "-n",
            default=None,
            help="Name of the sample (plot title). Base name of input file by default",
        ),
        click.option(
            "--blacklist",
            "-b",
            default=None,
            help="Exclude those chromosomes from the plot. List of comma-separated chromosome names.",
        ),
        click.option(
            "--whitelist",
            "-w",
            default=None,
            help="Only include those chromosomes in the plot. List of comma-separated chromosome names.",
        ),
        click.option(
            "--out",
            "-o",
            default=None,
            help="Output file where to write the plot. If not provided, the plot is shown interactively",
            type=click.Path(exists=False),
        ),
        click.option(
            "--max-depth",
            "-m",
            default=100000,
            help="Maximum read depth permitted. Position with higher coverage will set to this value",
            type=int,
        ),
        click.option(
            "--no-filter",
            "-N",
            help="Use all reads. By default, PCR duplicates and secondary alignments are excluded",
            is_flag=True,
        ),
    ]
    return functools.reduce(lambda x, opt: opt(x), options, f)


def split_commas(arg: str) -> List[str]:
    """Safely split input by commas if it is a string, return None otherwise."""
    try:
        return arg.split(",")
    except AttributeError:
        return None


@click.command()
@common_options
@click.version_option(version=__version__)
@click.argument("bam", type=click.Path(exists=True))
@click.help_option()
def covhist_cmd(
    bam, out, res, bins, skip, name, blacklist, whitelist, max_depth, no_filter
):
    """Visualise the histogram of coverage in rolling windows."""
    click.echo(
        "Visualise read coverage histogram in rolling windows from a bam file."
    )

    tc.covhist(
        bam,
        out=out,
        bins=bins,
        skip=skip,
        name=name,
        blacklist=split_commas(blacklist),
        whitelist=split_commas(whitelist),
        max_depth=max_depth,
        no_filter=no_filter,
    )


@click.command()
@common_options
@click.option(
    "-t",
    "--text",
    default=None,
    help="Output file where to write the raw data table.",
    type=click.Path(exists=False),
)
@click.option(
    "--ploidy",
    "-p",
    default=2,
    help="Ploidy of input sample, used to estimate coverage threshold for aneuploidies. Setting to 0 disables estimations.",
)
@click.version_option(version=__version__)
@click.argument("bam", type=click.Path(exists=True))
def covplot_cmd(
    bam,
    out,
    res,
    bins,
    skip,
    name,
    blacklist,
    whitelist,
    ploidy,
    text,
    max_depth,
    no_filter,
):
    """Visualise coverage in rolling windows, optionally save results to a bedgraph file."""
    click.echo("Visualise read coverage in rolling windows from a bam file.")
    tc.covplot(
        bam,
        out=out,
        res=res,
        bins=bins,
        skip=skip,
        name=name,
        blacklist=split_commas(blacklist),
        whitelist=split_commas(whitelist),
        ploidy=ploidy,
        text=text,
        max_depth=max_depth,
        no_filter=no_filter,
    )


@click.group()
@click.version_option(version=__version__)
def cli():
    """
    tinycov: visualisation of coverage from BAM files using rolling window averages.
    """
    ...


cli.add_command(covplot_cmd, name="covplot")
cli.add_command(covhist_cmd, name="covhist")
