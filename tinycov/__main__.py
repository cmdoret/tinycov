import click
from tinycov.version import __version__
import tinycov as tc

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
    "--bins",
    "-B",
    default=None,
    help="Tab-separated file of three columns (chromosome, start, end) without header containing a custom binning to use. Overrides --res and --skip, optional.",
    show_default=False,
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
@click.version_option(version=__version__)
@click.argument("bam", type=click.Path(exists=True))
@click.help_option(help="LOLOL")
def covhist_cmd(bam, out, res, bins, skip, name, blacklist, whitelist):
    """Visualise the histogram of coverage in rolling windows."""
    click.echo("Visualise read coverage histogram in rolling windows from a bam file.")
    tc.covhist(
        bam,
        out=out,
        bins=bins,
        skip=skip,
        name=name,
        blacklist=blacklist,
        whitelist=whitelist,
    )


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
    "--bins",
    "-B",
    default=None,
    help="Tab-separated file of three columns (chromosome, start, end) without header containing a custom binning to use. Overrides --res and --skip, optional.",
    show_default=False,
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
@click.option(
    "-t",
    "--text",
    default="",
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
def covplot_cmd(bam, out, res, bins, skip, name, blacklist, whitelist, ploidy, text):
    """Visualise coverage in rolling windows, optionally save results to a bedgraph file."""
    click.echo("Visualise read coverage in rolling windows from a bam file.")
    tc.covplot(
        bam,
        out=out,
        res=res,
        bins=bins,
        skip=skip,
        name=name,
        blacklist=blacklist,
        whitelist=whitelist,
        ploidy=ploidy,
        text=text,
    )

@click.group()
@click.version_option(version=__version__)
def cli():
    """
    tinycov: visualisation of coverage from BAM files using rolling window averages.
    """
    ...


cli.add_command(covplot_cmd, name='covplot')
cli.add_command(covhist_cmd, name='covhist')

