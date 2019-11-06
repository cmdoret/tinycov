# Compute and plot coverage in rolling windows along the genome of input bam file.
# cmdoret, 20190920

import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
import click
from tinycov.version import __version__


def check_gen_sort_index(bam, cores=4):
    """
    Index the input BAM file if needed. If the file is not coordinate-sorted,
    generate a sorted file. Returns the path to the sorted-indexed BAM file.
    Does nothing and return input file path if it was already indexed and
    coordinate-sorted.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Handle to the BAM file to index.
    cores : int
        Number of CPU cores to use for sorting.

    Returns
    -------
    sorted_bam : str
        The path to the sorted indexed BAM file
    """

    def check_bam_sorted(path):
        """Checks whether target BAM file is coordinate-sorted"""
        header = ps.view(path, "-H")
        if header.count('SO:coordinate') == 1:
            issorted = True
        else:
            issorted = False
        return issorted

    bam_path = bam.filename.decode()

    try:
        # If the file has an index, there is nothing to do
        bam.check_index()
        sorted_bam = bam_path
    except ValueError:
        # Make a new sorted BAM file and store name for indexing
        if not check_bam_sorted(bam_path):
            sorted_bam = os.path.splitext(bam_path)[0] + ".sorted.bam"
            print("Saving a coordinate-sorted BAM file as ", sorted_bam) 
            ps.sort(bam_path, "-O", "BAM", "-@", str(cores), "-o", sorted_bam)
        else:
            sorted_bam = bam_path
        # Index the sorted BAM file (input file if it was sorted)
        print("Indexing BAM file")
        ps.index("-@", str(cores), sorted_bam)

    return sorted_bam

def parse_bam(bam, chromlist, res):
    """
    Parse input indexed, coordinte-sorted bam file and yield chromosomes 
    one by one along with a rolling window mean of coverage.

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
    """
    Compute coverage thresholds for aneuploidies based on default ploidy.

    Parameters
    ----------
    depth : numpy.array of floats
        1D array of sequencing depths.
    ploidy : int
        The expected or known ploidy of the organism.

    Returns
    -------
    cn_cov : dict
        Map of ploidy to coverage thresholds. Also defines a line type
        for each threshold (for plotting). Has the form {ploidy: [lty, thresh], ...}.
    """

    med = np.nanmedian(depths)
    cn_values = np.array(range(1, (2 * ploidy) + 1), dtype=np.float)
    cov_mult = cn_values / ploidy
    ltypes = ["-", "--", "-.", ":"]
    cn_cov = {
        "{p}N".format(p=int(cn_values[i])): [med * mult, ltypes[i % len(ltypes)]]
        for i, mult in enumerate(cov_mult)
    }
    return cn_cov

def get_bp_scale(size):
    """
    Given a sequence length, compute the appropriate scale and associated suffix (bp, kb, Mb, Gb).

    Parameters
    ----------
    size : int
        The sequence size on which scale and suffix must be computed.

    Returns
    -------
    scale : int
        The number by which basepairs should be divided to achieve desired suffix.
    suffix : str
        The basepair unit associated with the scale.

    Examples
    --------
    >>> get_bp_scale(45000000)
    1000000, "Mb"
    >>> get_bp_scale(12)
    1, "bp"
    """
    # Define valid scales and associated suffixes
    pow_to_suffix = {0: 'bp', 3: 'kb', 6: 'Mb', 9: 'Gb', 12: "Tb"}
    # Compute power scale of genome
    genome_len_pow = int(np.log10(size))
    # Get the closest valid power below genome size
    sorted_pows = sorted(pow_to_suffix.keys())
    valid_pow_idx = max(0, np.searchsorted(sorted_pows, genome_len_pow) - 1)
    genome_valid_pow = sorted_pows[valid_pow_idx]
    # Convert power to scale and associated suffix
    suffix = pow_to_suffix[genome_valid_pow]
    scale = 10**genome_valid_pow
    return scale, suffix


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
    help="Ploidy of input sample, used to estimate coverage threshold for aneuploidies. Setting to 0 disables estimations.",
)
@click.version_option(
    version=__version__,
)
@click.argument("bam", type=click.Path(exists=True))
def covplot_cmd(bam, out, res, skip, name, blacklist, whitelist, ploidy, text):
    click.echo("Visualise read coverage in rolling windows from a bam file.")
    covplot(bam, out=out, res=res, skip=skip, name=name, blacklist=blacklist, whitelist=whitelist, ploidy=ploidy, text=text)

def covplot(bam, out, res=10000, skip=1000, name='', blacklist='', whitelist='', ploidy=2, text=''):
    sns.set_style("white")
    # Load BAM file, sort and index if needed
    bam_handle = ps.Samfile(bam)
    processed_bam = check_gen_sort_index(bam_handle)
    # Reload processed BAM file
    bam_handle = ps.Samfile(processed_bam)
    blacklist = blacklist.split(",")
    whitelist = whitelist.split(",")
    chromlist = []
    genome_len = sum(list(bam_handle.lengths))
    scale, suffix = get_bp_scale(genome_len)
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
                        text_out.write(
                            '{chrom}\t{start}\t{end}\t{cov}\n'.format(
                                chrom=chrom, start=bp-res//2, end=bp+res//2, cov=cov
                            )
                        )
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
    if ploidy > 0:
        aneuploidies = aneuploidy_thresh(all_depths, ploidy)
        for aneup, cov in aneuploidies.items():
            if aneup == str(ploidy) + "N":
                lw = 2
                color = "green"
            else:
                lw = 1
                color = "grey"
            plt.axhline(y=cov[0], label=aneup, ls=cov[1], lw=lw, color=color, alpha=0.5)
        plt.legend()
    plt.xlabel("Genomic position [%s]" % suffix)
    # plt.legend()
    plt.gca().set_ylim([min_count, 1.1 * max_count])
    res_scale, res_suffix = get_bp_scale(res)
    res_str = "%d %s" % (res / res_scale, res_suffix)
    plt.ylabel("coverage (%s averaged)" % res_str)
    plt.gca().set_xlim([0, offset[-1] / scale])
    if len(name) == 0:
        name = os.path.splitext(os.path.basename(bam))[0]
    plt.title(name)
    if len(out):
        plt.savefig(out)
    else:
        plt.show()

if __name__ == "__main__":
    covplot_cmd()
