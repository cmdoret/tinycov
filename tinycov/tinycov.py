# Compute and plot coverage in rolling windows along the genome of input bam file.
# cmdoret, 20190920

import os.path
from typing import Iterable, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
import tinycov.utils as tu

__all__ = ["covplot", "covhist"]


def covplot(
    bam: str,
    out: Optional[str],
    res: int = 10000,
    bins: Optional[str] = None,
    skip: int = 1000,
    name: Optional[str] = None,
    blacklist: Optional[Iterable[str]] = None,
    whitelist: Optional[Iterable[str]] = None,
    ploidy: Optional[int] = 2,
    text: Optional[str] = None,
    no_filter: bool = False,
    max_depth: int = 100000,
):
    """
    Compute read coverage from a BAM file in sliding windows and visualize the
    results.

    Parameters
    ----------
    bam : str
        Path to the input BAM file.
    out : str or None
        Output file for the figure. If none specified, show the figure without
        saving.
    res : int
        Size of the sliding window in basepair.
    skip : int
        Distance to skip between adjacent windows. If skip == res, windows will
        be non-overlapping.
    bins : str or None
        Define a custom genome segmentation in which to compute coverage, instead
        of regular sliding windows. Overrides res and skip. Set to None to ignore.
    name : str or None
        Display this name as the plot title.
    blacklist : list of strs or None
        Exclude contigs in this list.
    whitelist : list of strs or None
        Only include contigs in this list.
    ploidy : int or None
        Expected ploidy of the organism. Lines will be drawn at expected coverate
        ploidy for different copy number variations. Set to None to ignore ploidy.
    text : str or None
        Path to an optional output text file where the sliding window coverage
        will be written in bedgraph format.
    no_filter : bool
        Do not filter out PCR duplicates and secondary alignments (use all reads).
    max_depth : int
        Maximum read depth allowed. Positions above this value will be set to it.

    """
    sns.set_style("white")
    # Load BAM file, sort and index if needed
    bam_handle = ps.AlignmentFile(bam)
    processed_bam = tu.check_gen_sort_index(bam_handle)
    # Reload processed BAM file
    bam_handle = ps.AlignmentFile(processed_bam)
    chromlist = []
    genome_len = sum(list(bam_handle.lengths))
    scale, suffix = tu.get_bp_scale(genome_len)
    chromlist = tu.process_chromlist(
        bam_handle, black=blacklist, white=whitelist
    )
    all_depths = []
    if bins is not None:
        # Cannot skip windows if using custom binning
        skip = 1
        bins = pd.read_csv(
            bins, sep="\t", header=None, names=["chrom", "start", "end"]
        )
    if text:
        text_out = open(text, "w")
    with sns.color_palette("husl", bam_handle.nreferences):
        min_count, max_count = 0, 0
        offset, chrom_id = np.zeros(len(chromlist) + 1), 0
        for chrom, length, counts in tu.parse_bam(
            bam_handle,
            chromlist,
            res,
            bins,
            max_depth=max_depth,
            no_filter=no_filter,
        ):
            coverage = counts[counts.columns[0]].values[::skip]
            if bins is None:
                centers = counts.index.values[::skip]
            else:
                chrom_bins = bins.loc[bins.chrom == chrom, :]
                chrom_bins["depth"] = coverage
                centers = (chrom_bins.start + chrom_bins.end) / 2
            plt.scatter((centers + offset[chrom_id]) / scale, coverage, marker=".", edgecolors='none')
            # Write data as text, if requested
            if text:
                if bins is None:
                    for bp, cov in zip(centers, coverage):
                        if not np.isnan(cov):
                            text_out.write(
                                "{chrom}\t{start}\t{end}\t{cov}\n".format(
                                    chrom=chrom,
                                    start=bp - res // 2,
                                    end=bp + res // 2,
                                    cov=cov,
                                )
                            )
                else:
                    chrom_bins.to_csv(
                        text_out, mode="a", sep="\t", header=None, index=False
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
        aneuploidies = tu.aneuploidy_thresh(all_depths, ploidy)
        for aneup, cov in aneuploidies.items():
            if aneup == str(ploidy) + "N":
                lw = 2
                color = "green"
            else:
                lw = 1
                color = "grey"
            plt.axhline(
                y=cov[0], label=aneup, ls=cov[1], lw=lw, color=color, alpha=0.5
            )
        plt.legend()
    plt.xlabel("Genomic position [%s]" % suffix)
    # plt.legend()
    plt.gca().set_ylim([min_count, 1.1 * max_count])
    if bins is None:
        res_scale, res_suffix = tu.get_bp_scale(res)
        res_str = "%d %s" % (res / res_scale, res_suffix)
    else:
        res_str = "variable size windows"
    plt.ylabel("coverage (%s averaged)" % res_str)
    plt.gca().set_xlim([0, offset[-1] / scale])
    if name is None:
        name = os.path.splitext(os.path.basename(bam))[0]
    plt.title(name)
    if out is None:
        plt.show()
    else:
        plt.savefig(out)


def covhist(
    bam: str,
    out: Optional[str],
    res: int = 10000,
    bins: Optional[str] = None,
    skip: int = 1000,
    name: Optional[str] = None,
    blacklist: Optional[Iterable[str]] = None,
    whitelist: Optional[Iterable[str]] = None,
    no_filter: bool = False,
    max_depth: int = 100000,
):
    """
    Compute read coverage per sliding window from a BAM file and generate
    histograms of window coverage.

    Parameters
    ----------
    bam : str
        Path to the input BAM file.
    out : str or None
        Output file for the figure. If none specified, show the figure without
        saving.
    res : int
        Size of the sliding window in basepair.
    skip : int
        Distance to skip between adjacent windows. If skip == res, windows will
        be non-overlapping.
    bins : str or None
        Path to a BED file defining a custom genome segmentation in which to
        compute coverage, instead of regular sliding windows. Overrides res
        and skip. Set to None to ignore.
    name : str or None
        Display this name as the plot title.
    blacklist : list of strs or None
        Exclude contigs in this list.
    whitelist : list of strs or None
        Only include contigs in this list.
    no_filter : bool
        Do not filter out PCR duplicates and secondary alignments (use all reads).
    max_depth : int
        Maximum read depth allowed. Positions above this value will be set to it.
    """
    sns.set_style("white")
    # Load BAM file, sort and index if needed
    bam_handle = ps.AlignmentFile(bam)
    processed_bam = tu.check_gen_sort_index(bam_handle)
    # Reload processed BAM file
    bam_handle = ps.AlignmentFile(processed_bam)
    chromlist = []
    chromlist = tu.process_chromlist(
        bam_handle, black=blacklist, white=whitelist
    )
    all_depths = []
    all_chroms = []
    if bins is not None:
        # Cannot skip windows if using custom binning
        skip = 1
        bins = pd.read_csv(
            bins, sep="\t", header=None, names=["chrom", "start", "end"]
        )
    for chrom, length, counts in tu.parse_bam(
        bam_handle,
        chromlist,
        res,
        bins,
        max_depth=max_depth,
        no_filter=no_filter,
    ):
        coverage = counts[counts.columns[0]].values[::skip]
        all_depths.append(coverage)
        all_chroms += [chrom] * len(coverage)

    if bins is None:
        res_scale, res_suffix = tu.get_bp_scale(res)
        res_str = "%d %s" % (res / res_scale, res_suffix)
    else:
        res_str = "variable size windows"
    all_depths = np.concatenate(all_depths, axis=0).flatten()
    all_chroms = np.array(all_chroms)[~np.isnan(all_depths)]
    all_depths = all_depths[~np.isnan(all_depths)]
    hist_depths = pd.DataFrame({"chrom": all_chroms, "depth": all_depths})
    # Set x scale based on min and max coverage of whole genome
    upper_bound = np.percentile(all_depths, 95)
    lower_bound = np.percentile(all_depths, 5)
    hist_bins = np.linspace(lower_bound, upper_bound, 100)
    g = sns.FacetGrid(hist_depths, col="chrom", col_wrap=3)
    g = (
        g.map(plt.hist, "depth", color="c", bins=hist_bins)
        .set_titles("{col_name}")
        .set_axis_labels(
            "coverage (%s averaged)" % res_str, "Number of windows"
        )
    )
    if name is None:
        name = os.path.splitext(os.path.basename(bam))[0]
    g.fig.suptitle(name)
    if out is None:
        plt.show()
    else:
        plt.savefig(out)
