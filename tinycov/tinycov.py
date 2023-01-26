"""Compute and plot coverage in rolling windows along the genome of input bam file."""
# cmdoret, 20190920
# pylint: disable=E1101

import csv
from dataclasses import asdict, dataclass
import os.path
from typing import Iterable, Iterator, Optional, Tuple
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam as ps
from tinycov.io import Bin, CovIngester, get_ingester, write_bedgraph
import tinycov.utils as tu

__all__ = ["covplot", "covhist"]


@dataclass
class CovConfig:
    """Configuration for a CovIngester.

    Attributes
    ----------
    res:
        Size of the sliding window in basepair.
    skip:
        Distance to skip between adjacent windows. If skip == res, windows will
        be non-overlapping.
    bins:
        Define a custom genome segmentation in which to compute coverage, instead
        of regular sliding windows. Overrides res and skip. Set to None to ignore.
    no_filter:
        Do not filter out PCR duplicates and secondary alignments (use all reads).
    max_depth:
        Maximum read depth allowed. Positions above this value will be set to it.
    circular:
        Either the chromosomes are circular or not.
    """

    res: int = 10000
    step: int = 1000
    bins: Optional[Iterable[Bin]] = None
    pcr_filter: bool = True
    max_depth: int = 100000
    circular: bool = False


def covhist(
    cov_file: str,
    config: CovConfig = CovConfig(),
    blacklist: Optional[Iterable[str]] = None,
    whitelist: Optional[Iterable[str]] = None,
):
    ingester = get_ingester(cov_file)(cov_file, **asdict(config))
    chromlist = tu.filter_chroms(
        ingester.chrom_names, drop=blacklist, keep=whitelist
    )


def covtext(
    cov_file: str,
    out_file: Optional[str] = None,
    config: CovConfig = CovConfig(),
    blacklist: Optional[Iterable[str]] = None,
    whitelist: Optional[Iterable[str]] = None,
):
    ingester = get_ingester(cov_file)(cov_file, **asdict(config))
    chromlist = tu.filter_chroms(
        ingester.chrom_names, drop=blacklist, keep=whitelist
    )
    write_bedgraph(ingester.slide(chromlist), out_file)


def covplot(
    cov_file: str,
    out_file: Optional[str] = None,
    config: CovConfig = CovConfig(),
    blacklist: Optional[Iterable[str]] = None,
    whitelist: Optional[Iterable[str]] = None,
    name: Optional[str] = None,
    ploidy: Optional[int] = 2,
):
    """
    Compute read coverage from a BAM file in sliding windows and visualize the
    results.

    Parameters
    ----------
    cov_file : str
        Path to the input BAM or D4 file.
    out : str or None
        Output file for the figure. If none specified, show the figure without
        saving.
    """
    ingester = get_ingester(cov_file)(cov_file, **asdict(config))
    chromlist = tu.filter_chroms(
        ingester.chrom_names, drop=blacklist, keep=whitelist
    )
    if name is None:
        name = os.path.splitext(os.path.basename(cov_file))[0]
    plot_scatter(ingester.slide(chromlist), name, ploidy)
    if out_file is None:
        plt.show()
    else:
        plt.savefig(out_file)
