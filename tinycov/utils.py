"""Helper commands for bam files operations used in tinycov"""
# 20200206, cmdoret
# pylint: disable=E1101
import pathlib
from typing import (
    Optional,
    Iterable,
    Iterator,
    List,
    Tuple,
    Dict,
    Generator,
    Union,
)

from colorama import Fore
import numpy as np
import pandas as pd
import pysam as ps
from tqdm import tqdm


def filter_chroms(
    chromlist: Iterable[str],
    drop: Optional[Iterable[str]] = None,
    keep: Optional[Iterable[str]] = None,
) -> Iterator[str]:
    """Filter input chromosomes by removing those in drop or keeping only those in keep.

    Parameters
    ----------
    chroms:
        Chromosome names
    drop:
        Names of chromosomes to filter out

    """
    if drop and keep:
        raise ValueError("Either keep or drop must be set, not both.")
    if drop:
        return filter(lambda c: c not in drop, chromlist)
    if keep:
        return filter(lambda c: c in chromlist, keep)
    return (c for c in chromlist)


def check_gen_sort_index(bam: ps.AlignmentFile, cores: int = 4) -> str:
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

    def check_bam_sorted(path: str) -> bool:
        """Checks whether target BAM file is coordinate-sorted"""
        header = ps.view(str(path), "-H")
        issorted = bool(header.count("SO:coordinate") == 1)
        return issorted

    bam_path = pathlib.Path(bam.filename.decode())

    try:
        # If the file has an index, there is nothing to do
        bam.check_index()
        sorted_bam = bam_path
    except ValueError:
        # Make a new sorted BAM file and store name for indexing
        if not check_bam_sorted(bam_path):
            sorted_bam = str(bam_path.with_suffix(".sorted.bam"))
            print("Saving a coordinate-sorted BAM file as ", sorted_bam)
            ps.sort(
                str(bam_path), "-O", "BAM", "-@", str(cores), "-o", sorted_bam
            )
        else:
            sorted_bam = str(bam_path)
        # Index the sorted BAM file (input file if it was sorted)
        print("Indexing BAM file")
        ps.index("-@", str(cores), sorted_bam)

    return str(sorted_bam)


def get_bp_scale(size: int) -> Tuple[int, str]:
    """
    Given a sequence length, compute the appropriate scale and associated suffix
    (bp, kb, Mb, Gb).

    Parameters
    ----------
    size : int
        The sequence size on which scale and suffix must be computed.

    Returns
    -------
    scale : int
        The number by which basepairs should be divided to achieve desired
        suffix.
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
    pow_to_suffix = {0: "bp", 3: "kb", 6: "Mb", 9: "Gb", 12: "Tb"}
    # Compute power scale of genome
    genome_len_pow = int(np.log10(size))
    # Get the closest valid power below genome size
    sorted_pows = sorted(pow_to_suffix.keys())
    valid_pow_idx = max(0, np.searchsorted(sorted_pows, genome_len_pow) - 1)  # type: ignore
    genome_valid_pow = sorted_pows[valid_pow_idx]
    # Convert power to scale and associated suffix
    suffix = pow_to_suffix[genome_valid_pow]
    scale = 10**genome_valid_pow
    return scale, suffix


def aneuploidy_thresh(
    depths: "np.ndarray[float]", ploidy: int = 2
) -> Dict[str, List[Union[float, str]]]:
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
    cn_values = np.array(range(1, (2 * ploidy) + 1), dtype=float)
    cov_mult = cn_values / ploidy
    ltypes = ["-", "--", "-.", ":"]
    cn_cov = {
        "{p}N".format(p=int(cn_values[i])): [
            med * mult,
            ltypes[i % len(ltypes)],
        ]
        for i, mult in enumerate(cov_mult)
    }
    return cn_cov
