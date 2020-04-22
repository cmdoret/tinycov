# Helper commands for bam files operations used in tinycov
# 20200206, cmdoret
from typing import Optional, Iterable, List, Tuple, Dict, Generator, Union
import pysam as ps
import pathlib
import numpy as np
import pandas as pd
from colorama import Fore
from tqdm import tqdm


def process_chromlist(
    bam: "pysam.AlignmentFile",
    white: Optional[Iterable[str]] = None,
    black: Optional[Iterable[str]] = None,
) -> List[str]:
    """
    Given a bam handle, a chromosome blacklist and whitelist, define which
    chromosomes should be treated.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        The handle to the sam file containing the chromosomes information
    white : list of strs or None
        List of chromosomes to include in analyses. Unless None, all others
        will be dropped.
    black : list of strs or None
        List of chromosomes to exclude from analyses.

    Returns
    -------
    chroms : list of strs
        Chromosomes to be analysed.
    """
    # Require a list of strings
    for l in [white, black]:
        if not (isinstance(l, list) or l is None):
            raise ValueError("white and black must be lists of strings.")

    if white is None:
        chroms = list(bam.references)
        if black is not None:
            for chrom in black:
                chroms.remove(chrom)
    else:
        chroms = white

    return chroms


def check_gen_sort_index(bam: "pysam.AlignmentFile", cores: int = 4) -> str:
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
        if header.count("SO:coordinate") == 1:
            issorted = True
        else:
            issorted = False
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
    pow_to_suffix = {0: "bp", 3: "kb", 6: "Mb", 9: "Gb", 12: "Tb"}
    # Compute power scale of genome
    genome_len_pow = int(np.log10(size))
    # Get the closest valid power below genome size
    sorted_pows = sorted(pow_to_suffix.keys())
    valid_pow_idx = max(0, np.searchsorted(sorted_pows, genome_len_pow) - 1)
    genome_valid_pow = sorted_pows[valid_pow_idx]
    # Convert power to scale and associated suffix
    suffix = pow_to_suffix[genome_valid_pow]
    scale = 10 ** genome_valid_pow
    return scale, suffix


def parse_bam(
    bam: "pysam.AlignmentFile",
    chromlist: Iterable[str],
    res: int,
    bins: Optional["pandas.DataFrame"] = None,
    max_depth: int = 100000,
    no_filter: bool = False,
) -> Generator[Tuple[str, int, "pandas.DataFrame"], None, None]:
    """
    Parse input indexed, coordinte-sorted bam file and yield chromosomes 
    one by one along with a rolling window mean of coverage.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        A pysam file handle to an alignment file.
    chromlist : list of str
        A list of chromosome names to include in the analysis.
    res : int
        Number of bases to include in each window.
    bins : pandas.DataFrame, optional
        Predefined window segmentation in the case of variable-size windows.
        Will override res if used. The dataframe should have columns: "chrom",
        "start" and "end"
    no_filter : bool
        Do not filter out PCR duplicates and secondary alignments (use all reads).
    max_depth : int
        Maximum read depth allowed. Positions above this value will be set to it.

    Returns
    -------
    generator of str, int, pandas.DataFrame
        For each element in the generator, there are 3 values: The chromosome
        name, its length and an array of rolling coverage values.
    """
    stepper = "nofilter" if no_filter else "all"
    for chromo, length in tqdm(
        zip(bam.references, bam.lengths),
        total=len(bam.references),
        desc="chromosome",
        bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET),
    ):
        if chromo in chromlist:
            depths = np.zeros(length + 1)
            for base in tqdm(
                bam.pileup(chromo, stepper=stepper, max_depth=max_depth),
                total=len(depths),
                desc="basepair",
                bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET),
            ):
                try:
                    depths[base.reference_pos] += base.nsegments
                except AttributeError:
                    depths[base.pos] += len(base.pileups)
            df = pd.DataFrame(
                depths, index=np.arange(length + 1), columns=["depth"]
            )
            # If resolution is fixed, use pandas's method for rolling windows
            if bins is None:
                yield chromo, length, df.rolling(
                    window=res, center=True
                ).mean()
            # If using a custom binning, compute each window independently
            else:
                chrom_bins = bins.loc[bins.chrom == chromo, :]
                wins = np.zeros(chrom_bins.shape[0])
                # TODO: Use a faster trick (df.groupby ?)
                # TODO: Ensure windows are centered, test if results are identical to
                # fixed res
                for i, (s, e) in enumerate(
                    zip(chrom_bins.start, chrom_bins.end)
                ):
                    wins[i] = df.depth[s:e].mean()
                yield chromo, length, pd.DataFrame(wins)


def aneuploidy_thresh(
    depths: "numpy.array[float]", ploidy: int = 2
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
    cn_values = np.array(range(1, (2 * ploidy) + 1), dtype=np.float)
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
