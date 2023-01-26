"""Sliding window coverage ingesters."""
# TODO handle circular chromosomes
from abc import ABC, abstractmethod
import csv
from dataclasses import dataclass, field
from functools import cached_property
import os
import sys
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Tuple,
    Type,
)

from pyd4 import D4File
import pysam as ps
import tinycov.utils as tu


@dataclass(order=True)
class Bin:
    """A genomic region."""

    chrom: str
    start: int
    end: int

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def center(self):
        return (self.end + self.start) / 2


@dataclass(order=True)
class Chrom:
    """A chromosome."""

    idx: int
    name: str = field(compare=False)
    size: int = field(compare=False)

    def __len__(self) -> int:
        return self.size


class CovIngester(ABC):
    """Abstract class used to define sliding-window coverage
    ingesters for specific file formats.  Child classes need
    to define the cov and chroms methods.

    Attributes
    ----------
    path:
        Path to the input file.
    res:
        Size to use for sliding windows (in basepairs).
    step:
        Step between windows (in basepairs).
    bins:
        Predefined genome segmentation to use (overrides res and step).
    circular:
        Whether to assume a circular genome for sliding windows.
    """

    def __init__(
        self,
        path: str,
        res: int = 100000,
        step: int = 1000,
        bins: Optional[Iterable[Bin]] = None,
        circular: bool = False,
        pcr_filter: bool = True,
        max_depth: int = 100000,
    ):
        self.path = path
        self.res = res
        self.step = step
        self._bins = bins
        self.circular = circular
        self.pcr_filter = pcr_filter
        self.max_depth = max_depth

    @cached_property
    @abstractmethod
    def chroms(self) -> List[Chrom]:
        """List of chromosomes in the genome."""
        ...

    @abstractmethod
    def cov(self, region: Bin) -> float:
        """Compute mean coverage in target region."""
        ...

    def bins(self, chroms: Optional[Iterable[str]] = None) -> Iterator[Bin]:
        """Yields coordinates of genome segments in the form
        of bins (if provided), or using res and step to generate
        sliding windows.
        Parameters
        ----------
        chroms:
            Names of chromosomes to include.
        """
        if self._bins:
            for bin in self._bins:
                if chroms and bin.chrom in chroms:
                    yield bin
        else:
            for chrom in self.chroms:
                # Only include selected chromosomes
                if chroms and not (chrom.name in chroms):
                    continue
                start = 0
                for end in range(self.res, len(chrom), self.step):
                    yield Bin(chrom.name, start, end)
                    start = end

    def slide(
        self, chroms: Optional[Iterable[str]] = None
    ) -> Iterator[Tuple[Bin, float]]:
        """Iterates on bins (if set) or on sliding windows based on
        res and step, and returns the coverage in each window.

        Parameters
        ----------
        chroms:
            Names of chromosomes to process.
        """
        for bin in self.bins(chroms):
            cov = self.cov(bin)
            yield (bin, cov)

    @property
    def genome_len(self) -> int:
        """Returns the total genome length."""
        return sum([chrom.size for chrom in self.chroms])

    @property
    def chrom_names(self) -> List[str]:
        return [c.name for c in self.chroms]


class D4Ingester(CovIngester):
    """Sliding-window coverage ingester for the
    d4 binary file format."""

    def __post_init__(self):
        self._file = D4File(self.path)

    @cached_property
    def chroms(self) -> List[Chrom]:

        return [
            Chrom(idx, name, size)  # type: ignore
            for idx, (name, size) in enumerate(self._file.chroms())  # type: ignore
        ]

    def cov(self, region: Bin) -> float:
        return self._file[str(region)].mean()  # type: ignore


class BamIngester(CovIngester):
    """Sliding-window coverage ingester for the
    BAM file format."""

    def __post_init__(self):
        # Load BAM file, sort and index if needed
        self._file = ps.AlignmentFile(self.path)
        processed_bam = tu.check_gen_sort_index(self._file)
        self._file = ps.AlignmentFile(processed_bam)

    def cov(self, region: Bin) -> float:
        tot_cov: float = 0.0
        for pos in self._file.pileup(
            region.chrom,
            region.start,
            region.end,
            stepper="all" if self.pcr_filter else "nofilter",
            max_depth=self.max_depth,
        ):
            tot_cov += pos.nsegments
        return tot_cov / len(region)

    @cached_property
    def chroms(self) -> List[Chrom]:
        return [
            Chrom(idx, name, size)
            for idx, (name, size) in enumerate(
                zip(self._file.references, self._file.lengths)
            )
        ]


EXT2ING: Dict[str, Type[CovIngester]] = {
    "bam": BamIngester,
    "sam": BamIngester,
    "d4": D4Ingester,
}


def get_ingester(path: str) -> Type[CovIngester]:
    """Returns the correct ingester class based on input filename."""
    _, ext = os.path.splitext(path)
    try:
        return EXT2ING[ext]
    except IndexError:
        raise ValueError(
            f"Unrecognized file type: {ext}"
            f"Supported types are {list(EXT2ING.keys())}"
        )


def write_bedgraph(
    slider: Iterator[Tuple[Bin, float]],
    path: Optional[str] = None,
):
    if path:
        outf = open(path, "w")
    else:
        outf = sys.stdout
    bedg = csv.writer(outf, delimiter="\t")
    for bin, cov in slider:
        bedg.writerow(bin.chrom, bin.start, bin.end, cov)
