import os
from click.testing import CliRunner
import numpy as np
import pysam as ps
from tinycov import covplot, covhist
from tinycov.utils import (
    parse_bam,
    aneuploidy_thresh,
    get_bp_scale,
    process_chromlist,
)
from tinycov.__main__ import cli

TEST_BAM = "test_data/sorted.bam"
BAD_BAM = "test_data/unsorted.bam"
TEST_BINS = "test_data/custom_bins.bed"
OUT_TXT = "test_data/output.txt"
OUT_IMG = "test_data/output.png"


def test_process_chromlist():
    bam_handle = ps.AlignmentFile(TEST_BAM)
    assert len(process_chromlist(bam_handle)) == 3
    assert process_chromlist(bam_handle, white=["seq1"]) == ["seq1"]
    assert process_chromlist(bam_handle, black=["seq2"]) == ["seq1", "seq3"]


def test_parse_bam():
    """
    Check if parse_bam returns requested chromosomes and correct number of bins.
    """
    bam_handle = ps.AlignmentFile(TEST_BAM)
    res = 100
    # Test using only the first chromosome
    test_chrom = list(bam_handle.references)[0]
    chroms, lengths, counts = [], [], []
    for chrom, length, count in parse_bam(bam_handle, [test_chrom], res):
        chroms.append(chrom)
        lengths.append(length)
        counts.append(count)
    assert lengths[0] == bam_handle.lengths[0] == (counts[0].shape[0] - 1)
    assert len(chroms) == 1


def test_aneuploidy_thresh():
    """
    Test if the aneuploidy threshold returns expected values.
    """
    # Simulate normal distribution of mean=1000; sd=10
    mean = 1000
    depths = np.random.normal(size=10000, loc=mean, scale=10)

    # Test haploid case
    exp_haplo_thresh = np.array([x * mean for x in [1, 2]])
    obs_haplo_dict = aneuploidy_thresh(depths, ploidy=1)
    obs_haplo_thresh = np.sort(
        np.array([i[0] for i in obs_haplo_dict.values()])
    )
    assert np.all(np.isclose(exp_haplo_thresh, obs_haplo_thresh, rtol=1))

    # Test diploid case
    exp_diplo_thresh = np.array([x * mean for x in [0.5, 1, 1.5, 2]])
    obs_diplo_dict = aneuploidy_thresh(depths, ploidy=2)
    obs_diplo_thresh = np.sort(
        np.array([i[0] for i in obs_diplo_dict.values()])
    )
    assert np.all(np.isclose(exp_diplo_thresh, obs_diplo_thresh, rtol=1))

    # Test that output structure is {ploidy: [thresh, lty]}
    assert isinstance(obs_diplo_dict, dict)
    assert len(obs_diplo_dict["1N"]) == 2


def test_covplot():
    """
    Test whether the covplot function exits normally and handles unsorted
    BAM files as well.
    """
    # Remove index and output files if present from previous runs
    for f in [TEST_BAM + ".bai", OUT_TXT, OUT_IMG]:
        try:
            os.remove(f)
        except OSError:
            continue

    # Test sorted/not indexed and unsorted/not indexed cases
    for bam in [TEST_BAM, BAD_BAM]:
        covplot(
            bam,
            out=OUT_IMG,
            res=2000,
            skip=10,
            name="test_run",
            blacklist=None,
            whitelist=None,
            ploidy=2,
            text=OUT_TXT,
        )

    # Check if output files have been properly generated
    assert os.path.isfile(OUT_TXT) == True
    assert os.path.isfile(OUT_IMG) == True

    # Test sorted/indexed case (index generated in previous call)
    # Use a whitelist, custom bins and no ploidy thresholds
    covplot(
        bam,
        out="test_data/output.png",
        res=2000,
        skip=10,
        bins=TEST_BINS,
        name="test_run",
        blacklist=None,
        whitelist=["seq2"],
        ploidy=0,
        text="test_data/output.txt",
    )

    # Test sorted/indexed case (index generated in previous call)
    # Use a blacklist circular chromosomes and no output text or name
    covplot(
        bam,
        out=OUT_IMG,
        res=2000,
        skip=10,
        name=None,
        blacklist=["seq1"],
        whitelist=None,
        ploidy=0,
        text=None,
        circular=True,
    )


def test_covhist():
    """
    Test whether the covhist function exits normally and handles unsorted
    BAM files as well.
    """
    # Remove index and output files if present from previous runs
    for test_f in [TEST_BAM + ".bai", OUT_TXT, OUT_IMG]:
        if os.path.exists(test_f):
            os.remove(test_f)

    # Test sorted/not indexed and unsorted/not indexed cases
    for bam in [TEST_BAM, BAD_BAM]:
        covhist(
            bam,
            out=OUT_IMG,
            res=2000,
            skip=10,
            name="test_run",
            blacklist=None,
            whitelist=None,
        )

    # Check if output files have been properly generated
    assert os.path.isfile(OUT_IMG) == True

    # Test sorted/indexed case (index generated in previous call)
    # Use a whitelist, custom bins and no ploidy thresholds
    covhist(
        bam,
        out="test_data/output.png",
        res=2000,
        skip=10,
        bins=TEST_BINS,
        name="test_run",
        blacklist=None,
        whitelist=["seq2"],
    )

    # Test sorted/indexed case (index generated in previous call)
    # Use a blacklist, circular chromosomes and no output text or name
    covhist(
        bam,
        out=OUT_IMG,
        res=2000,
        skip=10,
        name=None,
        blacklist=["seq1"],
        whitelist=None,
        circular=True,
    )


def test_get_bp_scale():
    obs = [0] * 4
    for i, size in enumerate([17, 180, 14000, 150870320]):
        obs[i] = get_bp_scale(size)
    exp = [(1, "bp"), (1, "bp"), (1000, "kb"), (1000000, "Mb")]
    for o, e in zip(obs, exp):
        assert o == e


def test_cli():
    """Test exit codes of CLI"""
    runner = CliRunner()
    result_plot = runner.invoke(
        cli,
        ["covplot", TEST_BAM, "--out", OUT_IMG, "--res", 2000, "--skip", 10],
    )
    assert result_plot.exit_code == 0
    result_hist = runner.invoke(
        cli,
        ["covhist", TEST_BAM, "--out", OUT_IMG, "--res", 2000, "--skip", 10],
    )
    assert result_hist.exit_code == 0
