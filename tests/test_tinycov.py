import numpy as np
import pysam as ps
from tinycov import parse_bam, aneuploidy_thresh

TEST_BAM = "test_data/for.sorted.bam"

def test_parse_bam():
    """
    Check if parse_bam returns requested chromosomes and correct number of bins.
    """
    bam_handle = ps.Samfile(TEST_BAM)
    res = 100
    # Test using only the first chromosome
    test_chrom = list(bam_handle.references)[0]
    chroms, lengths, counts = [], [], []
    for chrom, length, count in parse_bam(bam_handle, [test_chrom], res):
        chroms.append(chrom)
        lengths.append(length)
        counts.append(count)
    assert(lengths[0] == bam_handle.lengths[0] == (counts[0].shape[0] - 1))
    assert(len(chroms) == 1)


def test_aneuploidy_thresh():
    """
    Test if the aneuploidy threshold returns expected values.
    """
    # Simulate normal distribution of mean=1000; sd=10
    mean = 1000
    depths= np.random.normal(size=10000, loc=mean, scale=10)

    # Test haploid case
    exp_haplo_thresh = np.array([x * mean for x in [1, 2]])
    obs_haplo_dict = aneuploidy_thresh(depths, ploidy=1)
    obs_haplo_thresh = np.sort(np.array([i[0] for i in obs_haplo_dict.values()]))
    assert(np.all(np.isclose(exp_haplo_thresh, obs_haplo_thresh, rtol=1)))

    # Test diploid case
    exp_diplo_thresh = np.array([x * mean for x in [0.5, 1, 1.5, 2]])
    obs_diplo_dict = aneuploidy_thresh(depths, ploidy=2)
    obs_diplo_thresh = np.sort(np.array([i[0] for i in obs_diplo_dict.values()]))
    assert(np.all(np.isclose(exp_diplo_thresh, obs_diplo_thresh, rtol=1)))

    # Test that output structure is {ploidy: [thresh, lty]}
    assert(isinstance(obs_diplo_dict, dict))
    assert(len(obs_diplo_dict['1N']) == 2)

