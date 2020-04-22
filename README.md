### tinycov
[![PyPI version](https://badge.fury.io/py/tinycov.svg)](https://badge.fury.io/py/tinycov) [![Build Status](https://travis-ci.com/cmdoret/tinycov.svg?branch=master)](https://travis-ci.com/cmdoret/tinycov) [![codecov](https://codecov.io/gh/cmdoret/tinycov/branch/master/graph/badge.svg)](https://codecov.io/gh/cmdoret/tinycov) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/cmdoret/tinycov.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cmdoret/tinycov/context:python) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Tinycov is a small standalone command line utility written in python to plot the coverage of a BAM file quickly.
This software was inspired by [Matt Edwards' genome coverage plotter](https://github.com/matted/genome_coverage_plotter).


### Installation

To install the stable version:
```pip3 install --user tinycov```

To install the development version:
```
git clone https://github.com/cmdoret/tinycov.git
cd tinycov
pip install .
```

### Input

Only a [BAM](https://www.htslib.org/) file alignment is required as input. If it is not coordinate-sorted, tinycov will make a sorted copy named `input.sorted.bam` if the file is named `input.bam`.

### Output

If no output is provided, the coverage plot will be displayed interactively using matplotlib. If `--out` is used, the plot will be saved in a format determined by the output file's extension.

Additionally, if `--text` is provided, an output text file will be saved in the [bedgraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format with overlapping windows (depending on the values of `--res` and `--skip`).

### Usage

The `tinycov` commands can be invoked from the command line to list subcommands.
```

Usage: tinycov [OPTIONS] COMMAND [ARGS]...

  tinycov: visualisation of coverage from BAM files using rolling window
  averages.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  covhist  Visualise the histogram of coverage in rolling windows.
  covplot  Visualise coverage in rolling windows, optionally save results
           to...
```

The covplot subcommand plots the coverage in sliding windows along the genome. A bedgraph file with the coverage of each window can be generated using the `--text` option.

```

Usage: tinycov covplot [OPTIONS] BAM

  Visualise coverage in rolling windows, optionally save results to a
  bedgraph file.

Options:
  -N, --no-filter          Use all reads. By default, PCR duplicates and
                           secondary alignments are excluded

  -m, --max-depth INTEGER  Maximum read depth permitted. Position with higher
                           coverage will set to this value

  -o, --out PATH           Output file where to write the plot. If not
                           provided, the plot is shown interactively

  -w, --whitelist TEXT     Only include those chromosomes in the plot. List of
                           comma-separated chromosome names.

  -b, --blacklist TEXT     Exclude those chromosomes from the plot. List of
                           comma-separated chromosome names.

  -n, --name TEXT          Name of the sample (plot title). Base name of input
                           file by default

  -B, --bins TEXT          Tab-separated file of three columns (chromosome,
                           start, end) without header containing a custom
                           binning to use. Overrides --res and --skip,
                           optional.

  -s, --skip INTEGER       Stride between windows, in basepairs.  [default:
                           1000]

  -r, --res INTEGER        Size of windows in which to compute coverage, in
                           basepairs.  [default: 10000]

  -t, --text PATH          Output file where to write the raw data table.
  -p, --ploidy INTEGER     Ploidy of input sample, used to estimate coverage
                           threshold for aneuploidies. Setting to 0 disables
                           estimations.

  --version                Show the version and exit.
  --help                   Show this message and exit.
```

The covhist subcommands generates a histogram of coverage values by window. To get a histogram of coverage by basepair, just set `--res` to 1.
```

Usage: tinycov covhist [OPTIONS] BAM

  Visualise the histogram of coverage in rolling windows.

Options:
  -N, --no-filter          Use all reads. By default, PCR duplicates and
                           secondary alignments are excluded

  -m, --max-depth INTEGER  Maximum read depth permitted. Position with higher
                           coverage will set to this value

  -o, --out PATH           Output file where to write the plot. If not
                           provided, the plot is shown interactively

  -w, --whitelist TEXT     Only include those chromosomes in the plot. List of
                           comma-separated chromosome names.

  -b, --blacklist TEXT     Exclude those chromosomes from the plot. List of
                           comma-separated chromosome names.

  -n, --name TEXT          Name of the sample (plot title). Base name of input
                           file by default

  -B, --bins TEXT          Tab-separated file of three columns (chromosome,
                           start, end) without header containing a custom
                           binning to use. Overrides --res and --skip,
                           optional.

  -s, --skip INTEGER       Stride between windows, in basepairs.  [default:
                           1000]

  -r, --res INTEGER        Size of windows in which to compute coverage, in
                           basepairs.  [default: 10000]

  --version                Show the version and exit.
  --help                   Show this message and exit.
```
