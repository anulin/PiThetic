# PiThetic

PiThetic is a command-line tool for calculating population-genetic statistics from a multisample BAM file or a set of BAM files.

The following statistics are supported:

* allele frequency
* nucleotide diversity (π)
* Watterson's theta (θ)
* D'

PiThetic uses `samtools mpileup` to read the BAM files. The BAM files must be aligned to the same reference sequence.

## Installation

### Requirements

* Conda
* BAM files suitable for `samtools mpileup`

The repository contains a Conda environment with the required Python packages and `samtools`.

Create the environment:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate PiThetic_env
```

Check that the program starts:

```bash
python PiThetic.py --h
```

## Usage

The general command is:

```bash
python PiThetic.py --[option] [samtools mpileup options] [BAM files]
```

For example:

```bash
python PiThetic.py --freq sample1.bam sample2.bam
```

PiThetic-specific options are listed below.

| Option             | Description                                                                                          |
| ------------------ | ---------------------------------------------------------------------------------------------------- |
| `--h`              | Show help                                                                                            |
| `--freq`           | Calculate allele frequency (default mode)                                                            |
| `--pi`             | Calculate nucleotide diversity (π)                                                                   |
| `--pi [window]`    | Calculate nucleotide diversity (π) in windows of the specified size                                  |
| `--theta [window]` | Calculate Watterson's theta (θ) in windows of the specified size                                     |
| `--D [window]`     | Calculate D' in windows of the specified size                                                        |
| `--accurate`       | Improve the accuracy of window-based π and D' calculations at the cost of increased computation time |
| `--t [threads]`    | Number of threads to use (default: `1`)                                                              |

### Examples

Calculate allele frequency:

```bash
python PiThetic.py --freq sample1.bam sample2.bam
```

Calculate nucleotide diversity:

```bash
python PiThetic.py --pi sample1.bam sample2.bam
```

Calculate nucleotide diversity in windows of 1000 positions:

```bash
python PiThetic.py --pi 1000 sample1.bam sample2.bam
```

Calculate Watterson's theta in windows of 1000 positions:

```bash
python PiThetic.py --theta 1000 sample1.bam sample2.bam
```

Calculate D' in windows of 1000 positions:

```bash
python PiThetic.py --D 1000 sample1.bam sample2.bam
```

Use multiple threads:

```bash
python PiThetic.py --pi 1000 --t 20 sample1.bam sample2.bam
```

The default output is written to standard output (`stdout`). Redirect it to a file when needed:

```bash
python PiThetic.py --pi 1000 --t 20 sample1.bam sample2.bam > pi.txt
```

## Using samtools mpileup options

Options not recognized by PiThetic are passed to `samtools mpileup`.

For example, a BED file or a position-list file can be supplied using the `-l` option:

```bash
python PiThetic.py --pi -l positions.txt sample1.bam sample2.bam
```

See the `samtools mpileup` documentation for the available options:

https://www.htslib.org/doc/samtools-mpileup.html

## Coordinate files

For the `samtools mpileup -l` option, PiThetic can be used with either:

* a position list containing chromosome and position
* a BED file containing chromosome, start, and end

Position-list coordinates are 1-based.

BED coordinates are 0-based and half-open.

Do not mix the two coordinate systems in the same file.

This can be useful, for example, when calculating statistics separately for synonymous and nonsynonymous positions.

## Output

PiThetic writes the calculated statistics to standard output together with the corresponding genomic position or window.

For window-based calculations, the output contains the chromosome, the window interval, and the calculated statistic.

For more details about the calculations and the input formats, see [Manual.md](Manual.md).

## Toy dataset

The `toydataset` directory contains a small example dataset based on approximately the first 200,000 bp of chromosome 20 from 10 individuals from the 1000 Genomes Project.

It includes:

* example BAM files
* example commands
* example output files

The commands used for the toy dataset are in `toydataset/command.sh`.


