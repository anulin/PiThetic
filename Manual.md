# PiThetic Manual

## 1. Overview

PiThetic is a command-line tool for calculating population-genetic statistics from a multisample BAM file or a set of BAM files.

The supported statistics are:

* allele frequency
* nucleotide diversity (π)
* Watterson's theta (θ)
* D'

PiThetic reads the BAM files through `samtools mpileup`.

All BAM files used in the same analysis must be aligned to the same reference sequence.

## 2. Requirements

PiThetic requires:

* Python 3.10
* NumPy
* SciPy
* Numba
* samtools

A Conda environment containing these dependencies is provided in `environment.yml`.

Create the environment with:

```bash
conda env create -f environment.yml
```

Activate it with:

```bash
conda activate PiThetic_env
```

## 3. Running PiThetic

The recommended way to run the program from the repository is:

```bash
python PiThetic.py --[option] [samtools mpileup options] [BAM files]
```

Display the built-in help:

```bash
python PiThetic.py --h
```

Alternatively, if `PiThetic.py` is made executable and placed in a directory in `PATH`, it can be run directly:

```bash
PiThetic.py --h
```

## 4. PiThetic options

### `--freq`

Calculate allele frequency.

This is the default mode, so the following two commands are equivalent:

```bash
python PiThetic.py --freq sample1.bam sample2.bam
```

```bash
python PiThetic.py sample1.bam sample2.bam
```

### `--pi`

Calculate nucleotide diversity (π).

```bash
python PiThetic.py --pi sample1.bam sample2.bam
```

When no window size is provided, the statistic is calculated for each position reported by `samtools mpileup`.

### `--pi [window]`

Calculate nucleotide diversity (π) using windows of the specified size.

For example, for windows of 1000 positions:

```bash
python PiThetic.py --pi 1000 sample1.bam sample2.bam
```

Without `--accurate`, the window calculation uses the default faster calculation.

With `--accurate`, a more accurate calculation is used:

```bash
python PiThetic.py --pi 1000 --accurate sample1.bam sample2.bam
```

The `--accurate` mode can substantially increase computation time, especially for large windows.

### `--theta [window]`

Calculate Watterson's theta (θ) in windows of the specified size.

Example:

```bash
python PiThetic.py --theta 1000 sample1.bam sample2.bam
```

### `--D [window]`

Calculate D' in windows of the specified size.

Example:

```bash
python PiThetic.py --D 1000 sample1.bam sample2.bam
```

### `--accurate`

Improve the accuracy of window-based π and D' calculations.

This option affects window-based calculations only and increases computation time.

Example:

```bash
python PiThetic.py --pi 1000 --accurate sample1.bam sample2.bam
```

### `--t [threads]`

Specify the number of threads.

The default is one thread.

Example:

```bash
python PiThetic.py --pi 1000 --t 20 sample1.bam sample2.bam
```

## 5. samtools mpileup options

Arguments that are not PiThetic-specific are passed to `samtools mpileup`.

This allows `samtools mpileup` options to be used together with PiThetic.

For example:

```bash
python PiThetic.py --pi -l positions.txt sample1.bam sample2.bam
```

The complete list of `samtools mpileup` options is available in the samtools documentation:

https://www.htslib.org/doc/samtools-mpileup.html

## 6. Using position lists and BED files

The `samtools mpileup -l` option can be used to restrict the analysis to selected genomic positions.

### Position list

A position-list file contains at least two columns:

```text
chromosome    position
```

Positions are 1-based.

Example:

```text
chr1    100
chr1    250
chr1    1000
```

Use it with:

```bash
python PiThetic.py --pi -l positions.txt sample1.bam sample2.bam
```

### BED file

A BED file contains at least three columns:

```text
chromosome    start    end
```

BED coordinates are 0-based and half-open.

Example:

```text
chr1    99    100
chr1    249   250
chr1    999   1000
```

Do not mix position-list coordinates and BED coordinates in the same file.

## 7. Example: synonymous and nonsynonymous positions

The `-l` option can be used to select a subset of genomic positions, for example synonymous or nonsynonymous sites.

For example:

```bash
python PiThetic.py --pi -l nonsynonymous.txt sample1.bam sample2.bam
```

and:

```bash
python PiThetic.py --pi -l synonymous.txt sample1.bam sample2.bam
```

The interpretation of the selected positions depends on how the input file was constructed.

## 8. Output

PiThetic writes results to standard output.

For example:

```bash
python PiThetic.py --pi 1000 sample1.bam sample2.bam
```

The output can be redirected to a file:

```bash
python PiThetic.py --pi 1000 sample1.bam sample2.bam > pi.txt
```

For window-based statistics, the output contains the chromosome, the corresponding window interval, and the calculated statistic.

For per-position calculations, the output contains the genomic position and the calculated value.

## 9. Toy dataset

The `toydataset` directory contains a small example dataset consisting of BAM files covering approximately the first 200,000 bp of chromosome 20 from 10 individuals from the 1000 Genomes Project.

The directory also contains example output files.

The commands used to generate the example results are stored in:

```text
toydataset/command.sh
```

The toy dataset can be used to verify that the installation is working and to understand the expected output format.

## 10. Example workflow

A typical workflow is:

### Create the environment

```bash
conda env create -f environment.yml
conda activate PiThetic_env
```

### Check the installation

```bash
python PiThetic.py --h
```

### Calculate π in 1000-position windows

```bash
python PiThetic.py --pi 1000 --t 20 sample1.bam sample2.bam > pi.txt
```

### Calculate Watterson's theta in 1000-position windows

```bash
python PiThetic.py --theta 1000 --t 20 sample1.bam sample2.bam > theta.txt
```

### Calculate D'

```bash
python PiThetic.py --D 1000 --t 20 sample1.bam sample2.bam > D.txt
```

## 11. Notes

PiThetic relies on the output generated by `samtools mpileup`. Therefore, the interpretation and availability of some input options depend on the installed samtools version.

For samtools-specific options and behavior, refer to the samtools documentation:

https://www.htslib.org/doc/samtools-mpileup.html
