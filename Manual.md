
# Manual
## Requirements
The tool is using [samtools](https://www.htslib.org/) 
## Installation
PiTetic.py is the only file that is needed. Most convenient is to place it in folder that is added to PATH, for example, ~/bin/
## Workflow
The tool can calculate allele frequency, nucleotide diversity and  Watterson theta from a multisample .bam or set of .bam files. 
The .bam files are assumed to be aligned to the same reference.
### Usage:

    PiThetic.py --[flag]  [in1.bam in2.bam ...  or other samtools mpileup input]
      --h - help
      mode flags: 
        --freq - frequency (default)
        --pi - nucleotide diversity
        --pi [window] - calculate nucleotide diversity in window (non-MLE)
        --theta [window] - Watterson theta caclulated within windows of provided size
        --D [window] -calculates D'
        --accurate - improves accuracy for window statistics (only) D' and Pi but increases computation time
        --t - number of threads (default 1)
      for samtools mpileup input you may call samtools mpileup --h 

The chosen statistics are printed in standard output with respect to corresponding samtools mpileup lines.

Samtools mpileup options can be used. They are described in detail on its [page](https://www.htslib.org/doc/samtools-mpileup.html). In order to use them,
they have to be placed between the PiThetic tool options and the .bam file paths. For example it is possible to calculate $\pi_n$ and $\pi_s$. 
For this purpose the list of nonsynonymous/synonymous positions has to be specified in form of BED file or position list via the samtools option **-l**. Position list files contain two columns 
(chromosome and position) and start counting from 1. BED files contain at least 3 columns (chromosome, start and end position) and are 0-based half-open.
While it is possible to mix both position-list and BED coordinates in the same file, this is strongly ill advised due to the differing coordinate systems.

$\pi_n$ usage example:

    PiThetic.py --pi -l nonsynSNP.txt alignment1.bam alignment2.bam alignment3.bam alignment4.bam 
__--accurate__ option increases the computation time for D' and $\pi$ proportionally to the window size specified. Thus a window of 1000 will increase computation time 1000 times

### Toy dataset
in the toydataset folder you may find a toy dataset - .bam files with first ~200000 bp from chr20 of 10 individuals from [1000genome project](www.internationalgenome.org). The example commands are provided in the command.sh file, while their outputs are provided in pi.txt and theta.txt files
