# PiThetic 
The tool can calculate allele frequency, nucleotide diversity and  Watterson theta from a multisample .bam or set of .bam files.

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

The chosen statistics are printed in standard output with respect to corresponding samtools mpileup lines

### Requirements:

[samtools](https://www.htslib.org/)

#### More detailed [manual](https://github.com/anulin/PiThetic/blob/master/Manual.md)
