## The difference between this directory and va_delta

va_delta has one code that does two operations:

1. computes entropy of each column of the 29k elements of genome.
2. given a genomic sequence, a `next' sequence is produced by perturbing one
or more of the 29k slots.

This directory splits these operations into two separate
executions of the code.


## Directory test.

There are tests cases.

To run a case, type:  _./run.XX.Y_

where _XX_ is _01_, _02_, etc., and _Y_ is either _a_ or _b_.

To compare the resulting output to the valid output, type:

_./run.diff.XX.Y_ where _XX_ and _Y_ are as above.

The file _epihyper_exp7_dendogram.csv_ was obtained from:

/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov-2-Biosurveillance-Simulation/data/dendogram

The file of the genomic sequence data (va_variant_BA.2.12.1_sequences.fasta, but also 
va_variant_BA.2.12.1_metadata.csv) was obtained from:

/project/biocomplexity/kuhlman/projects/vdh-andrew-2023/y2023/download-data-curl/real-01/variants-individ-plots


