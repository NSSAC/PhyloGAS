## The difference between this directory and va_delta

va_delta has one code that does two operations:

1. computes entropy of each column of the 29k elements of genome.
2. given a genomic sequence, a `next' sequence is produced by perturbing one
or more of the 29k slots.

This directory splits these operations into two separate
executions of the code.

---------------

## Directory src.

Contains the source code.

---------------

## Directory test.

There are tests cases.

For kuhlman, the conda venv is py39_andrew_entropy.

### How to run test cases.

To run a case, type:  _./run.XX.Y_

where _XX_ is _01_, _02_, etc., and _Y_ is either _a_ or _b_.

_a_ means compute the entropies based on existing sequences.

_b_ means compute next sequences from an existing sequence, using these entropies.

To compare the resulting output to the valid output, type:

_./run.diff.XX.Y_ where _XX_ and _Y_ are as above.

### Big data file for test cases.

To test the code with actual files, we use the following.

The file _epihyper_exp7_dendogram.csv_ was obtained from:

/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov-2-Biosurveillance-Simulation/data/dendogram

The file of the genomic sequence data (va_variant_BA.2.12.1_sequences.fasta, but also 
va_variant_BA.2.12.1_metadata.csv) was obtained from:

/project/biocomplexity/kuhlman/projects/vdh-andrew-2023/y2023/download-data-curl/real-01/variants-individ-plots

---------------

### Test case descriptions

**run.02 cases**

These run the two analyses with two separate code launches:
- one for computing entropy (_./run.02.a_)
- one for computing next sequences (_./run.02.b_)

In the firt code, a separate threshold DF is written for every
one of the 29k locations in the genomic sequence.
This is 29k files.

**run.03 cases**

These run the two analyses with two separate code launches:
- one for computing entropy (_./run.03.a_)
- one for computing next sequences (_./run.03.b_)

In the firt code, only --one-- threshold DF is written for all
of the 29k locations in the genomic sequence.

This is much better than run.02, in my opinion, because of the
well-known problems for unix to deal with a lot of files in one
directory.

**run.04 case**
This is just like the run,.02 cases, except both analyses
are done with one execution:  _./run.04_


**run.05 case**
This is just like the run,.03 cases, except both analyses
are done with one execution:  _./run.05_


----------

It has been checked that all of runs run.02 through run.05 give the
same results.


