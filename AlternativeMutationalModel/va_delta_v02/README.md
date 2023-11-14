## The difference between this directory and va_delta

va_delta has one code that does two operations:

1. computes entropy of each column of the 29k elements of genome.
2. given a genomic sequence, a `next' sequence is produced by perturbing one
or more of the 29k slots.

This directory splits these operations into two separate
executions of the code.
