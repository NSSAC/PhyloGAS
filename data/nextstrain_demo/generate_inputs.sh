#!/bin/bash

python ./freq.py ./delta.fa --reference ./reference.fasta --output_prefix syn_gen_neutral
python ./freq.py --proportional ./delta.fa --reference ./reference.fasta --proportional --output_prefix syn_gen_proportional
python ./freq.py --poor ./delta.fa --reference ./reference.fasta --poor --output_prefix syn_gen_poor
