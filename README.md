# SARS-Cov2-Biosurveillance-Simulation

This repository contains all scripts pertaining to the SARS Cov2 Bisurveillance Simulation Project. Large .fa and zipped .xz and .gz files were omitted from this repository and included in the Large Files release. There are notes in the Large Files Relase that describe what files are necessary for running which .py programs.


## AlternativeMutationalModel
Contains files for the model that propagates mutations based on an entropy threshold. The freq.py file in both the va_delta and va_omicron folders is the python script that generates the mutational model based on the propagation of mutations through a contact network. The output_abridges.py file is the file that defines an abridged contact network from t = 1 to t = 20, which the mutational mode is painted onto.

### Current status of the AlternativeMutationalModel
The most recent iteration of the code exists under AlternativeMutationalModel/va_delta_v02/src, with test files available under AlternativeMutationalModel/va_delta_v02/test. The best test files to look at for the latest functionality are run.03.a and run.03.b. Added features (all for the analysis_type=generate_sequence_analysis stage) include:
- compression: this allows the user to indicate whether to use no compression (the default) or "xz" compression, which significantly decreases the size of the output fasta file from 472M to 5M. NextStrain can read xz files.
- persontrait_file: This allows for the addition of a secondary input file based on the EpiHiper persontrait file that can be used to augment the metadata file (e.g., including the age, gender, age group, race, county, etc. in the metadata output. This argument is optional, but when present should be used with add_metadata argument. A sample fle for this is available at /project/biocomplexity/vdh_genomics/synthetic_biosurveillance/SARS-Cov2-Biosurveillance-Simulation/data/merged_population_files/va_merged_person.csv
- add_metadata: (should be used with persontrait_file parameter) indicates which additional columns should be pulled from the persontrait file for inclusion with in the metadata output.
- proportional: Although this was implemented in a previous version of the code, the introduction of the all_in_one_file option for the threshold_df_num_files parameter in a previous version of the code broke the format and hence the --proportial model. Unlike other parameters, in order to invoke the proportional model, the user needs only type --proportional (no argument) as the command-line argument. If the user does not specify the --proportional or --poor models on the command-line, the default is the "better" or "good" model.

Again, all of these features are highlighted in the test file run.03.b. 

## hmm
Contains files to generate hmm libraries based on defined genetic cutoffs. The hmm libaries could eventually be used to create an hmm-based mutaitonal model.

#### hmm/Global
Contains the files neccessary to build an hmm library based on all COVID sequences globally. Although msaprocessing_global.py should work in theory, the global multiple sequence alignment is too large to run even on a large memory partition. 

#### hmm/Virginia/hmm_library
Contains 27 hmms generated based on all Virginia COVID sequences. This libary was generated from the hmm/virginia/msaprocessing_va.py script.

#### hmm/Virginia_delta/hmm_library
Contains 27 hmms generated based on Virginia delta COVID sequences. The sequences defined as delta sequences are those that exised while the delta variant was in the majority. That is, from 7/13/2021 to 12/21/2021. This libarary was generated from the hmm/Virginia_delta/msaprocessing_delta.py script.

#### hmm/Virginia_omicron/hmm_library
Contains 27 hmms generated based on Virginia omicron COVID sequences. The sequences defined as omicron sequences are those that existed while the omicron variant was in the majority. That is, from 12/21/2021 to 3/7/2022. This library was generated from the hmm/Virginia_omicron/msaprocessing_omicron.py script. 

## network
Contains network.py which generates a visual contact network from t = 0 to 20 (colored by time) using the output_abridges.csv abridged contact network file. 

## trees
Contains newick tree files for all covid sequences globally, in Africa, Asia, Europe, North America, Oceania, and South America. These trees were downloaded from nextstrain.org. These trees can be compared using the robinson-foulds method. 

