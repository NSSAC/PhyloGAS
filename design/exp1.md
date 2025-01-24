# EpiHiper Simulation Design of exp1
Adapted based on [multi-variant surveilance exp7 design](https://github.com/NSSAC/2019-nCoV_variant-surveillance-EpiHiper/blob/master/design/exp7.md).

- **Disease model**: 2-variant COVID-19 model with age strtification, waning immunity, and asymmetric cross immunity.
- **Immunity waning**: Any node that has just recovered from infection (regardless of variant) becomes naively susceptible after W days, where W is sampled from an exponential distribution of mean 6 months.
- **Cross infection**: Individuals infected by and recovered from _current_ variant can be infected by _new_ variant again, but not the opposite. The susceptibility of an individual recovered from current variant is regulated by a cross infection parameter $\alpha \in$ {0.25, 0.50, 0.75, 1.00}. Larger $\alpha$ means more susceptible to the new variant.
- **Calibration**: We calibrate transmissibility of current variant against $R_0 = 4.0$. Transmissibility of new variant is 60% higher than that of current variant.
- **Initializations**: Randomly seed current variant and new variant in the whole state population based on state level time series of importations.
- **Vaccination**: We implement state level vaccine administration data from CDC. Vaccine efficacy is 95% protection on fully vaccinated individuals against infection.
- **NPIs**: We do not consider NPIs for now, but may add them if needed.
- **Network**: v2.4.0 networks of GA, MA, MN, VA, WA.
- **Duration**: 365 days.
- **Replicates**: 60.
- **Cells**: Each state has 4 cells for 4 levels of cross-infection parameter.
- Simulation path on Rivanna: `/project/bii_nssac/epihiper-simulations/pipeline-jc/run/20250120_1/output_root/proj/20250120_1/batch_1/`. The output files are `output*.csv.gz`.
