Code for "Genome-wide modelling of transcription kinetics reveals patterns of RNA processing delays"
====================================================================================================

This repository contains code used to produce the results in the paper
"Genome-wide modelling of transcription kinetics reveals patterns of
RNA processing delays".

All code is available under the BSD license (see LICENSE.txt).


Gaussian process (GP) modelling code
------------------------------------

The GP code is implemented in Matlab and uses the GPmat library
(https://github.com/SheffieldML/GPmat).

The main script for running the code is 'matlab/run_all_gpnddisim_hmc_convcheck_mixprior.m'
for the real data and 'matlab/run_all_gpnddisim_hmc_convcheck_mixprior_synthetic.m'
for the synthetic data.

The scripts make assumptions about file names and locations so they
would need editing before running.


HMC sample post-analysis code
-----------------------------

There is a number of Matlab scripts for post-processing the HMC samples:
'matlab/check_early_profiles.m'
'matlab/combine_profiles.m'
'matlab/analyse_hmc_new.m'
'matlab/get_hmc_medians_new.m'

The scripts make assumptions about file names and locations so they
would need editing before running.


Model plotting code (cf. Fig. 1)
--------------------------------

The main entry point is 'matlab/run_plot_sampled_predictions.m'.


Analysis of genomic features associated with delays
---------------------------------------------------

The analysis code is written in R. The main script is
'R/analyse_intron_lengths_meddelay.R'. The script run the analysis and
writes most of the figures in the paper.


pre-mRNA distribution analysis
------------------------------

The code is written in Python and is at 'python/plotwigs.py'.


Data pre-processing
-------------------

Pol-II pre-processing code is in 'matlab/pol2_processing'.
RNA-seq pre-processing code is in 'matlab/rnaseq_processing'.
