#Estimating undetected Ebola spillovers
Author: Emma Glennon
Contact: eeg31@cam.ac.uk
Last updated: 15 February 2019

This is the code used to generate the analysis for the paper "Estimating undetected Ebola spillovers," in press.

The primary file for this analysis is main.R, written for R version 3.4.4 ("Someone to Lean On"). This contains all steps of analysis, including generating secondary infection distributions, generating outbreak size distributions, generating figures, and performing sensitivity analyses. However, as many of these steps (particularly generating outbreak size distributions) are computationally intensive, .Rda files are provided so the user can pick up analysis from any point. To manually generate these distributions, please see the code chunk labelled "settings" in main.R. If generating new outbreak size distributions, we recommend parallelizing the code (in ./sim-by-name.R) and/or reducing the number of simulated outbreaks (nsim, adjustable in "settings" chunk of main.R), although this will slightly reduce the accuracy of the distribution.

Running main.R without editing the settings loads all secondary infection distributions, outbreak size distributions, and fit observation models and generates all figures (in ./figures/) as well as a file (./results/stats.txt) containing key statistics from the paper. 

This analysis requires the igraph package if generating outbreak distributions from scratch (but not if loading from the files provided); the numDeriv and foreach packages if fitting observation models; and the ggplot2, cowplot, reshape2, and scales packages for creating figures.
