# paper_MRTSampleSizeBinary

Reproducible code for paper "Sample Size Considerations for Micro-Randomized Trials with Binary Proximal Outcomes" by Eric Cohn, Tianchen Qian, Susan A. Murphy.

## File structure

- application code: code to reproduce results in Section 5 "Application".
- functions: code containing relevant user-defined R functions
- misc code: code to make plots to illustrate generative models in Section 6.1 "Generative Models".
- simulation code: code to reproduce simulation results in Section 6 "Detailed Simulation Results" (Sec. 6.2-6.7) and Supplementary Material D.

## How to reproduce the results

- To reproduce results in Section 5 "Application", run all R scripts in folder "application code". The R scripts do not depend on one another so no particular ordering in running.
- To reproduce simulation results in Section 6 "Detailed Simulation Results" (Sec. 6.2-6.7) and Supplementary Material D, do the following for each subfolder inside "simulation code":
    - First, run the R script(s) inside all subfolder(s) "simulationX.X". This conducts Monte Carlo simulations and saves result file. **Caution: each R script may take a long time (days) to finish. Also, X.X in simulationX.X does not correspond to Sections in the paper; these indices are for the authors' internal use.**
    - Second, run the R script(s) named "make figure X.R". This makes plots using the simulation result files. **Caution: X in the "make figure X.R" does not correspond to figure index in the paper. See table below for the figure index correspondence.**
    - For example, to reproduce everything in Section 6.3, go inside folder "2. WA-a violated (Sec 6.3)", then run the R scripts in subfolders "simulation2.2", "simulation2.3", "simulation2.4", "simulation2.5", "simulation4.1". Then, run the three R scripts "make figure 1.R", "make figure 2.R", "make figure 3.R".
- To reproduce Figure 4 in in Sec. 6.1 of the paper, run the three R scripts in folder "misc code".
