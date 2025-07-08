# Individual Variability
This repository allows to reproduce the results in [Assessing the impact of variance heterogeneity and misspecification in mixed-effects location-scale models](https://arxiv.org/abs/2505.18038). This analysis compares the importance of modeling variance heterogeneity and misspecification in mixed-effects location-scale models. 

## How to use ?
A detailed example is provided in `Tutorial_PBC.qmd`, comparing a traditional mixed-effects model with a model that accounts for variance heterogeneity. The example uses the `PBC` dataset from the `JMbayes2` library.

## Reproduce paper's results
To reproduce the paper's results:

0. Clone the repository: `git clone git@github.com:Jeanselme/IndividualVariability.git`
1. Create a R conda environment with all necessary packages `install.packages(c('brms', 'MASS', 'matrixStats', 'rstan', 'rlist', 'survival', 'dplyr', 'ggplot2', 'mvnfast')) `
3. Run `rscript simple_simulation.R` to run all simulations
5. Analysis the results using `Analysis.ipynb` to measure performance
