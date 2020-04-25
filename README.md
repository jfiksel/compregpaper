# Code for A Transformation-free Linear Regression for Compositional Outcomes and Predictors

### How to run
1. Run the following line of code to install the necessary packages: `install.packages(c("gtools", "here", "optimx", "SQUAREM", "boot"))`
2. You can see example code in `tests/example.R`
3. The R files in `simulations` can be used to conduct the simulation study. Note that we used `.sh` scripts to run the code on the JHPCE cluster, and this may need to be modified for your use
4. The two data analyses can be run in `data_analysis/educFMAnalysis.R` and `whiteCellsAnalysis.R`
