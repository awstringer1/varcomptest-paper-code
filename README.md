# Replication code for paper: Testing Linear Combinations of Multiple Variance Components

Paper: [arXiv](https://arxiv.org/abs/2604.25744)

Sections 1, 6, and 7 of the main text and Section B of the supplementary materials report empirical results.
This repository contains code and instructions for reproducing those results.

Reproduction of the results requires access to a number of readily available `julia` and `R` packages which
are listed in their respective scripts, as well as the `varcomptest-jl` `julia` package available [here](https://github.com/awstringer1/varcomptest-jl).

In all cases the reported p-values are random and will not be reproduced exactly. 
With high probability your results should be not statistically significantly different from the reported results.
This is more informative than reproducing a single realization of a random result as would be obtained by
setting the random seed. I deliberately do not set the seed, this isn't an error or omission.

For the loading of datasets that are not included in `RDatasets`, I have saved them in `.csv.` format
to the `data` folder in the root of this repository. It is assumed that your `julia` `home` directory
returned by `homedir()` is the root of this repository.

## Section 1

Table 1 in Section 1 reports results of seven data analyses.
You need the following packages:
```
using
  CSV,
  DataFrames,
  MixedModels,
  CategoricalArrays,
  varcomptest,
  RDatasets
```
The scripts to reproduce these analyses are the following:

1. `section1/01-pastes.jl`
2. `section1/02-oats.jl`
3. `section1/03-machines.jl`
4. `section1/04-penicillin.jl`
5. `section1/05-alfalfa.jl`
6. `section1/06-barley.jl`
7. `section1/07-oxide.jl` 


## Section 6 and Supplement B

These sections require the following packages:
```
using 
  varcomptest,
  DataFrames,
  MixedModels,
  Distributions,
  CategoricalArrays,
  JLD2,
  StatsPlots,
  CSV,
  Base.Threads
```

These results were obtained using [a high performance cluster](https://docs.alliancecan.ca/wiki/Trillium).
You may experience difficulty attempting to reproduce all of the simulations on smaller hardware.
However, any individual simulation should run easily on a modern laptop.

### Nested simulations

The reason for the extra naming arguments is to facilitate running these scripts in a slurm array.
The final three positional arguments are used for the slurm job ID, array node ID, and parallel task ID.
When running locally, set the last two to `1` and set the first to any value complicated enough to be parsed uniquely
out of the output file name, i.e. `12345` but not `1`.

- **Run the simulations**: `section6/01-nested.jl`. Call from the command line with the following positional arguments, e.g. like `julia 10 01-nested.jl 1000 20260427 1 300 1 12345 1 1`:
  - Number of threads (`10`)
  - File name
  - Number of simulations to do for each parameter combination (`1000`)
  - Date, used for naming output file (`20260427`)
  - Version number, used for naming output file (`1`)
  - Number of bootstrap samples to do per simulated dataset (`300`)
  - Indicator of sampling type, `1` is balanced, `2` is unbalanced (`1`)
  - ID, used for naming output file (`12345`)
  - ID, used for naming output file (`1`)
  - ID, used for naming output file (`1`) 
- **Summarize the results**: `section6/01-summarize-nested.R`. Run this from the command line with the job ID of the results from running the first file, e.g. `Rscript section6/01-summarize-nested.R 12345`.
  - The script assumes the results are saved in `getwd()` in `R`.
 


### Crossed simulations

- **Run the simulations**: `section6/02-crossed.jl`. Call from the command line with the following positional arguments, e.g. like `julia 10 01-nested.jl 1000 20260427 1 300 6789 1 1`:
  - Number of threads (`10`)
  - File name
  - Number of simulations to do for each parameter combination (`1000`)
  - Date, used for naming output file (`20260427`)
  - Version number, used for naming output file (`1`)
  - Number of bootstrap samples to do per simulated dataset (`300`)
  - ID, used for naming output file (`6789`)
  - ID, used for naming output file (`1`)
  - ID, used for naming output file (`1`) 
- **Summarize the results**: `section6/02-summarize-crossed.R`. Run this from the command line with the job ID of the results from running the first file, e.g. `Rscript section6/02-summarize-crossed.R 6789`.
  - The script assumes the results are saved in `getwd()` in `R`.

## Section 7

Section 7 contains slightly expanded analysis of three of the examples from Table 1 in Section 1.
In `section1` the code for these examples saves the relevant output (bootstrap samples) to the home directory.
To reproduce the histograms in Figures 3, 4, and 5 in the main text, run the following files:

1. `01-summarize-pastes.R`
2. `04-summarize-penicillin.R`
3. `07-summarize-oxide.R`
