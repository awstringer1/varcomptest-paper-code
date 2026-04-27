# Replication code for paper: Testing Linear Combinations of Multiple Variance Components

Paper: ... (to be updated when link available)

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
