# nemabiomeR

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ucvm/nemabiomeR.svg?branch=master)](https://travis-ci.org/ucvm/nemabiomeR)
<!-- badges: end -->

Helper functions to run the dada2 pipeline on nemabiome data.  This package provides some convience functions for running an analysis with dada2.  Primarily designed for analyzing nemabiome data but can be used for other dada2 analysis.  Includes a script to run the analysis

## Installation

You can install nemabiomeR from Github with:

``` r
devtools::install_github("ucvm/nemabiomeR")
```

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is required to use `clip_primers` function.  The recommended way to install cutadapt is from [bioconda](https://bioconda.github.io/).  

## Usage

Install the package to get access to the helper functions.  You can then download the `run_pipeline.R` script as a place to start your analysis.



