# SCclust R package

[![DOI](https://zenodo.org/badge/106557810.svg)](https://zenodo.org/badge/latestdoi/106557810)

The SCclust package implements feature selection based on
breakpoints, permutations for FDRs for Fisher test p-values and identification
of the clone structure in single cell copy number profiles.

Threre are two alternative ways to install SCclust package:

* installation from source;

* intallation as a conda package.

## Installation from source

* Download `SCclust` package source build from Github:
[https://github.com/KrasnitzLab/SCclust/archive/refs/tags/1.0.7.tar.gz](https://github.com/KrasnitzLab/SCclust/archive/refs/tags/1.0.7.tar.gz)
into your working directory;

* Run `R` and install the package from downloaded source:

    ```
    install.packages("SCclust-1.0.7.tar.gz")
    ```

## Conda installer

* The `SCclust` package can also be installed using KrasnitzLab Anaconda channel:

    ```
    conda install -c bioconda -c conda-forge -c krasnitzlab scclust
    ```
