# Phosphoproteomic identification of ERK1/2 targets in human neuromesodermal progenitors

This study used phospho-proteomics to identify phospho-dynamic proteins following acute inhibition of ERK1/2 with the MEK1/2 inhibitor PD0325901 in human embryonic stem cell-derived neuromesodermal progenitors.

**Maintainer**: Marek Gierlinski

**Collaborators**: Elisenda Raga Gil, Greg Findlay, Kate Storey

## Usage

The software in this repository is build as a [targets](https://books.ropensci.org/targets/) pipeline in R. We recommend using Positron or RStudio for running the pipeline.

Open an R console in Positron or RStudio and create the environment using `renv`:

```r
install.packages("renv")
renv::restore()
```

This will install all the required packages. If anything is missing, install additional packages locally using `renv::install()`. You can check if all packages are installed by running `source("_targets.R")`. If no errors are produced, you are good to go.

Once the environment is set up, run the `targets` pipeline:

```r
targets::tar_make()
```

This will perform all calculations, generate data objects, figures, and tables (as targets objects), and output TSV files in the `./tab` directory. An HTML report will be generated in the `./doc` directory.
