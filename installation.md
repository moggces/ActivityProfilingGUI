# Installation requirements

First clone the git repository:

```bash
git clone git@github.com:moggces/ActivityProfilingGUI.git
```

## Prerequisites in R

R must be installed to with the cairo support for PNG rendering, and some sort of graphical rendering package (tcltk, aqua, X11). R capabilities can be viewed using the [capabilities](https://stat.ethz.ch/R-manual/R-devel/library/base/html/capabilities.html) command in R.

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/) bundling system.  To install requirements, first ensure that packrat is installed globally in R:

```R
install.packages("packrat")
```

## Download datasets

Datasets are too large to be kept in the github repository. Request datasets from repository owner (latest version: 2017-03-21). The following files are required and expected to be in the `./data` path:

- activities_combined_170306.rds
- struct_mat.RData
- tox21_call_descriptions_v2.txt
- tox21_compound_id_v5a7.txt

## Starting the application (development)

To develop the application, use the [Rstudio](https://www.rstudio.com/) development environment. Create a new RStudio project in the root project path. RStudio should recognize the packrat module and begin installing requirements.

If for some reason modules are not automatically installed, can be used:

```R
library(packrat)
packrat::restore()
```

To run the application:

```R
library(shiny)
shiny::runApp()
```

## Starting the application (from the terminal)

When using a terminal, packrat does not bootstrap the initialization scripts and install dependencies. This can be done manually:

```bash
R -e "library(packrat); packrat::restore();"
```

After dependencies have been installed, we can start the application:

```bash
R -e "source(\"packrat/init.R\"); library(shiny); shiny::runApp(port=1234);"
```

### Development details

To contribute code, use [formatR](https://cran.r-project.org/web/packages/formatR/)
to automatically take a first pass at code-formatting.

```R
formatR::tidy_source('./server.R',indent=4, width.cutoff=78, recursive=F)
formatR::tidy_dir('./source',indent=4, width.cutoff=78)
```

Additional modification of lines is required, however, after using this tool,
it should only be used when creating a new file.
