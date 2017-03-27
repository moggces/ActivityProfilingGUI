# Installation requirements

First clone the git repository:

    git clone git@github.com:moggces/ActivityProfilingGUI.git

## Prerequisites in R

R must be installed to with the cairo support for PNG rendering, and some sort of graphical rendering package (tcltk, aqua, X11). R capabilities can be viewed using the [capabilities](https://stat.ethz.ch/R-manual/R-devel/library/base/html/capabilities.html) command in R.

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/) bundling system.  To install requirements, first ensure that packrat is installed globally in R:

    install.packages("packrat")

## Download datasets

Datasets are too large to be kept in the github repository. Request datasets
from repository owner (latest version: 2016-03-08).

## Starting the application (development)

To develop the application, use the [Rstudio](https://www.rstudio.com/) development enviornment. Create a new RStudio project in the root project path. RStudio should recognize the packrat module and begin installing requirements.

If for some reason modules are not automatically installed, can be used:

    library(packrat)
    packrat::restore()

To run the application:

    library(shiny)
    shiny::runApp()

## Starting the application (from the terminal)

When using a terminal, packrat does not bootstrap the initialization scripts and install dependencies. This can be done manually:

    R -e "library(packrat); packrat::restore();"

After dependencies have been installed, we can start the application:

    R -e "source(\"packrat/init.R\"); library(shiny); shiny::runApp(port=1234);"
