# Installation requirements

First clone the git repository: 

    git clone git@github.com:moggces/ActivityProfilingGUI.git

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/) bundling system.  To install requirements, first ensure that packrat is installed globally in R:

    install.packages("packrat")

## Starting the application (development)

To develop the application, use the [Rstudio](https://www.rstudio.com/) development enviornment. Create a new RStudio project in the root project path. RStudio should recognize the packrat module and begin installing requirements.

If for some reason modules are not automatically installed, can be used:
    
    library(packrat)
    packrat::restore()

To run the application:
    
    library(shiny)
    runApp()

## Starting the application (from the terminal)

When using a terminal, packrat does not bootstrap the initialization scripts and install dependencies. This can be done manually:

    R -e "library(packrat); packrat::restore();"

After dependencies have been installed, we can start the application:

    R -e "source(\"packrat/init.R\"); library(shiny); runApp(port=1234);"
