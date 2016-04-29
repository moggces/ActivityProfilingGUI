# Installation requirements

First clone the git repository: 

    git clone git@github.com:moggces/ActivityProfilingGUI.git

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/) bundling system. This requires use of the [Rstudio](https://www.rstudio.com/) IDE. To install requirements, first ensure that packrat is installed globally in R:

    install.packages("packrat")

## Starting the application (development)

If starting a development environment, create a new RStudio project in the root project path. RStudio should recognize the packrat module and begin installing requirements. Alternatively, the following commands could be given:
    
    library(packrat)
    packrat::restore()

To run the application:
    
    library(shiny)
    runApp()

## Starting the application (from the terminal)

When using a terminal, packrat does not bootstrap the initialization scripts and install dependencies. This can be done manually (again, from the root-path of the project):

    R -e "library(packrat); packrat::restore();"

After dependencies have been installed, we can star the application:

    R -e "source(\"packrat/init.R\"); library(shiny); runApp(port=1234);"
