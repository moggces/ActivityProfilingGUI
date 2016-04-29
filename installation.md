# Installation requirements

First clone the git repository: 

    git clone git@github.com:moggces/ActivityProfilingGUI.git

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/) bundling system. This requires use of the [Rstudio](https://www.rstudio.com/) IDE. To install requirements, first ensure that packrat is installed globally in R:

    install.packages("packrat")

Next, create a new Rstudio project for this application. Dependencies should install automatically for this application.

## Running the application

To run the application from the command-line, change directories to the root path of the applcation, then enter the following command in the terminal:

    R -e "library(shiny); runApp(port=1234)"
