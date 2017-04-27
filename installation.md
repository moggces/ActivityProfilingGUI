# Installation requirements

First clone the git repository:

```bash
git clone git@github.com:moggces/ActivityProfilingGUI.git
```

## Prerequisites in R

R must be installed to with the cairo support for PNG rendering, and some sort
of graphical rendering package (tcltk, aqua, X11). R capabilities can be viewed
using the [capabilities](https://stat.ethz.ch/R-manual/R-devel/library/base/html/capabilities.html)
command in R.

Library dependencies are managed using the [packrat](https://rstudio.github.io/packrat/)
bundling system.  To install requirements, first ensure that packrat is installed
globally in R:

```R
install.packages("packrat")
```

If this is a new project that currently doesn't have a `packrat` folder which
tracks the versions of libraries used, you'll need to initiate packrat:

```R
packrat::init()
```

## Download datasets

Datasets are too large to be kept in the github repository. Request datasets
from repository owner (latest version: 2017-03-21). The following files are
required and expected to be in the `./data` path:

- activities_combined_170306.rds
- struct_mat.RData
- tox21_call_descriptions_v2.txt
- tox21_compound_id_v5a7.txt

## Starting the application (R studio)

To develop the application, use the [Rstudio](https://www.rstudio.com/)
development environment. Create a new RStudio project in the root project path.
RStudio should recognize the packrat module and begin installing requirements.

If for some reason modules are not automatically installed, can be used:

```R
packrat::restore()
```

To run the application:

```R
shiny::runApp()
```

## Starting the application (from the terminal)

When using a terminal, packrat does not bootstrap the initialization scripts
and install dependencies. This can be done manually:

```bash
R -e "packrat::restore();"
```

After dependencies have been installed, we can start the application:

```bash
R -e "source('packrat/init.R'); shiny::runApp(port=1234);"
```

## Deploying to Shiny server

Create a bundle to ship to a [shiny server](https://www.rstudio.com/products/shiny/shiny-server/):

```bash
R -e "packrat::bundle(file='~/Desktop/myapp.tar.gz', overwrite=T);"
```

Copy this bundle from your development environment to your web server. In the
example below, we copied the file to `/srv/shiny-bundles/myapp.tar.gz`, and
shiny server is serving applications in the path `/srv/shiny-apps`,
but these could be changed as needed.

To deploy the bundle:

```bash
R -e "install.packages('packrat')"
R -e "packrat::unbundle(bundle='/srv/shiny-bundles/myapp.tar.gz', where='/srv/shiny-apps')"
R -e "setwd('/srv/shiny-apps/myapp'); packrat::on()"
```
