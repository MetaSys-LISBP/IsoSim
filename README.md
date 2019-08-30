# IsoSim: stoichiometric, kinetic and isotopic modeling of metabolic systems


## What is IsoSim?
**IsoSim is an R toolbox dedicated to simulation and analysis of metabolic systems under stationary or dynamic conditions.**

IsoSim is a highly versatile platform designed for integration of metabolomics, proteomics and isotopic data with kinetic, thermodynamic, regulatory and stoichiometric constraints.

IsoSim also implements ScalaFlux, a scalable <sup>13</sup>C-fluxomics approach to quantify fluxes in metabolic subnetworks (see [ScalaFlux publication](https://www.biorxiv.org/content/10.1101/735308v1) for details).

The code is open-source, and available under a GPLv3 license. Additional information can be found in [IsoSim](https://dx.doi.org/10.1186%2Fs12918-015-0213-8) and [ScalaFlux](https://www.biorxiv.org/content/10.1101/735308v1) publications.

## Key features

- construction and analysis of stoichiometric and kinetic models of metabolic systems
- simulation of the dynamics of metabolite concentration and/or isotope propagation
- <sup>13</sup>C-flux calculation under instationary and dynamic conditions (other isotopic tracers can also be used)
- <sup>13</sup>C-flux calculation using the ScalaFlux approach (e.g. functions to fit the labeling dynamics of metabolites used as local label inputs, see the ScalaFlux publication for details on this approach)
- estimation of kinetic parameters of enzymes
- non linear sensitivity analysis (Monte Carlo) to determine the precision on estimated parameters
- parallelization to reduce calculation times on multicore platforms


## Installation

R, a set of R development tools and additional R packages (nnls, numDeriv, Matrix, deSolve, stringr, RColorBrewer, pso, parallel) are required. This code was tested on Windows 7 and R 3.2, but should also run on Linux, MacOS and other platforms supporting R.

R can be downloaded online at http://cran.r-project.org/ and must be installed manually. You will also need a set of development tools. On Windows, 
download and install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) (Rtools should also be in the *PATH* variable of your system, see [Rtools documentation](https://cran.r-project.org/bin/windows/Rtools/) for details). On Mac, install the [Xcode command line tools](https://developer.apple.com/downloads). 
On Linux, install the R development package, usually called `r-devel` or `r-base-dev`.

To install all the packages required by IsoSim, run the following command in an R console:

```bash
> install.packages(c("nnls", "numDeriv", "Matrix", "deSolve", "stringr", "RColorBrewer", "pso", "parallel"))
```

Note: to analyse the network provided as [example 1](#1-example-network), you will also need the package vioplot, which can be installed with this command:

```bash
> install.packages("vioplot")
```

To verify that IsoSim works correctly on your system:
  
- go to IsoSim directory, e.g.:

```bash
$ cd /home/usr/isosim/isosim/
```

- open an R console:

```bash
$ R
```

- load IsoSim:

```bash
> source("isosim.R")
```

Note: A message will be displayed if some required packages are missing. In this case, follow the instructions to install the required package(s), then reload IsoSim.

- run the test function:

```bash
> isosim_test()
```

This function builds a test model (network shown in Figure 1A of IsoSim publication), simulates labeling kinetics, and recalculates some fluxes from this theoretical dataset. Most of 
the practical situations and features should be covered by these tests, have a look at this code for examples on network definition and IsoSim usage (we are working on a more complete documentation which should be available soon).
A folder `test` containing all the results should be created in the working directory (here `/home/usr/isosim/isosim/test/`), and no error should be displayed.
If an error is displayed at the compilation step, check that the set of R development tools is correctly installed. 
Please refers to IsoSim code or submit a new issue to our [GitHub issue tracker](https://github.com/MetaSys-LISBP/IsoSim/issues) for other problems.

## ScalaFlux: Scalable <sup>13</sup>C-metabolic flux analysis

We have implemented the [ScalaFlux approach](https://www.biorxiv.org/content/10.1101/735308v1) in IsoSim to quantity fluxes in metabolic subnetworks.

A documentation describing a typical workflow to quantify fluxes in metabolic subnetworks 
using the ScalaFlux approach implemented in IsoSim is provided in [`docs/scalaflux_workflow.md`](docs/scalaflux_workflow.md).

We have also analyzed two metabolic systems as illustrative examples:

- a theoretical network (shown in Figure 1A of [ScalaFlux](https://www.biorxiv.org/content/10.1101/735308v1) publication)

- the yeast prenyl pyrophosphate pathway (shown in Figure 5A of [ScalaFlux](https://www.biorxiv.org/content/10.1101/735308v1) publication)

The code used to perform the calculations detailed in the publication is provided in the `models` folder of this repository. To run the calculations, download 
the current repository and follow the instructions provided below.

Note: folders `isosim` and `models` should be in the same directory.

#### 1. Example network

Note: this script also requires the package vioplot (see [Installation](#installation)).

To run flux calculations on the **example network** and reproduce **Figure 3** (panels B and C) and **Figure 4** (panels B and C):

- go to the example directory, e.g.:

```bash
$ cd /home/usr/isosim/models/example_network/
```

- open an R console:

```bash
$ R
```

- run the example file:

```bash
> source("example_network.R")
```

All output files should be in a `res` folder of the example directory (here `/home/usr/isosim/models/example_network/res/`).

#### 2. Prenyl pyrophosphate biosynthetic pathway

To run flux calculations on the **prenyl pyrophosphate biosynthetic pathway** and reproduce **Figure 5** (panels B-E):

- go to the example directory, e.g.:

```bash
$ cd /home/usr/isosim/models/prenylpyrophosphate_pathway/
```

- open an R console:

```bash
$ R
```

- run the example file:

```bash
> source("example_prenylpyrophosphate.R")
```

All output files should be in a `res` folder of the example directory (here `/home/usr/isosim/models/prenylpyrophosphate_pathway/res/`).

## How to cite

Thank you for using IsoSim and citing us in your work! It means a lot to us and encourages us to continue its development.

Millard P., Schmidt U., Kiefer P., Vorholt J., Heux S., Portais J.C. (2019). ScalaFlux: a scalable approach to quantify fluxes in metabolic subnetworks. BioRxiv, DOI: [10.1101/735308](https://www.biorxiv.org/content/10.1101/735308v1).

## Authors

[Pierre Millard](https://orcid.org/0000-0002-8136-9963), [MetaSys Team](http://www.toulouse-biotechnology-institute.fr/en/research/molecular-physiology-and-metabolism/metasys.html), 
[Toulouse Biotechnology Institute](https://www.lisbp.fr/en/index.html), France

## Contact

:email: Pierre Millard, millard@insa-toulouse.fr
