2014-24-12 millard@insa-toulouse.fr

This archive contains the supplementary data of:

  Impact of kinetic isotope effects in isotope labeling experiments
  by P. Millard, J.C. Portais and P. Mendes

Copyright 2014, INRA, France
License: GNU General Public License v2

For the complete terms of the GNU General Public License, please see the file 'license.txt' distributed with this software or at this URL:
http://www.gnu.org/licenses/gpl-2.0.html


################
# Description  #
################

This archive contains all the scripts required to construct the model represented in figure 2 of the paper, perform the simulations and reproduce the figures 3 to 8.

  - simulations.r      : run this script to construct the models, perform the simulations and generate the figures
  - e_coli_13C_xch.r   : structure of the model shown in figure 2 of the paper, isotope exchange is considered
  - e_coli_13C_noxch.r : structure of the model shown in figure 2 of the paper, isotope exchange is not considered
  - IsoSim.r           : core functions, required to generate the equation system and perform the simulations
  - plot_fun.r         : functions used to plot the simulation results


################
# Requirements #
################

R, Rtools and additional R packages (Matrix, rootSolve, deSolve, stringr, RColorBrewer) are required. This code was tested on Windows 7 and R 3.0, but should run on Linux, MacOS and other platforms supporting R.

R and Rtools can be downloaded online at http://cran.r-project.org/ and must be installed manually (Rtools should also be in the PATH variable of your system, see Rtools documentation for details).
All the required packages can be installed automatically by running the following command in an R console:

    install.packages(c("Matrix", "rootSolve", "deSolve", "stringr", "RColorBrewer"))


################
# Usage        #
################

The figures can be created by running 'simulations.r'. 
All the simulations are performed by IsoSim. Therefore, one should check that IsoSim works correctly before running 'simulations.r'. In that aim:
  
  - open an R console
  - set the working directory to the supplementary data folder:

  setwd("C:/example/sup_data_folder")

  - load IsoSim:

  source("IsoSim.r")

A message will be displayed if some required packages are missing. In this case, please follow the instructions to install the required package(s), then reload IsoSim.

  - if no warning occurs when loading IsoSim, run the test function:

  isosim_test()

This function generates a model of the toy network shown in figure 1 of the paper and performs steady-state and time-course simulations (have a look at this function for examples on network definition and isosim usage). A folder 'test' containing the simulation results should be created in the working directory, and no error should be displayed. If an error is displayed at the compilation step, check that Rtools is correctly installed and is in the PATH variable of your system (update the PATH if needed, see Rtools documentation for help). Please refers to IsoSim code for other problems.
  
To construct the model represented in figure 2 of the paper, perform the simulations and generate the figures 3 to 8, run 'simulations.r' with the command:

  source("simulations.r")

The 'results' directory will be created, with the following files:

  - the figures in pdf format
  - the compiled libraries for each model (with and without taking KIEs into account, with and without isotope exchange)
  - the fortran code from which the librairies are compiled (*.f files)
