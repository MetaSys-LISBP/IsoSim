# ScalaFlux workflow

This documentation describes a typical workflow to quantify fluxes in metabolic subnetworks 
using the ScalaFlux approach.

IsoSim must be installed beforehand. Have a look [here](https://github.com/MetaSys-LISBP/IsoSim#installation) for installation instructions.

## Load IsoSim

- Go to IsoSim directory, e.g.:

  ```bash
  $ cd /home/usr/isosim/isosim/
  ```

- Open an R console:

  ```bash
  $ R
  ```
  
- Load IsoSim:

  ```bash
  > source("isosim.R")
  ```
  
A message will appear if some of the required packages are missing. In this case, follow the instructions to [install](https://github.com/MetaSys-LISBP/IsoSim#installation) the missing package(s), then reload IsoSim.
  
## Create results directory

By default, results will be saved in your current working directory. You can also gather all results in a new folder (here `res_wf`) of a specific directory with :

```bash
> wd <- "/path/to/working/directory/"
> res_folder <- "res_wf"
> setwd(wd)
> if (!file.exists(res_folder)) {dir.create(file.path(wd, res_folder))}
> setwd(file.path(wd, res_folder))
```

## Parallelization options

The overall computation time can be reduced by executing some calculations in parallel, i.e. simultaneously.

Parallelization options can be adapted with the variable `numCores` (int), which represents the maximal number of CPU cores 
which can be used in parallel by IsoSim.

To use only one core (i.e. no parallelization):

```bash
> numCores <- NULL
```

To use two core:

```bash
> numCores <- 2
```

The maximum number of processes that can run at a time is limited by the number of CPU cores on the current host. 
If want to use all of them, use the `detectCores()` function: 

```bash
> cat("IsoSim will use", detectCores(), "CPU cores.\n", sep=" ")
> numCores <- detectCores()
```

## Sensitivity analysis options

Number of Monte Carlo iterations for flux calculation:

```bash
> mc_iter <- 10
```

## Construct a metabolic model

- Define carefully the topology of the metabolic subsystem to model, 
using the following structure:

  ```bash
  network (list):
    $R (vector): reactant(s)
    $C (vector): stoichiometric coefficient(s)
    $E (vector): rate law (can be a constant, a variable or an expression)
    $T (vector): tracer atom transitions
  ```
  
  For instance, the model of the example network provided in figure 1A of ScalaFlux publication is:
    
  ```bash
  > rxn <- list(r1     = list("R"=c("Sout", "Sin"), "C"=c(-1, 1),     "E"="v1",   "T"=c("A", "A")),
                r2     = list("R"=c("Sin", "A"),    "C"=c(-1, 1),     "E"="v2",   "T"=c("A", "A")),
                r3     = list("R"=c("Sin", "D"),    "C"=c(-1, 1),     "E"="v3",   "T"=c("A", "A")),
                r4     = list("R"=c("A", "F"),      "C"=c(-1, 1),     "E"="v4",   "T"=c("A", "A")),
                r5     = list("R"=c("A", "B"),      "C"=c(-1, 1),     "E"="v5",   "T"=c("A", "A")),
                r6for  = list("R"=c("B", "C"),      "C"=c(-1, 1),     "E"="v6f",  "T"=c("A", "A")),
                r6rev  = list("R"=c("C", "B"),      "C"=c(-1, 1),     "E"="v6r",  "T"=c("A", "A")),
                r7     = list("R"=c("D", "C"),      "C"=c(-1, 1),     "E"="v7",   "T"=c("A", "A")),
                r8     = list("R"=c("C", "E"),      "C"=c(-1, 1),     "E"="v8",   "T"=c("A", "A")),
                r9     = list("R"=c("E", "F"),      "C"=c(-1, 1),     "E"="v9",   "T"=c("A", "A")),
                r10    = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                r11    = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v11",  "T"=c("A", "B", "BA")),
                r12    = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v12",  "T"=c("AB", "B", "A")),
                r13    = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v13",  "T"=c("B", "B")),
                r14    = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v14",  "T"=c("A", "A")),
                r15    = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v15",  "T"=c("A", "A")),
                r16    = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v16",  "T"=c("A", "A")),
                r17    = list("R"=c("N", "O"),      "C"=c(-1, 1),     "E"="v17",  "T"=c("A", "A")),
                r18    = list("R"=c("F", "O"),      "C"=c(-1, 1),     "E"="v18",  "T"=c("A", "A")),
                r19    = list("R"=c("O", "P"),      "C"=c(-1, 1),     "E"="v19",  "T"=c("A", "A")),
                r20    = list("R"=c("P"),           "C"=c(-1),        "E"="v20",  "T"=c("A"))
  )
  ```
  
- Define additional analytical equations (optional), e.g. to calculate determined fluxes from 
free fluxes or to estimate other variables used during calculations.

  ```bash
  > eq_det <- c("v3 = v1-v2",
                "v5 = v6f-v6r",
                "v4 = v2-v5",
                "v7 = v3",
                "v8 = v6f-v6r+v7",
                "v9 = v8-v10",
                "v11 = v10",
                "v12 = v11",
                "v13 = v12",
                "v14 = v12",
                "v15 = v14",
                "v16 = v15",
                "v17 = v16",
                "v18 = v4+v9",
                "v19 = v18+v17",
                "v20 = v19")
  ```
  
- Construct the model (isotopic and stoichiometric matrices, identification of sources, sinks, label inputs, etc) with:

  ```bash
  > net <- net2mat(rxn, add_eq=eq_det)
  ```
  
## Simulate labeling dynamics

Once you have constructed the model, we strongly encourage you to run simulations and to verify that the simulation results are consistent with the 
expected behaviour (e.g. assuming the investigated system operates at metabolic steady-state, metabolite concentrations and fluxes should be constant).

- Define the labeling dynamics of label input(s) (here S<sub>out</sub>) using analytical functions:

  For instance, to simulate a switch from unlabeled to fully-labeled nutrient:
  
  ```bash
  > anFun <- list("Sout_1-M0" = "0.0+0.0*t",
                  "Sout_1-M1" = "1.0+0.0*t")
  ```
  
  Each element of this list is a (time-dependent) analytical function of a given isotopologue of a label input EMU, with the appropriate name `X_Y-Mn` (`X_Y-Z-Mn`), 
  where `X` (str) is the metabolite name, `Y` (int) (and all `Z`s) is (are) the position(s) of atom(s) of tracer element present in the corresponding EMU, 
  and `n` (int) is the weight of the corresponding isotopologue. 
  The list of all label input EMUs isotopologues required to perform the simulations is automatically identified by IsoSim during model construction and can be found 
  in `net$min_meas`.
  
  Note: All the analytical functions *must* contain the time variable `t`, even if the corresponding label input(s) is (are) constant.
  			  
- Initialize fluxes and metabolite concentrations:

  ```bash
  > fluxes <- c("v1"=2,
                "v2"=1.5,
                "v6"=1.2,
                "v10"=1.0,
                "v6xch"=0.2)
  ```
  
  and
  
  ```bash
  > meta_conc <- c("Sout"=1,
                   "Sin"=0.5,
                   "A"=0.5,
                   "B"=0.5,
                   "C"=0.5,
                   "D"=0.5,
                   "E"=0.5,
                   "F"=0.5,
                   "G"=0.5,
                   "H"=0.5,
                   "I"=0.5,
                   "J"=0.5,
                   "K"=0.5,
                   "L"=0.5,
                   "M"=0.5,
                   "N"=0.5,
                   "O"=0.5,
                   "P"=20)
  ```
  
  Note: Concentration of all label input metabolites must be 1.
  
- Define simulation times:

  For instance, to simulate label propagation from t=0 to 100, with time steps of 1:
  
  ```bash
  > times <- seq(0, 100, by=1)
  ```
  
  As another example, to simulate an exponential sampling frequency from t=0 to 15:
  
  ```bash
  > times <- round(10**(seq(0, log10(16), length.out=30))-1, 2)
  ```
  
- Perform simulations:

  ```bash
  > res <- simulate(net       = net,
                    kp        = fluxes,
                    anFun     = anFun,
                    p         = 0.0,
                    trf       = xch2fb,
                    meta_conc = meta_conc,
                    method    = "FORTRAN",
                    unloadDll = TRUE,
                    times     = times)
  ```
  			  
Simulation results are saved in subfolder `sim`:

  `plot.pdf`: plot of time-course concentrations of metabolites, fluxes, isotopologue abundances and isotopic enrichments of all EMUs
  
  `lib_f.f`:  FORTRAN code of the model used to compile the dynamic library
 
  `res.txt`:  complete simulation results


## Fit label inputs

- Define the experimental labeling dynamics of label input(s):

  Experimental labeling dynamics of label input(s) EMU(s) should be provided as a matrix containing mean molecular enrichments of each label input(s) EMU(s) (column) at each time (row).
  
  Here, to define the time-course mean enrichment of C as label input:
  
  ```bash
  > enr_input <- res$res_dyn$enrichments[, "C_1", drop=FALSE]
  ```
  
  Several label inputs can also be selected. For instance, to process C and O:
  
  ```bash
  > enr_input <- res$res_dyn$enrichments[, c("C_1", "O_1")]
  ```
  
- Fit labeling dynamics using analytical functions:

  ```bash
  > enr_in <- fit_label_input(enr_input, t=times, file="res_fit_enr", mc.cores=numCores)
  ```
  
Results are saved in subfolder `res_fit_enr`:

  `res_fit_enr_X_Y[-Z].pdf`: plot of exp. vs fitted labeling dynamics of the EMU `Y[-Z]` of metabolite `X`
  
  `res_fit_enr_X_Y[-Z].txt`: detailed results, including analytical functions

## Calculate fluxes

- Define the subsystem(s) to analyze:

  ```bash
  subsystem (list):
    $name (str):                name of the subsystem
    $rxn_subnet (list):         definition of the subnetwork of interest (see #construct-a-metabolic-model)
    $meta_conc_subnet (vector): named vector of initial metabolite concentrations (see #simulate-labeling-dynamics)
    $kp_subnet (vector):        named vector of model parameters (see #simulate-labeling-dynamics)
    $te_subnet (vector):        names of free parameters to estimate (can be model parameters and metabolite concentrations)
    $te_upc_subnet (vector):    named vector of upper bound constraints on free parameters
    $te_loc_subnet (vector):    named vector of lower bound constraints on free parameters
    $data_meas_subnet (list):   experimental data to fit, using the same format as label input data (see #fit-label-inputs)
    $sd_meas (list):            standard deviations on experimental data to fit
    $times (vector):            simulation times (all measurement times *must* be included)
    $enr_in (list):             list of fitted label inputs returned by `fit_label_input()`, as detailed above (see #fit-label-inputs)
    $anFun (list):              analytical functions (see #simulate-labeling-dynamics), otherwise should be NULL and $enr_in is used
    $niter (int):               number of Monte Carlo iterations (see #sensitivity-analysis-options) for flux calculations
    $mc.cores (int):            number of cores for parallelization (see #parallelization-options)
  ```
  
  For instance, to estimate the flux through r8 based on the minimal subsystem S<sub>E</sub>:
  
  ```bash
  > subsystem_1 <- list(name             = "S_E",
                        rxn_subnet       = list(r8   = list("R"=c("C", "E"), "C"=c(-1, 1), "E"="v8", "T"=c("A", "A")),
                                                rout = list("R"=c("E"),      "C"=c(-1),    "E"="v8", "T"=c("A"))),
                        meta_conc_subnet = c("C"=1, "G"=1.0, "E"=1.0),
                        kp_subnet        = c("v8"=0.5),
                        te_subnet        = c("v8", "E"),
                        te_upc_subnet    = c("v8"=10, "E"=100),
                        te_loc_subnet    = c("v8"=1e-4, "E"=0.01),
                        sd_meas          = list(iso=0.02, conc=c("E"=0.05)),
                        times            = times,
                        enr_in           = enr_in,
                        anFun            = NULL,
                        niter            = mc_iter,
                        mc.cores         = NULL,
                        data_meas_subnet = list(conc=c("E"=0.5), iso=cbind(times, "E_1"=res$res_dyn$enrichments[, "E_1"])))
  ```
  
  If you want to calculate fluxes for several subsystems, conditions or replicates at once, just gather the different subsystems/datasets into a 
  list:
  
  ```bash
  > subsystems <- list(subsystem_1 = subsystem_1,
                       ...)
  ```
  
- Calculate fluxes:

  Calculate fluxes for a single subsystems with:

  ```bash
  > res_sub <- fit_subsystems(subsystem_1)
  ```

  Or pass directly the full list of subsystems/datasets to analyze:

  ```bash
  > res_sub <- fit_subsystems(subsystems, mc.cores=numCores)
  ```

  Flux calculation results are saved in subfolder `fit_subnet_n` (where `n` is the name of the subsystem, i.e. which is in `$name`):

  `results.pdf`: plot of exp. vs fitted labeling dynamics of all metabolites
  
  `lib_f_n.f`:  FORTRAN code of the model used to compile the dynamic library
  
  `res.txt`:  complete flux calculation results

- Save the complete flux calculation results:

  To save detailed results (containing the network structure, experimental data, optimization details, etc), run:

  ```bash
  > list2file(res_sub, file="results_minimal_subsystems.txt")
  ```

  To save only results summary (estimated fluxes with their statistics), run:

  ```bash
  > list2file(res_sub$summary, file="summary_minimal_subsystems.txt")
  ```

- Display results summary:
  
  ```bash
  > print(res_sub$summary)
  ```

## Return to the initial working directory

```bash
> setwd(wd)
```


