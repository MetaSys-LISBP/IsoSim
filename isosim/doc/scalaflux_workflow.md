# ScalaFlux workflow

This documentation describes a typical workflow to quantify fluxes in metabolic subnetworks 
using the ScalaFlux approach.

IsoSim must be installed beforehand. Have a look [here](https://github.com/MetaSys-LISBP/IsoSim) for installation instructions.

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

*Note:* A message will appear if some of the required packages are missing. In this case, follow the instructions to install the required package(s), then reload IsoSim.


## Parallelization options

Flux calculation can be fasten by parallelizing the calculation.

Parallelization options can be adapted with the variable `numCores` (int), which represents the number of CPU cores 
to use in parallel. To use a 
single-core version, set `numCores` to NULL.

To use all available CPU cores on the current host:

```bash
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

network (list):
  $R (vector): reactant(s)
  $C (vector): stoichiometric coefficient(s)
  $E (vector): rate law (can be a constant, a variable or an expression)
  $T (vector): tracer atom transitions

For instance, to construct the model of the example network provided in figure 1A of ScalaFlux publication:
  
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

Once you have constructed the model, we strongly encourage you to run some simulations and to verify that simulation results correspond to the 
expected behaviour (e.g. metabolite concentrations and fluxes should be constant, assuming the modeled system operates at metabolic steady-state).

- Define the labeling dynamics of label input(s) (here Sout) using analytical functions 
(here we simulate a switch from unlabeled to fully-labeled nutrient)

```bash
> anFun <- list("Sout_1-M0" = "0.0+0.0*t",
			  "Sout_1-M1" = "1.0+0.0*t")
```

Each element should contain the time-dependent analytical function of an isotopologue of all EMUs of label inputs, with the corresponding name. The required list 
of label input EMUs isotopologues is automatically identified when constructing the model (see `net$min_meas`).

Note: all the analytical functions *must* contain the time variable `t`, even if the label input(s) is (are) constant
			  
- Set initial values of fluxes and metabolite concentrations:

```bash
> fluxes <- c("v1"=2, "v2"=1.5, "v6"=1.2, "v10"=1.0, "v6xch"=0.2)
> meta_conc <- c("Sout"=1, "Sin"=0.5, "A"=0.5, "B"=0.5, "C"=0.5, "D"=0.5, "E"=0.5, "F"=0.5, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5, "O"=0.5, "P"=20)
```

- Define simulation times (here to simulate an exponential sampling frequency from 0 to 15 min):

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
			  
All simulation results (time course concentrations of metabolites, fluxes, isotopologue abundances and isotopic enrichments of all EMUs) are saved 
in a folder `sim`.


## Fit label inputs

- Experimental labeling dynamics of label input(s) to fit (here we fit the theoretical labeling dynamics of C):

```bash
> enr_input <- res$res_dyn$enrichments[, "C_1"]
```

- Fit labeling dynamics using analytical functions:

```bash
> enr_in <- fit_label_input(enr_input, t=times, file="res_fit_enr_C", mc.cores=numCores)
```

## Calculate fluxes

- Define the minimal subsystem(s) to analyze:

subsystem (list):
  $rxn_subnet (list): network definition (list), as detailed above
  $meta_conc_subnet (vector): named vector of initial metabolite concentrations
  $kp_subnet (vector): named vector of model parameters
  $te_subnet (vector): free parameters to estimate (can be model parameters and metabolite concentrations)
  $te_upc_subnet (vector): named vector of upper bound constraints on free parameters
  $te_loc_subnet (vector): named vector of lower bound constraints on free parameters
  $data_meas_subnet (list): experimental data to fit
  $sd_meas (list): standard deviations on experimental data to fit
  $times (vector): simulation times (all measurement times must be included)
  $enr_in (list): list of fitted label inputs returned by `fit_label_input()`
  $anFun (list): analytical functions (if any, otherwise should be NULL)
  $niter (int): number of Monte Carlo iterations for flux calculations
  $mc.cores (int): number of cores for parallelization

The flux through r8 can be estimated based on the minimal subsystem `S_E` of the example network:

```bash
> subsystem_1 <- list(name = "S_E",
		  rxn_subnet = list(r8     = list("R"=c("C", "E"), "C"=c(-1, 1), "E"="v8",      "T"=c("A", "A")),
							rout   = list("R"=c("E"),      "C"=c(-1),    "E"="v8", "T"=c("A"))),
		  meta_conc_subnet = c("C"=1, "G"=1.0, "E"=1.0),
		  kp_subnet = c("v8"=0.5),
		  te_subnet = c("v8", "E"),
		  te_upc_subnet = c("v8"=10, "E"=100),
		  te_loc_subnet = c("v8"=1e-4, "E"=0.01),
		  sd_meas = list(iso=0.02, conc=c("E"=0.01)),
		  times = times,
		  enr_in = enr_in,
		  anFun = NULL,
		  niter = mc_iter,
		  mc.cores = NULL,
		  data_meas_subnet = list(conc=c("E"=0.5), iso=cbind(times, "E_1"=res$res_dyn$enrichments[, "E_1"])))
```

If you want to calculate fluxes for different subsystems, conditions or replicates at once, just gather all the subsystems of interest into a 
list:

```bash
> subsystems <- list(subsystem_1=subsystem_1, ...)
```

- Calculate fluxes:

To calculate fluxes for a single subsystems with:

```bash
> res_sub <- fit_subsystem(subsystem_1)
```

if you want to calculate fluxes for all subsystems at once:

```bash
> res_subs <- fit_subsystems(subsystems, dirname="fit_subsystems", mc.cores=numCores)
```

- Save a summary of the flux calculation results:

```bash
> list2file(res_sub$summary, file="summary_minimal_subsystems.txt")
```

- Display the flux calculation results:

```bash
> print()
```
