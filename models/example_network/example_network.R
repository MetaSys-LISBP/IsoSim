# millard@insa-toulouse.fr
#
# IsoSim: stoichiometric, kinetic and isotopic modeling of metabolic systems
#
# https://github.com/MetaSys-LISBP/IsoSim
#
# https://doi.org/10.1101/735308
#
# Copyright 2019, INRA, France
# License: GNU General Public License v3 (see license.txt for details)

####################################
cat("
This code calculates fluxes through several subsystems of the
example network shown in figure 1A of the ScalaFlux publication.

It also generates figures 3 and 4.

")
####################################

####################################
cat("   ... initialize R environment ...\n\n")
####################################

# Clean workspace
rm(list= ls())

library("vioplot")

# Path of the example model
#setwd("C:/Users/millard/Documents/GIT/IsoSim/IsoSim/models/example_network")

# Get current directory
wd <- getwd()

# Load IsoSim
setwd(file.path(dirname(dirname(wd)), "isosim"))
source("isosim.R")

# Go back to working directory and create "res" folder to store the results
setwd(wd)
if (!file.exists("res")){
    dir.create(file.path(wd, "res"))
}
setwd(file.path(wd, "res"))


#########################
### GLOBAL PARAMETERS ###
#########################

# number of cores to use in parallel, i.e. at most how many child processes will be run simultaneously
# single-core version can be used by setting 'numCores' to NULL
numCores <- detectCores()

# number of Monte Carlo iterations for flux calculation
mc_iter <- 4


####################################
cat("\n   ... construct isotopic model of the example network ...\n\n")
####################################

# network definition
# reactions: Reactants, Coefficients, Rate law, Atom transitions
rxn <- list(r1     = list("R"=c("Sout", "Sin"), "C"=c(-1, 1),     "E"="v1",   "T"=c("A", "A")),
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

# equations for determined fluxes
eq_det <- c("v3 = v1-v2",
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

# construct the model
net <- net2mat(rxn, add_eq=eq_det)

####################################
cat("\n   ... simulate labeling dynamics ...\n\n")
####################################

# labeling dynamics of label input (here we simulate a switch from unlabeled to fully-labeled nutrient Sout)
anFun <- list("Sout_1-M0" = "0.0+0.0*t",
              "Sout_1-M1" = "1.0+0.0*t")

# fluxes and metabolite concentrations
fluxes <- c("v1"=2, "v2"=1.5, "v6"=1.2, "v10"=1.0, "v6xch"=0.2)
meta_conc <- c("Sout"=1, "Sin"=0.5, "A"=0.5, "B"=0.5, "C"=0.5, "D"=0.5, "E"=0.5, "F"=0.5, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5, "O"=0.5, "P"=20)

# simulation times (exponential sampling frequency from 0 to 15 min)
times <- round(10**(seq(0, log10(16), length.out=30))-1, 2)

# run simulations
res <- simulate(net       = net,
                kp        = fluxes,
                anFun     = anFun,
                p         = 0.0,
                trf       = xch2fb,
                meta_conc = meta_conc,
                method    = "FORTRAN",
                unloadDll = TRUE,
                times     = times)


####################################
cat("\n   ... fit labeling dynamics of all metabolites ...\n\n")
####################################

# fit analytical functions to the labeling dynamics of all metabolites
enr_input <- res$res_dyn$enrichments[,c("Sout_1", "Sin_1", "A_1", "B_1", "C_1", "D_1", "E_1", "F_1", "G_1", "H_1-2", "H_2", "H_1", "I_1", "J_1", "K_1", "L_1", "M_1", "N_1", "O_1")]

# fit analytical functions to the labeling dynamics of all metabolites
enr_in <- fit_label_input(enr_input, t=times, file="res_fit_enr", mc.cores=numCores)


####################################
cat("\n   ... calculate fluxes for all minimal subsystems ...\n\n")
####################################

# definition of all minimal subsystems to analyze
# here, fluxes and pools are estimated for all minimal subsystems of the network provided as example
subsystems <- list(
    s1 = list(name = "s1",
              rxn_subnet = list(r1   = list("R"=c("Sout", "Sin"), "C"=c(-1, 1), "E"="v1", "T"=c("A", "A")),
                                rout = list("R"=c("Sin"),        "C"=c(-1),    "E"="v1", "T"=c("A"))),
              meta_conc_subnet = c("Sin"=1, "Sout"=1.0),
              kp_subnet = c("v1"=0.5),
              te_subnet = c("v1", "Sin"),
              te_upc_subnet = c("v1"=10, "Sin"=100),
              te_loc_subnet = c("v1"=1e-4, "Sin"=0.01),
              sd_meas = list(iso=0.02, conc=c("Sin"=0.01)),
              times = times,
              enr_in = enr_in,
              anFun = NULL,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("Sin"=0.5), iso=cbind(times, "Sin_1"=res$res_dyn$enrichments[, "Sin_1"]))),
    
    s2 = list(name = "s2",
              rxn_subnet = list(r2   = list("R"=c("Sin", "A"), "C"=c(-1, 1), "E"="v2", "T"=c("A", "A")),
                                rout = list("R"=c("A"),        "C"=c(-1),    "E"="v2", "T"=c("A"))),
              meta_conc_subnet = c("Sin"=1, "A"=1.0),
              kp_subnet = c("v2"=0.5),
              te_subnet = c("v2", "A"),
              te_upc_subnet = c("v2"=10, "A"=100),
              te_loc_subnet = c("v2"=1e-4, "A"=0.01),
              sd_meas = list(iso=0.02, conc=c("A"=0.01)),
              times = times,
              enr_in = enr_in,
              anFun = NULL,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("A"=0.5), iso=cbind(times, "A_1"=res$res_dyn$enrichments[, "A_1"]))),
    
    s3 = list(name = "s3",
              rxn_subnet = list(r3   = list("R"=c("Sin", "D"), "C"=c(-1, 1), "E"="v3", "T"=c("A", "A")),
                                rout = list("R"=c("D"),        "C"=c(-1),    "E"="v3", "T"=c("A"))),
              meta_conc_subnet = c("Sin"=1, "D"=1.0),
              kp_subnet = c("v3"=1.2),
              te_subnet = c("v3", "D"),
              te_upc_subnet = c("v3"=10, "D"=100),
              te_loc_subnet = c("v3"=1e-4, "D"=0.01),
              sd_meas = list(iso=0.02, conc=c("D"=0.01)),
              times = times,
              enr_in = enr_in,
              anFun = NULL,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("D"=0.5), iso=cbind(times, "D_1"=res$res_dyn$enrichments[, "D_1"]))),
    
    s4_9 = list(name = "s4_9",
                rxn_subnet = list(r4   = list("R"=c("A", "F"), "C"=c(-1, 1), "E"="v4",    "T"=c("A", "A")),
                                  r9   = list("R"=c("E", "F"), "C"=c(-1, 1), "E"="v9",    "T"=c("A", "A")),
                                  rout = list("R"=c("F"),      "C"=c(-1),    "E"="v4+v9", "T"=c("A"))),
                meta_conc_subnet = c("A"=1, "E"=1, "F"=1.0),
                kp_subnet = c("v4"=0.5, "v9"=0.5),
                te_subnet = c("v4", "v9", "F"),
                te_upc_subnet = c("v4"=10, "v9"=10, "F"=100),
                te_loc_subnet = c("v4"=1e-4, "v9"=1e-4, "F"=0.01),
                sd_meas = list(iso=0.02, conc=c("F"=0.01)),
                times = times,
                enr_in = enr_in,
                anFun = NULL,
                niter = mc_iter,
                mc.cores = NULL,
                data_meas_subnet = list(conc=c("F"=0.5), iso=cbind(times, "F_1"=res$res_dyn$enrichments[, "F_1"]))),
    
    s5 = list(name = "s5",
              rxn_subnet = list(r5    = list("R"=c("A", "B"), "C"=c(-1, 1), "E"="v5",     "T"=c("A", "A")),
                                r6rev = list("R"=c("C", "B"), "C"=c(-1, 1), "E"="v6r",    "T"=c("A", "A")),
                                rout  = list("R"=c("B"),      "C"=c(-1),    "E"="v5+v6r", "T"=c("A"))),
              meta_conc_subnet = c("A"=1, "B"=1.0, "C"=1.0),
              kp_subnet = c("v5"=0.5, "v6r"=0.25),
              te_subnet = c("v5", "v6r", "B"),
              te_upc_subnet = c("v5"=10, "v6r"=10, "B"=100),
              te_loc_subnet = c("v5"=1e-4, "v6r"=1e-4, "B"=0.01),
              sd_meas = list(iso=0.02, conc=c("B"=0.01)),
              times = times,
              enr_in = enr_in,
              anFun = NULL,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("B"=0.5), iso=cbind(times, "B_1"=res$res_dyn$enrichments[, "B_1"]))),
    
    s6_7 = list(name = "s6_7",
                rxn_subnet = list(r6for = list("R"=c("B", "C"), "C"=c(-1, 1), "E"="v6f",    "T"=c("A", "A")),
                                  r7    = list("R"=c("D", "C"), "C"=c(-1, 1), "E"="v7",     "T"=c("A", "A")),
                                  rout  = list("R"=c("C"),      "C"=c(-1),    "E"="v6f+v7", "T"=c("A"))),
                meta_conc_subnet = c("B"=1, "C"=1.0, "D"=1.0),
                kp_subnet = c("v6f"=1.5, "v7"=0.7),
                te_subnet = c("v6f", "v7", "C"),
                te_upc_subnet = c("v6f"=10, "v7"=10, "C"=100),
                te_loc_subnet = c("v6f"=1e-4, "v7"=1e-4, "C"=0.01),
                sd_meas = list(iso=0.02, conc=c("C"=0.01)),
                times = times,
                enr_in = enr_in,
                anFun = NULL,
                niter = mc_iter,
                mc.cores = NULL,
                data_meas_subnet = list(conc=c("C"=0.5), iso=cbind(times, "C_1"=res$res_dyn$enrichments[, "C_1"]))),
    
    s8 = list(name = "s8",
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
              data_meas_subnet = list(conc=c("E"=0.5), iso=cbind(times, "E_1"=res$res_dyn$enrichments[, "E_1"]))),
    
    s10 = list(name = "s10",
               rxn_subnet = list(r10for = list("R"=c("E", "G"), "C"=c(-1, 1), "E"="v10", "T"=c("A", "A")),
                                 rout   = list("R"=c("G"),      "C"=c(-1),    "E"="v10", "T"=c("A"))),
               meta_conc_subnet = c("E"=1, "G"=1.0),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "G"),
               te_upc_subnet = c("v10"=10, "G"=100),
               te_loc_subnet = c("v10"=1e-4, "G"=0.01),
               sd_meas = list(iso=0.02, conc=c("G"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("G"=0.5), iso=cbind(times, "G_1"=res$res_dyn$enrichments[, "G_1"]))),
    
    s11 = list(name = "s11",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v11", "T"=c("A", "B", "AB")),
                                 rout = list("R"=c("H"),           "C"=c(-1),        "E"="v11", "T"=c("AB"))),
               meta_conc_subnet = c("G"=1, "K"=1.0, "H"=1.0),
               kp_subnet = c("v11"=0.5),
               te_subnet = c("v11", "H"),
               te_upc_subnet = c("v11"=10, "H"=100),
               te_loc_subnet = c("v11"=1e-4, "H"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5), iso=cbind(times, "H_1-2"=res$res_dyn$enrichments[, "H_1-2"]))),
    
    s12 = list(name = "s12",
               rxn_subnet = list(r12   = list("R"=c("H", "I", "J"),"C"=c(-1, 1, 1), "E"="v12", "T"=c("AB", "B", "A")),
                                 rout1 = list("R"=c("I"),          "C"=c(-1),       "E"="v12", "T"=c("A")),
                                 rout2 = list("R"=c("J"),          "C"=c(-1),       "E"="v12", "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=1.0, "J"=1.0),
               kp_subnet = c("v12"=0.5),
               te_subnet = c("v12", "I", "J"),
               te_upc_subnet = c("v12"=10, "I"=100, "J"=100),
               te_loc_subnet = c("v12"=1e-4, "I"=0.01, "J"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1")]))),
    
    s13 = list(name = "s13",
               rxn_subnet = list(r13  = list("R"=c("I", "K"), "C"=c(-1, 1), "E"="v13", "T"=c("A", "A")),
                                 rout = list("R"=c("K"),      "C"=c(-1),    "E"="v13", "T"=c("A"))),
               meta_conc_subnet = c("I"=1, "K"=1.0),
               kp_subnet = c("v13"=0.5),
               te_subnet = c("v13", "K"),
               te_upc_subnet = c("v13"=10, "K"=100),
               te_loc_subnet = c("v13"=1e-4, "K"=0.01),
               sd_meas = list(iso=0.02, conc=c("K"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("K"=0.5), iso=cbind(times, "K_1"=res$res_dyn$enrichments[, "K_1"]))),
    
    s14 = list(name = "s14",
               rxn_subnet = list(r14  = list("R"=c("J", "L"), "C"=c(-1, 1), "E"="v14", "T"=c("A", "A")),
                                 rout = list("R"=c("L"),      "C"=c(-1),    "E"="v14", "T"=c("A"))),
               meta_conc_subnet = c("J"=1, "L"=1.0),
               kp_subnet = c("v14"=0.5),
               te_subnet = c("v14", "L"),
               te_upc_subnet = c("v14"=10, "L"=100),
               te_loc_subnet = c("v14"=1e-4, "L"=0.01),
               sd_meas = list(iso=0.02, conc=c("L"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("L"=0.5), iso=cbind(times, "L_1"=res$res_dyn$enrichments[, "L_1"]))),
    
    s15 = list(name = "s15",
               rxn_subnet = list(r15  = list("R"=c("L", "M"), "C"=c(-1, 1), "E"="v15", "T"=c("A", "A")),
                                 rout = list("R"=c("M"),      "C"=c(-1),    "E"="v15", "T"=c("A"))),
               meta_conc_subnet = c("L"=1, "M"=1.0),
               kp_subnet = c("v15"=0.5),
               te_subnet = c("v15", "M"),
               te_upc_subnet = c("v15"=10, "M"=100),
               te_loc_subnet = c("v15"=1e-4, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("M"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("M"=0.5), iso=cbind(times, "M_1"=res$res_dyn$enrichments[, "M_1"]))),
    
    s16 = list(name = "s16",
               rxn_subnet = list(r16  = list("R"=c("M", "N"), "C"=c(-1, 1), "E"="v16", "T"=c("A", "A")),
                                 rout = list("R"=c("N"),      "C"=c(-1),    "E"="v16", "T"=c("A"))),
               meta_conc_subnet = c("M"=1, "N"=1.0),
               kp_subnet = c("v16"=0.5),
               te_subnet = c("v16", "N"),
               te_upc_subnet = c("v16"=10, "N"=100),
               te_loc_subnet = c("v16"=1e-4, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("N"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("N"=0.5), iso=cbind(times, "N_1"=res$res_dyn$enrichments[, "N_1"]))),
    
    s17_18 = list(name = "s17_18",
                  rxn_subnet = list(r17  = list("R"=c("N", "O"), "C"=c(-1, 1), "E"="v17",     "T"=c("A", "A")),
                                    r18  = list("R"=c("F", "O"), "C"=c(-1, 1), "E"="v18",     "T"=c("A", "A")),
                                    rout = list("R"=c("O"),      "C"=c(-1),    "E"="v17+v18", "T"=c("A"))),
                  meta_conc_subnet = c("N"=1, "F"=1.0, "O"=1.0),
                  kp_subnet = c("v17"=0.5, "v18"=0.5),
                  te_subnet = c("v17", "v18", "O"),
                  te_upc_subnet = c("v17"=10, "v18"=10, "O"=100),
                  te_loc_subnet = c("v17"=1e-4, "v18"=1e-4, "O"=0.01),
                  sd_meas = list(iso=0.02, conc=c("O"=0.01)),
                  times = times,
                  enr_in = enr_in,
                  anFun = NULL,
                  niter = mc_iter,
                  mc.cores = NULL,
                  data_meas_subnet = list(conc=c("O"=0.5), iso=cbind(times, "O_1"=res$res_dyn$enrichments[, "O_1"]))),
    
    s19 = list(name = "s19",
               rxn_subnet = list(r19  = list("R"=c("O", "P"), "C"=c(-1, 1), "E"="v19", "T"=c("A", "A")),
                                 rout = list("R"=c("P"),      "C"=c(-1),    "E"="v19", "T"=c("A"))),
               meta_conc_subnet = c("O"=1, "P"=1.0),
               kp_subnet = c("v19"=0.5),
               te_subnet = c("v19", "P"),
               te_upc_subnet = c("v19"=10, "P"=100),
               te_loc_subnet = c("v19"=1e-4, "P"=0.01),
               sd_meas = list(iso=0.02, conc=c("P"=0.01)),
               times = times,
               enr_in = enr_in,
               anFun = NULL,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("P"=20), iso=cbind(times, "P_1"=res$res_dyn$enrichments[, "P_1"])))
)

# calculate fluxes for all subsystems
res_sub <- fit_subsystems(subsystems, dirname="fit_minimal_subsystems", mc.cores=numCores)

# get summary
list2file(res_sub$summary, file="summary_minimal_subsystems.txt")

####################################
cat("\n      ... plot results ...\n\n")
####################################

all_cols <- fun_col(length(res_sub$res))
cols <- c()
data_plot <- list()
opt <- c()
ci95p <- c()
ci95m <- c()
for (i in seq(length(res_sub$res))){
    n <- names(res_sub$res)[i]
    a <- colnames(res_sub$res[[i]]$sens$par)
    ncol <- grep("v", a)
    # colors
    cols <- c(cols, rep(all_cols[i], length(ncol)))
    # get optimal value
    no <- names(opt)
    opt <- c(opt, res_sub$summary[[n]][ncol, "opt"])
    names(opt) <- c(no, a[ncol])
    # get higher bound of 95% CI
    nci95p <- names(ci95p)
    ci95p <- c(ci95p, res_sub$summary[[n]][ncol, "ci_97.5"])
    names(ci95p) <- c(no, a[ncol])
    # get lower bound of 95% CI
    nci95m <- names(ci95m)
    ci95m <- c(ci95m, res_sub$summary[[n]][ncol, "ci_2.5"])
    names(ci95m) <- c(no, a[ncol])
    if (length(ncol) == 1){
        data_plot[[a[1]]] <- res_sub$res[[i]]$sens$par[,ncol]
    }else if (length(ncol)>1){
        data_plot <- c(data_plot, as.list(as.data.frame(res_sub$res[[i]]$sens$par[,ncol])))
    }
}
odn <- names(data_plot)
true_values <- c(v1=2, v2=1.5, v3=0.5, v4=0.3, v5=1.2, v6r=0.25, v6f=1.45, v7=0.5, v8=1.7, v9=0.7, v10=1.0, v11=1, v12=1, v13=1, v14=1, v15=1, v16=1, v17=1, v18=1, v19=2)

# supplementary Figure S1
pdf(file="Fig_S1.pdf", width=10, height=8)
par(mfrow=c(4,5))
for (i in names(enr_in$all)){
    plot(times, enr_in$all[[i]]$exp, las=1, main=i, ylim=c(0,1), xlab="time", ylab="13C-enrichment", type="p", pch=20)
    lines(x=seq(from=min(times), to=max(times), length.out=100), enr_in$all[[i]]$sim)
}
dev.off()

# supplementary Figure S2
pdf(file="Fig_S2.pdf", width=10, height=8)
par(mfrow=c(4,5))
for (i in names(subsystems)){
    i_mes <- subsystems[[i]]$data_meas_subnet
    i_mes_e <- i_mes$iso[,-1]
    i_sim_all <- res_sub$res[[i]]$result$retres$sim
    if (is.matrix(i_mes_e)){
        matplot(times, i_mes_e, las=1, main=i, ylim=c(0,1), xlab="time", ylab="13C-enrichment", type="p", pch=20)
        matlines(times, i_sim_all)
    }else{
        plot(times, i_mes_e, las=1, main=i, ylim=c(0,1), xlab="time", ylab="13C-enrichment", type="p", pch=20)
        lines(times, i_sim_all)
    }
}
dev.off()

# plot Figure 3B
pdf(file="Fig_3B.pdf", width=7, height=6)
    plot(x=true_values[odn], y=opt[odn], xlim=c(0,2.5), col=cols, ylim=c(0,2.5), las=1, yaxs="i", xaxs="i", pch=16, xlab="true value", ylab="estimated value")
    abline(c(0,0), 1, lwd=2, lty=2)
    cor_res <- cor.test(x=true_values[odn], y=opt[odn], method = "pearson")
    points(x=true_values[odn], y=opt[odn], bg=cols, pch=21, cex=2)
    text(x=0.5, y=2, labels=paste("R2 = ", round(cor_res$estimate, 2), "\np = ", signif(cor_res$p.value, 2), sep=""))
dev.off()

# plot Figure 3C
pdf(file="Fig_3C.pdf", width=10, height=6)
    vioplot(data_plot[odn], names=odn, col=cols, las=2, ylim=c(0,2.5), outline=FALSE)
    points(true_values[odn], pch=16, col="red")
dev.off()


####################################
cat("\n   ... calculate the flux through the pathway r9-r16 based on different subsystems ...\n\n")
####################################

subsystems_l2 <- list(
    S1 = list(name = "S1",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   rout = list("R"=c("G"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G"),
                 te_upc_subnet = c("v10"=10, "G"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5), iso=cbind(times, "G_1"=res$res_dyn$enrichments[, "G_1"]))),
    S2 = list(name = "S2",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                   rout = list("R"=c("H"),           "C"=c(-1),        "E"="v10",  "T"=c("AB"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "K"=1),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G", "H"),
                 te_upc_subnet = c("v10"=10, "G"=100, "H"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2")]))),
    S3 = list(name = "S3",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                   r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                   rout1 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                   rout2 = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=1),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G", "H", "I", "J"),
                 te_upc_subnet = c("v10"=10, "G"=100, "H"=100, "I"=100, "J"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2", "I_1", "J_1")]))),
    S4 = list(name = "S4",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                   r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                   r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                   rout = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G", "H", "I", "J", "K"),
                 te_upc_subnet = c("v10"=10, "G"=100, "H"=100, "I"=100, "J"=100, "K"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2", "I_1", "J_1", "K_1")]))),
    S5 = list(name = "S5",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                   r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                   r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                   r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   rout = list("R"=c("L"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G", "H", "I", "J", "K", "L"),
                 te_upc_subnet = c("v10"=10, "G"=100, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2", "I_1", "J_1", "K_1", "L_1")]))),
    S6 = list(name = "S6",
                 rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                   r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                   r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                   r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                   rout = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
                 meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5),
                 kp_subnet = c("v10"=0.5),
                 te_subnet = c("v10", "G", "H", "I", "J", "K", "L", "M"),
                 te_upc_subnet = c("v10"=10, "G"=100, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100),
                 te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01),
                 sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01)),
                 times = times,
                 enr_in = enr_in,
                 niter = mc_iter,
                 mc.cores = NULL,
                 data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2", "I_1", "J_1", "K_1", "L_1", "M_1")]))),
    S7 = list(name = "S7",
               rxn_subnet = list(r10  = list("R"=c("E", "G"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("E"=1, "G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "G", "H", "I", "J", "K", "L", "M", "N"),
               te_upc_subnet = c("v10"=10, "G"=100, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("G"=0.01, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("G"=0.5, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("G_1", "H_1-2", "I_1", "J_1", "K_1", "L_1", "M_1", "N_1")]))),
    S8 = list(name = "S8",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 rout = list("R"=c("H"),           "C"=c(-1),        "E"="v10",  "T"=c("AB"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "K"=1),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H"),
               te_upc_subnet = c("v10"=10, "H"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5), iso=cbind(times, "H_1-2"=res$res_dyn$enrichments[, "H_1-2"]))),
    S9 = list(name = "S9",
              rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                rout1 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                rout2 = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
              meta_conc_subnet = c("G"=1, "H"=0.5, "I"=0.5, "J"=0.5, "K"=1),
              kp_subnet = c("v10"=0.5),
              te_subnet = c("v10", "H", "I", "J"),
              te_upc_subnet = c("v10"=10, "H"=100, "I"=100, "J"=100),
              te_loc_subnet = c("v10"=1e-4, "H"=0.01, "I"=0.01, "J"=0.01),
              sd_meas = list(iso=0.02, conc=c("H"=0.01, "I"=0.01, "J"=0.01)),
              times = times,
              enr_in = enr_in,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("H"=0.5, "I"=0.5, "J"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "I_1", "J_1")]))),
    S10 = list(name = "S10",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 rout = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H", "I", "J", "K"),
               te_upc_subnet = c("v10"=10, "H"=100, "I"=100, "J"=100, "K"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "I_1", "J_1", "K_1")]))),
    S11 = list(name = "S11",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("L"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H", "I", "J", "K", "L"),
               te_upc_subnet = c("v10"=10, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "I_1", "J_1", "K_1", "L_1")]))),
    S12 = list(name = "S12",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H", "I", "J", "K", "L", "M"),
               te_upc_subnet = c("v10"=10, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "I_1", "J_1", "K_1", "L_1", "M_1")]))),
    S13 = list(name = "S13",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H", "I", "J", "K", "L", "M", "N"),
               te_upc_subnet = c("v10"=10, "H"=100, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "I_1", "J_1", "K_1", "L_1", "M_1", "N_1")]))),
    S14 = list(name = "S14",
              rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                rout1 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                rout2 = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
              meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "K"=1),
              kp_subnet = c("v10"=0.5),
              te_subnet = c("v10", "I", "J"),
              te_upc_subnet = c("v10"=10, "I"=100, "J"=100),
              te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01),
              sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01)),
              times = times,
              enr_in = enr_in,
              niter = mc_iter,
              mc.cores = NULL,
              data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1")]))),
    S15 = list(name = "S15",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 rout1 = list("R"=c("J"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("K"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "K"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "K"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "K"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "K"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "K"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "K"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "K_1")]))),
    S16 = list(name = "S16",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("L"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("K"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "K", "L"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "K"=100, "L"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "K_1", "L_1")]))),
    S17 = list(name = "S17",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("B", "B")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("K"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "K", "L", "M"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "K_1", "L_1", "M_1")]))),
    S18 = list(name = "S18",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("K"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "K", "L", "M", "N"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "K"=100, "L"=100, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "K"=0.01, "L"=0.01, "M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "K"=0.5, "L"=0.5, "M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "K_1", "L_1", "M_1", "N_1")]))),
    S19 = list(name = "S19",
               rxn_subnet = list(r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("K"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("I"=1, "K"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "K"),
               te_upc_subnet = c("v10"=10, "K"=100),
               te_loc_subnet = c("v10"=1e-4, "K"=0.01),
               sd_meas = list(iso=0.02, conc=c("K"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("K"=0.5), iso=cbind(times, "K_1"=res$res_dyn$enrichments[, "K_1"]))),
    S20 = list(name = "S20",
               rxn_subnet = list(r11  = list("R"=c("G", "K", "H"), "C"=c(-1, -1, 1), "E"="v10",  "T"=c("A", "B", "BA")),
                                 r13  = list("R"=c("I", "K"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("H"),           "C"=c(-1),        "E"="v10",  "T"=c("AB"))),
               meta_conc_subnet = c("G"=1, "H"=0.5, "I"=1, "K"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "H", "K"),
               te_upc_subnet = c("v10"=10, "H"=100, "K"=100),
               te_loc_subnet = c("v10"=1e-4, "H"=0.01, "K"=0.01),
               sd_meas = list(iso=0.02, conc=c("H"=0.01, "K"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("H"=0.5, "K"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("H_1-2", "K_1")]))),
    S21 = list(name = "S21",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("L"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "L"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "L"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "L"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "L"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "L"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "L"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "L_1")]))),
    S22 = list(name = "S22",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "L"=0.5, "M"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "L", "M"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "L"=100, "M"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "L"=0.01, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "L"=0.01, "M"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "L"=0.5, "M"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "L_1", "M_1")]))),
    S23 = list(name = "S23",
               rxn_subnet = list(r12  = list("R"=c("H", "I", "J"), "C"=c(-1, 1, 1),  "E"="v10",  "T"=c("AB", "B", "A")),
                                 r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout1 = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A")),
                                 rout2 = list("R"=c("I"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("H"=1, "I"=0.5, "J"=0.5, "L"=0.5, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "I", "J", "L", "M", "N"),
               te_upc_subnet = c("v10"=10, "I"=100, "J"=100, "L"=100, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "I"=0.01, "J"=0.01, "L"=0.01, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("I"=0.01, "J"=0.01, "L"=0.01, "M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("I"=0.5, "J"=0.5, "L"=0.5, "M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("I_1", "J_1", "L_1", "M_1", "N_1")]))),
    S24 = list(name = "S24",
               rxn_subnet = list(r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("L"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("J"=1, "L"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "L"),
               te_upc_subnet = c("v10"=10, "L"=100),
               te_loc_subnet = c("v10"=1e-4, "L"=0.01),
               sd_meas = list(iso=0.02, conc=c("L"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("L"=0.5), iso=cbind(times, "L_1"=res$res_dyn$enrichments[, "L_1"]))),
    S25 = list(name = "S25",
               rxn_subnet = list(r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("J"=1, "L"=0.5, "M"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "L", "M"),
               te_upc_subnet = c("v10"=10, "L"=100, "M"=100),
               te_loc_subnet = c("v10"=1e-4, "L"=0.01, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("L"=0.01, "M"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("L"=0.5, "M"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("L_1", "M_1")]))),
    S26 = list(name = "S26",
               rxn_subnet = list(r14  = list("R"=c("J", "L"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("J"=1, "L"=0.5, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "L", "M", "N"),
               te_upc_subnet = c("v10"=10, "L"=100, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "L"=0.01, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("L"=0.01, "M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("L"=0.5, "M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("L_1", "M_1", "N_1")]))),
    S27 = list(name = "S27",
               rxn_subnet = list(r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("M"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("L"=1, "M"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "M"),
               te_upc_subnet = c("v10"=10, "M"=100),
               te_loc_subnet = c("v10"=1e-4, "M"=0.01),
               sd_meas = list(iso=0.02, conc=c("M"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("M"=0.5), iso=cbind(times, "M_1"=res$res_dyn$enrichments[, "M_1"]))),
    S28 = list(name = "S28",
               rxn_subnet = list(r15  = list("R"=c("L", "M"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("L"=1, "M"=0.5, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "M", "N"),
               te_upc_subnet = c("v10"=10, "M"=100, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "M"=0.01, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("M"=0.01, "N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("M"=0.5, "N"=0.5), iso=cbind(times, res$res_dyn$enrichments[, c("M_1", "N_1")]))),
    S29 = list(name = "S29",
               rxn_subnet = list(r16  = list("R"=c("M", "N"),      "C"=c(-1, 1),     "E"="v10",  "T"=c("A", "A")),
                                 rout = list("R"=c("N"),           "C"=c(-1),        "E"="v10",  "T"=c("A"))),
               meta_conc_subnet = c("M"=1, "N"=0.5),
               kp_subnet = c("v10"=0.5),
               te_subnet = c("v10", "N"),
               te_upc_subnet = c("v10"=10, "N"=100),
               te_loc_subnet = c("v10"=1e-4, "N"=0.01),
               sd_meas = list(iso=0.02, conc=c("N"=0.01)),
               times = times,
               enr_in = enr_in,
               niter = mc_iter,
               mc.cores = NULL,
               data_meas_subnet = list(conc=c("N"=0.5), iso=cbind(times, "N_1"=res$res_dyn$enrichments[, "N_1"])))
)


# calculate fluxes for all subsystems of the given flux module
res_sub_l2 <- fit_subsystems(subsystems_l2, dirname="fit_subsystems_r9-16", mc.cores=numCores)

# display results
res_sub_l2$summary

# get summary
list2file(res_sub_l2$summary, file="summary_subsystems_r9-16.txt")

####################################
cat("\n      ... plot results ...\n\n")
####################################

# create figure 4
pdf(file="Fig_4.pdf", width=10, height=10)
    par(mfrow=c(2,1))
    par(mar=c(5.1,4.1,4.1,2.1))
    # data used to plot Figure 4B
    msub <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1,
              2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 0, 3, 2, 2, 2, 0, 0, 0, 0, 0, 0,
              0, 0, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 0, 0, 0, 0, 0, 0,
              0, 0, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 3, 3, 3, 2, 2, 2, 0, 0, 0,
              0, 2, 2, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 3, 3, 3, 0, 0, 3, 3, 3, 0, 0, 3, 3, 3, 3, 3, 3, 2, 2, 0,
              0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 3, 3, 0, 0, 0, 3, 3, 0, 0, 0, 3, 3, 0, 3, 3, 3, 3, 2,
              0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 3, 0, 3, 3)
    msub_c <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29")
    msub_r <- c("r10", "r11", "r12", "r13", "r14", "r15", "r16", "E", "G", "H", "I", "J", "K", "L", "M", "N")
    subs <- matrix(msub, nrow=length(msub_r), ncol=length(msub_c), dimnames=list(row=msub_r, col=msub_c), byrow = TRUE)
    # gather flux calculation results to plot Figure 4C
    all_cols2 <- fun_col(length(res_sub_l2$res))
    cols2 <- c()
    data_plot2 <- list()
    opt2 <- c()
    for (i in seq(length(res_sub_l2$res))){
        ncol2 <- grep("v", colnames(res_sub_l2$res[[i]]$sens$par))
        cols2 <- c(cols2, rep(all_cols2[i], length(ncol2)))
        no <- names(opt2)
        opt2 <- c(opt2, res_sub_l2$summary[[i]][ncol2, "opt"])
        names(opt2) <- c(no, names(res_sub_l2$res)[i])
        if (length(ncol2) == 1){
            data_plot2[[names(res_sub_l2$res)[i]]] <- res_sub_l2$res[[i]]$sens$par[,ncol2]
        }else if (length(ncol2)>1){
            data_plot2 <- c(data_plot2, as.list(as.data.frame(res_sub_l2$res[[i]]$sens$par[,ncol2])))
        }
    }
    # colors
    col_breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
    col_col <- c("white", "#8FBBDA", "#EB9394", "#96D096")
    # plot figure 4B
    data_hm <- list(z=t(apply(subs, 2, rev)), x=seq(ncol(subs)+1), y=seq(nrow(subs)+1))
    image(data_hm, axes=F, col=col_col, breaks=col_breaks)#, ylab=rownames(subs))
    axis(side=2, at=seq(nrow(subs))+0.5, labels=rev(rownames(subs)), las=2)
    text(x=seq(ncol(subs))+0.5, y=par("usr")[4]+0.5, labels = colnames(subs), srt = 45, adj = c(0,0), xpd = TRUE, cex=1)
    for (i in seq(nrow(subs))){abline(h=i)}
    for (i in seq(ncol(subs))){abline(v=i)}
    box()
    # plot figure 4C
    par(mar=c(5.1,4.1,4.1,2.1))
    vioplot(data_plot2, names=names(data_plot2), col="#CAB3DE", las=2, ylim=c(0,1.2), yaxt="n", xaxt="n", xlim=c(0, length(data_plot2)), outline=FALSE, axes=FALSE, xaxs="i", at=-1+seq(length(data_plot2)))
    axis(side=2, at=seq(0, 1.2, by=0.2), las=2)
    mtext("estimated flux", side = 2, line=2.5, at=0.6)
    # true value: flux=1.0
    abline(h=1)
dev.off()

# go back to the initial working directory
setwd(wd)


####################################
cat("\nDone.\n")
####################################

