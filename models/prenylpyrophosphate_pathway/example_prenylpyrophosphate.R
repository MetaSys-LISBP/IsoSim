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
This code calculates the flux through the prenylpyrophosphate pathway of 
three S. cerevisiae strains, as detailed in the ScalaFlux publication.

It also generates figure 5.

")
####################################

####################################
cat("   ... initialize R environment ...\n\n")
####################################

# Clean workspace
rm(list= ls())

# Path of the prenyl pyrophosphate pathway model
#setwd("C:/Users/millard/Documents/GIT/IsoSim/IsoSim/models/prenylpyrophosphate_pathway")

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

# number of iterations for Monte Carlo sensitivity analysis
niter <- 4


####################################
cat("\n   ... construct isotopic model of the prenylpyrophosphate pathway ...\n\n")
####################################

rxn <- list(r1 = list("R"=c("IPP", "DMAPP"),        "C"=c(-1, 1),     "E"="v1", "T"=c("A", "A")),
            r2 = list("R"=c("IPP", "DMAPP", "GPP"), "C"=c(-1, -1, 1), "E"="v1", "T"=c("A", "B", "AB")),
            r3 = list("R"=c("GPP", "IPP", "FPP"),   "C"=c(-1, -1, 1), "E"="v1", "T"=c("AB", "C", "ABC")),
            r4 = list("R"=c("FPP", "IPP", "GGPP"),  "C"=c(-1, -1, 1), "E"="v1", "T"=c("ABC", "D", "ABCD")),
            r5 = list("R"=c("GGPP"),                "C"=c(-1),        "E"="v1", "T"=c("ABCD")))

net <- net2mat(rxn)


####################################
cat("\n   ... fit label inputs ...\n\n")
####################################

# measurement times (in min)
times <- c(0,1,2,5,10,15,20,30,45,60,90,120)

# 13C-enrichments of IPP (label input) in strains WT, S023 and S037
data_input <- cbind("wt_IPP_1"=c(0, 0.005794996, 0.04856902, 0.08116519, 0.55878448, 0.77030330, 0.85040355, 0.85153409, 0.91898236, 0.89379241, 0.89714548, 0.92250639),
                    "s23_IPP_1"=c(0, 0.005278256, 0.023006651, 0.123746028, 0.436174371, 0.760468662, 0.858349769, 0.926353215, 0.934527282, 0.946497311, 0.962170220, 0.950815889),
                    "s37_IPP_1"=c(0, 0.00000000, 0.02754242, 0.07194988, 0.52470797, 0.80775681, 0.86597886, 0.91210537, 0.92580177, 0.93348378, 0.94258971, 0.95688487))
                      
# fit labeling dynamics of IPP to define label input
enr_in <- fit_label_input(data_input, t=times, file="res_fit_enr", mc.cores=numCores)


####################################
cat("\n   ... calculate fluxes ...\n\n")
####################################


####################################
cat("\n      ... WT strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (GPP, FPP, GGPP)
enr_GPP <- c(0, 0.00000000, 0.00000000, 0.01420686, 0.17298781, 0.61105497, 0.72140360, 0.93089771, 0.96662547, 0.95718893, 0.89838470, 0.93903080)
enr_FPP <- c(0, 0.005118355, 0.057860687, 0.161687286, 0.430212533, 0.620583961, 0.695261466, 0.848903673, 0.874332750, 0.882266565, 0.898576273, 0.913389271)
enr_GGPP <- c(0, 0.005709117, 0.021023469, 0.039971326, 0.112270902, 0.255258056, 0.352550670, 0.485593768, 0.615767216, 0.722004558, 0.799543659, 0.764059945)
data_exp_wt <- cbind(times, enr_GPP, enr_FPP, enr_GGPP)
colnames(data_exp_wt) <- c("time", "GPP_1-2", "FPP_1-2-3", "GGPP_1-2-3-4")

# define analytical functions of label input
# here we use functions from the best fit
#anFun_fit <- list("IPP_1-M0" = enr_in$anFun[[1]],
#                  "IPP_1-M1" = enr_in$anFun[[2]])
anFun_fit <- list("IPP_1-M0" = "1-(1-((0.907600707560615)/(1+exp(-(-0.413807917938107)*(t-(9.00069668275438))))+(0.113875045170219)))",
                  "IPP_1-M1" = "1-((0.907600707560615)/(1+exp(-(-0.413807917938107)*(t-(9.00069668275438))))+(0.113875045170219))")

# initial parameters
#    metabolite concentrations
meta_conc <- c("IPP"=1, "DMAPP"=1, "GPP"=0.1687, "FPP"=0.804, "GGPP"=1.4196)
#    fluxes
kp <- c("v1"=1)
#    sd of measurements
sd_meas <- list(iso=0.02, conc=c("GPP"=0.01788, "FPP"=0.0606, "GGPP"=0.1292))
#    names of parameters (concentrations & fluxes) to estimate
te <- c("v1", "DMAPP", "GPP", "FPP", "GGPP")
#    upper & lower constraints on parameters
te_upc <- c("v1"=1e4, "DMAPP"=100, "GPP"=0.1687*1.3, "FPP"=0.804*1.3, "GGPP"=1.4196*1.3)
te_loc <- c("v1"=1e-4, "DMAPP"=1e-3, "GPP"=0.1687*0.7, "FPP"=0.804*0.7, "GGPP"=1.4196*0.7)

# run flux calculation
resF_wt <- fit(net        = net,
               times      = times,
               kp         = kp,
               to_est     = te,
               te_upc     = te_upc[te],
               te_loc     = te_loc[te],
               te_upc_det = NULL,
               te_loc_det = NULL,
               eq_det     = NULL,
               data_meas  = list(iso=data_exp_wt, conc=c("GPP"=0.1687, "FPP"=0.804, "GGPP"=1.4196)),
               sd_meas    = sd_meas,
               meta_conc  = meta_conc,
               iso_conc   = NULL,
               mode       = "enrichments",
               anFun      = anFun_fit,
               trf        = NULL,
               events     = NULL,
               p          = 0.0,
               subDir     = "fit_wt",
               nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
               method     = "FORTRAN",
               niter      = niter)

# display estimated parameters
resF_wt$sens$summary


####################################
cat("\n      ... S023 strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (GPP, FPP, GGPP)
enr_FPP_23 <- c(0, 0.002841308, 0.056693817, 0.145121747, 0.389343296, 0.562766707, 0.690064726, 0.791326222, 0.823459887, 0.849350967, 0.853482424, 0.794775448)
enr_GGPP_23 <- c(0, 0.0005139825, 0.0090328223, 0.0376496509, 0.1119685939, 0.2084228755, 0.2964691920, 0.5692766485, 0.7418035441, 0.7947631507, 0.8317840326, 0.7471978321)
enr_GPP_23 <- c(0.00000000, 0.00000000, 0.00000000, 0.03724387, 0.10475834, 0.53408860, 0.74880072, 0.79256680, 0.90912327, 1.00000000, 0.96958405, 0.95721379)
data_exp_23 <- cbind(times, enr_GPP_23, enr_FPP_23, enr_GGPP_23)
colnames(data_exp_23) <- c("time", "GPP_1-2", "FPP_1-2-3", "GGPP_1-2-3-4")

# define analytical functions of label input
# here we use functions from the best fit
#anFun_fit_23 <- list("IPP_1-M0" = enr_in$anFun[[3]],
#                     "IPP_1-M1" = enr_in$anFun[[4]])
anFun_fit_23 <- list("IPP_1-M0" = "1-(1-((0.98838994333855)/(1+exp(-(-0.293832173129716)*(t-(10.1882212485862))))+(0.0589972980325579)))",
                     "IPP_1-M1" = "1-((0.98838994333855)/(1+exp(-(-0.293832173129716)*(t-(10.1882212485862))))+(0.0589972980325579))")

# initial parameters
#    metabolite concentrations
meta_conc_23 <- c("IPP"=1, "DMAPP"=1, "GPP"=0.3258, "FPP"=3.492, "GGPP"=8.54)
#    fluxes
kp_23 <- c("v1"=1)
#    sd of measurements
sd_meas_23 <- list(iso=0.02)
#    names of parameters (concentrations & fluxes) to estimate
te_23 <- c("v1", "DMAPP", "GPP", "FPP", "GGPP")
#    upper & lower constraints on parameters
te_upc_23 <- c("v1"=1e4, "DMAPP"=100, "GPP"=0.3258*1.3, "FPP"=3.492*1.3, "GGPP"=8.54*1.3)
te_loc_23 <- c("v1"=1e-4, "DMAPP"=1e-3, "GPP"=0.3258*0.7, "FPP"=3.492*0.7, "GGPP"=8.54*0.7)

# run flux calculation
resF_23 <- fit(net        = net,
               times      = times,
               kp         = kp_23,
               to_est     = te_23,
               te_upc     = te_upc_23[te_23],
               te_loc     = te_loc_23[te_23],
               te_upc_det = NULL,
               te_loc_det = NULL,
               eq_det     = NULL,
               data_meas  = list(iso=data_exp_23, conc=c("GPP"=0.3258, "FPP"=3.492, "GGPP"=8.54)),
               sd_meas    = sd_meas_23,
               meta_conc  = meta_conc_23,
               iso_conc   = NULL,
               mode       = "enrichments",
               anFun      = anFun_fit_23,
               trf        = NULL,
               events     = NULL,
               p          = 0.0,
               subDir     = "fit_23",
               nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
               method     = "FORTRAN",
               niter      = niter)

# display estimated parameters
resF_23$sens$summary


####################################
cat("\n      ... S037 strain\n\n")
####################################

# 13C-enrichments of metabolic intermediates (GPP, FPP, GGPP)
enr_FPP_37 <- c(0, 0.005118355, 0.057860687, 0.161687286, 0.430212533, 0.620583961, 0.695261466, 0.848903673, 0.874332750, 0.882266565, 0.898576273, 0.913389271)
enr_GGPP_37 <- c(0, 0.0005139825, 0.0090328223, 0.0376496509, 0.1119685939, 0.2084228755, 0.2964691920, 0.5692766485, 0.7418035441, 0.7947631507, 0.8317840326, 0.8471978321)
enr_GPP_37 <- c(0.00000000, 0.0000000, 0.1617682, 0.1500000, 0.1896886, 0.6074216, 0.7330779, 0.9091959, 0.9366337, 0.9459965, 0.9581192, 0.8704253)
data_exp_37 <- cbind(times, enr_GPP_37, enr_FPP_37, enr_GGPP_37)
colnames(data_exp_37) <- c("time", "GPP_1-2", "FPP_1-2-3", "GGPP_1-2-3-4")

# define analytical functions of label input
# here we use functions from the best fit
anFun_fit_37 <- list("IPP_1-M0" = "1-(0.912176349774349+ (-0.999999999850669-(0.0535410795291259)*t) * ((1/(1 + exp((6.79315953738514-t)/(-3.07862487157989)))) -  (1/(1 + exp((625.083069497019-t)/(99.9999999394992))))))",
                     "IPP_1-M1" = "0.912176349774349+ (-0.999999999850669-(0.0535410795291259)*t) * ((1/(1 + exp((6.79315953738514-t)/(-3.07862487157989)))) -  (1/(1 + exp((625.083069497019-t)/(99.9999999394992)))))")

# initial parameters
#    metabolite concentrations
meta_conc_37 <- c("IPP"=1, "DMAPP"=1, "GPP"=0.2124, "FPP"=3.46, "GGPP"=12.835)
#    fluxes
kp_37 <- c("v1"=1)
#    sd of measurements
sd_meas_37 <- list(iso=0.02, conc=c("GPP"=0.0688, "FPP"=0.406, "GGPP"=0.20397))
#    names of parameters (concentrations & fluxes) to estimate
te_37 <- c("v1", "DMAPP", "GPP", "FPP", "GGPP")
#    upper & lower constraints on parameters
te_upc_37 <- c("v1"=1e4, "DMAPP"=100, "GPP"=0.2124*1.3, "FPP"=3.46*1.3, "GGPP"=12.835*1.3)
te_loc_37 <- c("v1"=1e-4, "DMAPP"=1e-3, "GPP"=0.2124*0.7, "FPP"=3.46*0.7, "GGPP"=12.835*0.7)

# run flux calculation
resF_37 <- fit(net        = net,
               times      = times,
               kp         = kp,
               to_est     = te_37,
               te_upc     = te_upc_37[te_37],
               te_loc     = te_loc_37[te_37],
               te_upc_det = NULL,
               te_loc_det = NULL,
               eq_det     = NULL,
               data_meas  = list(iso=data_exp_37, conc=c("GPP"=0.2124, "FPP"=3.46, "GGPP"=12.835)),
               sd_meas    = sd_meas_37,
               meta_conc  = meta_conc_37,
               iso_conc   = NULL,
               mode       = "enrichments",
               anFun      = anFun_fit_37,
               trf        = NULL,
               events     = NULL,
               p          = 0.0,
               subDir     = "fit_37",
               nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
               method     = "FORTRAN",
               niter      = niter)

# display estimated parameters
resF_37$sens$summary


####################################
cat("\n   ... plot results ...\n\n")
####################################


# colors for WT, S037 and S023 datasets
col_f <- c("#EB5726", "#8CAE29", "#5378B0")


pdf(paste("Fig_5B_37.pdf", sep=""), height=9, width=3)
    par(mfrow = c(4,1))
    sd_37 <- list("FPP"=c(0.02,0.02,0.023182386, 0.054210075, 0.021260870, 0.045893420, 0.028671184, 0.02, 0.02, 0.02, 0.025152512, 0.02),
                  "GPP"=c( 0.02, 0.02, 0.195489340, 0.212132034, 0.038874583, 0.078288615, 0.044612372, 0.02, 0.02, 0.02, 0.02, 0.121687604),
                  "GGPP"=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 3.297518e-02, 3.149120e-02, 4.130635e-02, 5.878546e-02, 6.262859e-02),
                  "IPP"=c(0.02, 0.02, 3.187939e-02, 2.450686e-02, 2.620640e-02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02))
    for (i in seq(ncol(data_exp_37[,-1]))){
        met <- strsplit(colnames(resF_37$result$retres$sim)[i], "_")[[1]][1]
        plot(times, data_exp_37[,i+1], type="p", ylim=c(0,1), pch=21, col=col_f[2], bg=col_f[2], main=met, las=1)
        segments(x0=times, y0=data_exp_37[,i+1]-sd_37[[met]], x1=times, y1=data_exp_37[,i+1]+sd_37[[met]]) 
        segments(x0=times-1.5, y0=data_exp_37[,i+1]+sd_37[[met]], x1=times+1.5, y1=data_exp_37[,i+1]+sd_37[[met]]) 
        segments(x0=times-1.5, y0=data_exp_37[,i+1]-sd_37[[met]], x1=times+1.5, y1=data_exp_37[,i+1]-sd_37[[met]]) 
        points(times, data_exp_37[,i+1], pch=21, col=col_f[2], bg=col_f[2])
        lines(times, resF_37$result$retres$sim[,i], col=col_f[2])
    }
    met="IPP"
    plot(times, data_input[,"s37_IPP_1"], type="p", ylim=c(0,1), pch=21, col=col_f[2], bg=col_f[2], main=met, las=1)
    segments(x0=times, y0=data_input[,"s37_IPP_1"]-sd_37[[met]], x1=times, y1=data_input[,"s37_IPP_1"]+sd_37[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"s37_IPP_1"]+sd_37[[met]], x1=times+1.5, y1=data_input[,"s37_IPP_1"]+sd_37[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"s37_IPP_1"]-sd_37[[met]], x1=times+1.5, y1=data_input[,"s37_IPP_1"]-sd_37[[met]]) 
    points(times, data_input[,"s37_IPP_1"], pch=21, col=col_f[2], bg=col_f[2])
    eval(parse(text=paste("fe=function(t){", anFun_fit_37[[2]], "}", sep="")))
    lines(times, fe(times), col=col_f[2])
dev.off()


pdf(paste("Fig_5B_23.pdf", sep=""), height=9, width=3)
    par(mfrow = c(4,1))
    sd_23 <- list("FPP"=c(0.02, 0.02, 0.02, 0.02, 0.02, 0.0250587627, 0.0548989016, 0.0235661496, 0.02, 0.02, 0.02, 0.0686234834),
                  "GPP"=c(0.02, 0.02, 0.02, 0.039713231, 0.134695378, 0.082516224, 0.116757632, 0.02142199, 0.07046334, 0.02201246, 0.07864324, 0.02606934),
                  "GGPP"=c(0.02,0.02,  0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
                  "IPP"=c(0.02, 0.02, 0.02, 0.02, 0.038416312, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02))
    for (i in seq(ncol(data_exp_23[,-1]))){
        met <- strsplit(colnames(resF_37$result$retres$sim)[i], "_")[[1]][1]
        plot(times, data_exp_23[,i+1], type="p", ylim=c(0,1), pch=21, col=col_f[i], bg=col_f[3], main=met, las=1)
        segments(x0=times, y0=data_exp_23[,i+1]-sd_23[[met]], x1=times, y1=data_exp_23[,i+1]+sd_23[[met]]) 
        segments(x0=times-1.5, y0=data_exp_23[,i+1]+sd_23[[met]], x1=times+1.5, y1=data_exp_23[,i+1]+sd_23[[met]]) 
        segments(x0=times-1.5, y0=data_exp_23[,i+1]-sd_23[[met]], x1=times+1.5, y1=data_exp_23[,i+1]-sd_23[[met]]) 
        points(times, data_exp_23[,i+1], pch=21, col=col_f[3], bg=col_f[3])
        lines(times, resF_23$result$retres$sim[,i], col=col_f[3])
    }
    met="IPP"
    plot(times, data_input[,"s23_IPP_1"], type="p", ylim=c(0,1), pch=21, col=col_f[3], bg=col_f[3], main=met, las=1)
    segments(x0=times, y0=data_input[,"s23_IPP_1"]-sd_23[[met]], x1=times, y1=data_input[,"s23_IPP_1"]+sd_23[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"s23_IPP_1"]+sd_23[[met]], x1=times+1.5, y1=data_input[,"s23_IPP_1"]+sd_23[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"s23_IPP_1"]-sd_23[[met]], x1=times+1.5, y1=data_input[,"s23_IPP_1"]-sd_23[[met]]) 
    points(times, data_input[,"s23_IPP_1"], pch=21, col=col_f[3], bg=col_f[3])
    eval(parse(text=paste("fe=function(t){", anFun_fit_23[[2]], "}", sep="")))
    lines(times, fe(times), col=col_f[3])
dev.off()


pdf(paste("Fig_5B_WT.pdf", sep=""), height=9, width=3)
    par(mfrow = c(4,1))
    sd_wt <- list("FPP"=c(0.02, 0.02, 0.078193767, 0.02, 0.02, 0.033831539, 0.050695759, 0.030982886, 0.092933280, 0.02, 0.020724313, 0.030005237),
                  "GPP"=c(0.02, 0.02, 0.02, 0.02009153, 0.08282480, 0.03178313, 0.08713216, 0.02142199, 0.03046334, 0.02201246, 0.07864324, 0.02606934),
                  "GGPP"=c(0.02, 0.02, 0.02, 0.02, 0.02, 3.343441e-02, 9.377963e-02, 1.213237e-01, 2.615769e-01, 1.993309e-02, 1.878597e-02, 1.645524e-01),
                  "IPP"=c(0.02, 0.081953621, 0.053290364, 0.02, 0.097212131, 0.034748075, 0.02, 0.031898181, 0.02, 0.020214841, 0.02, 0.02))
    for (i in seq(ncol(data_exp_wt[,-1]))){
        met <- strsplit(colnames(resF_37$result$retres$sim)[i], "_")[[1]][1]
        plot(times, data_exp_wt[,i+1], type="p", ylim=c(0,1), pch=21, col=col_f[1], bg=col_f[1], main=met, las=1)
        segments(x0=times, y0=data_exp_wt[,i+1]-sd_wt[[met]], x1=times, y1=data_exp_wt[,i+1]+sd_wt[[met]]) 
        segments(x0=times-1.5, y0=data_exp_wt[,i+1]+sd_wt[[met]], x1=times+1.5, y1=data_exp_wt[,i+1]+sd_wt[[met]]) 
        segments(x0=times-1.5, y0=data_exp_wt[,i+1]-sd_wt[[met]], x1=times+1.5, y1=data_exp_wt[,i+1]-sd_wt[[met]]) 
        points(times, data_exp_wt[,i+1], pch=21, col=col_f[1], bg=col_f[1])
        lines(times, resF_wt$result$retres$sim[,i], col=col_f[1])
    }
    met="IPP"
    plot(times, data_input[,"wt_IPP_1"], type="p", ylim=c(0,1), pch=21, col=col_f[1], bg=col_f[1], main=met, las=1)
    segments(x0=times, y0=data_input[,"wt_IPP_1"]-sd_wt[[met]], x1=times, y1=data_input[,"wt_IPP_1"]+sd_wt[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"wt_IPP_1"]+sd_wt[[met]], x1=times+1.5, y1=data_input[,"wt_IPP_1"]+sd_wt[[met]]) 
    segments(x0=times-1.5, y0=data_input[,"wt_IPP_1"]-sd_wt[[met]], x1=times+1.5, y1=data_input[,"wt_IPP_1"]-sd_wt[[met]]) 
    points(times, data_input[,"wt_IPP_1"], pch=21, col=col_f[1], bg=col_f[1])
    eval(parse(text=paste("fe=function(t){", anFun_fit[[2]], "}", sep="")))
    lines(times, fe(times), col=col_f[1])
dev.off()


pdf(paste("Fig_5C.pdf", sep=""), height=3, width=9)
    par(mfrow = c(1,3))
    # GGPP
    metab_fit_GGPP <- matrix(c(1.4196, 12.835, 8.54,
                               resF_wt$result$par["GGPP"], resF_37$result$par["GGPP"], resF_23$result$par["GGPP"]), nrow=2, ncol=3, byrow=TRUE,
                             dimnames=list(row=c("exp", "fit"), col=c("wt", "S037", "S023")))
    sd_GGPP <- c(0.03, 2.04, 0.68)
    barCenters <- barplot(metab_fit_GGPP, beside=TRUE, ylim=c(0,16), ylab="GGPP concentration (nmol/gDCW)", names.arg=colnames(metab_fit_GGPP), las=1, col=c(col_f[1], "white", col_f[2], "white", col_f[3], "white"))
    segments(barCenters[1,], metab_fit_GGPP[1,] - sd_GGPP, barCenters[1,],
             metab_fit_GGPP[1,] + sd_GGPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_GGPP[1,] + sd_GGPP, barCenters[1,]+0.1,
             metab_fit_GGPP[1,] + sd_GGPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_GGPP[1,] - sd_GGPP, barCenters[1,]+0.1,
             metab_fit_GGPP[1,] - sd_GGPP, lwd = 1.5)
    # GPP
    metab_fit_GPP <- matrix(c(0.17, 0.21, 0.33,
                              resF_wt$result$par["GPP"], resF_37$result$par["GPP"], resF_23$result$par["GPP"]), nrow=2, ncol=3, byrow=TRUE,
                            dimnames=list(row=c("exp", "fit"), col=c("wt", "S037", "S023")))
    sd_GPP <- c(0.02, 0.04, 0.07)
    barCenters <- barplot(metab_fit_GPP, beside=TRUE, ylim=c(0,0.4), ylab="GPP concentration (nmol/gDCW)", names.arg=colnames(metab_fit_GPP), las=1, col=c(col_f[1], "white", col_f[2], "white", col_f[3], "white"))
    segments(barCenters[1,], metab_fit_GPP[1,] - sd_GPP, barCenters[1,],
             metab_fit_GPP[1,] + sd_GPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_GPP[1,] + sd_GPP, barCenters[1,]+0.1,
             metab_fit_GPP[1,] + sd_GPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_GPP[1,] - sd_GPP, barCenters[1,]+0.1,
             metab_fit_GPP[1,] - sd_GPP, lwd = 1.5)
    # FPP
    metab_fit_FPP <- matrix(c(0.80, 3.46, 3.49,
                              resF_wt$result$par["FPP"], resF_37$result$par["FPP"], resF_23$result$par["FPP"]), nrow=2, ncol=3, byrow=TRUE,
                            dimnames=list(row=c("exp", "fit"), col=c("wt", "S037", "S023")))
    sd_FPP <- c(0.06, 0.91, 0.59)
    barCenters <- barplot(metab_fit_FPP, beside=TRUE, ylim=c(0,5), ylab="FPP concentration (nmol/gDCW)", names.arg=colnames(metab_fit_FPP), las=1, col=c(col_f[1], "white", col_f[2], "white", col_f[3], "white"))
    segments(barCenters[1,], metab_fit_FPP[1,] - sd_FPP, barCenters[1,],
             metab_fit_FPP[1,] + sd_FPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_FPP[1,] + sd_FPP, barCenters[1,]+0.1,
             metab_fit_FPP[1,] + sd_FPP, lwd = 1.5)
    segments(barCenters[1,]-0.1, metab_fit_FPP[1,] - sd_FPP, barCenters[1,]+0.1,
             metab_fit_FPP[1,] - sd_FPP, lwd = 1.5)
dev.off()


pdf(paste("Fig_5D.pdf", sep=""), height=3, width=3)
    flux_GGPP <- c(resF_wt$result$par["v1"],
                   resF_37$result$par["v1"],
                   resF_23$result$par["v1"],
                   1.332875)
    sd_flux_GGPP <- c(resF_wt$sens$summary["v1", "sd"],
                      resF_37$sens$summary["v1", "sd"],
                      resF_23$sens$summary["v1", "sd"],
                      0.1435506)
    
    barCenters <- barplot(flux_GGPP, beside=TRUE, ylim=c(0,1.5), ylab="GGPP production (nmol/gDCW/min)", names.arg=c("wt", "S037", "S023", "phytoene"), las=1, col=c(col_f[1:3], "grey"))
    segments(barCenters, flux_GGPP - sd_flux_GGPP, barCenters,
             flux_GGPP + sd_flux_GGPP, lwd = 1.5)
    segments(barCenters-0.1, flux_GGPP + sd_flux_GGPP, barCenters+0.1,
             flux_GGPP + sd_flux_GGPP, lwd = 1.5)
    segments(barCenters-0.1, flux_GGPP - sd_flux_GGPP, barCenters+0.1,
             flux_GGPP - sd_flux_GGPP, lwd = 1.5)
dev.off()


pdf(paste("Fig_5E.pdf", sep=""), height=3, width=3)

    turnover_GGPP <- c(resF_wt$result$par["GGPP"]/resF_wt$result$par["v1"],
                       resF_37$result$par["GGPP"]/resF_37$result$par["v1"],
                       resF_23$result$par["GGPP"]/resF_23$result$par["v1"])
    sd_turnover_GGPP <- c(sqrt(resF_wt$sens$summary["v1", "sd"]**2+resF_wt$sens$summary["GGPP", "sd"]**2),
                          sqrt(resF_37$sens$summary["v1", "sd"]**2+resF_37$sens$summary["GGPP", "sd"]**2),
                          sqrt(resF_23$sens$summary["v1", "sd"]**2+resF_23$sens$summary["GGPP", "sd"]**2))
    
    barCenters <- barplot(turnover_GGPP, beside=TRUE, ylim=c(0,16), ylab="GGPP turnover (min-1)", names.arg=c("wt", "S037", "S023"), las=1, col=c(col_f[1:3], "white"))
    segments(barCenters, turnover_GGPP - sd_turnover_GGPP, barCenters,
             turnover_GGPP + sd_turnover_GGPP, lwd = 1.5)
    segments(barCenters-0.1, turnover_GGPP + sd_turnover_GGPP, barCenters+0.1,
             turnover_GGPP + sd_turnover_GGPP, lwd = 1.5)
    segments(barCenters-0.1, turnover_GGPP - sd_turnover_GGPP, barCenters+0.1,
             turnover_GGPP - sd_turnover_GGPP, lwd = 1.5)
dev.off()

# go back to the initial working directory
setwd(wd)

####################################
cat("\nDone.\n")
####################################

