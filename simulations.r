# 2014-24-12 millard@insa-toulouse.fr
#
# This code reproduces the figures 3 to 8 of the folowing paper:
#
#   Impact of kinetic isotope effects in isotope labeling experiments
#   by P. Millard, S. Sokol, J.C. Portais and P. Mendes
#
# NOTE: Calculations are performed by IsoSim. To check that IsoSim works
#       correctly, load it (command 'source("IsoSim.r")') and run 'isosim_test()'.
#       No error should be displayed, and the 'test' folder created in the
#       working directory should contain simulation results.
#       Please refers to IsoSim code in case of problem.
#
# Copyright 2014, INRA, France
# License: GNU General Public License v2 (see license.txt for details)


# reinitialize environment
rm(list=ls(all=TRUE))

# set path
# setwd("D:/Users/millard/Desktop/Old_desktop/Drafts en cours/iso_effect/Sup_data")

cat("##########################################################################\n",
    "###            CONSTRUCT MODELS & SET KIES, CONDITIONS, ETC            ###\n",
    "##########################################################################\n", sep="")

cat("load IsoSim and plot functions...\n")

    source("IsoSim.r")
    source("plot_fun.r")

cat("construct models with and without isotope exchange...\n")

    source("e_coli_13C_xch.r")
    net <- net2mat(rxn)

    source("e_coli_13C_noxch.r")
    net_noxch <- net2mat(rxn)

cat("create 'results' directory...\n")

    mainDir <- getwd()
    subDir  <- "results"
    if (!file.exists(subDir)){
        dir.create(file.path(mainDir, subDir))
    }
    setwd(file.path(mainDir, subDir))

cat("generate fortran code and compile libraries (if needed)...\n")

    # check if the compiled models exist, build them if needed, and load them
    #     without KIEs, with isotope exchange
    check_lib("e_coli_13C_xch_noKIE", sys="i", net=net, kp=kp)
    #     with KIEs and isotope exchange
    check_lib("e_coli_13C_xch", sys="f", net=net, kp=kp)
    #     with KIEs, without isotope exchange
    check_lib("e_coli_13C_noxch", sys="f", net=net_noxch, kp=kp)

cat("set KIEs for 13C isotopes...\n")

    kie <- list("GND"    = c("PGN_1XXXXX" = 0.9905),  # ref: Hermes, Roeske, O'Leary and Clealand, Biochemistry 1982, 21(20):5106-5114
                "G6PDH"  = c("G6P_1XXXXX" = 0.9837),  # ref: Rendina, Hermes and Cleland, Biochemistry 1984, 23(25):6257-6262
                "ALD"    = c("FBP_XX1XXX" = 0.9843),  # ref: Gleixner and Schmidt, J. Biol. Chem. 1997, 272:5382-5387
                "ALDrev" = c("DHAP_1XX"   = 0.9843),
                "PDH"    = c("PYR_100"    = 0.9908,   # ref: Melzer and Schmidt, J. Biol. Chem. 1987, 262:8159-8164
                             "PYR_010"    = 0.9791,
                             "PYR_001"    = 0.9969,
                             "PYR_110"    = 0.9701,
                             "PYR_101"    = 0.9877,
                             "PYR_011"    = 0.9761,
                             "PYR_111"    = 0.9671),
                "RPE"    = c("RB5P_X1XXX" = 0.9931,   # ref: Lee, Vu and Cleland, Biochemistry 1999, 39:4808-4820
                             "RB5P_XX1XX" = 0.9818,
                             "RB5P_XXX1X" = 0.9852,
                             "RB5P_XXXX1" = 0.9980),
                "RPErev" = c("X5P_X1XXX"  = 0.9931,
                             "X5P_XX1XX"  = 0.9818,
                             "X5P_XXX1X"  = 0.9852,
                             "X5P_XXXX1"  = 0.9980))

cat("initialize conditions to simulate (label inputs, with or without KIEs)...\n")

    l_su <- list("NA"          = list("meta_conc"=c("GLC"=0.0563), iso_conc=NULL, kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "NA_noKIE"    = list("meta_conc"=c("GLC"=0.0563), iso_conc=NULL, kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "U"           = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("111111"=1)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "U_noKIE"     = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("111111"=1)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "80_20"       = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("100000"=0.8,"111111"=0.2)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "80_20_noKIE" = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("100000"=0.8,"111111"=0.2)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "50_50"       = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("000000"=0.5,"111111"=0.5)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "50_50_noKIE" = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("000000"=0.5,"111111"=0.5)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "2_6"         = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("011111"=1)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "2_6_noKIE"   = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("011111"=1)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "1_2"         = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("110000"=1)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "1_2_noKIE"   = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("110000"=1)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"),
                 "3_4"         = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("001100"=1)), kie=kie, lib_name="e_coli_13C_xch", sys="f"),
                 "3_4_noKIE"   = list("meta_conc"=c("GLC"=0.0563), iso_conc=list("GLC"=c("001100"=1)), kie=NULL, lib_name="e_coli_13C_xch_noKIE", sys="i"))


cat("##########################################################################\n",
    "###          SIMULATE STEADY-STATES FOR ALL THE CONDITIONS             ###\n",
    "##########################################################################\n", sep="")

    # simulate steady-states with isotope exchange (for all the label inputs, with & without KIEs)

    result <- list()
    for (su in names(l_su)){
        cat("condition : ", su, "\n", sep="")
        result[[su]] <- steady_state(net, kp, lib_mode="n", lib_name=l_su[[su]]$lib_name, sys=l_su[[su]]$sys, kie=l_su[[su]]$kie, meta_conc=l_su[[su]]$meta_conc, iso_conc=l_su[[su]]$iso_conc)
    }

    # simulate steady-states without isotope exchange (only for 80% 1-13C- + 20% U-13C-glucose as label input)

    cat("condition : 80_20_noKIE_noXCH\n", sep="")
    result[["80_20_noKIE_noxch"]] <- steady_state(net_noxch, kp, lib_mode="n", lib_name="e_coli_13C_noxch", sys="f", kie=NULL, meta_conc=l_su[["80_20_noKIE"]]$meta_conc, iso_conc=l_su[["80_20_noKIE"]]$iso_conc)

    cat("condition : 80_20_noXCH\n", sep="")
    result[["80_20_noxch"]] <- steady_state(net_noxch, kp, lib_mode="n", lib_name="e_coli_13C_noxch", sys="f", kie=kie[c("GND","G6PDH","ALD","PDH","RPE")], meta_conc=l_su[["80_20"]]$meta_conc, iso_conc=l_su[["80_20"]]$iso_conc)

    # calculate the impact of KIEs on metabolite concentrations and fluxes
    # at metabolic steady-state for each condition

    diff_meta <- matrix(nrow=0, ncol=length(result[[1]]$metabolites))
    diff_flx  <- matrix(nrow=0, ncol=length(result[[1]]$fluxes))

    for (i in seq(1, length(l_su), 2)){
        rn <- c(rownames(diff_meta), names(result)[i])
        diff_meta <- rbind(diff_meta, 100*(result[[i]]$metabolites/result[[i+1]]$metabolites - 1))
        diff_flx  <- rbind(diff_flx,  100*(result[[i]]$fluxes/result[[i+1]]$fluxes - 1))
        rownames(diff_meta) <- rn
        rownames(diff_flx)  <- rn
    }

    # calculate the impact of KIEs on isotopologue abundances (with and without isotope exchange),
    # isotopomer abundances, and molecular enrichments at metabolic steady-state, for three label inputs

    diff_id <- rbind(result[["NA"]]$isotopologues    - result[["NA_noKIE"]]$isotopologues,
                     result[["80_20"]]$isotopologues - result[["80_20_noKIE"]]$isotopologues,
                     result[["50_50"]]$isotopologues - result[["50_50_noKIE"]]$isotopologues)

    diff_en <- rbind(result[["NA"]]$enrichments    - result[["NA_noKIE"]]$enrichments,
                     result[["80_20"]]$enrichments - result[["80_20_noKIE"]]$enrichments,
                     result[["50_50"]]$enrichments - result[["50_50_noKIE"]]$enrichments)

    diff_ip <- rbind(result[["NA"]]$isotopomers    - result[["NA_noKIE"]]$isotopomers,
                     result[["80_20"]]$isotopomers - result[["80_20_noKIE"]]$isotopomers,
                     result[["50_50"]]$isotopomers - result[["50_50_noKIE"]]$isotopomers)

    diff_id_no_xch <- rbind(result[["80_20_noxch"]]$isotopologues - result[["80_20_noKIE_noxch"]]$isotopologues,
                            result[["80_20"]]$isotopologues       - result[["80_20_noKIE"]]$isotopologues)


cat("##########################################################################\n",
    "###                   GLOBAL SENSITIVITY ANALYSIS                      ###\n",
    "##########################################################################\n",
    "notes: by default only 10 random sets are generated (vs 10,000 in the\n",
    "         paper) to keep computation time low. To reproduce the figures,\n",
    "         please increase the number of sets (variable 'niter' in the code)\n",
    "       errors are sometimes returned by rootSolve if a steady-state cannot\n",
    "         be reached, you can safely ignore these messages\n", sep="")

    # number of sets of random enzyme levels to generate (10,000 sets were generated in the publication)
    niter <- 10

    # initialize parameters
    res_sensitivity <- list("KIE"=list(), "noKIE"=list(), "params"=list())
    kp_tmp <- kp
    lrates <- grep("rmax", names(kp_tmp))
    ss <- 0
	
    # calculate steady-states for niter sets of random Vmaxes
    for (i in seq(niter)){
        cat("iteration ", i, "/", niter, "\n", sep="")
		
        # log uniform sampling of Vmaxes between 0.1 and 10 times their initial value
        kp_tmp[lrates] <- kp[lrates] * 10 ** runif(length(lrates), min=-1, max=1)
		
        # calculate steady-state without KIEs (library generated by R2fortran using 'sys=i' argument => faster)
        sim_nokie <- steady_state(net, kp_tmp, lib_mode="n", norm="PTS", lib_name="e_coli_13C_xch_noKIE", sys="i", meta_conc=c("GLC"=0.0563), iso_conc=list("GLC"=c("100000"=0.8,"111111"=0.2)), stol=1e-8, atol=1e-6, rtol=1e-6, times=c(0, 1e+6), hmin=1e-10)
		
        # if a steady-state is reached
        if (sim_nokie$steady == TRUE){
            # save simulation results
            ss <- ss + 1
            res_sensitivity$params[[ss]] <- kp_tmp
            res_sensitivity$noKIE[[ss]]  <- sim_nokie
            # calculate the corresponding steady-state vith KIEs
            res_sensitivity$KIE[[ss]] <- steady_state(net, kp_tmp, lib_mode="n", norm="PTS", lib_name="e_coli_13C_xch", sys="f", kie=kie, meta_conc=sim_nokie$metabolites, iso_conc=iso2list(sim_nokie$isotopomers), stol=1e-8, atol=1e-6, rtol=1e-6, times=c(0, 1e+6), hmin=1e-10)
        }
        cat("  steady-states reached: ", ss, "/", i, " (", round(ss/i * 100, 1), "%)", "\n", sep="")
    }
    # save and load the results
    # save(res_sensitivity, file="res_sens_10k.RData")
    # load("res_sens_10k.RData")

    # glucose uptake and relative flux through PGI (for all the steady-states)
    vec_PGI_rel <- sapply(res_sensitivity$noKIE, "[[", "fluxes_n")["PGI",]
    vec_PTS     <- sapply(res_sensitivity$noKIE, "[[", "fluxes")["PTS",]

    # errors of each isotopic dataset (for all the steady-states)
    vec_id <- unlist(lapply(res_sensitivity$noKIE, "[[", "isotopologues")) - unlist(lapply(res_sensitivity$KIE, "[[", "isotopologues"))
    vec_ip <- unlist(lapply(res_sensitivity$noKIE, "[[", "isotopomers"))   - unlist(lapply(res_sensitivity$KIE, "[[", "isotopomers"))
    vec_en <- unlist(lapply(res_sensitivity$noKIE, "[[", "enrichments"))   - unlist(lapply(res_sensitivity$KIE, "[[", "enrichments"))



cat("##########################################################################\n",
    "###                        PLOT THE RESULTS                            ###\n",
    "##########################################################################\n", sep="")

    # set color palette
    fun_col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

cat("  figure 3...\n")
    
    # impact of KIEs on metabolite concentrations and fluxes

    pdf(file="Figure_3.pdf", width=8, height=8)

      # set layout
      layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE))

      # figure A: impact on fluxes
      ord <- c("PTS","PGI","PFK","ALD","TPI","GAPDH","PGK","PGM","ENO","PYK","PPC","PDH","G6PDH","GND","RPE","RPI","TA","TK1","TK2", "G1PAT","GPM","MURSYN","GDH","SERSYN","SYN1","SYN2","RPP","DAHPSYN","G6Pdilution","F6Pdilution","FDPdilution","GAPdilution","DHAPdilution","BPGdilution","PG3dilution","PG2dilution","PEPdilution","RB5Pdilution","R5Pdilution", "X5Pdilution","S7Pdilution","PYRdilution","PGNdilution","E4Pdilution","G1Pdilution")
      ndim <- nrow(diff_flx)
      plot(diff_flx[1, ord], col=fun_col(ndim)[1], ylim=c(-2.5, 2.5), xaxt="n", las=2, xlab="reaction", ylab="flux change (%)", main="Figure 3A", cex.axis=0.8, pch=15)
      for (i in seq(2, ndim)){
        par(new=TRUE)
        plot(diff_flx[i, ord], col=fun_col(ndim)[i], ylim=c(-2.5, 2.5), xaxt="n", yaxt="n", las=2, xlab="", ylab="", cex.axis=0.8, pch=15)
      }
      abline(h=0)
      axis(1, at=seq(length(ord)), labels=ord, las=2)

      # figure B: impact on metabolites
      col_plot <- colnames(diff_meta) %ni% c("GLC")
      plot(diff_meta[1, col_plot], col=fun_col(ndim)[1], ylim=c(-1.4, 1), xaxt="n", las=2, xlab="metabolite", ylab="concentration change (%)", main="Figure 3B", cex.axis=0.8, pch=15)
      for (i in seq(2, ndim)){
        par(new=TRUE)
        plot(diff_meta[i, col_plot], col=fun_col(ndim)[i], ylim=c(-1.4, 1), yaxt="n", xaxt="n", las=2, xlab="", ylab="", cex.axis=0.8, pch=15)
      }
      abline(h=0)
      axis(1, at=seq(sum(col_plot)), labels=colnames(diff_meta)[col_plot], las=2)

      # legend
      plot(0, type='n', axes=FALSE, ann=FALSE)
      legend("topleft", legend=rownames(diff_flx), fill=fun_col(ndim))

    dev.off()

cat("  figure 4...\n")

    # impact of KIEs on isotopic data

    pdf(file="Figure_4.pdf", width=7, height=7)

      # figure A-C: impact on isotopologues
      nam <- sapply(strsplit(colnames(diff_id), "_"), "[[", 1) %ni% c("GLC")
      hist_detailed(diff_id[1,nam], ylab="# isotopologues", xlab="error", main="Fig. 4A - nat. ab.", nbin=20, norm=FALSE)
      hist_detailed(diff_id[2,nam], ylab="# isotopologues", xlab="error", main="Fig. 4B - 80% 1-13C + 20% U-13C", nbin=20, norm=FALSE)
      hist_detailed(diff_id[3,nam], ylab="# isotopologues", xlab="error", main="Fig. 4C - 50% 12C + 50% U-13C", nbin=20, plot_leg=TRUE, norm=FALSE)

      # figure D-F: impact on isotopomers
      nam <- sapply(strsplit(colnames(diff_ip), "_"), "[[", 1) %ni% c("GLC")
      hist_detailed(diff_ip[1,nam], ylab="# isotopomers", xlab="error", main="Fig. 4D - nat. ab.", nbin=20, gap=c(19,99), norm=FALSE)
      hist_detailed(diff_ip[2,nam], ylab="# isotopomers", xlab="error", main="Fig. 4E - 80% 1-13C + 20% U-13C", nbin=20, gap=c(19,99), norm=FALSE)
      hist_detailed(diff_ip[3,nam], ylab="# isotopomers", xlab="error", main="Fig. 4F - 50% 12C + 50% U-13C", nbin=20, gap=c(19,99), plot_leg=TRUE, norm=FALSE)

      # figure G-I: impact on enrichments
      nam <- colnames(diff_en) %ni% c("GLC")
      hist_detailed(diff_en[1,nam], ylab="# enrichments", xlab="error", main="Fig. 4G - nat. ab.", nbin=20, norm=FALSE)
      hist_detailed(diff_en[2,nam], ylab="# enrichments", xlab="error", main="Fig. 4H - 80% 1-13C + 20% U-13C", nbin=20, norm=FALSE)
      hist_detailed(diff_en[3,nam], ylab="# enrichments", xlab="error", main="Fig. 4I - 50% 12C + 50% U-13C", nbin=20, plot_leg=TRUE, norm=FALSE)

    dev.off()

cat("  figure 5...\n")

    # distribution of the glucose uptake rate and the glycolytic flux (relatve glc uptake) for the
    # steady-states reached at random enzyme levels

    pdf(file="Figure_5.pdf", width=7, height=5)

      hist_detailed(vec_PTS, ylab="% of steady-states", xlab="Glucose uptake (mM/s)", main="Figure 5A", xlim=c(0,1.6), nbin=16, las=1)
      hist_detailed(vec_PGI_rel, ylab="% of steady-states", xlab="PGI flux (% glucose uptake)", main="Figure 5B", xlim=c(0,1), nbin=20, las=1)

    dev.off()

cat("  figure 6...\n")

    # impact of KIEs on isotopic data for the steady-states calculated from random enzyme levels

    # get metabolite names for isotopologues and isotopomers
    nid <- sapply(strsplit(names(vec_id), "_"), "[[", 1)
    nip <- sapply(strsplit(names(vec_ip), "_"), "[[", 1)

    # plot distributions (exclude glucose since it is fixed)
    pdf(file="Figure_6.pdf", width=7, height=7)

      hist_detailed(vec_id[nid != "GLC"], ylab="% of isotopologues", xlab="error", main="Figure 6A", nbin=20, gap=c(1.1,10), xlim=c(-0.01,0.014))
      hist_detailed(vec_en[names(vec_en) != "GLC"], ylab="% of enrichments", xlab="error", main="Figure 6B", nbin=20, gap=c(1.5,2))
      hist_detailed(vec_ip[nip != "GLC"], ylab="% of isotopomers", xlab="error", main="Figure 6C", nbin=20, gap=c(1.1,10), xlim=c(-0.008,0.009), plot_leg=TRUE)

    dev.off()

cat("  figure 7...\n")

    # impact of KIEs on isotopic data measured on pyruvate

    # extract errors for pyruvate (and change names to get detailed histograms)
    pyr_id <- vec_id[nid == "PYR"]
    names(pyr_id) <- paste(sapply(strsplit(names(pyr_id), "_"), "[[", 2), "PYR", sep="_")
    pyr_ip <- vec_ip[nip == "PYR"]
    names(pyr_ip) <- paste(sapply(strsplit(names(pyr_ip), "_"), "[[", 2), "PYR", sep="_")

    # plot distributions
    pdf(file="Figure_7.pdf", width=7, height=5)

      hist_detailed(pyr_id, ylab="% of pyruvate isotopologues", xlab="error", main="Figure 7A", nbin=20)
      hist_detailed(pyr_ip, ylab="% of pyruvate isotopomers", xlab="error", main="Figure 7B", nbin=20, plot_leg=TRUE)

    dev.off()

cat("  figure 8...\n")

    pdf(file="Figure_8.pdf", width=12, height=4)

      # plot impact without isotope exchange
      plot(diff_id_no_xch[1,-(1:7)], col=c("steelblue"), ylim=c(-0.008,0.008), xaxt = "n", yaxt = "n", las=2, xlab="metabolite/isotopologue", main="Figure 8", ylab="error", cex.axis=0.8, pch=15)

      # plot impact with isotope exchange
      par(new=TRUE)
      plot(diff_id_no_xch[2,-(1:7)], col=c("firebrick"), ylim=c(-0.008,0.008), xaxt = "n", yaxt = "n", las=2, xlab="", ylab="", cex.axis=0.8, pch=15)
      abline(h=0)

      # add legend and axis
      legend("topright", legend=c("with isotope exchange", "without isotope exchange"), fill=c("firebrick", "steelblue"), bty="n")
      axis(1, at=seq(length(colnames(diff_id_no_xch[,-(1:7)]))), labels=colnames(diff_id_no_xch[,-(1:7)]), las=2, cex.axis=0.6)
      axis(2, at=c(seq(-0.008, 0.008, 0.002)), las=1, cex.axis=0.7)

    dev.off()

setwd(mainDir)

cat("Done.\n", sep="")

