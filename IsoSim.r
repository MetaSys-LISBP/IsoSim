# 2014-20-09 millard@insa-toulouse.fr
#
# IsoSim v0.1
#
# Modelling metabolic and isotopic dynamics (kinetic isotope effects can be considered).
#
# Rtools and some R packages (Matrix, rootSolve, deSolve, stringr, and RColorBrewer)
# are required by IsoSim. Rtools can be downloaded online and must be installed manually.
# Other packages can be installed automatically by running the following command in R:
#
#     install.packages(c("Matrix", "rootSolve", "deSolve", "stringr", "RColorBrewer"))
#
# A test function containing the network shown in Figure 1 of the paper and performing
# steady-state and dynamic simulations can be run with the command 'isosim_test()', please
# refer to this function for examples on IsoSim usage.
#
# To check that IsoSim works correctly, run the test function. A folder 'test' containing
# simulation results should be created in the working directory, and no error should be displayed.
#
# Copyright 2014, INRA, France
# License: GNU General Public License v2 (see license.txt for details)


####################
# LOAD R LIBRARIES #
####################

load_pack <- function(packages){
    for (i in packages){
        if (!suppressWarnings(require(i, character.only=TRUE, quietly=TRUE))){
            cat("The package '", i, "' is required by IsoSim. To install this package, please run the following command:\n   install.packages('", i, "')\n", sep="")
            # install.packages(i, dependencies=TRUE, verbose=FALSE)
            # suppressWarnings(require(i, character.only=TRUE, quietly=TRUE))
        }
    }
    return(invisible(NULL))
}

load_pack(c("Matrix", "rootSolve", "deSolve", "stringr", "RColorBrewer"))

"%ni%" <- Negate("%in%")


####################
# ISOSIM FUNCTIONS #
####################

check_lib <- function(lib_n, sys=sys, net=net, kp=kp, force=FALSE){
    # generate the fortran code and compile it if the library is not found,
    # or just load the library it exists already

    lib_n_ext <- paste(lib_n, .Platform$dynlib.ext, sep="")
    # check if the library exists
    if (!file.exists(lib_n_ext) | force){
        # generate the fortran code and compile it
        cat("  library '", lib_n_ext, "' not found, compilation in progress...\n", sep="")
        R2Fortran(net, kp, sys=sys, lib_name=lib_n)
    }else{
        # (re)load the library
        suppressWarnings(try(dyn.unload(lib_n_ext), silent=TRUE))
        dyn.load(lib_n_ext)
    }

    return(invisible(NULL))
}

net2mat <- function(rxn, fixed=NULL, add_eq=NULL){
    # generate stoichiometric and isotopic matrices
    #   'rxn'    input network for which the equation system is built, 'rxn' must be a list
    #            where each key is a reaction name that contains a list of the following objects:
    #      $su     names of the substrates of the reaction.
    #      $pr     names of the products of the reaction (can be empty)
    #      $tr     carbon atom transitions; if none are given, only metabolite dynamics can be simulated.
    #              Atom transitions can be omitted only for metabolites not involved in carbon scrambling.
    #              The carbon transitions of a given metabolite must be consistent through the network.
    #              A parallel network can be decoupled from the carbon network only if it does not produce
    #              metabolites involved in isotope scrambling.
    #      $eq     rate equation; can be a constant or a kinetic rate law (with names of kinetic
    #              parameters matching those of the corresponding vector)
    #   'add_eq' vector containing additional equations to calculate before rate equations
    #            (can be useful to assign the value of some variables, e.g. 'p_1=2*v1' and
    #            a rate equation can refers to the parameter 'p_1'); can also be some fortran code
    #            to include in the subroutine or some R code (depending if ODEs are solved using
    #            fortran or R).
    #   'fixed'  names of metabolites with a fixed concentrations (sinks are automatically identified
    #            and their concentration is fixed, internal metabolites can also be fixed using this
    #            argument)
    # notes:
    #    an example network is provided in isosim_test()
    #    in the current version, stoichiometric coefficients are -1 for substrates and +1 for products
    #
    # return a list containing the following objects:
    #   $s_mat          stoichiometric matrix
    #   $i_mat          isotopic matrix (rows=isotopomers, columns=fluxomers, 0 * 0 if no carbon
    #                   transition are given in rxn)
    #   $i_rows_info    vector of metabolite names (for rows of i_mat)
    #   $i_cols_info    information on each column of i_mat:
    #                      $flx: flux name
    #                      $m_s1, $m_s2, $m_p1 and $m_p2: name(s) of the substrate(s)/product(s)
    #                      $i_s1, $i_s2, $i_p1 and $i_p2: isotopomer(s) of the substrate(s)/product(s)
    #   $m_map          'isotopomers -> metabolites' mapping matrix
    #   $e_map          'isotopomers -> enrichments' mapping matrix
    #   $id_map         'isotopomers -> isotopologues' mapping matrix
    #   $eq             kinetic equations used to calculate fluxes (one per reaction)
    #   $add_eq         additional equations (calculated before $eq, so can be used e.g. to assign
    #                   some variables); can also be fortran or R code
    #   $fixed          sink metabolites ($in and $out for substrates and products, respectively, and
    #                   $user for those constrained by user - names passed in the vector 'fixed' -)

    s_mat_app_m <- c()
    s_mat_app_r <- c()
    s_mat_app_c <- c()
    i_mat_app_i <- c()
    i_mat_app_r <- c()
    i_mat_app_c <- c()
    i_cols_info <- list()
    eq          <- paste(names(rxn), sapply(rxn, "[[", "eq"), sep="=")
    isIso       <- !(identical(unique(as.vector(sapply(rxn, "[[", "tr"))), c("")) | is.null(unlist(unique(as.vector(sapply(rxn, "[[", "tr"))))))
    meta_no_iso <- c()

    # for each reaction
    for (r in names(rxn)){

        # get substrates and products
        subs  <- rxn[[r]]$"su"
        produ <- rxn[[r]]$"pr"

        # append metabolites, reactions and stoichiometric coefficients to construct s_mat
        # (note: in this version, -1 for substrates and +1 for products)
        s_red <- (subs != "")
        p_red <- (produ != "")
        s_mat_app_m <- c(s_mat_app_m, subs[s_red], produ[p_red])
        s_mat_app_r <- c(s_mat_app_r, rep(r, sum(s_red)+sum(p_red)))
        s_mat_app_c <- c(s_mat_app_c, rep(-1, sum(s_red)), rep(1, sum(p_red)))

        # if carbon transitions are given, generate the isotopic vectors to construct i_mat
        if (isIso){

            # atom mapping vector & number of fluxomers
            tot <- unlist(strsplit(paste(rxn[[r]]$"tr"[1:length(subs)], collapse=""), ""))
            n_fluxomers <- 2 ** length(tot)

            # generate isotopomers and stoichiometric coefficients for substrates
            s1 <- add_substrate(nchar(rxn[[r]]$"tr"[1]), subs[1], n_fluxomers)
            s2 <- add_substrate(nchar(rxn[[r]]$"tr"[2]), subs[2], n_fluxomers, nrep=length(s1$isn))

            # create template for mapping atom transitions
            cb <- as.vector(outer(s1$lsn, s2$lsn, paste, sep=""))

            # generate isotopomers and stoichiometric coefficients for products
            p1 <- add_product(rxn[[r]]$"tr"[3], produ[1], cb, tot, n_fluxomers)
            p2 <- add_product(rxn[[r]]$"tr"[4], produ[2], cb, tot, n_fluxomers)

            # append substrate(s) and product(s) isotopomers, stoichiometric coefficients, and fluxomers
            i_mat_app_i <- c(i_mat_app_i, s1$app_i, s2$app_i, p1$app_i, p2$app_i)
            i_mat_app_c <- c(i_mat_app_c, s1$app_c, s2$app_c, p1$app_c, p2$app_c)
            i_mat_app_r <- c(i_mat_app_r, rep(paste(r, as.vector(outer(s1$isn, s2$isn, paste, sep="_")), p1$pnl, p2$pnl, sep="_"), sum(s_red)+sum(p_red)))

            # vector of flux names
            i_cols_info$flx <- c(i_cols_info$flx, rep(r, n_fluxomers))

            # vector of isotopomer names
            i_cols_info$i_s1 <- c(i_cols_info$i_s1, rep(s1$isn, n_fluxomers/length(s1$isn)))
            i_cols_info$i_s2 <- c(i_cols_info$i_s2, rep(s2$isn, each=length(s1$isn)))
            i_cols_info$i_p1 <- c(i_cols_info$i_p1, p1$pnl)
            i_cols_info$i_p2 <- c(i_cols_info$i_p2, p2$pnl)

            # vector of metabolite names
            i_cols_info$m_s1 <- c(i_cols_info$m_s1, rep(subs[1], n_fluxomers))
            i_cols_info$m_s2 <- c(i_cols_info$m_s2, rep(subs[2], n_fluxomers))
            i_cols_info$m_p1 <- c(i_cols_info$m_p1, rep(produ[1], n_fluxomers))
            i_cols_info$m_p2 <- c(i_cols_info$m_p2, rep(produ[2], n_fluxomers))

            # metabolites without atom transitions
            meta_no_iso <- c(meta_no_iso, s1$no_iso, s2$no_iso, p1$no_iso, p2$no_iso)
        }
    }

    # construct the stoichiometric matrix
    metabolites <- unique(s_mat_app_m)
    s_mat <- matrix(0, nrow=length(metabolites), ncol=length(rxn), dimnames=list("meta"=metabolites, "flux"=names(rxn)))
    s_mat[cbind(s_mat_app_m, s_mat_app_r)] <- s_mat_app_c

    # identify sink metabolites
    ext_in  <- metabolites[rowSums(s_mat > 0) == 0]
    ext_out <- metabolites[rowSums(s_mat < 0) == 0]

    # set stoichiometric coefficients of sink and fixed metabolites to 0 (i.e. constant concentrations)
    sinks <- unique(c(ext_in, ext_out, fixed))
    s_mat[sinks, ] <- 0

    if (isIso){

        # generate vectors of isotopomers, fluxomers, metabolites involved (or not) in
        # the isotopic transitions network, etc.
        meta_no_iso <- unique(meta_no_iso)
        isotopomers <- unique(i_mat_app_i[i_mat_app_i %ni% meta_no_iso])
        fluxomers   <- unique(i_mat_app_r)
        meta_iso    <- unlist(sapply(strsplit(isotopomers, "_"), "[[", 1))
        i_rows_info <- c(meta_iso, meta_no_iso)
        var_num     <- length(isotopomers) + length(meta_no_iso)
        var_name    <- c(isotopomers, meta_no_iso)
    
        # construct the fluxomer matrix
        i_mat <- matrix(0, nrow = var_num, ncol = length(fluxomers), dimnames = list("isotopomer" = var_name, "fluxomer" = fluxomers))
        i_mat[cbind(i_mat_app_i, i_mat_app_r)] <- i_mat_app_c
    
        # set stoichiometric coefficients of sink and fixed metabolites
        i_mat[i_rows_info %in% sinks, ] <- 0
    
        # calculate the weight of each isotopomer
        weight_iso <- sapply(sapply(strsplit(isotopomers, "_"), "[[", 2), FUN=function(x) str_count(x, "1"))
    
        # construct the 'isotopomers -> isotopologues' mapping matrix
        id_name <- paste(meta_iso, weight_iso, sep="_M")
        isotopologues <- unique(id_name)
        id_map <- matrix(0, nrow = length(isotopologues), ncol = var_num, dimnames = list("isotopologue" = isotopologues, "isotopomer" = var_name))
        id_map[cbind(id_name, isotopomers)] <- 1
    
        # construct the 'isotopomers -> enrichments' mapping matrix
        rows_e_map <- unique(meta_iso)
        e_map <- matrix(0, nrow = length(rows_e_map), ncol = var_num, dimnames = list("meta" = rows_e_map, "isotopomer" = var_name))
        e_map[cbind(meta_iso, isotopomers)] <- weight_iso/nchar(names(weight_iso))
    
        # construct the 'isotopomers -> metabolites' mapping matrix
        rows_m_map <- unique(i_rows_info)
        m_map <- matrix(0, nrow = length(rows_m_map), ncol = var_num, dimnames = list("meta" = rows_m_map, "isotopomer" = var_name))
        m_map[cbind(c(meta_iso, meta_no_iso), c(isotopomers, meta_no_iso))] <- 1
        
    }else{ 
        i_mat  <- matrix(0, nrow=0, ncol=0)
        m_map  <- matrix(0, nrow=0, ncol=0)
        e_map  <- matrix(0, nrow=0, ncol=0)
        id_map <- matrix(0, nrow=0, ncol=0)
        i_rows_info <- c()
        i_cols_info <- c()
    }

    # return objects
    return(list(i_mat       = Matrix(i_mat, sparse=TRUE),
                s_mat       = Matrix(s_mat, sparse=TRUE),
                i_rows_info = i_rows_info,
                i_cols_info = i_cols_info,
                m_map       = Matrix(m_map, sparse=TRUE),
                e_map       = Matrix(e_map, sparse=TRUE),
                id_map      = Matrix(id_map, sparse=TRUE),
                eq          = eq,
                add_eq      = add_eq,
                fixed       = list("in"   = ext_in,
                                   "out"  = ext_out,
                                   "user" = fixed)))
}

add_substrate <- function(len_meta, s_name, n_fluxomers, nrep=NULL, coeff=-1){
    # generate the vector of isotopomers for the substrate 's_name',
    # and returns information required to construct the isotopic
    # matrices as a list.
    #    'len_meta'       total number of positions that may be labeled
    #    's_name'         name of the metabolite
    #    'n_fluxomers'    number of fluxomers
    #    'nrep'           output order of the isotopic vector (refers
    #                     to the code for details)
    #    'coeff'          stoichiometric coefficient (-1 by default)

    if (len_meta == 0){
        lsn <- c("")
        isn <- s_name
        if (isn == ""){
            no_iso <- c()
            app_i  <- c()
            app_c  <- c()
        }else{
            no_iso <- s_name
            app_i  <- rep(isn, n_fluxomers)
            app_c  <- rep(coeff, n_fluxomers)
        }
    }else{
        no_iso <- c()
        lsn    <- as.vector(apply(format(expand.grid(rep(list(0:1), len_meta))), 1, paste, collapse=""))
        isn    <- paste(s_name, lsn, sep="_")
        if (is.null(nrep)){
            app_i <- rep(isn, n_fluxomers/length(isn))
        }else{
            app_i <- rep(isn, each=nrep)
        }
        app_c <- rep(coeff, n_fluxomers)
    }
    substrate <- list(no_iso=no_iso, lsn=lsn, isn=isn, app_i=app_i, app_c=app_c)

    return(substrate)
}

add_product <- function(tr, p_name, cb, tot, n_fluxomers, coeff=1){
    # generate the vector of isotopomers for the product 'pname',
    # and returns information required to construct the isotopic
    # matrices as a list.
    #    'tr'            atom transitions
    #    'pname'         name of the metabolite
    #    'cb'            template for atom transition mapping
    #    'tot'           merged substrate(s) isotopomers
    #    'n_fluxomers'   number of fluxomers
    #    'coeff'         stoichiometric coefficient (1 by default)

    if (nchar(tr) == 0){
        pnl <- rep(p_name, n_fluxomers)
        if (p_name == ""){
            no_iso <- c()
            app_i  <- c()
            app_c  <- c()
        }else{
            no_iso <- p_name
            app_i  <- pnl
            app_c  <- rep(coeff, n_fluxomers)
        }
    }else{
        no_iso <- c()
        pnl    <- paste(p_name, sapply(cb, FUN=function(x) paste(unlist(strsplit(x,""))[match(unlist(strsplit(tr,"")), tot)], collapse="")), sep="_")
        app_i  <- pnl
        app_c  <- rep(coeff, n_fluxomers)
    }
    product <- list(no_iso=no_iso, pnl=pnl, app_i=app_i, app_c=app_c)

    return(product)
}

R2Fortran <- function(net, parms, sys="f", lib_name="lib_f"){
    # generate the fortran code of the ODE system, compile and load the library
    #   'net'         network generated by net2mat()
    #   'parms'       named vector of kinetic parameters
    #   'lib_name'    name of the fortran library (optional, 'lib_f' by default)
    #   'sys'         equation system to generate:
    #                    'm'  metabolite dynamics only, no isotopic equations
    #                    'i'  all the isotopomers are simulated, enzyme kinetic effects are not considered
    #                    'f'  full system, i.e. all the isotopomers are simulated, and one enzyme
    #                         isotope kinetic factor by substrate(s)-flux pair must be given in 'parms'
    #                    to be implemented: 'e' (emu), 'c' (cumomers)

    # calculate system-dependent parameters (headers, variable names, parameters, etc.)
    rnm  <- rownames(net$s_mat)
    if (sys == "m"){
        lp   <- length(parms)
        nr   <- length(rnm)
        pnm  <- names(parms)
        mc_h <- "C   set METABOLITE concentrations"
        mc_v <- paste(rnm, "=y(", seq(1, nr), ")", sep="")
        d_h  <- "C   calculate derivatives of METABOLITE concentrations"
    }else if (sys == "i"){
        lp   <- length(parms)
        nr   <- length(net$i_rows_info)
        pnm  <- names(parms)
        mc_h <- "C   calculate METABOLITE concentrations by summing ISOTOPOMER concentrations"
        mc_v <- unlist(lapply(rnm, FUN=function(x) paste(x, "=", paste(paste("y(", which(net$i_rows_info == x), ")", sep=""), collapse="+"), sep="")))
        d_h  <- "C   calculate derivatives of ISOTOPOMER concentrations"
    }else if (sys == "f"){
        lk   <- length(net$i_cols_info$flx)
        lp   <- length(parms) + lk
        nr   <- length(net$i_rows_info)
        pnm  <- c(names(parms), paste("kie_", seq(1, lk), sep=""))
        mc_h <- "C   calculate METABOLITE concentrations by summing ISOTOPOMER concentrations"
        mc_v <- unlist(lapply(rnm, FUN=function(x) paste(x, "=", paste(paste("y(", which(net$i_rows_info == x), ")", sep=""), collapse="+"), sep="")))
        d_h  <- "C   calculate derivatives of ISOTOPOMER concentrations"
    }

    # generate equations for metabolite or isotopomer derivatives
    d_eq <- mat2eq(net, sys)

    # regular expression to split lines after 64 characters
    rex <- paste("(?<=", paste(rep(".", 64), collapse=""), ")", sep="")

    # generate the code
    fcode <- c("C---------------------------------------------------------------",
               paste("C file ", lib_name, ".f", sep=""),
               "C Equation system to simulate metabolic and/or isotopic dynamics.",
               "C This code is automatically generated by IsoSim, do not edit.",
               "C---------------------------------------------------------------",
               "",
               "C initialiser for parameter common block",
               "      subroutine initmod(odeparms)",
               "       external odeparms",
               paste("       double precision parms(", lp, ")", sep=""),
               "       common /myparms/parms",
               paste("       call odeparms(", lp, ", parms)", sep=""),
               "       return",
               "      end",
               "",
               "C rate of change and output variables",
               "      subroutine derivs (neq, t, y, ydot, yout, ip)",
               "       implicit none",
               "       double precision t",
               paste("       double precision y(", nr, ")", sep=""),
               paste("       double precision ydot(", nr, ")", sep=""),
               "       double precision yout(*)",
               "       integer neq, ip(*)",
               "C   declare METABOLITES",
               paste("       double precision ", rnm, sep=""))

    if (!is.null(net$add_eq)){
        fcode <- c(fcode, "C   declare VARIABLES",
                   paste("       double precision ", sapply(strsplit(net$add_eq,"="), "[[", 1), sep=""))
    }

    fcode <- c(fcode, "C   declare FLUXES",
               paste("       double precision ", colnames(net$s_mat), sep=""),
               "C   declare KINETIC PARAMETERS",
               paste("       double precision ", pnm, sep=""),
               "C   set values of KINETIC PARAMETERS",
               paste("       common /myparms/ ", paste(pnm, collapse=",\n     * "), sep=""),
               mc_h,
               split_lines(mc_v, rex))

    if (!is.null(net$add_eq)){
        fcode <- c(fcode, "C   calculate VARIABLES",
                   split_lines(net$add_eq, rex))
    }

    fcode <- c(fcode, "C   calculate FLUXES",
               split_lines(net$eq, rex),
               "C   ---------------",
               "C   equation system",
               "C   ---------------",
               d_h,
               split_lines(d_eq, rex),
               "       return",
               "      end",
               "C---------------------------------------------------------------",
               paste("C end of file ", lib_name, ".f", sep=""),
               "C---------------------------------------------------------------")

    # write the code in 'lib_name'.f
    f_file <- file(paste(lib_name, ".f", sep=""))
    writeLines(fcode, f_file)
    close(f_file)

    # compile and load the library (Rtools is required!)
    lib_f <- paste(lib_name, .Platform$dynlib.ext, sep="")
    suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
    system(paste("R CMD SHLIB --preclean --clean ", lib_name, ".f", sep=""))
    dyn.load(lib_f)

    return(invisible(NULL))
}

mat2eq <- function(net, sys="f"){
    # generate (fortran) equations used to calculate derivative of isotopomers or metabolites
    # 'net'   network generated by net2mat()
    # 'sys'   equation system to generate:
    #           'm'  metabolite dynamics only, no isotopic equations
    #           'i'  all the isotopomers are simulated, enzyme KIEs are not considered
    #           'f'  full system, i.e. all the isotopomers are simulated KIEs are considered
    #           to be implemented: 'e' (emu), 'c' (cumomers)

    nf <- unlist(net$fixed)
    d_eq <- c()

    # looping n times, where n is the number of cells different from 0 in the isotopic or
    # stoichiometric matrix (usually n is low since matrices are sparse)
    if (sys == "m"){
        fluxes <- colnames(net$s_mat)
        rowS   <- rownames(net$s_mat)
        for (meta in seq(length(rowS))){
            if (rowS[meta] %ni% nf){
                sub_s_mat <- net$s_mat[meta,]
                f_prod <- paste(fluxes[sub_s_mat > 0], collapse="+")
                f_cons <- paste(fluxes[sub_s_mat < 0], collapse="-")
                l      <- paste("ydot(", meta, ")=", paste(f_prod, "-", f_cons, sep=""), sep="")
                d_eq   <- c(d_eq, l)
            }else{
                d_eq <- c(d_eq, paste("ydot(", meta, ")=0", sep=""))
            }
        }

    }else if (sys == "i"){
        isotopomers <- rownames(net$i_mat)
        for (isot in seq(length(isotopomers))){
            if (net$i_rows_info[isot] %ni% nf){
                sub_i_mat <- net$i_mat[isot,]
                f_prod <- c()
                rowS   <- which(sub_i_mat > 0)
                for (v in unique(net$i_cols_info$flx[rowS])){
                    index  <- rowS[net$i_cols_info$flx[rowS] == v]
                    m1     <- paste("y(", match(net$i_cols_info$i_s1[index], isotopomers), ")", sep="")
                    m2     <- paste("y(", match(net$i_cols_info$i_s2[index], isotopomers), ")", sep="")
                    f_prod <- c(f_prod, paste(paste(v, unique(net$i_cols_info$m_s1[index]), unique(net$i_cols_info$m_s2[index]), sep="/"),
                                              paste("(", paste(paste(m1, m2, sep="*"), collapse="+"), ")", sep=""), sep="*"))
                }
                f_cons <- paste("(", paste(unique(net$i_cols_info$flx[sub_i_mat < 0]), collapse="+"), ")*", paste("y(", isot, ")", sep=""), "/", net$i_rows_info[isot], sep="")
                l      <- paste("ydot(", isot, ")=", gsub("/*", "*", gsub("*y(NA)", "", paste(paste(f_prod, collapse="+"), f_cons, sep="-"), fixed=TRUE), fixed=TRUE), sep="")
                d_eq   <- c(d_eq, l)
            }else{
                d_eq <- c(d_eq, paste("ydot(", isot, ")=0", sep=""))
            }
        }

    }else if (sys == "f"){
        isotopomers <- rownames(net$i_mat)
        for (isot in seq(length(isotopomers))){
            if (net$i_rows_info[isot] %ni% nf){
                f_prod <- c()
                rowS   <- which(net$i_mat[isot,] > 0)
                for (v in unique(net$i_cols_info$flx[rowS])){
                    index  <- rowS[net$i_cols_info$flx[rowS] == v]
                    m1     <- paste("y(", match(net$i_cols_info$i_s1[index], isotopomers), ")", sep="")
                    m2     <- paste("y(", match(net$i_cols_info$i_s2[index], isotopomers), ")", sep="")
                    f_prod <- c(f_prod, paste(paste(v, unique(net$i_cols_info$m_s1[index]), unique(net$i_cols_info$m_s2[index]), sep="/"),
                                              paste("(", paste(paste(paste("kie_", index, sep=""), m1, m2, sep="*"), collapse="+"), ")", sep=""), sep="*"))
                }
                f_cons <- c()
                rowP   <- which(net$i_mat[isot,] < 0)
                for (v in unique(net$i_cols_info$flx[rowP])){
                    index  <- rowP[net$i_cols_info$flx[rowP] == v]
                    m1     <- paste("y(", match(net$i_cols_info$i_s1[index], isotopomers), ")", sep="")
                    m2     <- paste("y(", match(net$i_cols_info$i_s2[index], isotopomers), ")", sep="")
                    f_cons <- c(f_cons, paste(paste(v, unique(net$i_cols_info$m_s1[index]), unique(net$i_cols_info$m_s2[index]), sep="/"),
                                              paste("(", paste(paste(paste("kie_", index, sep=""), m1, m2, sep="*"), collapse="+"), ")", sep=""), sep="*"))
                }
                l    <- paste("ydot(", isot, ")=", gsub("/*", "*", gsub("*y(NA)", "", paste(paste(f_prod, collapse="+"), paste(f_cons, collapse="-"), sep="-"), fixed=TRUE), fixed=TRUE), sep="")
                d_eq <- c(d_eq, l)
            }else{
                d_eq <- c(d_eq, paste("ydot(", isot, ")=0", sep=""))
            }
        }
    }

    return(d_eq)
}

split_lines <- function(v_str, rex){
    # split each line of the vector 'v_str' according to fortran format
    return(unlist(lapply(paste("       ", v_str, sep=""), FUN=function(x) paste(unlist(strsplit(x, rex, perl=T)), collapse="\n     &"))))
}

calc_fluxes <- function(net, kp, conc_m){
    with(as.list(c(net, kp, conc_m)),{
        # calculate fluxes from metabolite concentrations and kinetic parameters

        flx <- colnames(net$s_mat)
        eval(parse(text = paste(net$add_eq, collapse=";\n")))
        eval(parse(text = paste(net$eq, collapse=";\n")))
        eval(parse(text = paste("rates=c(", paste(paste(flx, flx, sep="="), collapse=","), ")", sep="")))

        return(rates)
    })
}

calc_fluxomers <- function(net, kp, iso_conc, iso_eff=NULL){
    # calculate fluxomers from isotopomer concentrations and kinetic parameters

    # calculate metabolite concentrations and fluxes
    conc_m <- (net$m_map %*% iso_conc)[,1]
    fluxes <- calc_fluxes(net, kp, conc_m)

    # abundance of each substrate isotopomer
    ff1 <- iso_conc[net$i_cols_info$i_s1]/conc_m[net$i_cols_info$m_s1]
    ff2 <- iso_conc[net$i_cols_info$i_s2]/conc_m[net$i_cols_info$m_s2]
    ff1[is.na(ff1)] <- 1
    ff2[is.na(ff2)] <- 1

    # calculate fluxomers
    if (is.null(iso_eff)){
        fluxomers <- fluxes[net$i_cols_info$flx] * ff1 * ff2
    }else{
        fluxomers <- fluxes[net$i_cols_info$flx] * ff1 * ff2 * iso_eff
    }
    names(fluxomers) <- colnames(net$i_mat)

    return(fluxomers)
}

dyn_metaR <- function(t, conc_m, params){
    with(as.list(c(conc_m, params)),{
        # calculate time derivatives of metabolite concentrations, this function
        # can be called by rootSolve and deSolve.
        # use preferentially the fortran library to solve ODEs for large systems

        flx_v <- calc_fluxes(params$net, params$kp, conc_m)
        dM_dt <- params$net$s_mat %*% flx_v
        list(dM_dt@x)
    })
}

dyn_isoR <- function(t, iso_conc, params){
    with(as.list(c(iso_conc, params)),{
        # calculate time derivatives of isotopomer concentrations, this function
        # can be called by rootSolve and deSolve.
        # use preferentially the fortran library to solve ODEs for large systems

        flx_v <- calc_fluxomers(params$net, params$kp, iso_conc, iso_eff=params$iso_eff)
        dI_dt <- params$net$i_mat %*% flx_v
        list(dI_dt@x)
    })
}

kie_ini <- function(net, kie){
    # generate the kinetic isotope effects vector from the list 'kie'
    # where each key is a flux name and the corresponding object is a
    # vector c('isotopomer_i' = v_isotopomer_i / v_unlabeled, ...)
    # notes:
    #    cumomers can be given instead of - or in addition to - isotopomers
    #      (e.g. '1X0', where '0' represents the light isotope, '1' denotes
    #      the heavy isotope, and 'X' denotes that the isotope at the
    #      corresponding position can be either a light or heavy isotope)
    #    if a given isotopic specie is represented more than once for a
    #      particular reaction, KIEs are assumed to be cumulative
    #      (i.e. they are multiplied one each others, see Wasylenko and
    #      Stephanopoulos 2013 for details)

    # initialize kinetic isotope effects to 1
    iso_eff <- rep(1, length(net$i_cols_info$flx))
    names(iso_eff) <- colnames(net$i_mat)

    # update 'iso_eff' according to 'kie'
    for (v in names(kie)){
        idf <- which(net$i_cols_info$flx == v)
        for (k in names(kie[[v]])){
            # if cumomer {...} else isotopomer {...}
            if (str_detect(k, "X")){
                idn1 <- net$i_cols_info$i_s1[idf]
                idn2 <- net$i_cols_info$i_s2[idf]
                pos  <- paste("^(.{", which(unlist(strsplit(k, "")) == "X") - 1, "}).", sep="")
                for (i in pos){
                    idn1 <- gsub(i, "\\1X", idn1)
                    idn2 <- gsub(i, "\\1X", idn2)
                }
                ids <- idf[idn1 == k | idn2 == k]
            }else{
                ids <- idf[net$i_cols_info$i_s1[idf] == k | net$i_cols_info$i_s2[idf] == k]
            }
            iso_eff[ids] <- iso_eff[ids] * kie[[v]][k]
        }
    }

    return(iso_eff)
}

iso2list <- function(iso_conc){
    # build a list from the vector 'iso_conc', where each key
    # is a metabolite and the corresponding object is a named
    # vector of isotopologue abundances

    nmi  <- names(iso_conc)
    mn   <- sapply(strsplit(nmi[str_detect(nmi, "_")], "_"), "[[", 1)
    isol <- list()
    for (i in unique(mn)){
        isol[[i]][sapply(strsplit(nmi[mn==i], "_"), "[[", 2)] <- iso_conc[mn==i]
    }

    return(isol)
}

iso_ini <- function(net, p=0.0107, meta_conc=NULL, iso_conc=NULL){
    # generate isotopomer concentrations vector with a given isotopic enrichment
    #   'p'           isotopic enrichment
    #   'meta_conc'   metabolite concentrations (optional, set at 1 by default)
    #   'iso_conc'    isotopomer concentrations (optional, normalized to 1 by
    #                 default, or to the initial concentration of the
    #                 corresponding metabolite when it is given in 'meta_conc')

    # initialize the isotopomers vector at 0
    cnm <- rownames(net$i_mat)
    isotopomers <- rep(0, length(cnm))
    names(isotopomers) <- cnm

    # calculate the (absolute) isotopomer concentrations for each metabolite
    rnm <- rownames(net$e_map)
    rnm_iso <- names(iso_conc)
    rnm_meta <- names(meta_conc)
    for (i in unique(net$i_rows_info)){
        if (i %in% rnm){
            if (i %in% rnm_iso){
                if (i %in% rnm_meta){
                    isotopomers[paste(i, names(iso_conc[[i]]), sep="_")] <- iso_conc[[i]]/sum(iso_conc[[i]]) * meta_conc[i]
                }else{
                    isotopomers[paste(i, names(iso_conc[[i]]), sep="_")] <- iso_conc[[i]]/sum(iso_conc[[i]])
                }
            }else{
                ip <- sapply(strsplit(cnm[net$i_rows_info == i], "_"), "[[", 2)
                weight <- str_count(ip, "1")
                if (i %in% rnm_meta){
                    isotopomers[paste(i, ip, sep="_")] <- (p**weight * (1-p)**(max(weight)-weight)) * meta_conc[i]
                }else{
                    isotopomers[paste(i, ip, sep="_")] <-  p**weight * (1-p)**(max(weight)-weight)
                }
            }
        }else{
            if (i %in% rnm_meta){
                isotopomers[i] <- meta_conc[i]
            }else{
                isotopomers[i] <- 1
            }
        }
    }

    return(isotopomers)
}

xch2net <- function(fxch){
    # calculate net fluxes from forward and backward fluxes,
    # names of backward fluxes must end by 'rev'

    nm   <- names(fxch)
    nmb  <- nm[str_detect(nm, "rev")]
    nmf  <- gsub("rev", "", nmb)
    fnet <- c(fxch[nm %ni% c(nmf, nmb)], fxch[nmf]-fxch[nmb])

    return(fnet)
}

steady_state <- function(net, kp, sys="f", lib_mode="b", lib_name="lib_f", meta_conc=NULL, iso_conc=NULL, norm=NULL, p=0.0107, kie=NULL, stol=1e-9, atol=1e-6, rtol=1e-6, times=c(0, 1e+9), hmin=0, hmax = NULL){
    # simulate the metabolic and isotopic steady-state of the system, with the following arguments:
    #   'net'           network generated by net2mat()
    #   'kp'            named vector of kinetic parameters
    #   'sys'           equation system to generate:
    #                    'm'  metabolite dynamics only, no isotopic equations
    #                    'i'  all the isotopomers are simulated, kinetic isotope effects are not considered
    #                    'f'  full system, i.e. all the isotopomers are simulated and kinetic isotope effects
    #                         are considered
    #   'lib_mode'      options to build and load a compiled version of the model
    #                    'b'  generate the fortran code, compile and load the library
    #                    'c'  check if the library exists, build it if needed, and load it
    #                    'l'  just load the library
    #                    'n'  do nothing, the library must be already loaded
    #   'lib_name'      name of the fortran library to compile and/or load ('lib_f' by default)
    #   'meta_conc'     named vector of initial concentrations of metabolites (optional, by default set at 1)
    #   'iso_conc'      named vector of initial concentrations of isotopomers (optional, by default set at
    #                   natural abundance, with the total concentration of the corresponding metabolite at 1)
    #   'p'             initial enrichment of metabolites (0.0107 by default - natural abundance of 13C -)
    #   'kie'           list defining the kinetic isotope effects, please refers to kie_ini() for details on
    #                   its structure
    #   'norm'          if not NULL, also return the flux distribution normalized to the corresponding flux
    #   'stol'          steady-state tolerance (1.e-9 by default); it is assumed that steady-state is reached
    #                   if the average of absolute values of the derivatives drops below this number
    #   'rtol'          relative error tolerance of integrator (1.e-6 by default); refers to the documentation of
    #                   runsteady() for details
    #   'atol'          absolute error tolerance of integrator (1.e-6 by default); refers to the documentation of
    #                   runsteady() for details
    #   'times'         2-valued vector containing the initial time and the end time of simulations; the last
    #                   time value should be large enough to make sure that steady-state is effectively
    #                   reached in this period (by default c(0, 1e+9))
    #
    # return a list containing the following objects:
    #   $metabolites    steady-state metabolite concentrations
    #   $isotopomers    steady-state isotopomer abundances (normalized to 1)
    #   $isotopologues  steady-state isotopologue abundances (normalized to 1)
    #   $enrichments    steady-state metabolite enrichments (normalized to 1)
    #   $fluxes_xch     steady-state exchange fluxes (i.e. contains forward and reverse fluxes)
    #   $fluxes         steady-state net fluxes (i.e. only forward-backward)
    #   $fluxomers      steady-state fluxomers
    #   $fluxes_xch_n   steady-state exchange fluxes normalized to the flux 'norm'
    #   $fluxes_n       steady-state net fluxes normalized to the flux 'norm'
    #   $steady         TRUE if a steady-state is reached, FALSE if no steady-state is reached
    #   $run            calculation time
    #
    # notes: arguments 'iso_conc' and 'p' are used only if 'sys' = 'i' or 'f',
    #        argument 'kie' is used only if 'sys' = 'f'

    # filename of the compiled model
    lib_f <- paste(lib_name, .Platform$dynlib.ext, sep="")
    
    # compilation and loading options
    if (lib_mode == "b"){
        cat("  build and load '", lib_f, "'...\n", sep="")
        R2Fortran(net, kp, sys=sys, lib_name=lib_name)
    }else if (lib_mode == "l"){
        cat("  load '", lib_f, "'...\n", sep="")
        suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
        dyn.load(lib_f)
    }else if (lib_mode == "c"){
        cat("  build '", lib_f, "' (only if it does not exit) and load it...\n", sep="")
        check_lib(lib_name, sys=sys, net=net, kp=kp)
    }else if (lib_mode == "n"){
        cat("  '", lib_f, "' is assumed to be already loaded...\n", sep="")
    }else{
        stop("The argument 'lib_name' can be only 'b', 'c', 'l' or 'n'.")
    }

    cat("  set initial conditions...\n")
    if (sys == "m"){
        yini                   <- rep(1, nrow(net$s_mat))
        names(yini)            <- rownames(net$s_mat)
        yini[names(meta_conc)] <- meta_conc
        parms                  <- kp
        msg                    <- "metabolic steady-state"
        iso_eff                <- NULL
    }else if (sys == "i"){
        yini    <- iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc)
        iso_eff <- NULL
        parms   <- kp
        msg     <- "metabolic and isotopic steady-states"
    }else if (sys == "f"){
        yini    <- iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc)
        iso_eff <- kie_ini(net, kie)
        parms   <- c(kp, iso_eff)
        msg     <- "metabolic and isotopic steady-states (with KIEs)"
    }else{
        stop("The argument 'sys' can be only 'm', 'i', or 'f'.")
    }

    cat("  calculate ", msg, "...\n", sep="")
    sim  <- c()
    ptmi <- proc.time()
    withRestarts(
        tryCatch(
            sim <- runsteady(y        = yini,
                             func     = "derivs",
                             times    = times,
                             parms    = parms,
                             stol     = stol,
                             atol     = atol,
                             rtol     = rtol,
                             dllname  = lib_name,
                             initfunc = "initmod",
                             nout     = 0,
                             hmin     = hmin,
                             hmax     = hmax),
            finally = cat("")),
        abort = function(){})
    ptmt <- proc.time() - ptmi

    cat("  process and return the results...\n")
    res <- list('run' = ptmt["elapsed"], 'iso_eff'=iso_eff)
    
    # if a steady-state has been found
    if (!is.null(head(sim))){
        if (attr(sim, "steady")){
            
            res$steady <- TRUE
            
            # calculate steady-state variables (metabolite concentrations, fluxes, etc)
            if (sys == "m"){
                res$metabolites   <- sim$y
                res$fluxes_xch    <- calc_fluxes(net, kp, sim$y)
                res$fluxes        <- xch2net(res$fluxes_xch)
            }else{
                m_iso <- rownames(net$e_map)
                res$metabolites   <- (net$m_map %*% sim$y)[,1]
                res$isotopomers   <- sim$y[net$i_rows_info %in% m_iso]/res$metabolites[net$i_rows_info[net$i_rows_info %in% m_iso]]
                res$isotopologues <- (net$id_map %*% sim$y)[,1]/res$metabolites[sapply(strsplit(rownames(net$id_map), "_"), "[[", 1)]
                res$enrichments   <- (net$e_map %*% sim$y)[,1]/res$metabolites[m_iso]
                res$fluxomers     <- calc_fluxomers(net, kp, sim$y, iso_eff)
                res$fluxes_xch    <- sapply(colnames(net$s_mat), FUN=function(x) sum(res$fluxomers[net$i_cols_info$flx == x]))
                res$fluxes        <- xch2net(res$fluxes_xch)
            }
            
            # normalize fluxes
            if (!is.null(norm)){
                if (norm %in% names(res$fluxes)){
                    res$fluxes_xch_n <- res$fluxes_xch/res$fluxes_xch[norm]
                    res$fluxes_n     <- res$fluxes/res$fluxes_xch[norm]
                }else{
                    warning("The argument norm='", norm, "' is not a reaction, fluxes cannot be normalized.")
                    res$fluxes_xch_n <- NULL
                    res$fluxes_n     <- NULL
                }
            }
        }else{
            warning("No steady-state found.")
            res$steady <- FALSE
        }
    }else{
        warning("No steady-state found.")
        res$steady <- FALSE
    }
    return(res)
}

dynamic <- function(net, kp, sys, times, lib_mode="n", lib_name="lib_f", meta_conc=NULL, iso_conc=NULL, p=0.0107, norm=NULL, kie=NULL, events=NULL, stol=1e-9, atol=1e-6, rtol=1e-6){
    # simulate the metabolic and isotopic dynamics of the system, with the following parameters:
    #   'net'           network generated by net2mat()
    #   'kp'            named vector of kinetic parameters
    #   'sys'           equation system:
    #                      'm'  metabolite dynamics only, no isotopic equations
    #                      'i'  all the isotopomers are simulated, kinetic isotope effects are not considered
    #                      'f'  full system, i.e. all the isotopomers are simulated and kinetic isotope effects
    #                           are considered
    #   'times'         vector containing the times at which explicit estimates for system variables are desired
    #   'lib_mode'      options to build and load a compiled version of the model
    #                      'b'  generate the fortran code, compile and load the library
    #                      'c'  check if the library exists, build it if needed, and load it
    #                      'l'  just load the library
    #                      'n'  do nothing, the library must be already loaded
    #   'lib_name'      name of the fortran library to compile and/or load ('lib_f' by default)
    #   'meta_conc'     named vector of initial concentrations of metabolites (optional, by default set at 1)
    #   'iso_conc'      named vector of initial concentrations of isotopomers (optional, by default set at
    #                   natural abundance, with the total concentration of the corresponding metabolite at 1)
    #   'p'             initial enrichment of metabolites (0.0107 by default - natural abundance of 13C -)
    #   'kie'           list defining the kinetic isotope effects, please refers to kie_ini() for details on
    #                   its structure
    #   'norm'          if not NULL, also return the flux distribution normalized to the flux named 'norm'
    #   'events'        events, please refer to the documentation of the 'deSolve' package for details
    #   'stol'          steady-state tolerance (1.e-9 by default); it is assumed that steady-state is reached
    #                   if the average of absolute values of the derivatives drops below this number
    #   'rtol'          relative error tolerance of integrator (1.e-6 by default); refers to the documentation of
    #                   runsteady() for details
    #   'atol'          absolute error tolerance of integrator (1.e-6 by default); refers to the documentation of
    #                   runsteady() for details
    #
    # return a list containing the following objects (for matrices: columns=variables, rows=times):
    #   $times          vector of times at which system variables are simulated
    #   $metabolites    matrix of metabolite concentrations
    #   $isotopomers    matrix of isotopomer abundances (normalized to 1)
    #   $isotopologues  matrix of isotopologue abundances (normalized to 1)
    #   $enrichments    matrix of metabolite enrichments (normalized to 1)
    #   $fluxes_xch     matrix of exchange fluxes (i.e. includes forward and reverse fluxes)
    #   $fluxes         matrix of net fluxes (i.e. only forward-backward)
    #   $fluxomers      matrix of fluxomers
    #   $fluxes_xch_n   matrix of exchange fluxes normalized to the flux 'norm'
    #   $fluxes_n       matrix of net fluxes normalized to the flux 'norm'
    #   $run            calculation time
    #
    # notes: arguments 'iso_conc' and 'p' are used only if 'sys' = 'i' or 'f',
    #        argument 'kie' is used only if 'sys' = 'f'

    # filename of the compiled model
    lib_f <- paste(lib_name, .Platform$dynlib.ext, sep="")

    # compilation and loading options
    if (lib_mode == "b"){
	   cat("  build and load '", lib_f, "'...\n", sep="")
	   R2Fortran(net, kp, sys=sys, lib_name=lib_name)
    }else if (lib_mode == "l"){
	   cat("  load '", lib_f, "'...\n", sep="")
	   suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
	   dyn.load(lib_f)
    }else if (lib_mode == "c"){
	   cat("  build '", lib_f, "' (only if it does not exit) and load it...\n", sep="")
	   check_lib(lib_name, sys=sys, net=net, kp=kp)
    }else if (lib_mode == "n"){
	   cat("  '", lib_f, "' is assumed to be already loaded...\n", sep="")
    }else{
	   stop("The argument 'lib_name' can be only 'b', 'c', 'l' or 'n'.")
    }
	
    cat("  set initial conditions...\n")
    if (sys == "m"){
        yini                   <- rep(1, nrow(net$s_mat))
        names(yini)            <- rownames(net$s_mat)
        yini[names(meta_conc)] <- meta_conc
        parms                  <- kp
        iso_eff                <- NULL
        msg                    <- "metabolic dynamic"
    }else if (sys == "i"){
        yini    <- iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc)
        iso_eff <- NULL
        parms   <- kp
        msg     <- "metabolic and isotopic dynamics"
    }else if (sys == "f"){
        yini    <- iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc)
        iso_eff <- kie_ini(net, kie)
        parms   <- c(kp, iso_eff)
        msg     <- "metabolic and isotopic dynamics (with KIEs)"
    }else{
        stop("The argument 'sys' can be only 'm', 'i', or 'f'.")
    }


    cat("  simulate ", msg, "...\n", sep="")
    ptmi <- proc.time()
    sim <- lsoda(y        = yini,
                 func     = "derivs",
                 parms    = parms,
                 dllname  = lib_name,
                 initfunc = "initmod",
                 events   = list(data=events),
                 times    = times,
                 stol     = stol,
                 atol     = atol,
                 rtol     = rtol)
    ptmt <- proc.time() - ptmi

    cat("  process and return the results...\n")
    
    res <- list('times' = times, 'run' = ptmt["elapsed"], 'iso_eff'=iso_eff)
    
    # calculate time-course variables (metabolite concentrations, fluxes, etc)
    if (sys == "m"){
        res$metabolites   <- sim[, colnames(sim) != "time"]
        res$fluxes_xch    <- t(apply(res$metabolites, 1, FUN=function(x) calc_fluxes(net, kp, x)))
        res$fluxes        <- t(apply(res$fluxes_xch, 1, FUN=function(x) xch2net(x)))
    }else{
        m_iso             <- rownames(net$e_map)
        sim_red           <- sim[, colnames(sim) != "time"]
        res$metabolites   <- t(apply(sim_red, 1, FUN=function(x) (net$m_map %*% x)[,1]))
        res$isotopomers   <- sim_red[,net$i_rows_info %in% m_iso]/res$metabolites[, net$i_rows_info[net$i_rows_info %in% m_iso]]
        res$isotopologues <- t(apply(sim_red, 1, FUN=function(x) (net$id_map %*% x)[,1])/apply(res$metabolites, 1, FUN=function(x) x[sapply(strsplit(rownames(net$id_map), "_"), "[[", 1)]))
        res$enrichments   <- t(apply(sim_red, 1, FUN=function(x) (net$e_map %*% x)[,1]))/res$metabolites[, m_iso]
        res$fluxomers     <- t(apply(sim_red, 1, FUN=function(x) calc_fluxomers(net, kp, x, iso_eff)))
        res$fluxes_xch    <- t(apply(res$fluxomers, 1, FUN=function(x) sapply(colnames(net$s_mat), FUN=function(y) sum(x[net$i_cols_info$flx == y]))))
        res$fluxes        <- t(apply(res$fluxes_xch, 1, FUN=function(x) xch2net(x)))
    }
    
    # normalize fluxes
    if (!is.null(norm)){
        if (norm %in% colnames(res$fluxes)){
            res$fluxes_xch_n <- t(apply(res$fluxes_xch, 1, FUN=function(x) x/x[norm]))
            res$fluxes_n     <- t(apply(res$fluxes, 1, FUN=function(x) x/x[norm]))
        }else{
            warning("The argument norm='", norm, "' is not a reaction, fluxes cannot be normalized.")
            res$fluxes_xch_n <- NULL
            res$fluxes_n     <- NULL
        }
    }

    return(res)
}

list2file <- function(l, file = paste(deparse(substitute(l)), ".txt", sep = "")) {
   # save the list 'l' in a file

    tmp <- getOption("width")
    options(width=10000)
    sink(file)
    print(l)
    sink()
    options(width=tmp)
    return(invisible(NULL))
}



#####################
### TEST FUNCTION ###
#####################

isosim_test <- function(){
    # Build the equation system for the toy network given in Figure 1
    # of the paper, simulate i) metabolic and isotopic steady-states
    # and ii) metabolic and isotopic dynamics, save the results in
    # text files ('model_example_ss.txt' and 'model_example_dyn.txt')
    # and plot the results in pdf ('model_example_ss.pdf' and
    # 'model_example_dyn.pdf').
    # Return a list containing the network, matrices and simulation results
    # (also saved in 'model_example_res.txt').

    # create 'test' directory
    mainDir <- getwd()
    subDir  <- "test"
    if (!file.exists(subDir)){
        dir.create(file.path(mainDir, subDir))
    }
    setwd(file.path(mainDir, subDir))
    
    # function used to create color palettes
    # usage: fun_col(X), where X is the number of colors to generate
    fun_col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

    # start timer
    ptmi <- proc.time()

    ###########################
    ### NETWORK DEFINITION ####

    cat("load the example network and construct matrices...\n")

    #            reaction      substrate(s)     product(s)       atom transition(s)         rate law
    rxn <- list("v1"=    list("su"=c("S",""),  "pr"=c("A",""),  "tr"=c("abc","","abc",""), "eq"="e1*S/(S+p1)"),
                "v2"=    list("su"=c("A",""),  "pr"=c("B",""),  "tr"=c("abc","","acb",""), "eq"="e2*A/(A+p2)"),
                "v3"=    list("su"=c("B",""),  "pr"=c("C","D"), "tr"=c("abc","","a","bc"), "eq"="e3*B"),
                "v3rev"= list("su"=c("C","D"), "pr"=c("B",""),  "tr"=c("a","bc","abc",""), "eq"="e3*C*D/p3"),
                "v4"=    list("su"=c("B",""),  "pr"=c("E",""),  "tr"=c("ghi","","ghi",""), "eq"="e4*B/(p4+B)"),
                "v5"=    list("su"=c("C",""),  "pr"=c("",""),   "tr"=c("a","","",""),      "eq"="e5*C/(p5+C)"),
                "v6"=    list("su"=c("D",""),  "pr"=c("",""),   "tr"=c("ab","","",""),     "eq"="e6*D/(p6+D)"),
                "v7"=    list("su"=c("E",""),  "pr"=c("",""),   "tr"=c("abc","","",""),    "eq"="e7*E/(p7+E)"))

    # kinetic parameters
    kp <- c(e1=1, e2=1, e3=1, e4=1, e5=3, e6=3, e7=3, p1=10, p2=10, p3=10, p4=10, p5=1, p6=1, p7=1)

    # kinetic isotope effects
    kie <- list("v1" = c("S_1XX" = 0.9),
                "v2" = c("A_100" = 0.8,
                         "A_101" = 0.8,
                         "A_110" = 0.8,
                         "A_111" = 0.7))

    # construct matrices
    net <- net2mat(rxn)

    ############################################
    ### ISOTOPIC AND METABOLIC STEADY-STATES ###

    cat("steady-state simulations...\n")

    # calculate steady-state metabolite concentrations, isotopomers, isotopologues, enrichments, fluxes, etc.
    # when the concentration of the substrate S is fixed at 0.5 mM, and its labeling is 80% 1-13C / 20% 3-13C
    result_ss <- steady_state(net       = net,
                              kp        = kp,
                              lib_mode  = "b",
                              lib_name  = "model_example",
                              sys       = "f",
                              meta_conc = c("S" = 0.5),
                              iso_conc  = list("S" = c("100" = 0.8, "001" = 0.2)),
                              kie       = kie,
                              norm      = "v1")

    cat("save and plot the data (in 'model_example_ss.txt' and 'model_example_ss.pdf')...\n")
    list2file(result_ss, file="model_example_ss.txt")
    pdf(file="model_example_ss.pdf", width=10, height=7)
      par(mar=c(5.1, 5.1, 4.1, 1.8))
      layout(matrix(seq(1,6), 2, 3, byrow = TRUE))
      barplot(result_ss$metabolites, col=fun_col(length(result_ss$metabolites)), xlab="metabolite", ylab="concentration (mM)", log="y", ylim=c(1e-3, 1))
      barplot(result_ss$fluxes, col=fun_col(length(result_ss$fluxes)), xlab="reaction", ylab="flux (mM/s)", las=1)
      barplot(result_ss$fluxes_n, col=fun_col(length(result_ss$fluxes_n)), xlab="reaction", ylab="relative flux (% uptake)", las=1)
      barplot(result_ss$isotopomers, col=fun_col(length(result_ss$isotopomers)), xlab="isotopomer", ylab="abundance", las=2, cex.names=0.6)
      barplot(result_ss$isotopologues, col=fun_col(length(result_ss$isotopologues)), xlab="isotopologue", ylab="abundance", las=2, cex.names=0.85)
      barplot(result_ss$enrichments, col=fun_col(length(result_ss$enrichments)), xlab="metabolite", ylab="enrichment", las=1)
    dev.off()

    #######################################
    ### ISOTOPIC AND METABOLIC DYNAMICS ###

    cat("time-course simulations...\n")

    # initial condition: steady-state (result_ss)
    # event:             at t=50, the concentration of S is increased from 0.5 to 2 mM
    #                    and S is switched to 50% 1,2-13C / 30% 1-13C / 20% U-13C
    events <- data.frame(var    = c("S_001", "S_100", "S_111", "S_110"),
                         time   = c(50, 50, 50, 50),
                         value  = c(0, 0.6, 0.4, 1.0),
                         method = c("rep", "rep", "rep", "rep"))

    result_dyn <- dynamic(net       = net,
                          kp        = kp,
                          times     = seq(0, 200, 1),
                          lib_mode  = "n",
                          lib_name  = "model_example",
                          meta_conc = result_ss$metabolites,
                          iso_conc  = iso2list(result_ss$isotopomers),
                          sys       = "f",
                          kie       = kie,
                          events    = events,
                          norm      = "v1")

    cat("save and plot the data (in 'model_example_dyn.txt' and 'model_example_dyn.pdf')...\n")
    list2file(result_dyn, file="model_example_dyn.txt")
    plot_iso <- function(x, y, xlab, ylab, type="l", lty=1, lwd=2, lpos="topright", inset=c(-0.3, 0), pch=20, cex=0.7, bty="n", ncolumns=1, las=1){
        yn <- colnames(y)
        colpal <- fun_col(length(yn))
        matplot(x, y, type=type, lty=lty, lwd=lwd, col=colpal, xlab=xlab, ylab=ylab, las=las)
        legend(lpos, inset=inset, legend=yn, col=colpal, pch=pch, cex=cex, bty=bty, ncol=ncolumns)
    }
    pdf(file="model_example_dyn.pdf", width=10, height=7)
      par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
      layout(matrix(seq(1, 6), 2, 3, byrow = TRUE))
      plot_iso(result_dyn$times, result_dyn$metabolites, "time (s)", "metabolite concentration (mM)")
      plot_iso(result_dyn$times, result_dyn$fluxes, "time (s)", "flux (mM/s)")
      plot_iso(result_dyn$times, result_dyn$fluxes_n, "time (s)", "relative flux (% uptake)")
      plot_iso(result_dyn$times, result_dyn$isotopomers, "time (s)", "isotopomer abundance", inset=c(-0.5, 0), ncolumns=2)
      plot_iso(result_dyn$times, result_dyn$isotopologues, "time (s)", "isotopologue abundance")
      plot_iso(result_dyn$times, result_dyn$enrichments, "time (s)", "metabolite enrichment")
    dev.off()

    cat("save all the R objects in 'model_example_res.txt')...\n")
    res <- list(rxn=rxn, kp=kp, kie=kie, net=net, events=events, res_steady_state=result_ss, res_dyn=result_dyn)
    list2file(res, file="model_example_res.txt")

    ptmt <- proc.time() - ptmi
    cat("completed (total running time = ", ptmt["elapsed"], " seconds)\n", sep="")

    # unload the library, set the original working directory and return the results
    suppressWarnings(try(dyn.unload(paste("model_example", .Platform$dynlib.ext, sep="")), silent=TRUE))
    setwd(mainDir)
    return(invisible(res))
}



## ####################################################
## ### ODEs CAN BE SOLVED USING EITHER R OR FORTRAN ###
## ####################################################
##
##  # NOTE: calculations are much slower when R is used instead of fortran,
##  #       use R only for small models and tests
##
##  # initialize isotopomer/metabolite concentrations
##  iso_conc <- iso_ini(net)
##
##  ##########################
##  ### SOLVE ODEs USING R ###
##  #   sys="m" metabolite dynamics only
##  #   sys="i" metabolite and isotope dynamics without enzyme isotope effect
##  #   sys="f" metabolite and isotope dynamics with enzyme isotope effect
##  sys <- "f"
##  if (sys == "m"){
##      yini   <- (net$m_map %*% iso_conc)[,1]
##      parms  <- list(kp=kp, net=net)
##      fun    <- dyn_metaR
##  }else if (sys == "i"){
##      yini    <- iso_conc
##      parms   <- list(kp=kp, net=net)
##        iso_eff <- NULL
##      fun     <- dyn_isoR
##  }else if (sys == "f"){
##      yini    <- iso_conc
##        iso_eff <- kie_ini(net, kie)
##      parms   <- list(kp=kp, net=net, iso_eff=iso_eff)
##      fun     <- dyn_isoR
##  }
##
##  system.time(
##      sim <- lsoda(y=yini, times=seq(0, 100), func=fun, parms=parms)
##  )
##  matplot(x=sim[,"time"], y=sim[,colnames(sim)!="time"], type = "l", lwd = 2, ylab = "concentration", xlab="time", main = "metabolite dynamics")
##
##  # calculate the derivative of isotopomer/metabolite concentrations, here at t=0
##  if (sys == "m"){
##    dX_dt <- dyn_metaR(NULL, yini, list(net=net, kp=kp))[[1]]
##  }else{
##    dX_dt <- dyn_isoR(NULL, yini, list(net=net, kp=kp, iso_eff=iso_eff))[[1]]
##  }
##  names(dX_dt) <- names(yini)
##
##  ################################
##  ### SOLVE ODEs USING FORTRAN ###
##
##  # build the fortran library initialize variables and parameters
##  sys <- "m"
##  R2Fortran(net, kp, sys=sys)
##  if (sys == "m"){
##      parms <- kp
##      yini  <- (net$m_map %*% iso_conc)[,1]
##  }else if (sys == "i"){
##      yini    <- iso_conc
##        iso_eff <- NULL
##      parms   <- kp
##  }else if (sys == "f"){
##      yini    <- iso_conc
##        iso_eff <- kie_ini(net, kie)
##      parms   <- c(kp, iso_eff)
##  }
##  # calculate the derivative of isotopomer/metabolite concentrations, here at t=0
##  # dX_dt <- DLLfunc(func="derivs", times=times, y=yini, parms=parms, dllname=lib_name, initfunc = "initmod")$dy
##  # steady-state simulations
##  #    sim <- runsteady(y=iso_conc, func="derivs", times=c(0,1e+8), parms=parms, dllname="lib_f", initfunc="initmod", nout=0)
##  # time-course simulations
##  #    sim <- lsoda(y=yini, func="derivs", parms=parms, dllname="lib_f", initfunc="initmod", times=times)
##  system.time(
##      sim <- lsoda(y=yini, func="derivs", parms=parms, dllname="lib_f", initfunc="initmod", times=seq(0,100))
##  )
##  matplot(x=sim[,"time"], y=sim[,names(yini)], type = "l", lwd = 2, ylab = "concentration", xlab="time", main = "isotopomer dynamics")

