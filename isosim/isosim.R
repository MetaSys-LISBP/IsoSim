# millard@insa-toulouse.fr
#
# Toolbox for modelling of isotopic and/or metabolic dynamics.
#
# Some R packages (nnls, numDeriv, Matrix, deSolve, stringr, pso, parallel and RColorBrewer) are
# required. They can be installed automatically by running the following command in an R console:
#
#     install.packages(c("nnls", "numDeriv", "Matrix", "deSolve", "stringr", "pso", "RColorBrewer", "parallel"))
#
# Rtools must be installed manually.
#
# Copyright 2019, INRA, France
# License: GNU General Public License v3 (see license.txt for details)


isosimVersion = "2.0.0"

###################
# LOAD R PACKAGES #
###################

load_pack <- function(packages){
  #
  # Load packages
  #
  # Args:
  #   packages (vector): names of packages to load
  #
  # Returns (vector):
  #   names of packages that could not be loaded
  #
  err <- c()
  for (i in packages){
    if (!suppressWarnings(require(i, character.only=TRUE, quietly=TRUE))){
      err <- c(err, i)
      # cat("The package '", i, "' is required. To install this package, please run the following command:\n   install.packages('", i, "')\n", sep="")
      # install.packages(i, dependencies=TRUE, verbose=FALSE)
      # suppressWarnings(require(i, character.only=TRUE, quietly=TRUE))
    }
  }
  return(invisible(list(err=err)))
}

res_load <- load_pack(c("nnls", "numDeriv", "Matrix", "deSolve", "stringr", "RColorBrewer", "pso", "parallel"))

if (length(res_load$err)){
  lpack <- paste(res_load$err, collapse="', '")
  msg <- c("The following R package(s) are required:", paste("   ", res_load$err, sep=""))
  msg <- c(msg, "To install missing package(s), please run the following command:", paste("   install.packages(c('", lpack, "'))\n", sep=""))
  stop(paste(msg, collapse="\n"))
}

source("nlsic.R")

#########################
# LOAD ISOSIM FUNCTIONS #
#########################

# Opposite of %in% to identify elements which are *not* in an object
"%ni%" <- Negate("%in%")

# rBind and cBind are deprecated since R version 3.2.0,  base's rbind() and cbind() now work fine with S4 objects
rBind <- rbind
cBind <- cbind

# Create color palettes
# usage: fun_col(X), where X is the number of colors to generate
fun_col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

net2mat <- function(rxn, meas=NULL, fixed=NULL, add_eq=NULL){
  #
  # Construct a model (generate stoichiometric and isotopic matrices, identify label inputs, etc)
  #
  # Args:
  #   rxn (list): input network for which the equation system is built, 'rxn' must be a list
  #               where each key is a reaction name that contains the following objects:
  #      $R      names of reactants
  #      $C      stoichiometric coefficients
  #      $T      carbon atom transitions (optional)
  #              Atom transitions can be omitted only for metabolites not involved in carbon scrambling
  #              The carbon transitions of a given metabolite must be consistent through the network
  #              A parallel network can be decoupled from the carbon network only if it does not produce
  #              metabolites involved in isotope scrambling
  #      $E      rate equation; can be a constant or a kinetic rate law (with names of kinetic
  #              parameters matching those of the corresponding vector)
  #   meas:   names of EMUs which must be simulated
  #   fixed:  names of metabolites with a fixed concentrations (sinks are automatically identified
  #            and their concentration is fixed, so this is mainly for internal metabolites)
  #   add_eq: vector containing additional equations to calculate before rate equations
  #            (can be useful to assign the value of some variables, e.g. 'p_1=2*v1' and
  #            a rate equation can refers to the parameter 'p_1'); can also be some code
  #            to run before flux calculation
  #
  # Returns (list):
  #   $s_mat          stoichiometric matrix
  #   $i_mat          isotopic matrix (rows=EMUs, columns=fluxomers, 0 * 0 if no carbon
  #                   transition are given in rxn)
  #   $i_rows_info    vector of metabolite names (for rows of i_mat)
  #   $i_cols_info    information on each column of i_mat:
  #                      $flx: flux name
  #                      $m: substrate(s)
  #                      $i: substrate(s) EMUs
  #   $m_map          'EMUs -> metabolites' mapping matrix
  #   $e_map          'EMUs -> enrichments' mapping matrix
  #   $eq             kinetic equations used to calculate fluxes (one per reaction)
  #   $add_eq         additional equations (calculated before $eq, so can be used e.g. to assign
  #                   some variables); can also be fortran or R code
  #   $fixed          sink metabolites ($in and $out for substrates and products, respectively, and
  #                   $user for those constrained by user - names passed in the vector 'fixed' -)
  #   $min_meas       minimal set of EMUs of source metabolites required to simulate measurements
  #
  
  # create stoichiometric model
  Smodel <- createSmodel(rxn, fixed=fixed)
  
  # create measurement dataset (i.e. all metabolites) if not provided
  if (is.null(meas)){
    meas <- create_meas_set(Smodel)
  }
  
  # identify minimal set of EMUs to simulate measurements
  lEMUs <- EMU_decomp(rxn, meas, Smodel$s_mat)
  
  # create isotopic matrix
  colS <- which(colnames(lEMUs)=="substrate")
  colP <- which(colnames(lEMUs)=="product")
  
  nSmax <- length(colS)
  i_cols_info <- list(flx = c(),
                      m = matrix(0, nrow=nSmax, ncol=0, dimnames = list("num" = paste(seq(nSmax)), "fluxomer" = NULL)),
                      i = matrix(0, nrow=nSmax, ncol=0, dimnames = list("num" = paste(seq(nSmax)), "fluxomer" = NULL)))
  
  i_app_c <- c()
  i_app_flm <- c()
  i_app_ip <- c()
  
  for (i in seq(nrow(lEMUs))){
    rxnn <- lEMUs$reaction[i]
    r <- rxn[[rxnn]]
    
    # identify substrates EMUs
    colSf <- colS[!is.na(lEMUs[i, colS])]
    lS <- unlist(lEMUs[i, colSf])
    lC <- unlist(lEMUs[i, colSf+1])
    lco <- unlist(lEMUs[i, colSf+3])
    lnC <- gsub(",", "-", lC)
    nIsot <- as.numeric(lEMUs[i, colSf+2])
    isoMap <- expand.grid(lapply(nIsot, FUN=function(x) seq(0,x)))
    nFluxomer <- nrow(isoMap)
    lIsot <- sapply(seq(length(lC)), FUN=function(x){ paste(paste(lS[x], lnC[x], sep="_"), isoMap[, x], sep="-M") })
    flm <- paste(lEMUs$reaction[i], apply(lIsot, 1, paste, collapse="_"), sep="_")
    
    # add substrates for the corresponding fluxomers
    i_app_c <- c(i_app_c, rep(lco, each=nFluxomer))
    i_app_flm <- c(i_app_flm, rep(flm, length(lS)))
    i_app_ip <- c(i_app_ip, as.vector(lIsot))
    
    # get reaction, EMUs, and substrates for each fluxomer
    i_cols_info$m <- mergevMat(lS, i_cols_info$m, nFluxomer)
    i_cols_info$i <- mergemMat(t(lIsot), i_cols_info$i)
    i_cols_info$flx <- c(i_cols_info$flx, rep(rxnn, nFluxomer))
    
    # add product
    lP <- unlist(lEMUs[i, colP])
    if (lP != "none"){
      lC_p <- unlist(lEMUs[i, colP+1])
      lco_p <- unlist(lEMUs[i, colP+3])
      lnC_p <- gsub(",", "-", lC_p)
      isoMap_p <- rowSums(isoMap)
      lIsot_p <- paste(lP, paste(lnC_p, "-M", isoMap_p, sep=""), sep="_")
      # append fluxomers, coefficients, EMUs
      i_app_c <- c(i_app_c, rep(lco_p, each=nFluxomer))
      i_app_flm <- c(i_app_flm, flm)
      i_app_ip <- c(i_app_ip, as.vector(lIsot_p))  
    }
  }
  
  i_app_m_m <- c()
  i_app_m_c <- c()
  i_app_m_flx <- c()
  
  for (r in names(rxn)){
    i_cols_info$m <- mergevMat(rxn[[r]]$R, i_cols_info$m, 1)
    i_cols_info$flx <- c(i_cols_info$flx, r)
    i_cols_info$i <- mergevMat(rxn[[r]]$R, i_cols_info$i, 1)
    i_app_m_m <- c(i_app_m_m, rxn[[r]]$R)
    i_app_m_c <- c(i_app_m_c, rxn[[r]]$C)
    i_app_m_flx <- c(i_app_m_flx, rep(r, length(rxn[[r]]$R)))
  }
  
  fluxomers <- unique(i_app_flm)
  EMUs <- unique(i_app_ip)
  i_mat <- matrix(0, nrow = length(EMUs)+nrow(Smodel$s_mat), ncol = length(fluxomers)+length(Smodel$flx), dimnames = list("EMU" = c(EMUs, rownames(Smodel$s_mat)), "fluxomer" = c(fluxomers, Smodel$flx)))
  i_mat <- fillMat(i_mat, i_app_ip, i_app_flm, i_app_c)
  i_mat <- fillMat(i_mat, i_app_m_m, i_app_m_flx, i_app_m_c)
  
  i_cols_info$i_filled   <- (i_cols_info$i != "0")
  i_rows_info <- unlist(sapply(strsplit(rownames(i_mat), "_"), "[[", 1))
  i_mat[i_rows_info %in% Smodel$sinks, ] <- 0
  
  m_map <- matrix(0, nrow = nrow(Smodel$s_mat), ncol = nrow(i_mat), dimnames = list("meta" = rownames(Smodel$s_mat), "EMU" = rownames(i_mat)))
  m_map[cbind(rownames(Smodel$s_mat), rownames(Smodel$s_mat))] <- 1
  
  nCiso <- str_count(EMUs, "-")
  weight_iso <- as.numeric(unlist(sapply(strsplit(EMUs, "-M"), "[[", 2)))
  rn_emap <- unlist(sapply(strsplit(EMUs, "-M"), "[[", 1))
  urn_emap <- unique(rn_emap)
  e_map <- matrix(0, nrow = length(urn_emap), ncol = nrow(i_mat), dimnames = list("meta" = urn_emap, "EMU" = rownames(i_mat)))
  e_map[cbind(rn_emap, EMUs)] <- weight_iso/nCiso
  
  f_map <- matrix(0, nrow = length(urn_emap), ncol = nrow(i_mat), dimnames = list("meta" = urn_emap, "EMU" = rownames(i_mat)))
  f_map[cbind(rn_emap, EMUs)] <- 1
  
  min_meas <- rownames(i_mat)[(i_rows_info %in% Smodel$ext_in) & (rownames(i_mat) %ni% Smodel$ext_in)]
  
  eq <- paste(names(rxn), sapply(rxn, "[[", "E"), sep="=")
  
  # return objects
  return(list(rxn         = rxn,
              i_mat       = Matrix(i_mat, sparse=TRUE),
              s_mat       = Matrix(Smodel$s_mat, sparse=TRUE),
              m_map       = Matrix(m_map, sparse=TRUE),
              e_map       = Matrix(e_map, sparse=TRUE),
              f_map       = Matrix(f_map, sparse=TRUE),
              meas        = meas,
              i_rows_info = i_rows_info,
              i_cols_info = i_cols_info,
              eq          = eq,
              add_eq      = add_eq,
              min_meas    = min_meas,
              fixed       = list("in"   = Smodel$ext_in,
                                 "out"  = Smodel$ext_out,
                                 "user" = Smodel$sinks)))
}


fillMat <- function(mat, a, b, d){
  #
  # Fill isotopic matrix.
  #
  tmp <- cbind(a, b)
  rDuplicated <- (duplicated(tmp) | duplicated(tmp, fromLast=TRUE))
  mat[tmp[!rDuplicated,]] <- as.numeric(d[!rDuplicated])
  if (sum(!rDuplicated) != length(a)){
    dup <- unique(tmp[rDuplicated,])
    for (j in seq(nrow(dup))){
      mat[t(dup[j,])] <- sum(as.numeric(d[apply(tmp, 1, FUN=function(x) identical(x, dup[j,]))]))
    }
  }
  return(mat)
}

create_meas_set <- function(Smodel){
  #
  # Returns the set of internal metabolites of the stoichiometric model.
  #
  meas <- list()
  for(i in rownames(Smodel$s_mat)){
    if (i %ni% Smodel$ext_in){
      meas[[length(meas)+1]] <- list(name=i)
    }
  }
  return(meas)
}

mergevMat <- function(lS, m, nFluxomer){
  #
  # Create an EMU matrix block.
  #
  toap <- matrix(lS, ncol=nFluxomer, nrow=length(lS))
  res <- mergemMat(toap, m)
  return(res)
}

mergemMat <- function(toap, m){
  #
  # Merge two EMU matrix blocks.
  #
  if (nrow(toap) < nrow(m)){
    toap <- rbind(toap, matrix("0", nrow=nrow(m)-nrow(toap), ncol=ncol(toap)))
  } else if (nrow(toap) > nrow(m)){
    m <- rbind(m, matrix("0", nrow=nrow(toap)-nrow(m), ncol=ncol(m)))
  }
  return(cbind(m, toap))
}

createSmodel <- function(rxn, fixed=c()){
  #
  # Construct the stoichiometric model.
  #
  rn <- names(rxn)
  mn <- unique(unlist(lapply(rxn, "[[", 1)))
  s_mat <- matrix(0, nrow=length(mn), ncol=length(rn), dimnames=list("meta"=mn, "flux"=rn))
  
  for (r in rn){
    rr <- rxn[[r]]$R
    for (m in seq(length(rr))){
      s_mat[rr[m], r] <- s_mat[rr[m], r] + rxn[[r]]$C[m]
    }
  }
  
  # identify sinks
  metabolites <- rownames(s_mat)
  ext_in  <- metabolites[rowSums(s_mat > 0) == 0]
  ext_out <- metabolites[rowSums(s_mat < 0) == 0]
  
  # set stoichiometric coefficients of sink and fixed metabolites at 0 (i.e. fixed or assigned concentrations)
  sinks <- unique(c(ext_in, ext_out, fixed))
  s_mat[sinks, ] <- 0
  
  return(list(s_mat=s_mat, sinks=sinks, ext_in=ext_in, ext_out=ext_out, metabolites=metabolites, flx=rn))
}

EMU_decomp <- function(rxn, meas, s_mat){
  #
  # Identify all EMUs required to simulate measurements.
  #
  
  # initialize the minimal set of EMUs
  lres <- unique(getCTsin(meas, rxn, s_mat))
  
  # add all EMUs producing reactions
  end <- FALSE
  processed <- 0
  while (end == FALSE){
    ns <- getS(lres, processed, rxn, s_mat)
    processed <- nrow(lres)+1
    lres <- unique(mergedf(lres, ns))
    if (nrow(lres)+1 == processed){
      end <- TRUE
    }
  }
  
  # add missing EMUs consumming reactions
  id <- which(colnames(lres) == "substrate")
  for (r in seq(nrow(lres))){
    lre <- lres[r,]
    m <- lre$product
    mC <- lre$Catoms
    rcp <- colnames(s_mat)[s_mat[m,]<0]
    for (ra in rcp){
      lrf <- lres[lres$reaction == ra,]
      tt <- as.vector(sapply(id, FUN=function(x) apply(lrf[,c(x, x+1)], 1, paste, collapse="")))
      if (paste(m, mC, sep="") %ni% tt){
        lres <- appendRcp(lres, m, mC, ra, s_mat[m, ra])
      }
    }
  }
  
  # remove rownames
  rownames(lres) <- c()
  
  # remove duplicated EMUs for reactions involving several EMUs
  # of the same metabolites
  nS <- (ncol(lres)-5)/4
  if (nS >= 2){
    nE <- nrow(lres)
    lres_lim <- lres[,seq(1,5)]
    # linearize blocks of lres
    for (i in seq(1, nS)){
      tmp <- lres[,c(1, seq(2+i*4, 5+i*4))]
      names(tmp) <- names(lres_lim)
      lres_lim <- rbind(lres_lim, tmp)
    }
    # identify & remove duplicated coefficients
    irow <- ((duplicated(lres_lim) | duplicated(lres_lim[nrow(lres_lim):1, ])[nrow(lres_lim):1]) & !is.na(lres_lim[,3]))
    tmp <- unique(lres_lim[irow,])
    for (i in seq(nrow(tmp))){
      dup <- apply(lres_lim, 1, FUN=function(x, j) {identical(paste(x,collapse=";"), j)}, j=paste(tmp[i,],collapse=";"))
      id <- which(dup)[-1]
      lres[cbind(id %% nE, 5+4*(id %/% nE))] <- 0
    }
  }
  
  # return results
  return(lres)
}

appendRcp <- function(lres, m, mC, ra, coeff){
  #
  # Append reactants.
  #
  lresn <- data.frame(reaction = ra,
                      product = "none", 
                      Catoms = NA, 
                      size = NA,
                      coeff = NA,
                      stringsAsFactors = FALSE)
  lresn2 <- data.frame(substrate = m,
                       Catoms = mC, 
                       size = str_count(mC,",")+1,
                       coeff = coeff,
                       stringsAsFactors = FALSE)
  lresn <- cbind(lresn, lresn2)
  lresn <- mergedf(lres, lresn)
  return(lresn)
}

getS <- function(lres, processed, rxn, s_mat){
  #
  # Identify substrate EMUs.
  #
  lresn <- data.frame(reaction = character(),
                      product = character(), 
                      Catoms = character(), 
                      size = integer(),
                      coeff = numeric(),
                      stringsAsFactors = FALSE)
  lresr <- lres[seq(from=processed, to=nrow(lres), by=1),]
  id <- which(colnames(lresr) == "substrate")
  for (i in id){
    for (j in seq(nrow(lresr))){
      if (!is.na(lresr[j,i])){
        m <- list(name = lresr[j,i],
                  Catoms = eval(parse(text=paste("c(", lresr[j,i+1], ")", sep=""))))
        lrm <- getCTin(m, rxn, s_mat)
        lresn <- mergedf(lresn, lrm)
      }
    }
  }
  return(lresn)
}

mergedf <- function(df1, df2){
  #
  # Merge dataframes.
  #
  nco1 <- ncol(df1)
  nco2 <- ncol(df2)
  if (nco2>nco1){
    NAnro <- rep(NA, nrow(df1))
    lapp <- data.frame(substrate = NAnro, Catoms = NAnro, size = NAnro, coeff=NAnro)
    for (z in seq((nco2-nco1)/4)){
      df1 <- cbind(df1, lapp)
    }
  }else if (nco2<nco1){
    NAnro <- rep(NA, nrow(df2))
    lapp <- data.frame(substrate = NAnro, Catoms = NAnro, size = NAnro, coeff=NAnro)
    for (z in seq((nco1-nco2)/4)){
      df2 <- cbind(df2, lapp)
    }
  }
  return(rbind(df1, df2))
}

getCTin <- function(meta, rxn, s_mat){
  #
  # Get carbon transitions of a given metabolite.
  #
  m <- meta$name
  Cskel <- meta$Catoms
  lres <- data.frame(reaction = character(),
                     product = character(), 
                     Catoms = character(), 
                     size = integer(),
                     coeff = numeric(),
                     stringsAsFactors = FALSE)
  r_in <- colnames(s_mat)[s_mat[m,] > 0]
  for (r in r_in){
    id <- which(rxn[[r]][["R"]] == m & rxn[[r]][["T"]] != "")
    if (length(id)){
      for (i in id){
        if (is.null(Cskel)){
          Cskel <- seq(nchar(rxn[[r]][["T"]][i]))
        }
        lp <- rxn[[r]][["R"]][which(rxn[[r]][["C"]] < 0)]
        ta <- c(r, m, paste(Cskel, collapse=","), length(Cskel), rxn[[r]][["C"]][i])
        tmp <- unlist(strsplit(rxn[[r]][["T"]][i], ""))[Cskel]
        for (p in unique(lp)){
          idp <- which(rxn[[r]][["R"]] == p)
          for (j in idp){
            ct_p <- unlist(strsplit(rxn[[r]][["T"]][j], ""))
            ct_match <- which(!is.na(match(ct_p, tmp)))
            if (length(ct_match) > 0){
              ta <- c(ta, p, paste(ct_match, collapse=","), length(ct_match), rxn[[r]][["C"]][j])
            }
          }
        }
        lta <- length(ta)
        nco <- ncol(lres)
        if (nco < lta){
          NAnro <- rep(NA, nrow(lres))
          lapp <- data.frame(substrate = NAnro, Catoms = NAnro, size = NAnro, coeff=NAnro)
          for (z in seq((lta-nco)/3)){
            lres <- cbind(lres, lapp)
          }
        }else if (lta < nco){
          ta <- c(ta, rep(NA, nco-lta))
        }
        lres[nrow(lres)+1,] <- ta
      }
    }
  }
  return(lres)
}

getCTsin <- function(products, rxn, s_mat){
  #
  # Get all carbon transitions.
  #
  lres <- data.frame(reaction = character(),
                     product = character(), 
                     Catoms = character(), 
                     size = integer(),
                     coeff = numeric(),
                     stringsAsFactors = FALSE)
  for(m in products){
    lrm <- getCTin(m, rxn, s_mat)
    lres <- mergedf(lres, lrm)
  }
  return(lres)
}

mat2eq <- function(net, sys="i", anFun=NULL){
  #
  # Generate ODEs (FORTRAN code) from the stoichiometric and isotopic matrices
  # to calculate derivative of EMU and metabolite concentrations
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #
  # Returns (vector):
  #   FORTRAN code, each element is a line of the ODEs system
  #
  
  nf <- unlist(net$fixed)
  d_eq <- c()
  
  if (sys == "m"){
    fluxes <- colnames(net$s_mat)
    rowS   <- rownames(net$s_mat)
    for (meta in seq(length(rowS))){
      if (rowS[meta] %in% names(anFun) | rowS[meta] %in% nf){
        d_eq <- c(d_eq, paste("ydot(", meta, ")=0", sep=""))
      }else{
        sub_s_mat <- net$s_mat[meta,]
        f_prod <- paste(paste(sub_s_mat[sub_s_mat > 0], fluxes[sub_s_mat > 0], sep="*"), collapse="+")
        f_cons <- paste(paste(-sub_s_mat[sub_s_mat < 0], fluxes[sub_s_mat < 0], sep="*"), collapse="-")
        l      <- paste("ydot(", meta, ")=", paste(f_prod, "-", f_cons, sep=""), sep="")
        d_eq   <- c(d_eq, l)
      }
    }
  }else if (sys == "i"){
    EMUs <- rownames(net$i_mat)
    for (isot in seq(length(EMUs))){
      if (EMUs[isot] %in% names(anFun) | net$i_rows_info[isot] %in% nf){
        d_eq <- c(d_eq, paste("ydot(", isot, ")=0", sep=""))
      }else{
        id    <- which(net$i_mat[isot,] != 0)
        uname <- unique(net$i_cols_info$flx[id])
        fmat  <- paste(paste("y(", match(net$i_cols_info$i[,id], EMUs), ")", sep=""), net$i_cols_info$m[,id], sep="/")
        subs  <- rbind(net$i_mat[isot, id],
                       matrix(fmat, ncol=length(id), nrow=nrow(net$i_cols_info$i)),
                       rep("null", length(id)))
        iff   <- sapply(uname, function(x){paste(paste(subs[,which(net$i_cols_info$flx[id] == x)], collapse="*"), collapse="+")})
        flxv  <- paste(uname, paste("(", iff, ")", sep=""), sep="*")
        d_eq  <- c(d_eq, paste("ydot(", isot, ")=", paste(flxv, collapse="+"), sep=""))
      }
    }
  }
  
  # replace some terms
  lrep <- list("*null*"="+", "*null"="", "+1*"="+", "-1*"="-", "*1*"="*", "/1*"="*", "(1*"="(", "=1*"="=", "+-"="-", "*y(NA)/0"="", "y(NA)/0"="")
  for (i in names(lrep)){
    d_eq <- gsub(i, lrep[i], d_eq, fixed=TRUE)
  }
  
  return(d_eq)
}

R2Fortran <- function(net, parms, sys="i", lib_name="lib_f", anFun=NULL){
  #
  # Generate FORTRAN code of the model, compile and load the library
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   parms (vector): named vector of kinetic parameters
  #   lib_name (str): name of the fortran library (optional, 'lib_f' by default)
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #
  # Returns (vector):
  #   FORTRAN code, each element is a line of the model library
  #
  
  # calculate system-dependent parameters (headers, variable names, parameters, etc.)
  rnm <- rownames(net$s_mat)
  anF_eq <- NULL
  if (sys == "m"){
    lp   <- length(parms)
    nr   <- length(rnm)
    pnm  <- names(parms)
    # update metabolite concentrations for which an analytical function is provided
    notf <- which(rnm %ni% names(anFun))
    if (!is.null(anFun)){
      anF_eq <- unlist(lapply(names(anFun), FUN=function(x) paste("y(", which(rownames(net$s_mat) == x), ")=", anFun[[x]], sep="")))
    }
    mc_h <- "C   set METABOLITE concentrations"
    mc_v <- paste(rnm[notf], "=y(", notf, ")", sep="")
    d_h  <- "C   calculate derivatives of METABOLITE concentrations"
  }else if (sys == "i"){
    lp   <- length(parms)
    nr   <- nrow(net$i_mat)
    pnm  <- names(parms)
    # update EMU concentrations for which an analytical function is provided
    if (!is.null(anFun)){
      anF_eq <- unlist(lapply(names(anFun), FUN=function(x) paste("y(", which(rownames(net$i_mat) == x), ")=", anFun[[x]], sep="")))
    }
    mc_h <- "C   calculate METABOLITE concentrations by summing EMU concentrations"
    mc_v <- unlist(lapply(rnm, FUN=function(x) paste(x, "=", paste(paste("y(", which(rownames(net$i_mat) == x), ")", sep=""), collapse="+"), sep="")))
    d_h  <- "C   calculate derivatives of EMU concentrations"
  }
  
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
             "C calculate derivatives",
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
             "C   set values of MODEL PARAMETERS",
             paste("       common /myparms/ ", paste(pnm, collapse=",\n     * "), sep=""))
  
  # add analytical functions for metabolite/EMU concentrations
  if (!is.null(anFun)){
    fcode <- c(fcode,
               "C   set variables with analytical functions",
               split_lines(anF_eq))
  }
  
  # calculate metabolite concentrations
  fcode <- c(fcode,
             mc_h,
             split_lines(mc_v))
  
  if (!is.null(net$add_eq)){
    fcode <- c(fcode, "C   calculate VARIABLES",
               split_lines(net$add_eq))
  }
  
  # generate equations for metabolite or EMU derivatives
  d_eq <- mat2eq(net, sys, anFun=NULL)
  
  fcode <- c(fcode, "C   calculate FLUXES",
             split_lines(net$eq),
             "C   ---------------",
             "C   equation system",
             "C   ---------------",
             d_h,
             split_lines(d_eq),
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
  
  return(invisible(fcode))
}

split_lines <- function(v_str){
  #
  # Split a FORTRAN code line according to FORTRAN format (max 64 characters).
  #
  # Args:
  #   v_str (str): code line to split
  #
  # Returns (vector):
  #   FORTRAN code, where each element is a splitted code line
  #
  
  # regular expression to split lines after 64 characters
  rex <- paste("(?<=", paste(rep(".", 64), collapse=""), ")", sep="")
  # split line
  return(unlist(lapply(paste("       ", v_str, sep=""), FUN=function(x) paste(unlist(strsplit(x, rex, perl=T)), collapse="\n     &"))))
}

calc_fluxes <- function(net, kp, conc_m){
  #
  # Calculate fluxes from concentrations and kinetic parameters
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   conc_m (vector): named vector of concentrations
  #
  # Returns (vector):
  #   flux values
  #
  
  with(as.list(c(net, kp, conc_m)),{
    flx <- colnames(net$s_mat)
    eval(parse(text = paste(net$add_eq, collapse=";\n")))
    eval(parse(text = paste(net$eq, collapse=";\n")))
    eval(parse(text = paste("rates=c(", paste(paste(flx, flx, sep="="), collapse=","), ")", sep="")))
    
    return(rates)
  })
}

calc_fluxomers <- function(params, iso_conc, t){
  #
  # Calculate fluxomers from EMU concentrations and kinetic parameters
  #
  # Args:
  #   params (list): internal model object
  #   iso_conc (vector): name vector of EMUs concentrations
  #   t (float): time
  #
  # Returns (vector):
  #   fluxomer values
  #
  
  # update EMUs for which an analytical function is provided
  if (!is.null(params$anFun)){
    #iso_conc[names(params$anFun)] <- sapply(params$anFun, FUN=function(x,t) x(t), t=t)
    iso_conc[names(params$anFun)] <- sapply(params$anFun, FUN=function(x,t) {eval(parse(text=paste("fe=function(t){", x, "}", sep=""))); fe(t)}, t=t)
  }
  
  # calculate metabolite concentrations and fluxes
  #conc_m <- (params$net$m_map %*% iso_conc)[,1]
  conc_m <- iso_conc[grep("-M",names(iso_conc), invert=TRUE)]
  fluxes <- calc_fluxes(params$net, params$kp, conc_m)
  
  # calculate fluxomers
  fluxomers <- fluxes[params$net$i_cols_info$flx]
  for (ns in seq(nrow(params$net$i_cols_info$i))){
    tmp <- params$net$i_cols_info$i_filled[ns,]
    fluxomers[tmp] <- fluxomers[tmp] * iso_conc[params$net$i_cols_info$i[ns, tmp]]/conc_m[params$net$i_cols_info$m[ns, tmp]]
  }
  
  return(fluxomers)
}

dyn_isoR <- function(t, iso_conc, params){
  #
  # Calculate time derivatives of EMU concentrations at time t.
  #
  # This function must be callable by rootSolve and deSolve.
  # Note: use preferentially the fortran library to solve ODEs of large systems
  #
  # Args:
  #   params (list): internal model object
  #   iso_conc (vector): name vector of EMUs concentrations
  #   t (float): time
  #
  # Returns (list):
  #   time-course derivatives of EMUs
  # 
  
  with(as.list(c(iso_conc, params, t)),{
    
    flx_v <- calc_fluxomers(params, iso_conc, t)
    dI_dt <- params$net$i_mat %*% flx_v
    
    # return derivatives
    list(dI_dt@x)
  })
}

iso2list <- function(iso_conc){
  #
  # Build a list of EMUs concentrations
  #
  # Args:
  #   iso_conc (vector): name vector of EMUs concentrations
  #
  # Returns (list):
  #   concentrations of EMUs, each key is a metabolite and the corresponding object is
  #   a named vector of EMUs abundances)
  # 
  
  nmi  <- names(iso_conc)
  mn   <- sapply(strsplit(nmi[str_detect(nmi, "_")], "_"), "[[", 1)
  isol <- list()
  for (i in unique(mn)){
    isol[[i]][sapply(strsplit(nmi[mn==i], "_"), "[[", 2)] <- iso_conc[mn==i]
  }
  
  return(isol)
}

iso_ini <- function(net, p=0.0, meta_conc=NULL, iso_conc=NULL, anFun=NULL){
  #
  # Initialize EMUs concentration vector at a given isotopic enrichment.
  #
  # Args:
  #   p (float): isotopic enrichment, between 0 and 1
  #   meta_conc (vector): metabolite concentrations (optional, set at 1 by default)
  #   iso_conc (vector): EMU concentrations (optional, normalized to 1 by
  #                      default, or to the initial concentration of the
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #
  # Returns (list):
  #   initial concentrations of all EMUs
  # 
  
  # initialize EMUs concentrations at 0
  EMUs <- setNames(rep(0, nrow(net$i_mat)), rownames(net$i_mat))
  
  # initialize metabolite concentrations at 1, and update those given in 'meta_conc'
  metabolites <- setNames(rep(1, nrow(net$s_mat)), rownames(net$s_mat))
  metabolites[names(meta_conc)] <- meta_conc
  
  # get EMUs, without the metabolite prefix
  rne <- sapply(strsplit(rownames(net$e_map), "_"), "[[", 1)
  tmp <- rownames(net$i_mat)[net$i_rows_info %in% rne]
  rnm_isoall <- sapply(strsplit(tmp[tmp %ni% rne], "_"), "[[", 2)
  
  # get names of metabolites for which an analytical time-dependent function is provided
  if (is.null(anFun)){
    rnm_anFun <- NULL
  }else{
    rnm_anFun <- unique(sapply(strsplit(names(anFun), "_"), "[[", 1))
  }
  
  # calculate the (absolute) EMU concentrations for each metabolite
  for (i in unique(net$i_rows_info)){
    # if metabolite i is in the isotope scrambling network
    if (i %in% rne){
      # if the concentrations of some of its EMU are given in 'iso_conc', use these values (other EMUs are set at 0)
      if (i %in% names(iso_conc)){
        if (!is.null(names(iso_conc[[i]]))){
          EMUs[paste(i, names(iso_conc[[i]]), sep="_")] <- iso_conc[[i]]/sum(iso_conc[[i]]) * metabolites[i]
        }
        # if no analytical function is provided for its EMUs, initialize at the abundance 'p'
      }else{
        if (i %ni% rnm_anFun){
          ip <- rnm_isoall[net$i_rows_info == i]
          ip <- ip[!is.na(ip)]
          nC <- str_count(ip, "-")
          inc <- as.numeric(sapply(strsplit(ip, "-M"), "[[", 2))
          abundance <- choose(nC, inc)*(1-p)**(nC-inc)*p**inc
          EMUs[paste(i, ip, sep="_")] <- abundance * metabolites[i]
        }
      }
    }
    # just set metabolite concentration
    #}else{
    if (i %in% names(meta_conc)){
      EMUs[i] <- metabolites[i]
    }else{
      EMUs[i] <- 1
    }
    #}
  }
  
  return(EMUs)
}

check_lib <- function(lib_n, sys, net, kp, anFun=NULL, force=FALSE){
  #
  # Just load the library if it already exists, otherwise generate the FORTRAN code, compile
  # and load the library.
  #
  # Args:
  #   lib_n (str): name of the fortran library (optional, 'lib_f' by default)
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #   force (bool): of True, recreate and recompile the library even if it is already found
  # 
  
  lib_n_ext <- paste(lib_n, .Platform$dynlib.ext, sep="")
  # check if the library exists
  if (!file.exists(lib_n_ext) | force){
    # generate the fortran code and compile it
    cat("  library '", lib_n_ext, "' not found, compilation in progress...\n", sep="")
    R2Fortran(net, kp, sys=sys, lib_name=lib_n, anFun=anFun)
  }else{
    # (re)load the library
    suppressWarnings(try(dyn.unload(lib_n_ext), silent=TRUE))
    dyn.load(lib_n_ext)
  }
  
  return(invisible(NULL))
}

fb2net <- function(fxch){
  #
  # Calculate net fluxes from forward and backward fluxes.
  #
  # Args:
  #   fxch (vector): named vector of fluxes (for reversible reactions, names of backward
  #                  fluxes must end by 'rev')
  #
  # Returns (vector):
  #   flux vector updated with all net fluxes
  # 
  
  nm   <- names(fxch)
  nmb  <- nm[str_detect(nm, "rev")]
  nmf  <- gsub("rev", "", nmb)
  fnet <- c(fxch[nm %ni% c(nmf, nmb)], fxch[nmf]-fxch[nmb])
  
  return(fnet)
}

xch2fb <- function(par){
  #
  # Calculate forward and backward fluxes from net fluxes and exchange coefficients.
  #
  # Args:
  #   par (vector): named vector of fluxes (for reversible reactions, names of xch
  #                 coefficients must end by 'xch')
  #
  # Returns (vector):
  #   flux vector updated with all forward and backward fluxes
  # 
  
  nm   <- names(par)
  nmb  <- nm[str_detect(nm, "xch")]
  for (i in nmb){
    nmf  <- gsub("xch", "", i)
    rev <- par[i]/(1-par[i])
    par[paste(nmf, "f", sep="")] <- rev + par[nmf]
    par[paste(nmf, "r", sep="")] <- rev
  }
  return(par)
}

dynamic <- function(net, kp, times, meta_conc=NULL, iso_conc=NULL, p=0.0, events=NULL, atol=1e-6, rtol=1e-6, anFun=NULL, method="R", dllname="lib_f"){
  #
  # Simulate dynamics of metabolite concentrations and EMU abundances.
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   times (vector): simulation times
  #   meta_conc (vector): metabolite concentrations (optional, set at 1 by default)
  #   iso_conc (vector): EMU concentrations (optional, normalized to 1 by
  #                      default, or to the initial concentration of the
  #   p (float): isotopic enrichment, between 0 and 1
  #   events (dataframe): events to pass to deSolve
  #   atol (float): absolute tolerance of the ODEs solver (1.e-6 by default)
  #   rtol (float): relative tolerance of the ODEs solver (1.e-6 by default)
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #   method (str): solver engine ('R' or 'FORTRAN', default: 'FORTRAN')
  #   dllname (str): name of the fortran library (optional, 'lib_f' by default)
  #
  # Returns (list):
  #   $times          vector of times at which system variables are simulated
  #   $metabolites    matrix of metabolite concentrations (columns=variables, rows=times)
  #   $isotopologues  matrix of isotopologue abundances (normalized to 1; columns=variables, rows=times)
  #   $enrichments    matrix of metabolite enrichments (normalized to 1; columns=variables, rows=times)
  #   $run            calculation time
  # 
  
  # solve ODEs with lsoda solver
  if (method == "R"){
    # run simulations
    sim <- lsoda(y        = iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc, anFun=anFun),
                 func     = dyn_isoR,
                 parms    = list(kp=kp, net=net, anFun=anFun),
                 events   = list(data=events),
                 times    = times,
                 atol     = atol,
                 rtol     = rtol)
    # remove time
    sim <- sim[, -1]
    # update variables for which an analytical function is provided
    if (!is.null(anFun)){
      sim[, names(anFun)] <- sapply(anFun, FUN=function(x,t) {eval(parse(text=paste("fe=function(t){", x, "}", sep=""))); fe(t)}, t=times)
      # sapply(anFun, FUN=function(x,t) x(t), t=times)
    }
  }else if (method == "FORTRAN"){
    # run simulations
    sim <- lsoda(y        = iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc, anFun=anFun),
                 func     = "derivs",
                 parms    = kp,
                 dllname  = dllname,
                 initfunc = "initmod",
                 events   = list(data=events),
                 times    = times,
                 atol     = atol,
                 rtol     = rtol)
    # update variables for which an analytical function is provided
    if (!is.null(anFun)){
      sim[, names(anFun)] <- sapply(anFun, FUN=function(x,t) {eval(parse(text=paste("fe=function(t){", x, "}", sep=""))); fe(t)}, t=times)
    }
    # remove time
    sim <- sim[, -1]
  }
  
  # calculate time-course variables (metabolite concentrations, EMUs/isotopologues abundance and enrichments)
  res <- list("an_fun"=anFun, "times"=times, "sim"=sim)
  res$metabolites <- sim[, grep("-M",colnames(sim), invert=TRUE)]
  res$isotopologues <- sim[, grep("-M",colnames(sim))]/res$metabolites[, net$i_rows_info[grep("-M",colnames(sim))]]
  res$enrichments <- t(apply(res$isotopologues, 1, FUN=function(x) (net$e_map[,colnames(res$isotopologues)]  %*% x)[,1]))
  res$flx <- t(apply(res$metabolites, 1, FUN=function(x) calc_fluxes(net, kp, x)))
  
  return(res)
}

calc_deriv <- function(net, kp, anFun=NULL, times=seq(0,100,1), meta_conc=NULL, iso_conc=NULL, events=NULL, trf=NULL, p=0.0, rtol=1e-6, atol=1e-6, subDir="deriv", sys="i", method="R", compile="b", dllname="lib_f", unloadDll=FALSE, plot=TRUE, write=TRUE, plot_opt=list(inset=c(-0.35, 0))){
  #
  # Calculate derivatives of metabolite concentrations and EMU abundances. Usefull for debugging models.
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   times (vector): simulation times
  #   meta_conc (vector): metabolite concentrations (optional, set at 1 by default)
  #   iso_conc (vector): EMU concentrations (optional, normalized to 1 by
  #                      default, or to the initial concentration of the
  #   p (float): isotopic enrichment, between 0 and 1
  #   events (dataframe): events to pass to deSolve
  #   trf (function): transformation function applied on parameters (kp) before running simulations (e.g. to
  #                   convert net/xch fluxes to forw/back)
  #   atol (float): absolute tolerance of the ODEs solver (1.e-6 by default)
  #   rtol (float): relative tolerance of the ODEs solver (1.e-6 by default)
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   method (str): solver engine ('R' or 'FORTRAN', default: 'FORTRAN')
  #   compile (str): model library options (only used if 'method'='FORTRAN')
  #                    'b': generate FORTRAN code, compile and load the library
  #                    'l': only load the library (it is assumed to be already compiled)
  #                    'n': do nothing (the library is assumed to be already compiled and loaded)
  #   unloadDll (bool): if True, unload the library after runing simulations (default: False)
  #   plot (bool): if True, plot simulated data (default: True)
  #   plot_opt (list): plotting options
  #   write (bool): if True, save simulated data (default: True)
  #   dllname (str): name of the fortran library (optional, 'lib_f' by default)
  #   subDir (str): name of the folder to save simulation results
  #
  # Returns (list):
  #   $net                    network and matrices
  #   $kp                     vector of parameters
  #   $res_dyn$times          vector of times at which system variables are simulated
  #   $res_dyn$metabolites    matrix of metabolite concentrations (columns=variables, rows=times)
  #   $res_dyn$isotopologues  matrix of isotopologue abundances (normalized to 1; columns=variables, rows=times)
  #   $res_dyn$enrichments    matrix of metabolite enrichments (normalized to 1; columns=variables, rows=times)
  #   $res_dyn$run            calculation time
  #
  
  # create results directory
  mainDir <- getwd()
  if (write | plot | compile != "n"){
    if (!file.exists(subDir)){
      dir.create(file.path(mainDir, subDir))
    }
    setwd(file.path(mainDir, subDir))
  }
  
  # if a transformation function is provided (e.g. to convert net/xch to forw/rev,
  # or to calculate determined fluxes from free fluxes), apply this transformation
  if (!is.null(trf)){
    kp <- trf(kp)
  }
  
  # compile library if needed
  if (method == "FORTRAN"){
    if (compile == "b"){
      R2Fortran(net, kp, sys=sys, lib_name=dllname, anFun=anFun)
    }else if (compile == "l"){
      lib_f <- paste(dllname, .Platform$dynlib.ext, sep="")
      suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
      dyn.load(lib_f)
    }else if (compile == "c"){
      check_lib(dllname, sys=sys, net=net, kp=kp, anFun=anFun)
    }
  }
  
  # initialize concentrations
  yini <- iso_ini(net, p=p, meta_conc=meta_conc, iso_conc=iso_conc, anFun=anFun)
  
  # calculate derivatives
  if (method == "R"){
    stop("Not implemented yet for 'method'='R'.")
  }else if (method == "FORTRAN"){
    sim <- DLLfunc(func     = "derivs",
                   times    = times,
                   y        = yini,
                   parms    = kp,
                   dllname  = dllname,
                   initfunc = "initmod")
  }
  
  # go back to working directory
  setwd(mainDir)
  
  # initial concentrations & parameters
  sim$yini <- yini
  sim$par <- kp
  
  # identify non-null derivatives and the reactions involved
  sim$nonnull <- list()
  sim$nonnull$all <- sim$dy[sim$dy!=0]
  for (i in names(sim$nonnull$all)){
    sim$nonnull[[paste("fluxomers_", i, sep="")]] <- net$i_mat[i, net$i_mat[i,]!=0]
  }
  
  return(sim)
}

list2file <- function(l, file = paste(deparse(substitute(l)), ".txt", sep = "")) {
  #
  # Save a list in a text file.
  #
  # Args:
  #   l (list): list to save
  #   file (str): filename
  #  
  
  tmpw <- getOption("width")
  tmpp <- getOption("max.print")
  options(width=10000, max.print=99999999)
  sink(file)
  print(l)
  sink()
  options(width=tmpw, max.print=tmpp)
  return(invisible(NULL))
}

uplo2uco <- function(param, upper=NULL, lower=NULL) {
  #
  # Return a list with a matrix u and a vector co such that u %*% param - co >= 0
  # to translate the inequalities param <= upper and param >= lower
  #
  # Args:
  #   param (vector): named vector of parameters
  #   upper (vector): named vector of upper bound constraints on parameters (must be in 'param')
  #   lower (vector): named vector of lower bound constraints on parameters (must be in 'param')
  #
  # Returns (list):
  #   $u (matrix)
  #   $co (vector)
  #
  
  u=matrix(0., nrow=length(upper)+length(lower), ncol=length(param))
  co=numeric(nrow(u))
  colnames(u)=names(param)
  rownames(u)=c(paste(names(upper), " <= ", upper, sep=""), paste(names(lower), " >= ", lower, sep=""))
  names(co)=rownames(u)
  # fill u and co
  u[iseq(length(upper)), names(upper)]=diag(-1, length(upper))
  co[iseq(length(upper))]=-upper
  u[length(upper)+iseq(length(lower)), names(lower)]=diag(1, length(lower))
  co[length(upper)+iseq(length(lower))]=lower
  return(list(u=u, co=co))
}

iseq=function(n) {
  #
  # Positive sequence of integer numbers 1:n.
  # if n=0 then integer(0) is returned
  #
  
  seq.int(from=1, to=n, length=n)
}

plot_iso <- function(x, y, xlab, ylab, type="l", lty=1, lwd=2, lpos="topright", inset=c(-0.14, 0), pch=20, cex=0.45, bty="n", ncolumns=1, las=1){
  #
  # Plot isotopic dynamics.
  #
  
  colpal <- fun_col(ncol(y))
  #colnames(y) <- gsub("\\_[0-9]", "_", gsub("\\-[0-9]*", "", colnames(y)))
  matplot(x, y, type=type, lty=lty, lwd=lwd, col=colpal, xlab=xlab, ylab=ylab, las=las)
  legend(lpos, inset=inset, legend=colnames(y), col=colpal, pch=pch, cex=cex, bty=bty, ncol=ncolumns)
}

simulate <- function(net, kp, anFun=NULL, times=seq(0,100,1), meta_conc=NULL, iso_conc=NULL, events=NULL, trf=NULL, p=0.0, rtol=1e-6, atol=1e-6, subDir="sim", sys="i", method="R", compile="b", dllname="lib_f", unloadDll=FALSE, plot=TRUE, write=TRUE, plot_opt=list(inset=c(-0.35, 0))){
  #
  # Run time-course simulations of metabolite concentrations and EMU abundances.
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   times (vector): simulation times
  #   meta_conc (vector): metabolite concentrations (optional, set at 1 by default)
  #   iso_conc (vector): EMU concentrations (optional, normalized to 1 by
  #                      default, or to the initial concentration of the
  #   p (float): isotopic enrichment, between 0 and 1
  #   events (dataframe): events to pass to deSolve
  #   trf (function): transformation function applied on parameters (kp) before running simulations (e.g. to
  #                   convert net/xch fluxes to forw/back)
  #   atol (float): absolute tolerance of the ODEs solver (1.e-6 by default)
  #   rtol (float): relative tolerance of the ODEs solver (1.e-6 by default)
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   method (str): solver engine ('R' or 'FORTRAN', default: 'FORTRAN')
  #   compile (str): model library options (only used if 'method'='FORTRAN')
  #                    'b': generate FORTRAN code, compile and load the library
  #                    'l': only load the library (it is assumed to be already compiled)
  #                    'n': do nothing (the library is assumed to be already compiled and loaded)
  #   unloadDll (bool): if True, unload the library after runing simulations (default: False)
  #   plot (bool): if True, plot simulated data (default: True)
  #   plot_opt (list): plotting options
  #   write (bool): if True, save simulated data (default: True)
  #   dllname (str): name of the fortran library (optional, 'lib_f' by default)
  #   subDir (str): name of the folder to save simulation results
  #
  # Returns (list):
  #   $net                    network and matrices
  #   $kp                     vector of parameters
  #   $res_dyn$times          vector of times at which system variables are simulated
  #   $res_dyn$metabolites    matrix of metabolite concentrations (columns=variables, rows=times)
  #   $res_dyn$isotopologues  matrix of isotopologue abundances (normalized to 1; columns=variables, rows=times)
  #   $res_dyn$enrichments    matrix of metabolite enrichments (normalized to 1; columns=variables, rows=times)
  #   $res_dyn$run            calculation time
  #
  
  # perform some checks
  if (!setequal(net$min_meas, names(anFun))){
    stop("Analytical functions do not match the minimal set of EMUs required to perform simulations.")
  }
  
  # create results directory
  mainDir <- getwd()
  if (write | plot | compile != "n"){
    if (!file.exists(subDir)){
      dir.create(file.path(mainDir, subDir))
    }
    setwd(file.path(mainDir, subDir))
  }
  
  # if a transformation function is provided (e.g. to convert net/xch to forw/rev,
  # or to calculate determined fluxes from free fluxes), apply this transformation
  if (!is.null(trf)){
    kp <- trf(kp)
  }
  
  # compile library if needed
  if (method == "FORTRAN"){
    if (compile == "b"){
      R2Fortran(net, kp, sys=sys, lib_name=dllname, anFun=anFun)
    }else if (compile == "l"){
      lib_f <- paste(dllname, .Platform$dynlib.ext, sep="")
      suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
      dyn.load(lib_f)
    }else if (compile == "c"){
      check_lib(dllname, sys=sys, net=net, kp=kp, anFun=anFun)
    }
  }
  
  # run simulations
  result_dyn <- dynamic(net       = net,
                        kp        = kp,
                        times     = times,
                        meta_conc = meta_conc,
                        iso_conc  = iso_conc,
                        anFun     = anFun,
                        events    = events,
                        rtol      = rtol,
                        atol      = atol,
                        p         = p,
                        method    = method,
                        dllname   = dllname)
  
  # unload library
  if (method == "FORTRAN" & unloadDll){
    suppressWarnings(try(dyn.unload(paste(dllname, .Platform$dynlib.ext, sep="")), silent=TRUE))
  }
  
  # save simulation results in 'res.txt'
  if (write){
    res <- list(net=net, kp=kp, res_dyn=result_dyn)
    list2file(res, file="res.txt")
  }
  
  # plot simulated and measured data
  if (plot){
    pdf(file="plot.pdf", width=10, height=7)
    par(mar=c(5.1, 5.1, 4.1, 7.1), xpd=TRUE)
    layout(matrix(seq(1, 4), 2, 2, byrow = TRUE))
    plot_iso(result_dyn$times, result_dyn$metabolites, "time", "metabolite conc.", inset=plot_opt$inset, ncolumns=2)
    plot_iso(result_dyn$times, result_dyn$flx, "time", "flux", inset=plot_opt$inset, ncolumns=2)
    plot_iso(result_dyn$times, result_dyn$isotopologues, "time", "isotopologue", inset=plot_opt$inset, ncolumns=2)
    plot_iso(result_dyn$times, result_dyn$enrichments, "time", "enrichment", inset=plot_opt$inset, ncolumns=2)
    dev.off()
  }
  
  # go back to the main working directory and return the results
  setwd(mainDir)
  return(invisible(res))
}

saveEnv <- function(subDir, sname="Renv.RData"){
  #
  # Save the R environment.
  #
  # Args:
  #   subDir (str): name of the subfolder in which the R environment should be saved
  #   sname (str): filename of the RData file
  #  
  
  # create results directory
  mainDir <- getwd()
  if (!file.exists(subDir)){
    dir.create(file.path(mainDir, subDir))
  }
  setwd(file.path(mainDir, subDir))
  # save R session
  save.image(file = sname)
  # go back to the original working directory
  setwd(mainDir)
}

addConstraints <- function(luco, te_loc_add, te_upc_add, eq){
  #
  # Return a list with a matrix u and a vector co such that u %*% param - co >= 0
  # to translate the inequalities param <= upper and param >= lower
  #
  # This functions only adds non-linear constraints (on determined fluxes) to linear
  # constraints constructed by uplo2uco(), which should thus be called first.
  #
  # Args:
  #   luco (list): set of constraints constructed by uplo2uco()
  #   eq (list): set of equations of determined fluxes
  #   te_upc_add (vector): named vector of upper bound constraints on determined fluxes (must be in 'eq')
  #   te_loc_add (vector): named vector of lower bound constraints on determined fluxes (must be in 'eq')
  #
  # Returns (list):
  #   $u (matrix)
  #   $co (vector)
  #
  
  u=luco$u
  co=luco$co
  for (i in seq(length(te_loc_add))){
    nam=names(te_loc_add)[i]
    # parse equation
    l <- c("+",unlist(strsplit(gsub("[-]"," - ",gsub("[+]"," + ",eq[nam])), " ")))
    # get indices of fluxes
    idx <- which(l %in% colnames(u))
    if (length(idx)){
      # name of dependent flux
      n <- paste(nam, ">=", te_loc_add[i])
      u_add <- matrix(0, nrow=1, ncol=ncol(u), dimnames=list(n, colnames(u)))
      u_add[,l[idx]] <- eval(parse(text=paste("c(", paste(paste(l[idx-1], rep(1,length(idx)),sep=""),collapse=","), ")",sep="")))
      u <- rbind(u, u_add)
      co_add <- eval(parse(text=paste("c('", n, "'=", as.character(te_loc_add[i]), ")", sep="")))
      co <- c(co, co_add)
    }
  }
  for (i in seq(length(te_upc_add))){
    nam <- names(te_upc_add)[i]
    # parse equation
    l <- c("+",unlist(strsplit(gsub("[-]"," - ",gsub("[+]"," + ",eq[nam])), " ")))
    idx <- which(l %in% colnames(u))
    if (length(idx)){
      n <- paste(nam, "<=", te_upc_add[i])
      u_add <- matrix(0, nrow=1, ncol=ncol(u), dimnames=list(n, colnames(u)))
      u_add[,l[idx]] <- -eval(parse(text=paste("c(", paste(paste(l[idx-1], rep(1,length(idx)),sep=""),collapse=","), ")",sep="")))
      u <- rbind(u, u_add)
      co_add <- -eval(parse(text=paste("c('", n, "'=", as.character(te_upc_add[i]), ")", sep="")))
      co <- c(co, co_add)
    }
  }
  return(list(u=u, co=co))
}

run_mc <- function(resF, niter=10, net, kp, times, to_est, te_upc, te_loc, te_upc_det, te_loc_det, eq_det, data_meas, data_meas_conc, sd_meas, meta_conc=NULL, iso_conc=NULL, mode="EMUs", anFun=NULL, events=NULL, trf=NULL, p=0.0, atol=1e-6, rtol=1e-6, subDir="fit", nlsic_ctrl=list(errx=1.e-4, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T), method="R", sys="i",compile="c", dllname="lib_f", unloadDll=TRUE, mc.cores=NULL, save_res=TRUE){
  #
  # Runs Monte-Carlo sensitivity analysis.
  #
  # For details on arguments and usage, see fit().
  #
  
  res_mc <- matrix(NA, ncol=length(resF$result$par)+1, nrow=niter, dimnames=list(r=NULL, par=c(names(resF$result$par), "res")))
  res_mc_sim <- array(NA, dim=c(dim(resF$result$retres$sim_all[[mode]]), niter), dimnames=list(time=times, isotopologue=colnames(resF$result$retres$sim_all[[mode]]), iter=NULL))
  
  data_enr_ori <- resF$result$retres$sim_all[[mode]][,colnames(data_meas)[-1]]
  
  kp_mc <- kp
  ntr <- names(resF$result$par)[names(resF$result$par) %in% names(kp)]
  kp_mc[ntr] <- resF$result$par[ntr]
  
  meta_conc_mc <- meta_conc
  ntr <- names(meta_conc)[names(meta_conc) %in% to_est]
  meta_conc_mc[ntr] <- resF$result$par[ntr]
  data_meas_conc <- resF$result$par[names(data_meas_conc)]
  
  if (ncol(data_meas) < 3){
    nv <- 1
  }else{
    nv <- ncol(data_enr_ori)
  }
  
  nlsic_ctrl$trace = 0
  
  cat("run monte carlo analysis (", niter, " iterations)", "\n", sep="")
  # run mc iterations
  if (is.null(mc.cores)){
    resmc <- lapply(seq(niter), mc_internal,
                    data_enr_ori=data_enr_ori,
                    data_meas=data_meas,
                    nv=nv,
                    data_meas_conc=data_meas_conc,
                    net = net,
                    times = times,
                    kp = kp_mc,
                    to_est = to_est,
                    te_upc = te_upc,
                    te_loc = te_loc,
                    te_upc_det = te_upc_det,
                    te_loc_det = te_loc_det,
                    eq_det = eq_det,
                    sd_meas = sd_meas,
                    meta_conc = meta_conc_mc,
                    iso_conc = iso_conc,
                    mode = mode,
                    anFun = anFun,
                    trf = trf,
                    events = events,
                    p = p,
                    nlsic_ctrl = nlsic_ctrl,
                    method = method,
                    compile = compile,
                    unloadDll = unloadDll,
                    subDir=subDir,
                    save_res=save_res,
                    dllname=dllname)
  }else{
    resmc <- mclapply(seq(niter), mc_internal,
                      data_enr_ori=data_enr_ori,
                      data_meas=data_meas,
                      nv=nv,
                      data_meas_conc=data_meas_conc,
                      net = net,
                      times = times,
                      kp = kp_mc,
                      to_est = to_est,
                      te_upc = te_upc,
                      te_loc = te_loc,
                      te_upc_det = te_upc_det,
                      te_loc_det = te_loc_det,
                      eq_det = eq_det,
                      sd_meas = sd_meas,
                      meta_conc = meta_conc_mc,
                      iso_conc = iso_conc,
                      mode = mode,
                      anFun = anFun,
                      trf = trf,
                      events = events,
                      p = p,
                      nlsic_ctrl = nlsic_ctrl,
                      method = method,
                      compile = compile,
                      unloadDll = unloadDll,
                      subDir=subDir,
                      save_res=save_res,
                      mc.cores=mc.cores,
                      dllname=dllname)
  }
  # extract summary
  for (i in seq(niter)){
    res_mc[i,] <- c(resmc[[i]]$result$par, sum(resmc[[i]]$result$res**2))
  }
  #res_mc[i,] <- c(restmp$result$par, sum(restmp$result$res**2))
  #res_mc_sim[,,i] <- restmp$result$retres$sim_all[[mode]]
  
  if (save_res){
    pdf(file="res_mc.pdf")
    for (i in colnames(res_mc)){
      hist(res_mc[, i], xlab="value", main=i, breaks=10)
    }
    dev.off()
  }
  
  mc_summary <- matrix(NA, nrow=length(to_est)+1, ncol=6, dimnames=list(row=c(to_est, "res"), col=c("opt", "mean", "median", "sd", "ci_2.5", "ci_97.5")))
  mc_summary[, "opt"] <- c(resF$result$par, sum(resF$result$res**2))
  mc_summary[, "mean"] <- apply(res_mc, 2, mean)
  mc_summary[, "median"] <- apply(res_mc, 2, median)
  mc_summary[, "sd"] <- apply(res_mc, 2, sd)
  mc_summary[, "ci_2.5"] <- apply(res_mc, 2, FUN=function(x) {quantile(x,0.025)})
  mc_summary[, "ci_97.5"] <- apply(res_mc, 2, FUN=function(x) {quantile(x,0.975)})
  
  return(list(summary=mc_summary, sim=res_mc_sim, par=res_mc))
}

mc_internal <- function(i, data_enr_ori, data_meas, nv, data_meas_conc, net , times , kp , to_est , te_upc , te_loc , te_upc_det , te_loc_det , eq_det , sd_meas , meta_conc , iso_conc , mode , anFun , trf , events , p , nlsic_ctrl , method , compile , unloadDll, subDir, save_res, dllname){
  # 
  # Generate and fit noisy datasets for each Monte Carlo iteration. 
  #
  
  # create noisy isotopic dataset
  noise <- matrix(rnorm(length(as.numeric(data_enr_ori)), mean=0.0, sd=0.02), ncol=nv, dimnames=list(r=NULL, c=colnames(data_meas)[-1]))
  data_enr_noise <- data_enr_ori + noise
  data_enr_noise[data_enr_noise<0] <- 0
  data_enr_noise[data_enr_noise>1] <- 1
  data_enr_noise <- cbind(times, data_enr_noise)
  # create noisy concentration dataset
  data_meas_conc_noise <-  data_meas_conc + rnorm(length(data_meas_conc), mean=0, sd=data_meas_conc*0.1)
  # estimate parameters
  restmp <- fit(net        = net,
                times      = times,
                kp         = kp,
                to_est     = to_est,
                te_upc     = te_upc[to_est],
                te_loc     = te_loc[to_est],
                te_upc_det = te_upc_det,
                te_loc_det = te_loc_det,
                eq_det     = eq_det,
                data_meas  = list(iso=data_enr_noise, conc=data_meas_conc_noise),
                sd_meas    = sd_meas,
                meta_conc  = meta_conc,
                iso_conc   = iso_conc,
                mode       = mode,
                anFun      = anFun,
                trf        = trf,
                events     = events,
                p          = p,
                subDir     = paste(subDir, "_", i, sep=""),
                nlsic_ctrl = nlsic_ctrl,
                method     = method, 
                niter      = NULL,
                compile    = "n",
                unloadDll  = unloadDll,
                save_res   = save_res,
                dllname    = dllname)
  # save results
  #res_mc[i,] <- c(restmp$result$par, sum(restmp$result$res**2))
  #res_mc_sim[,,i] <- restmp$result$retres$sim_all[[mode]]
  return(restmp)
}

test_khi2 <- function(nb_points, k_val, nb_par){
  #
  # Perform Khi2 statistical test.
  #
  # Args:
  #   nb_points (int): number of data points
  #   k_val (float): khi2 value (cost)
  #   nb_par (int): number of free parameters
  #
  # Returns (list):
  #   $'khi2 value' (float): khi2 value (cost)
  #   $'data points' (int): number of data points
  #   $'fitted parameters' (int): number of free parameters
  #   $'degrees of freedom' (int): degrees of freedom
  #   $'khi2 reduced value' (float): chi2 reduced value
  #   $'p-value, i.e. P(X^2<=value)' (float): p value
  #   $conclusion (str): message indicating whether the models fits (or not) the data at 95% confidence interval
  #
  
  df <- nb_points - nb_par
  p_val <- pchisq(k_val, df=df)
  khi2test <- list("khi2 value"                  = k_val,
                   "data points"                 = nb_points,
                   "fitted parameters"           = nb_par,
                   "degrees of freedom"          = df,
                   "khi2 reduced value"          = k_val/df,
                   "p-value, i.e. P(X^2<=value)" = p_val)
  if (p_val > 0.95){
    khi2test$conclusion <- "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD."
  }else{
    khi2test$conclusion <- "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD."
  }
  return(khi2test)
}

fit <- function(net, kp, times, to_est, te_upc, te_loc, te_upc_det, te_loc_det, eq_det, data_meas, sd_meas, meta_conc=NULL, iso_conc=NULL, mode="EMUs", anFun=NULL, events=NULL, trf=NULL, p=0.0, atol=1e-6, rtol=1e-6, subDir="fit", nlsic_ctrl=list(errx=1.e-4, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T), method="R", sys="i", compile="b", dllname="lib_f", unloadDll=TRUE, niter=10, mc.cores=NULL, save_res=TRUE){
  #
  # Estimate model parameters by fitting experimental data.
  #
  # All the matrices and the results are saved in a text file ('res.txt'), and simulated/experimental data are plotted in pdf ('fit.pdf').
  #
  # Args:
  #   net (list): network object constructed by net2mat()
  #   kp (vector): named vector of kinetic parameters
  #   times (vector): vector containing the times at which explicit estimates for system variables are desired,
  #                   measurement times MUST be included!
  #   to_est (vector): names of parameters/metabolite concentrations to estimate
  #   te_upc (vector): named vector of upper constraints of free parameters (optional, but recommended)
  #   te_loc (vector): named vector of lower constraints of free parameters (optional, but recommended)
  #   data_meas (matrix): experimental data, the first column must be named 'time' and contain measurement times
  #   sd_meas (matrix): standard deviation on measurements (used to calculate the AIC)
  #   meta_conc (vector): metabolite concentrations (optional, set at 1 by default)
  #   iso_conc (vector): EMU concentrations (optional, normalized to 1 by
  #                      default, or to the initial concentration of the
  #   p (float): isotopic enrichment, between 0 and 1
  #   events (dataframe): events to pass to deSolve
  #   trf (function): transformation function applied on parameters (kp) before running simulations (e.g. to
  #                   convert net/xch fluxes to forw/back)
  #   atol (float): absolute tolerance of the ODEs solver (1.e-6 by default)
  #   rtol (float): relative tolerance of the ODEs solver (1.e-6 by default)
  #   anFun (list): analytical functions of some metabolites/EMUs, where each key contains
  #                 the analytical functions of its time-course labeling (MUST be function of t)
  #   sys (str): ODEs to generate:
  #                'm': metabolite dynamics only, no isotopic equations
  #                'i': metabolite and EMUs dynamics
  #   method (str): solver engine ('R' or 'FORTRAN', default: 'FORTRAN')
  #   compile (str): model library options (only used if 'method'='FORTRAN')
  #                    'b': generate FORTRAN code, compile and load the library
  #                    'l': only load the library (it is assumed to be already compiled)
  #                    'n': do nothing (the library is assumed to be already compiled and loaded)
  #   unloadDll (bool): if True, unload the library after runing simulations (default: False)
  #   plot (bool): if True, plot simulated data (default: True)
  #   plot_opt (list): plotting options
  #   write (bool): if True, save simulated data (default: True)
  #   dllname (str): name of the fortran library (optional, 'lib_f' by default)
  #   subDir (str): name of the folder to save simulation results
  #
  # Returns (list):
  #   $net (list): network and matrices
  #   $params_ini (vector): named vector of initial parameters and metabolite concentrations
  #   $AIC (float): Akaike's information criterion, calculated with the following formula:
  #           AIC = 2k + n*log(RSS/n), where k is the number of free parameters, n is the
  #           number of measurements, and RSS is the residuum sum of square
  #   $chi2 (list): results of the khi2 statistical test:
  #       $'khi2 value' (float): khi2 value (cost)
  #       $'data points' (int): number of data points
  #       $'fitted parameters' (int): number of free parameters
  #       $'degrees of freedom' (int): degrees of freedom
  #       $'khi2 reduced value' (float): chi2 reduced value
  #       $'p-value, i.e. P(X^2<=value)' (float): p value
  #       $conclusion (str): message indicating whether the models fits (or not) the data at 95% confidence interval
  #   $result (list): optimization results:
  #       $times: vector of times at which system variables are simulated
  #       $metabolites: matrix of metabolite concentrations (columns=variables, rows=times)
  #       $isotopologues: matrix of isotopologue abundances (normalized to 1; columns=variables, rows=times)
  #       $enrichments: matrix of metabolite enrichments (normalized to 1; columns=variables, rows=times)
  #       $run: calculation time
  #
  
  # create results directory
  mainDir <- getwd()
  if (save_res | (compile %in% c("b", "c"))){
    if (!file.exists(subDir)){
      dir.create(file.path(mainDir, subDir))
    }
    setwd(file.path(mainDir, subDir))
  }
  
  if (!is.null(trf)){
    kp <- trf(kp)
  }
  # compile library if needed
  if (method == "FORTRAN"){
    if (compile == "b"){
      R2Fortran(net, kp, sys=sys, lib_name=dllname, anFun=anFun)
    }else if (compile == "l"){
      lib_f <- paste(dllname, .Platform$dynlib.ext, sep="")
      suppressWarnings(try(dyn.unload(lib_f), silent=TRUE))
      dyn.load(lib_f)
    }else if (compile == "c"){
      check_lib(dllname, sys=sys, net=net, kp=kp, anFun=anFun)
    }
  }
  
  data_meas_conc <- data_meas$conc
  data_meas <- data_meas$iso
  
  # create free parameters vector
  nconc <- intersect(names(meta_conc), to_est)
  npar <- intersect(names(kp), to_est)
  fpar <- c(kp[npar], meta_conc[nconc]) 
  # set constraints on free & determined fluxes
  # construct matrix u and vector co such that u%*%param-co>=0
  luco=uplo2uco(fpar, te_upc, te_loc)
  if (length(eq_det)){
    luco=addConstraints(luco, te_loc_det, te_upc_det, eq_det)
  }
  u=luco$u
  co=luco$co
  # fit experimental data using nlsic
  result <- nlsic(fpar,
                  cost,
                  u=u, co=co,
                  control=nlsic_ctrl,
                  e=NULL, eco=NULL, flsi=lsi_ln,
                  data_meas,
                  data_meas_conc,
                  sd_meas,
                  times,
                  kp,
                  net,
                  nconc,
                  npar,
                  meta_conc,
                  iso_conc,
                  mode,
                  anFun,
                  events,
                  trf,
                  p,
                  atol,
                  rtol,
                  method,
                  dllname)
  
  # calculate Akaike's information criterion
  #AIC <- length(result$res) * log10(sum(result$res**2)/(as.numeric(sd_meas)**2)/length(result$res)) + 2*length(to_est)
  
  # calculate chi2
  chi2 <- test_khi2(nb_points=length(result$res), k_val=sum(result$res**2), nb_par=length(to_est))
  
  ret <- list(params_ini=list(kp=kp, meta_conc=meta_conc), chi2=chi2, net=net, result=result, constraints=luco)
  
  # run monte carlo sensitivity analysis
  if (!is.null(niter)){
    ret[["sens"]] <- run_mc(ret, niter, net, kp, times, to_est, te_upc, te_loc, te_upc_det, te_loc_det, eq_det, data_meas, data_meas_conc, sd_meas, meta_conc=meta_conc, iso_conc=iso_conc,
                            mode=mode, anFun=anFun, events=events, trf=trf, p=p, atol=atol, rtol=rtol,
                            subDir="mc", nlsic_ctrl=nlsic_ctrl, method=method, sys=sys, compile="c", dllname=dllname, unloadDll=FALSE, mc.cores=mc.cores, save_res=FALSE)
  }
  
  # unload library
  if (method == "FORTRAN" & unloadDll){
    suppressWarnings(try(dyn.unload(paste(dllname, .Platform$dynlib.ext, sep="")), silent=TRUE))
  }
  
  if (save_res){
    # save simulation results in 'res.txt'
    list2file(ret, file="res.txt")
    
    # plot simulated and measured data
    pdf("results.pdf", width=6, height=4.5)
    metab <- sapply(strsplit(colnames(data_meas)[-1], "_"), "[[", 1)
    for (m in unique(metab)){
      tp <- colnames(data_meas)[which(metab==m)+1]
      cols <- fun_col(sum(metab==m))
      matplot(x=times, y=result$retres$sim_all[[mode]][, tp], ylim=c(0,1), lty = 1, type = "l", lwd = 1, ylab = mode, xlab="time", main = "fitting results", col=cols)
      matplot(x=times, y=data_meas[, tp], type="p", lwd=2, ylim=c(0,1), add=T, pch="+", col=cols)
      legend("topright", 4, tp, pch="-", col=cols, cex=0.6, ncol=2)
    }
    dev.off()
  }
  
  # go back to the main working directory and return the results
  setwd(mainDir)
  return(invisible(ret))
}

cost <- function(fp, cjac=F, data, data_meas_conc, sd_meas, times, kp, net, nconc, npar, meta_conc=NULL, iso_conc=NULL, mode="EMUs", anFun=NULL, events=NULL, trf=NULL, p=0.0, atol=1e-6, rtol=1e-6, method="R", dllname="lib_f") {
  #
  # Cost function used for optimization.
  #
  
  # update parameters and concentrations
  kp[npar] <- fp[npar]
  meta_conc[nconc] <- fp[nconc]
  
  # if a transformation function is provided (e.g. to convert net/xch to forw/rev,
  # or to calculate determined fluxes from free fluxes), apply this transformation
  if (!is.null(trf)){
    kp <- trf(kp)
  }
  
  # run simulations
  result_dyn <- dynamic(net       = net,
                        kp        = kp,
                        times     = times,
                        meta_conc = meta_conc,
                        iso_conc  = iso_conc,
                        events    = events,
                        p         = p,
                        anFun     = anFun,
                        atol      = atol,
                        rtol      = rtol,
                        method    = method,
                        dllname   = dllname)
  
  # calculate residuum vector (difference between measured and simulated data)
  # isotopic data
  datasim <- result_dyn[[mode]][, colnames(data)[-1]]
  res <- as.numeric((datasim - data[, -1])/sd_meas$iso)
  # metabolite concentrations
  if (!is.null(data_meas_conc)){
    res <- c(res, as.numeric((meta_conc[names(data_meas_conc)] - data_meas_conc)/sd_meas$conc))
  }
  
  return(list(res=res, sim=datasim, sim_all=result_dyn))
}

fit_subsystem <- function(subsystem){
  #
  # Estimate model parameters (e.g. fluxes or kinetic parameters) for each subsystem.
  #
  
  cat("\nEstimate fluxes and concentrations for subsystem ", subsystem$name, "\n\n", sep="")
  res_subsystem <- with (subsystem, {
    # create the equation system
    #subnet <- net2mat(rxn_subnet, meas=colnames(data_meas_subnet)[-1])
    subnet <- net2mat(rxn_subnet)
    # fitted label input
    anFun_fit_subnet <- anFun4net(enr_in$all, subnet)
    # run optimization
    resF_subnet <- fit(net        = subnet,
                       times      = times,
                       kp         = kp_subnet,
                       to_est     = te_subnet,
                       te_upc     = te_upc_subnet[te_subnet],
                       te_loc     = te_loc_subnet[te_subnet],
                       te_upc_det = NULL,
                       te_loc_det = NULL,
                       eq_det     = NULL,
                       data_meas  = data_meas_subnet,
                       sd_meas    = sd_meas,
                       meta_conc  = meta_conc_subnet,
                       iso_conc   = NULL,
                       mode       = "enrichments",
                       anFun      = anFun_fit_subnet,
                       trf        = NULL,
                       events     = NULL,
                       p          = 0.0,
                       subDir     = paste("fit_subnet_", subsystem$name, sep=""),
                       nlsic_ctrl = list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                       method     = "FORTRAN",
                       niter      = niter,
                       mc.cores   = mc.cores,
                       dllname    = paste("lib_f_", subsystem$name, sep=""))
    return(resF_subnet)
  })
  return(res_subsystem)
}

cleanDLL <- function(){
  #
  # Unload all the libraries loaded by IsoSim (their names
  # contains *lib_f*).
  #
  getLoadedDLLs()
  a <- getLoadedDLLs()
  for (i in seq(length(a))){
    if (length(grep("lib_f", as.vector(unlist(a[i][[1]][2])))>0)){
      dyn.unload(as.vector(unlist(a[i][[1]][2])))
    }
  }
  getLoadedDLLs()
}

fit_subsystems <- function(subsystems, dirname="fit_subsystems", mc.cores=NULL){
  #
  # Calculate fluxes for all subsystems.
  #
  # Args:
  #   subsystems (list): subsystems to process
  #   dirname (str): subdirectory where results will be saved
  #   mc.cores (int): number of cores for parallel flux calculations (default = NULL, i.e. no parallelization).
  #
  # Returns (list):
  #   Flux calculation results for all subsystems.
  #
  
  # if there is only one subsystem, use the corresponding function
  if (("name" %in% names(subsystems)) & ("rxn_subnet" %in% names(subsystems))){
    
    res_fluxes <- fit_subsystem(subsystems)
    res_subsystems <- res_fluxes$sens$summary
    
  }else{
  
    # create results directory
    mainDir <- getwd()
    if (!file.exists(dirname)){
      dir.create(file.path(mainDir, dirname))
    }
    setwd(file.path(mainDir, dirname))
    
    # estimate fluxes for all subsystems
    if (is.null(mc.cores)){
      res_fluxes <- lapply(subsystems, fit_subsystem)
    }else{
      res_fluxes <- mclapply(subsystems, fit_subsystem, mc.cores=mc.cores)
    }
    
    # get summary
    res_subsystems <- lapply(lapply(res_fluxes, "[[", 6), "[[", 1)
    
    # go back to the initial working directory
    setwd(mainDir)
    
  }
  
  # return results
  return(list(res=res_fluxes, summary=res_subsystems))
}

#################
# TEST FUNCTION #
#################


isosim_test <- function(){
  #
  # IsoSim test function.
  #
  # Have a look to this code for examples of IsoSim usage.
  #
  
  cat("
Welcome to IsoSim!
      
This function tests the different features of IsoSim
(construction of isotopic models, simulation of labeling
dynamics, flux calculation using the ScalaFlux approach).

Have a look to the code (by running 'isosim_test' in R) for
detailed examples of IsoSim usage.

")
  
  ####################################
  cat("   ... initialize R environment ...\n")
  ####################################
  
  # create results directory
  
  wd <- getwd()
  if (!file.exists("test")){
    dir.create(file.path(wd, "test"))
  }
  setwd(file.path(wd, "test"))
  
  # parallelization options

  # single-core version can be used by setting 'numCores' to NULL,
  # here we use all CPU cores on the current host
  numCores <- detectCores()
  
  # number of Monte Carlo iterations for flux calculation
  mc_iter <- 4

  ####################################
  cat("   ... construct isotopic model ...\n")
  ####################################
  
  # network definition (list):
  #   $R (vector): reactant(s)
  #   $C (vector): stoichiometric coefficient(s)
  #   $E (vector): rate law (can be a constant, a variable or an expression)
  #   $T (vector): tracer atom transitions
  
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
  
  # construct the model (isotopic and stoichiometric matrices, identify label inputs, etc)
  net <- net2mat(rxn, add_eq=eq_det)
  
  ####################################
  cat("   ... simulate labeling dynamics ...\n")
  ####################################

  # labeling dynamics of label input (here we simulate a switch from unlabeled to fully-labeled nutrient Sout)
  # note: time variable (t) MUST be present in the analytical functions, even if the label input(s) is (are) constant
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
  cat("   ... fit label inputs ...\n")
  ####################################
  
  # fit labeling dynamics of some metabolites considered as label input (C and O) using analytical functions
  enr_input <- res$res_dyn$enrichments[, c("C_1", "O_1")]
  enr_in <- fit_label_input(enr_input, t=times, file="res_fit_enr", mc.cores=numCores)
  
  ####################################
  cat("   ... calculate fluxes ...\n\n")
  ####################################
  
  # definition of each minimal subsystem to analyze
  # (here, fluxes and pools are estimated for two minimal subsystems - S_E and S_P - of the example network)
  # subsystems (list):
  #   $rxn_subnet (list): network definition (list), as detailed above
  #   $meta_conc_subnet (vector): named vector of initial metabolite concentrations
  #   $kp_subnet (vector): named vector of model parameters
  #   $te_subnet (vector): free parameters to estimate (can be model parameters and metabolite concentrations)
  #   $te_upc_subnet (vector): named vector of upper bound constraints on free parameters
  #   $te_loc_subnet (vector): named vector of lower bound constraints on free parameters
  #   $data_meas_subnet (list): experimental data to fit
  #   $sd_meas (list): standard deviations on experimental data to fit
  #   $times (vector): simulation times (all measurement times must be included)
  #   $enr_in (list): list of fitted label inputs returned by 'fit_label_input()'
  #   $anFun (list): analytical functions (if any, otherwise should be NULL)
  #   $niter (int): number of Monte Carlo iterations for flux calculations
  #   $mc.cores (int): number of cores for parallelization
  
  subsystems <- list(
    
    S_E = list(name = "S_E",
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
    
    S_P = list(name = "S_P",
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
  
  # save summary of flux calculation results
  list2file(res_sub$summary, file="summary_minimal_subsystems.txt")
  
  ####################################
  cat("IsoSim test results:\n\n")
  ####################################
  
  cat("Estimation of fluxes (v8 and v19) and metabolite concentrations (E and P) in subsystems S_E and S_P.\n\n")

  res_test <- rbind(res_sub$summary$S_E[c("v8", "E"),], res_sub$summary$S_P[c("v19", "P"),])
  variables <- c("v8", "E", "v19", "P")
  true_value <- c(1.7, 0.5, 2.0, 20)
  error <- res_test[,"opt"] - true_value
  rel_error <- error/true_value
  subsystem <- c("S_E", "S_E", "S_P", "S_P")
  df_test <- data.frame(subsystem, variables, res_test, true_value, error, rel_error, row.names = NULL)
  
  print(df_test)
  
  cat("\n95% confidence intervals of estimated parameters include the true values ? ")
  if (all((true_value < res_test[,"ci_97.5"]) & (true_value > res_test[,"ci_2.5"]))){
    cat("Yes\n")
    T1 <- TRUE
  }else{
    cat("No\n")
    T1 <- FALSE
  }
  
  cat("Relative error of estimated parameters < 10% ? ")
  if (all(rel_error < 0.1)){
    cat("Yes\n\n")
    T2 <- TRUE
  }else{
    cat("No\n\n")
    T2 <- FALSE
  }
  
  if (T1 & T2){
    cat("All tests passed successfully.\n")
  }else{
    cat("Some tests did not passed, details provided above.\n")
  }
  
  # go back to initial directory
  setwd(wd)
  
}



########################
# LABEL INPUT ANALYSIS #
########################


# double logistic model
fn_doublelogistic <- function(par, times){
  # Simulate labeling dynamics for a given isotopologue using a double logistic model
  xpred <- par[1] + (par[2] - par[7] * times) * ((1/(1 + exp((par[3] - times)/par[4]))) - (1/(1 + exp(((par[5]-times)/par[6])))))
  return(xpred)
}

fn_doublelogistic_parinit <- function(x){
  # Initialize parameters
  mx <- max(x, na.rm=TRUE)
  mn <- min(x, na.rm=TRUE)
  par <- c(mn, mx - mn, 200, 1.5, 300, 1.5, 0.002)
  return(list(par=par, upper=c(1,1, 1000, 100, 1000, 100, 10), lower=-c(1,1, 1000, 100, 1000, 100, 10)))
}

fn_doublelogistic_eq <- function(m){
  # Returns analytical function as a string (can be evaluated by R)
  eq <- paste(m[1], "+ (", m[2],"-(", m[7], ")*t) * ((1/(1 + exp((", m[3] ,"-t)/(", m[4], ")))) -  (1/(1 + exp((", m[5] ,"-t)/(", m[6], ")))))", sep="")
  return(eq)
}

# logistic model
fn_M0_logistic <- function(par, times){
  # Simulate labeling dynamics for a given isotopologue using a double logistic
  xpred <- par[1] / (1+exp(-par[2]*(times-par[3]))) + par[4]
  return(xpred)
}

fn_M0_logistic_parinit <- function(x, times){
  # Initialize parameters
  par <- c(max(x)-min(x), 1, mean(times), min(x))
  lower <- c(0, -100, -1000, 0)
  upper <- c(1, 100, 1000, 1)
  return(list(par=par, upper=upper, lower=lower))
}

fn_M0_logistic_eq <- function(m){
  eq <- paste(m[1], ")/(1+exp(-(", m[2], ")*(t*(", m[4], ")-(", m[3], ")))", sep="")
  return(eq)
}

# logistic model enr
fn_enr_logistic <- function(par, times){
  # Simulate labeling dynamics for a given isotopologue using a double logistic
  xpred <- 1-(par[1] / (1+exp(-par[2]*(times-par[3]))) + par[4])
  return(xpred)
}

fn_enr_logistic_parinit <- function(x, times){
  # Initialize parameters
  par <- c(1+(min(x)-max(x)), 1, mean(times), 1-max(x))
  lower <- c(0, -100, -1000, 0)
  upper <- c(1, 100, 1000, 1)
  return(list(par=par, upper=upper, lower=lower))
}

fn_enr_logistic_eq <- function(m){
  eq <- paste("1-((", m[1], ")/(1+exp(-(", m[2], ")*(t-(", m[3], "))))+(", m[4], "))", sep="")
  return(eq)
}

cost_labin_single <- function(par, id, times, fun_dyn){
  # Simulate labeling dynamics for all isotopologues of the CID
  return(sum((fun_dyn(par, times) - id)**2))
}

cost_labin_enr <- function(par, id, times, fun_dyn){
  # Simulate labeling dynamics for all isotopologues of the CID
  sim <- fun_dyn(par, times)
  d <- sim - id
  d[sim<0] <- d[sim<0]*10
  d[sim>1] <- d[sim>1]*10
  return(sum(d**2))
}

fit_M0 <- function(id, times, niter){
  ctrl <- list(maxit=2000, trace=0, vectorize=TRUE, hybrid="improved", hybrid.control=list(trace=0))
  resopt <- NULL
  opt_resi <- Inf
  for (i in seq(niter)){
    par=fn_M0_logistic_parinit(id, times)
    res_i <- psoptim(par=par$par, fn=cost_labin_single, lower=par$lower, upper = par$upper, control=ctrl, id=id, times=times, fun_dyn=fn_M0_logistic)
    if (res_i$value < opt_resi){
      opt_resi <- res_i$value
      resopt <- res_i$par
    }
  }
  return(resopt)
}

fit_Mx_doublelog <- function(id, times, niter){
  ctrl <- list(maxit=2000, trace=0, vectorize=TRUE, hybrid="off", hybrid.control=list(trace=0))
  resopt <- NULL
  opt_resi <- Inf
  for (i in seq(niter)){
    par <- fn_doublelogistic_parinit(id)
    res_i <- psoptim(par=par$par, fn=cost_labin_single, lower=par$lower, upper = par$upper, control=ctrl, id=id, times=times, fun_dyn=fn_doublelogistic)
    if (res_i$value < opt_resi){
      opt_resi <- res_i$value
      resopt <- res_i$par
    }
  }
  return(resopt)
}

fit_CID_ind <- function(x, times){
  if (ncol(x) < 2){
    stop("at least 2 isotopologues should be provided\n")
  }
  cat("fit individual isotopologues to determine optimal starting points\n")
  name <- c(rep("M0", 4), unlist(lapply(seq(1, ncol(x)-1), FUN=function(x){rep(paste("M", x, sep=""), 7)})))
  id <- c(seq(4), rep(seq(1,7), ncol(x)-1))
  val <- rep(NA, length(id))
  lower <- c(0, -100, -1000, 0, rep(-1000, length(id)-4))
  upper <- c(1, 100, 1000, 1, rep(1000, length(id)-4))
  
  params <- data.frame(name, id, val, lower, upper)
  
  cat("  M0...\n")
  res_M0 <- fit_M0(id=x[,1], times=times, niter=10)
  params[params$name == "M0", "val"] <- res_M0
  plot(x=times,x[,1], type="p")
  lines(x=seq(from=min(times), to=max(times), length.out=100), y=fn_M0_logistic(par = res_M0, times = seq(from=min(times), to=max(times), length.out=100)))
  
  for (i in seq(2, ncol(x))){
    cat("  M", i-1, "...\n", sep="")
    res_Mx_dlog <- fit_Mx_doublelog(id=x[,i], times=times, niter=10)
    params[params$name == paste("M", i-1, sep=""), "val"] <- res_Mx_dlog
    plot(x=times,x[,i], type="p")
    lines(x=seq(from=min(times), to=max(times), length.out=100), y=fn_doublelogistic(par = res_Mx_dlog, times = seq(from=min(times), to=max(times), length.out=100)))
  }
  return(params)
}

fit_enr_log_doublelog <- function(id, times, niter){
  ctrl <- list(maxit=3000, trace=0, vectorize=TRUE, hybrid="improved", hybrid.control=list(trace=0, pgtol=1e-8, factr=1e-8))
  resopt <- list(value=Inf)
  for (i in seq(niter)){
    par <- fn_doublelogistic_parinit(id)
    res_i <- psoptim(par=par$par, fn=cost_labin_enr, lower=par$lower, upper = par$upper, control=ctrl, id=id, times=times, fun_dyn=fn_doublelogistic)
    if (res_i$value < resopt$value){
      resopt <- res_i
      funSet <- list(sim=fn_doublelogistic, eq=fn_doublelogistic_eq)
    }
  }
  for (i in seq(niter)){
    par <- fn_enr_logistic_parinit(id, times=times)
    res_i <- psoptim(par=par$par, fn=cost_labin_enr, lower=par$lower, upper = par$upper, control=ctrl, id=id, times=times, fun_dyn=fn_enr_logistic)
    if (res_i$value < resopt$value){
      resopt <- res_i
      funSet <- list(sim=fn_enr_logistic, eq=fn_enr_logistic_eq)
    }
  }
  return(list(resopt=resopt, funSet=funSet))
}

fit_ENR_ind <- function(x, times, fun_dyn=fit_enr_log_doublelog){
  cat("fit enrichment\n")
  params <- fun_dyn(id=x, times=times, niter=10)
  
  plot(x=times,x, type="p")
  lines(x=seq(from=min(times), to=max(times), length.out=100), y=params$funSet$sim(par = params$resopt$par, times = seq(from=min(times), to=max(times), length.out=100)))
  
  return(params)
}

fit_enr_custom <- function(id, times, ctrl=list(maxit=2000, trace=0, vectorize=TRUE, hybrid="improved", hybrid.control=list(trace=0, pgtol=1e-8, factr=1e-8)), niter=1, funs=list(parinit=NULL, sim=NULL, eq=NULL)){
  resopt <- list(value=Inf)
  for (i in seq(niter)){
    par <- funs$parinit(id, times)
    res_i <- psoptim(par=par$par, fn=cost_labin_enr, lower=par$lower, upper = par$upper, control=ctrl, id=id, times=times, fun_dyn=funs$sim)
    if (res_i$value < resopt$value){
      resopt <- res_i
      funSet <- list(sim=funs$sim, eq=funs$eq)
    }
  }
  return(list(resopt=resopt, funSet=funSet))
}

sim_all <- function(par, x, times){
  sim <- matrix(NA, nrow=length(times), ncol=ncol(x))
  sim[,1] <- fn_M0_logistic(par[seq(4)], times)
  for (i in seq(2,ncol(x))){
    sim[,i] <- fn_doublelogistic(par[4+seq((i-2)*7+1, (i-1)*7)], times)
  }
  return(sim)
}

sim_all_1 <- function(par, x, times){
  nc <- ncol(x)
  sim <- matrix(NA, nrow=length(times), ncol=nc)
  sim[,1] <- fn_M0_logistic(par[seq(4)], times)
  if (nc>2){
    for (i in seq(2,nc-1)){
      sim[,i] <- fn_doublelogistic(par[4+seq((i-2)*7+1, (i-1)*7)], times)
    }
  }
  sim[,nc] <- 1 - apply(sim, 1, sum)
  return(sim)
}

cost_labin_all <- function(par, x, times, n){
  # Simulate labeling dynamics for all isotopologues of the CID
  sim <- sim_all(par, x, times)
  resi <- sim-x
  resi[sim<0] <- resi[sim<0]*(2**n)
  resi[sim>1] <- resi[sim>1]*(2**n)
  st <- (apply(sim, 1, sum)-1)*(2**n)
  
  return(sum(resi**2) + sum(st**2))
}

fit_CID_penalty <- function(params_ind, x, times, nseq){
  cat("fit CID\n")
  cat("  n=1, PSO\n", sep="")
  ctrl <- list(maxit=2000, trace=0, vectorize=TRUE, hybrid="off", hybrid.control=list(trace=1))
  res_opt <- psoptim(par=params_ind$val, fn=cost_labin_all, lower=params_ind$lower, upper = params_ind$upper, control=ctrl, x=x, times=times, n=1)
  
  res <- list("1"=res_opt)
  #optim(params_ind$val, cost_labin_all, lower=params_ind$lower, upper=params_ind$upper, method=c("L-BFGS-B"),control=list(trace=1,maxit=100), x=x, times=times, n=1)
  #psoptim(par=params_ind$val, fn=cost_labin_all, lower=params_ind$lower, upper = params_ind$upper, control=ctrl, x=x, times=times, n=1)
  if (nseq<1){
    stop("nseq should be > 0\n")
  }
  # create color palette
  col <- fun_col(ncol(x))
  for (n in seq(1,nseq)){
    cat("  n=", n, ", L-BFGS-B\n", sep="")
    params_ind$val <- res_opt$par
    # if (method=="BFGS"){
    res_n <- optim(params_ind$val, cost_labin_all, lower=params_ind$lower, upper=params_ind$upper, method=c("L-BFGS-B"), control=list(trace=0, maxit=500), x=x, times=times, n=n)
    res_n$times <- times
    res_n$meas <- x
    res_n$sim <- sim_all(res_n$par, x, times=times)
    # }else{
    #     cat("   ...PSO\n")
    #     res_n <- psoptim(par=params_ind$val, fn=cost_labin_all, lower=params_ind$lower, upper = params_ind$upper, control=ctrl, x=x, times=times, n=n)
    # }
    res[[n]] <- res_n
    #optim(params_ind$val, cost_labin_all, lower=params_ind$lower, upper=params_ind$upper, method=c("L-BFGS-B"),control=list(trace=1,maxit=100), x=x, times=times, n=n)
    #psoptim(par=params_ind$val, fn=cost_labin_all, lower=params_ind$lower, upper = params_ind$upper, control=ctrl, x=x, times=times, n=n)
    
    # plot simulations vs measurements
    matplot(x=res_n$times, y=res_n$meas,type="p", pch=20, main=paste("penalty:", n, sep=""), col=col, xlab = "time", ylab = "isotopologue abundance")
    matplot(x=res_n$times, y=res_n$sim, type = "l", lty=1, col=col, add = TRUE)
    for (j in seq(ncol(res_n$meas))){
      polygon(x=c(times, rev(times)),
              y=c(pmin(res_n$meas[,j], res_n$sim[,j]), rev(pmax(res_n$meas[,j], res_n$sim[,j]))),
              col=paste(col[j], "55", sep=""),
              border=NA)
    }
  }
  return(res)
}

get_anFun <- function(par, x){
  if (ncol(x)>1){
    cn <- colnames(x)
    anFun <- list()
    anFun[cn[1]] <- fn_M0_logistic_eq(par[seq(1,4)])
    for (i in seq(2,ncol(x))){
      #parx <- par[4+seq((i-2)*7+1, (i-1)*7)]
      anFun[cn[i]] <- fn_doublelogistic_eq(par[4+seq((i-2)*7+1, (i-1)*7)])
    }
  }else{
    anFun <- fn_doublelogistic_eq(par)
  }
  return(anFun)
}

createAnFun_ENR <- function(x, times, name, fun_dyn=fit_enr_log_doublelog){
  # transform to CIDs
  
  params_ind <- fit_ENR_ind(x, times, fun_dyn=fun_dyn)
  
  sim <- params_ind$funSet$sim(par = params_ind$resopt$par, times = seq(from=min(times), to=max(times), length.out=100))
  
  anFun <- list()
  anFunMx <- params_ind$funSet$eq(params_ind$resopt$par)
  lc <- str_count(name, "-") + 1
  anFun[[paste(name, "-M0", sep="")]] <- paste("1-(", anFunMx, ")", sep="")
  if (lc>1){
    for (i in seq(1, lc-1)){
      anFun[[paste(name, "-M", i, sep="")]] <- "0+0*t"
    }
  }
  anFun[[paste(name, "-M", lc, sep="")]] <- anFunMx
  
  return(list(anFun=anFun, sim=sim, exp=x))
}

createAnFun_CID <- function(x, times, nseq=2){
  params_ind <- fit_CID_ind(x, times)
  if (ncol(x)>1){
    params_all <- fit_CID_penalty(params_ind, x, times, nseq=nseq)
  }else{
    params_all <- list('1'=list(par=params_ind$val))
  }
  anFun <- get_anFun(params_all[[length(params_all)]]$par, x)
  sim <- sim_all(params_all[[length(params_all)]]$par, x, times=seq(0,100))
  
  
  return(list(anFun=anFun, sim=sim, exp=x))
}

fit_label_input <- function(x, t=1:nrow(x), net=NULL, file="res_label_input", subDir="fit_input", nseq=2, mc.cores=NULL, fun_dyn=fit_enr_log_doublelog){
  # times should not be NaNs
  if (any(is.na(t))){
    stop("some time values are NaNs")
  }
  # identify metabolites
  metab <- sapply(strsplit(colnames(x), "-M"), "[[", 1)
  # check measurements are sufficient to simulate label propagation through the network provided
  if (!is.null(net)){
    # if enrichments, else CIDs
    if (length(unique(metab)) == length(metab)){
      min_meas <- unique(sapply(strsplit(net$min_meas, "-M"), "[[", 1))
    }else{
      min_meas <- net$min_meas
    }
    ind <- (min_meas %in% colnames(x))
    if (!all(ind)){
      warning(paste(c("missing measurement to simulate label propagation:", min_meas[!ind]), collapse="\n"))
    }
  }
  # Fit CIDs of several metabolites, and plot meas. vs sim.)
  mainDir <- getwd()
  if (!file.exists(subDir)){
    dir.create(file.path(mainDir, subDir))
  }
  setwd(file.path(mainDir, subDir))
  # fit CID of each metabolite
  if (is.null(mc.cores)){
    lres <- lapply(as.list(unique(metab)), fit_input_job, data=x, t=t, metab=metab, file=file, nseq=nseq, fun_dyn=fun_dyn)
  }else{
    lres <- mclapply(as.list(unique(metab)), fit_input_job, data=x, t=t, metab=metab, file=file, nseq=nseq, fun_dyn=fun_dyn, mc.cores=mc.cores)
  }
  names(lres) <- unique(metab)
  # get analytical functions
  anFun <- sapply(lres, "[[", "anFun")
  # if 'net' is provided, automatically set the adequate analytical functions
  anFun_input <- list()
  if (!is.null(net)){
    anFun_input <- anFun4net(lres, net)
  }
  # return results
  setwd(mainDir)
  return(list(all=lres, anFun=anFun, anFun_input=anFun_input))
}

fit_input_job <- function(i, data, t, metab, file, nseq=NULL, fun_dyn=fit_enr_log_doublelog){
  cat(paste("fitting labeling dynamics of: ", i, sep=""), sep="\n")
  # create pdf
  pdf(paste(file, "_", i, ".pdf", sep=""))
  # fit experimental CID and get analytical functions
  if (sum(metab==i)>1){
    filt <- !(is.na(apply(data[,metab==i], 1, sum)))
    resin <- createAnFun_CID(data[filt, metab==i], times=t[filt], nseq=nseq)
  }else{
    filt <- !(is.na(data[,i]))
    resin <- createAnFun_ENR(t(data[filt, i]), times=t[filt], name=i, fun_dyn=fun_dyn)
  }
  # close pdf
  dev.off() 
  # save all results in a text file
  list2file(resin, file=paste(file, "_", i, ".txt", sep=""))
  return(resin)
}

anFun4net <- function(lres, net){
  # create the set of analytical functions
  anFun <- list()
  for (i in net$min_meas){
    metab <- unlist(strsplit(i, "-M"))[1]
    anFun[[i]] <- lres[[metab]][["anFun"]][[i]]
  }
  return(anFun)
}


# costScalar <- function(fp, optn){
#     fpar <- c(optn$kp,optn$meta_conc)[c(optn$npar,optn$nconc)]
#     fpar[c(optn$npar,optn$nconc)] <- fp
#     cos=cost(fpar,
#              data      = optn$data,
#              times     = optn$times,
#              kp        = optn$kp,
#              net       = optn$net,
#              nconc     = optn$nconc,
#              npar      = optn$npar,
#              meta_conc = optn$meta_conc,
#              iso_conc  = optn$iso_conc,
#              mode      = optn$mode,
#              anFun     = optn$anFun,
#              events    = optn$events,
#              trf       = optn$trf,
#              p         = optn$p,
#              atol      = optn$atol,
#              rtol      = optn$rtol)
#     return(sum(cos$res**2))  
# }

# costVector <- function(fp, optn){
#     fpar <- c(optn$kp,optn$meta_conc)[c(optn$npar,optn$nconc)]
#     fpar[c(optn$npar,optn$nconc)] <- fp
#     cos=cost(fpar,
#              data      = optn$data,
#              times     = optn$times,
#              kp        = optn$kp,
#              net       = optn$net,
#              nconc     = optn$nconc,
#              npar      = optn$npar,
#              meta_conc = optn$meta_conc,
#              iso_conc  = optn$iso_conc,
#              mode      = optn$mode,
#              anFun     = optn$anFun,
#              events    = optn$events,
#              trf       = optn$trf,
#              p         = optn$p,
#              atol      = optn$atol,
#              rtol      = optn$rtol)
#     return(cos$res)  
# }

# ### OTHER ALGORITHMS
# optn <- list(data      = data_meas,
#              times     = times,
#              kp        = kp,
#              net       = net,
#              to_est    = to_est,
#              meta_conc = meta_conc,
#              iso_conc  = NULL,
#              mode      = mode,
#              anFun     = anFun,
#              events    = NULL,
#              trf       = NULL,
#              p         = 0,
#              atol      = 1e-6,
#              rtol      = 1e-6,
#              nconc     = intersect(names(meta_conc), to_est),
#              npar      = intersect(names(kp), to_est))

# # Particle swarm
# #install.packages("pso")
# library("pso")

# tPSO <- system.time(
#     resPSO <- psoptim(par     = rep(1,length(to_est)),
#                       fn      = costScalar,
#                       lower   = lower,
#                       upper   = upper,
#                       control = list(maxit=10, s=5, abstol=1e-8, REPORT=1, trace=1, trace.stats=TRUE),
#                       optn    = optn)
# )

# # Evolution
# #install.packages("DEoptim")
# library("DEoptim")
# tDE <- system.time(
#     resDE <- DEoptim(fn      = costScalar,
#                      lower   = lower,
#                      upper   = upper,
#                      control = DEoptim.control(storepopfrom=1, itermax = 10, trace=TRUE),
#                      optn    = optn)
# )

# # Levenberg-Marquardt
# #install.packages("minpack.lm")
# library(minpack.lm)
# tLM <- system.time(
#     resLM <- nls.lm(c(1),
#                     lower   = lower,
#                     upper   = upper,
#                     fn      = costVector,
#                     jac     = NULL,
#                     control = nls.lm.control(nprint=1),
#                     optn)
# )

# # here L-BFGS
# # Note: optimx wraps several optimization methods:
# # "Nelder-Mead","BFGS","CG","L-BFGS-B","nlm","nlminb","spg","ucminf","newuoa","bobyqa","nmkb","hjkb","Rcgmin","Rvmmin"
# install.packages("optimx")
# library("optimx")
# tLBFGS <- system.time(
#     resLBFGS <- optimx(par     = rep(1,length(to_est)),
#                        fn      = costScalar,
#                        gr      = NULL,
#                        hess    = NULL,
#                        lower   = lower,
#                        upper   = upper,
#                        method  = c("L-BFGS-B"),
#                        itnmax  = NULL,
#                        hessian = FALSE,
#                        control = list(trace=1, maxit=10),
#                        optn    = optn)
# )

# # display loaded DLLs
# getLoadedDLLs()

####################################################
### CODE FOR PARALLELIZATION ON WINDOWS MACHINES ### 
####################################################

# from http://www.stat.cmu.edu/~nmv/2014/07/14/implementing-mclapply-on-windows
#
#' Define a sockets version of mclapply
#'
#' An implementation of \code{\link[parallel]{mclapply}} using \code{parallel::parLapply}
#'
#' Windows does not support forking. This makes it impossible to use mclapply on Windows to
#' farm out work to additional cores.
#'
#' @param ... What you pass to mclapply
#' @return mclapply like list
#' @import parallel
#' @export
mclapply_socket <- function(
  X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = NULL,
  mc.cleanup = TRUE, mc.allow.recursive = TRUE
) {
  ## Create a cluster
  if (is.null(mc.cores)) {
    mc.cores <- min(length(X), detectCores())
  }else{
    mc.cores <- min(length(X), detectCores(), mc.cores)
  }
  cl <- parallel::makeCluster( mc.cores )
  
  tryCatch( {
    ## Find out the names of the loaded packages
    loaded.package.names <- c(
      ## Base packages
      sessionInfo()$basePkgs,
      ## Additional packages
      names( sessionInfo()$otherPkgs ))
    
    ### Ship it to the clusters
    parallel::clusterExport(cl,
                            'loaded.package.names',
                            envir=environment())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parallel::parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        require(yy , character.only=TRUE)})
    })
    
    clusterExport_function(cl, FUN)
    
    ## Run the lapply in parallel, with a special case for the ... arguments
    if( length( list(...) ) == 0 ) {
      return(parallel::parLapply( cl = cl, X=X, fun=FUN) )
    } else {
      return(parallel::parLapply( cl = cl, X=X, fun=FUN, ...) )
    }
  }, finally = {
    ## Stop the cluster
    parallel::stopCluster(cl)
  })
}

#' Overwrite the serial version of mclapply on Windows only
#'
#' @param empty it takes nothing
#' @return mclapply like list
#' @export
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {mclapply_socket},
                    Linux   = {mclapply},
                    Darwin  = {mclapply})

clusterExport_function <- function(cl, FUN ) {
  ## We want the enclosing environment, not the calling environment
  env <- environment(FUN)
  while(!identical(env, globalenv())) {
    env <- parent.env(env)
    parallel::clusterExport(cl, ls(all.names=TRUE, envir = env), envir = env)
  }
  parallel::clusterExport(cl, ls(all.names=TRUE, envir = env), envir = env)
}
