# 2014-20-12 millard@insa-toulouse.fr
#
# hist_detailed() generates detailed histograms based on names vectors.
#
# Copyright 2014, INRA, France
# License: GNU General Public License v2 (see license.txt for details)

require("RColorBrewer")
fun_col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

cnvrt.coords <-function(x, y=NULL){
# Stolen from the teachingDemos library, simplified for this use case
    xy    <- xy.coords(x, y, recycle=TRUE)
    cusr  <- par('usr')
    cplt  <- par('plt')
    plt   <- list()
    plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
    plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
    fig   <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
    return(list(fig=fig))
}

subplot <- function(fun, x, y=NULL, ax=NULL, tickPretty=FALSE){
# Stolen from the teachingDemos library, simplified for this use case
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    xy <- xy.coords(x, y)
    xy <- cnvrt.coords(xy)$fig
    par(plt=c(xy$x, xy$y), new=TRUE)
    fun
    if (!is.null(ax)){
        if (tickPretty){
            axis(1, at=axTicks(1), labels=ax$bins[axTicks(1) + 1], las=ax$las)
        }else{
            axis(1, at=ax$at, labels=ax$bins, las=ax$las)
        }
    }
    tmp.par <- par(no.readonly=TRUE)
    return(invisible(tmp.par))
}

barplot_gap <- function(r, gap, lum, las, xlab, ylab, main, npbins, bins, tickPretty=FALSE){
    lower   <- c(0, gap[1])
    upper   <- c(gap[2], max(colSums(r)) * 1.1)
    y_outer <- 0.21 * (lower[2] - lower[1] + upper[2] - upper[1])
    lowspan <- c(0, 0.1 * (lower[2] - lower[1] + upper[2] - upper[1]) + 1)
    topspan <- c(lowspan[2] + 1, y_outer)
    # initialize the plot area
    plot(c(0, 1), c(0, y_outer), type='n', axes=FALSE, ylab=ylab, xlab=xlab)
    subplot(barplot(r, space=0, col=fun_col(lum), ylim=lower, xpd=FALSE, las=las), x=c(0,1), y=lowspan, ax=list("at"=seq(0, npbins), "bins"=bins, "las"=las), tickPretty=tickPretty)
    subplot(barplot(r, space=0, las=las, main=main, col=fun_col(lum), ylim=upper, xpd=FALSE), x=c(0, 1), y=topspan, tickPretty=tickPretty)
    # all the following plots are in units of the outer coordinate system
    lowertop     <- lowspan[2] + 0.1        # where to end the lower axis
    breakheight  <- 0.5                    # height of the break
    upperbot     <- lowertop + breakheight # Where to start the upper axes
    markerheight <- 0.4                    # height difference for the break markers
    markerwidth  <- 0.04                # with of the break markers
    # draw the break markers
    lines(c(0, 0), c(1, lowertop))
    lines(c(markerwidth/-2, markerwidth/2), c(lowertop-markerheight/2, lowertop+markerheight/2))
    lines(c(0, 0), c(upperbot, 14))
    lines(c(markerwidth/-2, markerwidth/2), c(upperbot-markerheight/2, upperbot+markerheight/2))
}

hist_detailed <- function(x, nbin=20, las=2, plot_leg=FALSE, norm=TRUE, xlab=NULL, ylab=NULL, main=NULL, gap=NULL, leg_title=NULL, cex.legend=1, draw=TRUE, xlim=NULL, tickPretty=FALSE, btn=TRUE){
    if (is.null(xlim)){
        rang <- pretty(x, nbin)
        bins <- seq(min(rang), max(rang), length.out=nbin+1)
    }else{
        bins <- seq(xlim[1], xlim[2], length.out=nbin+1)
    }
    npbins <- length(bins) - 1
    if (!is.null(names(x))){
        meta <- sapply(strsplit(names(x), "_"), "[[", 1)
    }else{
        meta <- rep("unknown", length(x))
    }
    um  <- unique(meta)
    lum <- length(um)
    r   <- matrix(0, nrow=lum, ncol=npbins, dimnames=list("metabolite"=um, "bin"=NULL))
    for (m in um){
        xm <- x[meta == m]
        for (b in seq(npbins)){
            r[m, b] <- sum(xm >= bins[b] & xm < bins[b+1])
        }
    }
    if (norm){
        r <- r/sum(r) * 100
    }
    ymax <- max(colSums(r))
    if (draw){
        if (is.null(gap)){
            barplot(r, space=0, las=las, col=fun_col(lum), xlab=xlab, ylab=ylab, main=main, ylim=c(0, ymax * 1.1))
            if (tickPretty){
                axis(1, at=axTicks(1), labels=bins[axTicks(1) + 1], las=las)
            }else{
                axis(1, at=seq(0, npbins), labels=bins, las=las)
            }
            if (btn){
                box()
            }
        }else{
            if (ymax > gap[2]){
                barplot_gap(r=r, gap=gap, lum=lum, las=las, xlab=xlab, ylab=ylab, main=main, npbins=npbins, bins=bins, tickPretty=tickPretty)
            }else{
                stop(paste("'gap' out of range (max(gap) > max(y), with max(gap)=", gap[2], " and max(y)=", max(colSums(r)), ")", sep=""))
            }
        }
        if (plot_leg){
            plot(0, type='n', axes=FALSE, ann=FALSE)
            legend("topleft", legend=um, fill=fun_col(lum), cex=cex.legend, title=leg_title)
        }
    }
    return(invisible(list(frequency=r, bins=bins)))
}


hist_fixed <- function(x, nbin=20, las=2, xlab=NULL, ylab=NULL, main=NULL, xlim=NULL, btn=TRUE, labels=NULL, at=NULL, col="steelblue",border="steelblue4"){
    if (is.null(xlim)){
        rang <- pretty(x, nbin)
        bins <- seq(min(rang), max(rang), length.out=nbin+1)
    }else{
        bins <- seq(xlim[1], xlim[2], length.out=nbin+1)
    }
    npbins <- length(bins) - 1
    meta <- rep("unknown", length(x))
    um  <- unique(meta)
    lum <- length(um)
    r   <- matrix(0, nrow=lum, ncol=npbins, dimnames=list("metabolite"=um, "bin"=NULL))
    for (m in um){
        xm <- x[meta == m]
        for (b in seq(npbins)){
            r[m, b] <- sum(xm >= bins[b] & xm < bins[b+1])
        }
    }
    ymax <- max(colSums(r))
    barplot(r, space=0, las=las, col=col, xlab=xlab, ylab=ylab, main=main, border=border, ylim=c(0, ymax * 1.1))
    if (btn){
        box()
    }
    axis(1, at=at, labels=labels)
    return(invisible(list(frequency=r, bins=bins)))
}
