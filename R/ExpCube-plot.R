##
## ExpCube analysis
## Author: Tomasz Konopka
##
## Some plot functions (uses Rcssplot package)
##





##' Draw a scatter plot, including line of best fit, correlation measures and labels
##' (using Rcssplot)
##'
##' @param xx - numeric vector with names. x coordinates for points.
##' @param yy - numeric vector with names. y coordinates for points
##' @param items - character vector. ids of data points to display
##' @param outliers - character vector. ids of data points drawn as outliers.
##' @param col -
##' @param onecol -
##' @param items.highlight - character vector. Names of items x,y to highlight.
##' @param xlim - numeric vector with two elements. range for x axis
##' @param ylim - numeric vector with two elements. range for y axis
##' @param xlimwiden - numeric. Determines how much wider the actual x range will be relative
##' to xlim. Use this to add some padding around the xlim. 
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title, above heatmap.
##' @param labelq - numeric vector of length two. Quantiles that determine which
##' items should be highlihted.
##' @param cex.rescale - numeric. Rescaling factor for items marked for highlighting.
##' @param show.correlation - logical. Toggle display of Spearman correlation info.
##' @param show.glm - logical. Toggle display of best fit line.
##' @param show.labels - character vector. Names of items to label on the plot.
##' @param correlation.threshold - numeric. Correlations with p-value below this threshold
##' will be highlighted in bold.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
E3PlotScatter = function(xx, yy, items=NULL, outliers=NULL, col=NULL, onecol=NULL,
    items.highlight=c(),
    xlim=NULL, ylim=NULL, xlimwiden=0.1,    
    xlab=NULL, ylab=NULL, main="",
    labelq=c(0.02, 0.98), cex.rescale=NULL,
    show.correlation = TRUE, show.glm = FALSE, show.labels = c(),
    correlation.threshold=1e-4,
    Rcss="default", Rcssclass=c()) {
    
    if (is.null(items)) {
        items = names(xx)
    }
    
    ## just look at requested items 
    xx = xx[items]
    yy = yy[items]

    ## for colors, either request a vector wiht one color per item,
    ## or force a single color on all the points
    if (!is.null(onecol)) {
        col = rep(onecol, length(xx))
        names(col) = names(xx)
    }
    col = col[items]

    ## allow user to specify a rescaling factor for each item, to make
    ## some points appear larger than others
    if (!is.null(cex.rescale)) {
        cex.rescale = cex.rescale[items]
    } else {
        cex.rescale = 1;
    }
    
    ## shorthand for coding
    if (identical(Rcss, "default")) {
        RC = RcssGetDefaultStyle()
    } else {
        RC=Rcss
    }
    RCC=Rcssclass
    
    ## Estimate ranges for axis (if this doesn't work, compute the outside this plot function)
    if (is.null(xlim)) {
        xlim = range(xx[is.finite(xx)])
        if (!is.finite(xlim[1]) | !is.finite(xlim[2])) {
            xlim = c(0,1)
        }
        xlimwiden = (xlim[2]-xlim[1])*xlimwiden/2;
        xlim = xlim+c(-xlimwiden, xlimwiden)
    }    
    if (is.null(ylim)) {
        ylim = c(0, max(yy[is.finite(yy)])*1.2);
    }

    ## start a new scatter plot
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xlab="", ylab="", 
         xaxs="i", yaxs="i", type="n", frame=F, axes=F)
    
    ## display some axes
    E3Axes(xlim, ylim, xlab=xlab, ylab=ylab, RC=RC, RCC=RCC)
    
    ## evaluate and display a loess fit
    if (show.glm) {
        if (length(xx[is.finite(xx)])>1 & length(yy[is.finite(yy)])>1) {
            xydata = data.frame(x=xx, y=yy)[order(xx),]
            loessy = xydata[,"y"]
            loessx = xydata[,"x"]
            xyglm = glm(loessy~loessx)
            plotglmy = predict(xyglm, data.frame(loessx))
            plotxy = cbind(loessx, plotglmy)
            plotxy = plotxy[plotxy[,2]>ylim[1] & plotxy[,2]<ylim[2],]
            Rcsslines(plotxy[,1], plotxy[,2], Rcss=RC, Rcssclass=c(Rcssclass, "glm"));
        }
    }
    
    basecex = RcssGetPropertyValueOrDefault(RC, "points", "cex", default=1, Rcssclass=RCC)
    highcex = RcssGetPropertyValueOrDefault(RC, "points", "cex", default=1,
        Rcssclass=c(RCC, "highlight"))
    
    items.basic = items[!(items %in% items.highlight)]
    
    Rcsspoints(xx[items.basic], yy[items.basic], cex=basecex*cex.rescale, col=col[items.basic],
               Rcss=RC, Rcssclass=RCC)
    if (length(items.highlight)>0) {
        Rcsspoints(xx[items.highlight], yy[items.highlight],
                   cex=highcex*cex.rescale, col=col[items.highlight],                   
                   Rcss=RC, Rcssclass=c(RCC, "highlight"))
        Rcsstext(xx[items.highlight], yy[items.highlight], labels=items.highlight,
                 Rcss=RC, Rcssclass=c(RCC, "highlight"))
    }
    
    
    if (length(show.labels)>0) {
        Rcsstext(xx[names(show.labels)], yy[names(show.labels)], labels=show.labels,
                 col=col[names(show.labels)], Rcss=RC, Rcssclass=RCC)
    }
    
    ## evaluate Spearman correlation and display in the corner
    if (show.correlation) {
        if (length(xx[is.finite(xx)])>1 & length(yy[is.finite(yy)])>1) {
            temp = cor.test(xx, yy, method="spearman");
            spval = signif(temp$p.value,2);
            srho = signif(temp$estimate,2);        
        } else {
            temp = list(p.value="NA", estimate="NA");
            spval = "NA";
            srho = "NA";        
        }
        legtext = paste0("Spearman rho: ", srho,"\np: ", spval)
        if (spval<correlation.threshold) {
            Rcssmtext(legtext, side=3, Rcss=RC, font=2, col="#ff0000",
                      Rcssclass=c(RCC, "legend"));
        } else {
            Rcssmtext(legtext, side=3, Rcss=RC, Rcssclass=c(RCC, "legend"));
        }
    }
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    
    return(plotxy)
    
    ## and that's all!
}







##' Draw a boxplot (using Rcssplot)
##' 
##' @param vlist - names list of numeric vectors. Analogous to input to boxplot()
##' @param quantiles - numeric vector of five elements. The five elements will
##' be treated as quantiles for lower-whisker, lowe-box, middle, upper-box, upper-whisker
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title, above heatmap.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotBoxes = function(vlist, quantiles=c(0, 0.25, 0.5, 0.75, 1),
    xlab="", ylab="", main="", 
    Rcss="default", Rcssclass="plotbox") {
    
    if (length(quantiles)!=5) {
        stop("vector quantiles must have five numers\n")
    }
    RC = Rcss
    RCC = Rcssclass

    ## find the space between the boxes
    spacer = RcssGetPropertyValueOrDefault(RC, "plotbox", "spacer",
        Rcssclass=RCC, default=0.2)
    HS = spacer/2;
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)
    
    xlim = c(-HS, length(vlist)+HS)
    ylim = c(0, max(unlist(vlist))*1.15)

    plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n", 
         frame=F, axes=F, xlab="", ylab="")

    ## display some axes
    E3Axes(xlim, ylim, xlab=xlab, ylab=ylab, RC=RC, RCC=RCC, where="y")    
    Rcssaxis(1, at=xlim, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(1, at=seq(0.5, length(vlist), 1), labels=names(vlist), tck=0, 
             Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))    
    
    for (i in 1:length(vlist)) {
        nowqs = quantile(vlist[[i]], p=quantiles)
        Rcssrect(i-1+HS, nowqs[2], i-HS, nowqs[4], Rcss=RC, Rcssclass=RCC)
        Rcsslines(rep(i-0.5, 5), c(nowqs[1:2], NA, nowqs[4:5]), Rcss=RC, Rcssclass=RCC)
        Rcsslines(c(i-1+HS, i-HS), rep(nowqs[3],2), Rcss=RC, Rcssclass=c(RCC, "median"))
    }
    
}





##' Draw a histogram (using Rcssplot)
##' 
##' @param histobj - histogram object
##' @param xlim - numeric vector of length 2. Determines range to plot on x axis
##' @param ylim - numeric vector of length 2. Determines range to plot on y axis
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title, above heatmap.
##' @param fill - logical. Set TRUE to obtain a fill color under the histogram line.
##' @param xat - numeric vector. Positions of labels on the x axis.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotHist = function(histobj, xlim=NULL, ylim=NULL, 
    xlab="", ylab="", main="",  fill=TRUE, xat=NULL, 
    Rcss="default", Rcssclass="hist") {
    
    if (is.null(xlim)) {
        xlim = range(histobj$breaks)
    }
    if (is.null(ylim)) {
        ylim = c(0, max(histobj$density))
    }
    RC = Rcss
    RCC = Rcssclass
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n", 
         frame=F, axes=F, xlab="", ylab="")
    
    ## 
    nowxvals = rep(histobj$breaks, each=2)
    nowxvals = nowxvals[2:(length(nowxvals)-1)]
    nowyvals = rep(histobj$density, each=2)    
    ## perhaps draw a filling polygon
    if (fill) {
        Rcsspolygon(nowxvals, nowyvals, xpd=FALSE, Rcss=RC, Rcssclass=RCC)
    }
    Rcsslines(nowxvals, nowyvals, xpd=FALSE, Rcss=RC, Rcssclass=RCC)
    
    ## display some axes
    Rcssaxis(2, at=ylim, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "y"), line=0)    
    Rcssaxis(1, at=xlim, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "x"))
    if (is.null(xat)) {
        Rcssaxis(1, Rcss=RC, Rcssclass=c(RCC, "x"))
    } else {
        Rcssaxis(1, at=xat, Rcss=RC, Rcssclass=c(RCC, "x"))
    }
    Rcssmtext(xlab, side=1, Rcss=RC, Rcssclass=c(RCC, "x"));
    Rcssmtext(ylab, side=2, Rcss=RC, Rcssclass=c(RCC, "y"));
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
}





##' Draw series data as points with smoothed lines (using Rcssplot)
##'
##' @param vdat - data.frame with numbers
##' @param xcolumn - character/integer. Column containing x values 
##' @param ycolumns - vector of of columns to show as lines
##' @param lineclass - character vector. Classes used to style each lines using Rcssplot. 
##' @param lineloess - logical vector. Determines whether for each line, the line will
##' be plotted smoothed (loess) or not.
##' @param xlim - numeric vector with two elements. range for x axis
##' @param ylim - numeric vector with two elements. range for y axis
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title.
##' @param legend - character vector. Text to display in legend
##' @param legend.at - numeric vector of length equal to ycolumns. Determines where
##' to place the legends.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotLines = function(vdat, xcolumn=1, ycolumns=c(2,3),
    lineclass=rep("plain", length(ycolumns)), lineloess=rep(TRUE, length(ycolumns)),
    xlim = NULL, ylim=NULL,
    xlab="", ylab="", main="", 
    legend=rep("", length(ycolumns)), legend.at = rep(1, length(ycolumns)),
    Rcss="default", Rcssclass="plotline") {
    
    RC = Rcss
    RCC = Rcssclass
    
    ## find the space between the boxes
    if (is.null(xlim)) {
        xlim = c(0, max(vdat[,xcolumn]))
    }
    if (is.null(ylim)) {
        ylim = c(0, max(vdat[,ycolumns])*1.15)
    }
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)
    
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n", 
         frame=F, axes=F, xlab="", ylab="")
    
    ## display some axes
    E3Axes(xlim, ylim, xlab=xlab, ylab=ylab, RC=RC, RCC=RCC, where=c("x", "y"))    
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## draw the individual data series
    for (i in 1:length(ycolumns)) {
        ## check if draw a loess fit for the points
        if (lineloess[i]) {
            tempy = vdat[,ycolumns[i]]
            tempx = vdat[,xcolumn]
            temploess = loess(tempy~tempx, span=0.5)
            predx = seq(xlim[1], xlim[2], length=64)
            loessy = as.numeric(predict(temploess, data.frame(tempx=predx)))
            Rcsslines(predx, loessy, Rcss=RC, Rcssclass=c(RCC, lineclass[i]))
        } else {
            Rcsslines(vdat[,xcolumn], vdat[,ycolumns[i]], Rcss=RC, Rcssclass=c(RCC, lineclass[i]))
        }
        ## beside the fit line, also draw the raw data points
        Rcsspoints(vdat[,xcolumn], vdat[,ycolumns[i]], Rcss=RC, Rcssclass=c(RCC, lineclass[i]))
        ## draw a label for each series
        Rcsstext(vdat[legend.at[i], xcolumn],
                 vdat[legend.at[i], ycolumns[i]], labels=legend[i],
                 Rcss=RC, Rcssclass=c(RCC, lineclass[i]))
    }
    
}












##' Draws a heatmap representation of a genomic segmentation (using Rcssplot)
##' 
##' @param seglist - list of cbs segmentation data frames
##' @param chrlens - named vector with lengths of chromosomes (chromosomes
##' will appear in this order)
##' @param categories - a list splitting the groups into categories (eg. plates)
##' groups on the vertical axis will be grouped accordingly
##' @param heatcols -
##' @param label.chr - names of chromosomes to label on horizontal axis
##' @param seg.rescale - numeric. Used to modulate color intensity 
##' @param seg.hlab - numeric. Determines spacing between labels near x axis.
##' @param seg.vlab - numeric. Determines spacing between labels near y axis.
##' @param legend - numeric. Values to use on the legend box.
##' @param axes.labels - character vector of length two. Text to write near
##' y axis to the left and to the right of the heatmap.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotSegmentationOverview = function(seglist, chrlens=NULL,
    categories=list(All=names(seglist)), 
    heatcols=c("#00ffff", "#ff00ff"), label.chr=c(), 
    seg.rescale = 2, seg.hlab=0.3, seg.vlab=1e6,
    legend=NULL, axes.labels=NULL,
    Rcss="default", Rcssclass="segmentation") {
    
    RC=Rcss;
    RCC=Rcssclass;

    ## get length of genome
    genome.length = sum(chrlens)
    xlim = c(0, sum(chrlens))
    ylim = c(0, length(seglist))
    ## learn about the chromosome coordinates
    chrends = cumsum(chrlens)
    chrstarts = c(0, chrends)
    names(chrends) = names(chrlens)
    names(chrstarts) = names(chrlens)
    
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
         xlab="", ylab="", frame=F, axes=F)

    if (seg.hlab<0) {
        seg.hlab = (ylim[2]-ylim[1])*(-seg.hlab)
    }
    
    ## plot chromosome legends on horizontal axis    
    for (i in 1:length(chrlens)) {
        chrname = names(chrlens)[i]
        if (i%%2==0) {
            xdown = -2*seg.hlab;
        } else {
            xdown = -3*seg.hlab;
        }
        Rcssrect(chrstarts[i], xdown, chrends[i], xdown-seg.hlab,
                 Rcss=RC, Rcssclass=c(RCC, "chroms"))
        if (chrname %in% label.chr) {
            Rcsstext(mean(chrstarts[c(i, i+1)]), -4*seg.hlab, chrname,
                     pos=1, Rcss=RC, Rcssclass=RCC)
        }                
    }

    ## plot plate legends on vertical axis (right hand side)
    platedividers = c(0, cumsum(sapply(categories, length)))
    temp = rev(1:length(categories))
    if (is.null(axes.labels))  {
        for (i in 1:length(categories)) {
            catname = names(categories)[i]
            j = temp[i]
            Rcsstext(genome.length, mean(platedividers[c(i, i+1)]), catname,
                     pos=4, srt=0, Rcss=RC, Rcssclass=RCC)                         
        }
    }
    
    ## for each segmented sample...
    for (i in 1:length(seglist)) {
        nowseg = seglist[[i]]
        ## extract only those segments that are abnormal (segmentation mean > threshold)
        nowseg[,"col"] = E3Val2Col(nowseg[,"seg.mean"]/2, col=heatcols)
        for (j in 1:nrow(nowseg)) {
            nowchr = as.character(nowseg[j, "chrom"])
            nowstart = chrstarts[nowchr]
            ##cat(nowchr, " ", nowseg[j,"loc.start"], " ", nowseg[j, "loc.end"], "\n")
            Rcssrect(nowstart + nowseg[j, "loc.start"], i-1,
                     nowstart + nowseg[j, "loc.end"], i, col=nowseg[j, "col"],
                     Rcss=RC, Rcssclass=c(RCC, "heatmap"))                                        
        }

        if (is.null(axes.labels)) {
            Rcsstext(0, i-0.6, names(seglist)[i], pos=2, Rcss=RC,  Rcssclass=c(RCC, "cells"))
        }
    }
    
    ## Draw a color-code legend
    legend = E3Val2Col(log2(legend), col=heatcols)
    E3DrawScaleLegend(legend=legend, xlim=xlim, ylim=ylim, RC=RC, RCC=RCC)
    
    ## draw some dividing lines?
    temp = unique(c(0, cumsum(sapply(categories, length))))
    for (i in temp) {
        Rcsslines(c(0, genome.length), y=rep(i, 2),
                  Rcss=RC, Rcssclass=c(RCC, "divider"))
    }
    for (i in 1:length(chrstarts)) {
        Rcsslines(rep(chrstarts[i], 2), c(0, length(seglist)),
                  Rcss=RC, Rcssclass=c(RCC, "divider"))
    }
    
    ## draw a box around the whole plot
    Rcsslines(c(xlim, rev(xlim), xlim[1]), c(rep(ylim, each=2), ylim[1]),
              Rcss=RC, Rcssclass=c(RCC, "box"))
        
    Rcssmtext("Regional expression changes relative to WT", side=3,
              Rcss=RC, Rcssclass=c(RCC, "main"))
    if (!is.null(axes.labels)) {
        Rcssmtext(axes.labels[1], side=2, Rcss=RC, Rcssclass=c(RCC, "y"))
        Rcssmtext(axes.labels[2], side=4, Rcss=RC, Rcssclass=c(RCC, "y4"))
    }  
}





##' Draws a heatmap representation of regional genomic patterns (using Rcssplot)
##'
##' The function takes numeric data associated with genes, then performs a smoothing using
##' bins and plots the resulting smoothed signal.
##' 
##' @param fcdata - numeric matrix with gene/feature names in rownames and samples in colnames.
##' This data will be log2 transformed before the plotting takes place
##' @param gene.positions - data frame linking gene names to gene positions
##' @param categories - named list of vectors, splitting samples into groups
##' (e.g. Plates or batches). This object also determines the order of the samples on
##' the y axis
##' @param label.chr - names of chromosomes to label on horizontal axis (
##' @param chrlens - named vector with lengths of chromosomes (chromosomes will appear
##' in this order)
##' @param heatcols - numeric vector of lengths two. Determines the color scale of the plot.
##' @param seg.rescale - numeric. Used to modulate color intensity 
##' @param seg.resolution - numeric. Width of bins to use for smoothing
##' @param seg.hlab - numeric. Determines spacing between labels near x axis.
##' @param seg.vlab - numeric. Determines spacing between labels near y axis.
##' @param legend - numeric. Values to use on the legend box.
##' @param axes.labels - character vector of length two. Text to write near
##' y axis to the left and to the right of the heatmap.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
E3PlotSmoothGenomeOverview = function(fcdata, gene.positions=NULL,
    categories=list(All=colnames(fcdata)),  label.chr=c(), 
    chrlens=NULL, heatcols=c("#00ffff", "#ff00ff"), 
    seg.rescale = 2, seg.resolution=5e6, seg.hlab=0.3, seg.vlab=1e6,
    legend=NULL, axes.labels=NULL, 
    Rcss="default", Rcssclass="diamond") {
    
    RC=Rcss;
    RCC=Rcssclass;
    
    ## get length of genome
    genome.length = sum(chrlens)
    xlim = c(0, sum(chrlens))
    ylim = c(0, ncol(fcdata))
    ## learn about the chromosome coordinates
    chrends = cumsum(chrlens)
    chrstarts = c(0, chrends)
    names(chrends) = names(chrlens)
    names(chrstarts) = names(chrlens)
    ## split up each chromosome into windows
    chrwindows = list()
    for (nowchr in names(chrlens)) {
        chrwindows[[nowchr]] = c(seq(0, chrlens[nowchr], by=seg.resolution), chrlens[nowchr])
    }
    
    ## adjust fcdata so that it is log transformed
    gene.positions = gene.positions[!is.na(gene.positions[,"pos"]),]    
    fcdata = log2(fcdata[rownames(gene.positions),]);
    fcdata[!is.finite(fcdata)] = 1;
    ## get positions of genes per chromosome
    gene.positions = split(gene.positions, gene.positions[,"chr"])[names(chrlens)]
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
         xlab="", ylab="", frame=F, axes=F)

    if (seg.hlab<0) {
        seg.hlab = (ylim[2]-ylim[1])*(-seg.hlab)
    }
    
    ## plot chromosome legends on horizontal axis    
    for (i in 1:length(chrlens)) {
        chrname = names(chrlens)[i]
        if (i%%2==0) {
            xdown = -2*seg.hlab;
        } else {
            xdown = -3*seg.hlab;
        }
        Rcssrect(chrstarts[i], xdown, chrends[i], xdown-seg.hlab,
                 Rcss=RC, Rcssclass=c(RCC, "chroms"))
        if (chrname %in% label.chr) {
            Rcsstext(mean(chrstarts[c(i, i+1)]), -4*seg.hlab, chrname,
                     pos=1, Rcss=RC, Rcssclass=RCC)
        }          
    }

    ## plot plate legends on vertical axis (right hand side)
    platedividers = c(0, cumsum(sapply(categories, length)))
    temp = rev(1:length(categories))
    if (is.null(axes.labels)) {
        for (i in 1:length(categories)) {
            catname = names(categories)[i]
            j = temp[i]
            Rcsstext(genome.length, mean(platedividers[c(i, i+1)]), catname,
                 pos=4, srt=0, Rcss=RC, Rcssclass=RCC)                       
        }
    }
        
    ## for each segmented sample, for each chromosome, compute and draw smooth FC
    for (i in 1:ncol(fcdata)) {        
        for (nowchr in names(chrlens)) {
            
            ## create a matrix of fold changes in one sample and in one
            nowdata = data.frame(FC=fcdata[gene.positions[[nowchr]][,"Gene"], i],
                pos= gene.positions[[nowchr]][,"pos"])
            nowstart = chrstarts[nowchr]
            
            ## plot one box per window
            nowbreaks = chrwindows[[nowchr]];
            kkmin = 1;            
            for (kk in 1:(length(nowbreaks)-1)) {
                nowhits = nowdata[,"pos"]>= nowbreaks[kkmin] & nowdata[,"pos"]<nowbreaks[kk+1];
                if (sum(nowhits)>10) {
                    nowfc = median(nowdata[nowhits, "FC"])                    
                    if (nowfc>1) {
                        nowfc=1;
                    } else if (nowfc<(-1)) {
                        nowfc=-1;
                    }
                    nowcol = E3Val2Col(nowfc*nowfc*sign(nowfc)/seg.rescale, col=heatcols)
                    
                    Rcssrect(nowstart+nowbreaks[kkmin], i-1, nowstart+nowbreaks[kk+1], i,
                             col=nowcol, Rcss=RC, Rcssclass=c(RCC, "heatmap"))
                    
                    kkmin = kk+1;
                } 
            }                        
        }

        if (is.null(axes.labels)) {
            Rcsstext(0, i-0.6, colnames(fcdata)[i], pos=2, Rcss=RC, Rcssclass=c(RCC, "cells"))
        }
    }
    
    ## Draw a color-code legend
    legend = E3Val2Col(log2(legend), col=heatcols)
    E3DrawScaleLegend(legend=legend, xlim=xlim, ylim=ylim, RC=RC, RCC=RCC)
    
    ## draw some dividing lines?
    temp = unique(c(0, cumsum(sapply(categories, length))))
    for (i in temp) {
        Rcsslines(c(0, genome.length), y=rep(i, 2),
                  Rcss=RC, Rcssclass=c(RCC, "divider"))
    }
    for (i in 1:length(chrstarts)) {
        Rcsslines(rep(chrstarts[i], 2), c(0, ncol(fcdata)),
                  Rcss=RC, Rcssclass=c(RCC, "divider"))
    }

    ## draw a box around the whole plot
    Rcsslines(c(xlim, rev(xlim), xlim[1]), c(rep(ylim, each=2), ylim[1]),
              Rcss=RC, Rcssclass=c(RCC, "box"))
    
    Rcssmtext("Regional expression changes relative to WT", side=3,
              Rcss=RC, Rcssclass=c(RCC, "main"))
    if (!is.null(axes.labels)) {
        Rcssmtext(axes.labels[1], side=2, Rcss=RC, Rcssclass=c(RCC, "y"))
        Rcssmtext(axes.labels[2], side=4, Rcss=RC, Rcssclass=c(RCC, "y4"))
    }    
}









##' Draw a complex diagram containing a tree/diamond clustering and a barplot of signature size
##' (using Rcssplot)
##'
##' This function produces a figure in two panels. On the left is a diamond clustering.
##' On the right is barplot showing the signature size for each element in the clustering.
##' 
##' @param gjh - hclust object. Used to produce diamond clustering on LHS.
##' @param colmat - matrix of colors (or values in range [0,1])
##' @param gsize - numeric vector with names. Used to produce barplot on RHS. 
##' @param n - integer. Sets vertical scale. (Set equal or more than number of elements
##' in clustering. Fixing this value helps making multiple clusterings with the same
##' box sizes)
##' @param spacer - numeric. Determines empty space between bars on RHS.
##' @param lab - character string. Text to write next to clustering on LHS.
##' @param col.group - named vector of colors. Associated a color with each element in the
##' clustering. 
##' @param highlight - character vector. Names of items to highlight on RHS.
##' @param mai - numeric vector of length four. Sets margin. 
##' @param maxlim - numeric. Determines horizontal scale of barplot. 
##' @param xticks - numeric vector. Values to label on x-axis in the barplot.
##' @param lineext - numeric. determines lenght of line labeling highlighted items.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
E3PlotStimuliDiamond = function(gjh, colmat, gsize, n=10, spacer=0.1, lab="All", 
    col.group=NULL, highlight=c(), mai=rep(0.2, 4), maxlim=NULL, xticks=seq(0, 300, 100),
    lineext=1.2,
    Rcss="default", Rcssclass="diamond") {
    
    RC=Rcss;
    RCC=Rcssclass;

    rotate45 = function(xy, aa=-pi/4) {
        ans = xy;
        ans[1] = (xy[1]*cos(aa)) - (xy[2]*sin(aa))
        ans[2] = (xy[1]*sin(aa)) + (xy[2]*cos(aa))
        return(ans)
    }
    
    Rcsspar(mai=mai, mfrow=c(1,2))
    
    ## draw the diamon on the left hand side
    xlim = c(-n/2, 0)
    xlim = xlim + 0.5;
    ylim=c(0, n)-0.25;
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, axes=F, frame=F, 
         xaxs="i", yaxs="i", xlab="", ylab="")
    
    ## plot the diamond
    ngroups = length(col.group);
    
    ## plot the individual boxes based on gjh
    for (i in 1:ngroups) {
        groupi = gjh$labels[gjh$order[i]]
        for (j in 1:ngroups) {
            groupj = gjh$labels[gjh$order[j]]
            ## create a color based on overlap and group colors
            nowval = min(1, 2*(1-colmat[groupi, groupj]))
            nowcol = paste0(avgcol(col.group[c(groupi, groupj)]), val2hex(nowval))
            ## plot a marker on the heatmap
            if (j<=i & nowcol!="#ffffff") {
                ##cat(groupi, " ", groupj, " ", nowcol, "\n")
                nowpolygon = rbind(c(-i, j-1), c(-i+1, j-1), c(-i+1, j), c(-i, j))
                
                for (kk in 1:4) {
                    nowpolygon[kk, ] = rotate45(nowpolygon[kk,])
                }
                Rcsspolygon(nowpolygon[,1]/sqrt(2), nowpolygon[,2]/sqrt(2),
                            Rcss=RC, Rcssclass=c(RCC), col=nowcol,
                            border=NA)
            }
        }
    }

    ## draw the axes around the outside
    Rcsslines(c(0-0.05, -ngroups/2, -0.05), c(0.05, ngroups/2, ngroups-0.05),
              Rcss=RC, Rcssclass=RCC)

    ## draw the stimulus set label on the plot
    Rcsstext((-ngroups/4)-0.25, (ngroups*3/4)-0.25, lab, pos=3, srt=45, Rcss=RC, Rcssclass=RCC)
    
    
    ## draw a barplot on the right hand side
    maxsize = max(gsize)*lineext;
    if (is.null(maxlim)) {
        xlim = c(0, maxsize*3)
    } else {
        xlim = c(0, maxlim)
    }
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", #
         frame=F, axes=F, xlab="", ylab="")
    
    ## draw the bars and labels
    for (i in 1:length(gjh$order)) {
        nowgroup = gjh$labels[gjh$order[i]]
        nowsize = gsize[nowgroup]
        if (nowsize>0) {
            if (nowgroup %in% highlight) {
                Rcsslines(c(nowsize, maxsize), rep(i-0.5, 2), col = col.group[nowgroup],
                          Rcssclass=c(RCC, "highlight"))
                Rcssrect(0, i-1+spacer, nowsize, i-spacer, col=col.group[nowgroup],
                         Rcss=RC, Rcssclass=c(RCC, "highlight"))
                Rcsstext(maxsize, i-0.6, nowgroup, col=col.group[nowgroup],
                         Rcssclass=c(RCC, "highlight"))
            } else {
                Rcssrect(0, i-1+spacer, nowsize, i-spacer,
                         col=paste0(col.group[nowgroup],"99"), border=NA, 
                         Rcss=RC, Rcssclass=c(RCC, "plain"))
                Rcsstext(nowsize+2, i-0.6, nowgroup, col=col.group[nowgroup],
                         Rcssclass=c(RCC, "plain"))
            }
            
        }
    }
    
    Rcssaxis(1, at=c(0, signif(maxsize/1.2, 1)), labels=c("", ""),
         Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(1, at=xticks, labels=xticks, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssmtext("Signature size", side=1, Rcss=RC, Rcssclass=c(RCC, "x"))
    
}







##' Draw a legend with boxes and albels (using Rcssplot)
##'
##' @param col.stim - vector of colors with names.
##' @param boxwidth - numeric
##' @param file - filename
##' @param spacer - numeric
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotStimulusLegend = function(col.stim, boxwidth=0.1, file=NULL, spacer=0.1,
    Rcss="default", Rcssclass="diamond") {
    
    RC = Rcss;
    RCC = Rcssclass;

    temp = boxwidth*(length(col.stim)+length(col.stim)+2)
    Rcsspdf(file=file, width=temp, height=temp)
    Rcsspar(mai=rep(boxwidth, 4))
    
    xlim = c(0, 2*length(col.stim))
    ylim = xlim;
    
    plot(xlim, ylim, type="n", xaxs="i", yaxs="i", xlab="", ylab="",
         frame=F, axes=F)

    col.stim = rev(col.stim)
    
    ## draw one box per stimulus
    for (i in 1:length(col.stim)) {
        
        j = 2*i;
        Rcssrect(0, j-1+spacer, 1.62, j-spacer, col=col.stim[i],
                 Rcss=RC, Rcssclass=c(RCC, "highlight"))
        Rcsstext(1.62, j-0.6, names(col.stim)[i], 
                 Rcssclass=c(RCC, "highlight"))
        
    }
    
    dev.off();
}








##' Draw a bar-plot with two bars per label (using Rcssplot)
##'
##' @param before - numeric vector with names. First set of values to plot
##' @param after - numeric vector with names. Second set of values to plot.
##' @param rhs.lab - character vector with names. Labels to display next to the bars.
##' @param col.stim - vector with colors.
##' @param pixel - numeric
##' @param spacer - numeric. Determines empty space between boxes.
##' @param file - filename
##' @param xticks - numeric vector. Determines values to place on x-axis
##' @param maxlim - numeric. rescaling factor
##' @param legend.x - numeric vector of length two. Determines
##' the position and size (left and right extemes) of boxes in the legend
##' @param legend.y - numberic in range [0,1]. Use to position the legend vertically.
##' @param density - numeric. Density of shading lines on one set of boxes
##' @param legend - character vector of length two. Text labels describing the
##' two series of numbers.
##' @param main - character string. Title for the plot.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotBeforeAfter = function(before, after, rhs.lab = NULL, col.stim=NULL,
    pixel=0.1, spacer=0.04, file=NULL, xticks=seq(0, 250, 50), maxlim=NULL,
    legend.x=c(200, 220), legend.y=0.8, density=50,
    legend = c("Preliminary", "Final"), main="", 
    Rcss="default", Rcssclass="beforeafter") {
    
    RC = Rcss;
    RCC = Rcssclass;
    
    nowmai = RcssGetPropertyValueOrDefault(RC, "par", "mai", default=rep(0.2, 4),
        Rcssclass=RCC)
    nowheight = nowmai[1]+nowmai[3]+((1*pixel+pixel+spacer)*length(col.stim))+spacer;
    
    Rcsspdf(file=file, height=nowheight, Rcss=RC, Rcssclass=RCC)            
    Rcsspar(mai=nowmai, Rcss=RC, Rcssclass=RCC)

    maxsize = max(c(before, after))*1.4;
    if (is.null(maxlim)) {
        xlim = c(0, maxsize*3)
    } else {
        xlim = c(0, maxlim)
    }
    ylim = c(-spacer, ((1*pixel+pixel+spacer)*length(col.stim)));
    
    plot(xlim, ylim, xlim=xlim, ylim=ylim, type="n", xaxs="i", yaxs="i", xlab="", ylab="",
         frame=F, axes=F)
    
    Rcssaxis(1, at=range(xticks), labels=c("", ""),
             Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(1, at=xticks, labels=xticks, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(2, at=xlim+c(-2,2), labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssmtext("Signature size", side=1, Rcss=RC, Rcssclass=c(RCC, "x"))
            
    nowy  = 0;
    ## draw one box per stimulus
    for (i in 1:length(col.stim)) {        
        nowstim = names(col.stim)[i];        
        Rcssrect(0, nowy, after[nowstim], nowy+pixel, col=col.stim[i], 
                 Rcss=RC, Rcssclass=c(RCC, "highlight"))
        Rcssrect(0, nowy+1*pixel, before[nowstim], nowy+(2*pixel),
                 col=col.stim[i], density=density,
                 Rcss=RC, Rcssclass=c(RCC, "highlight", "before"))
        
        Rcsstext(0, nowy+pixel, nowstim, pos=2, Rcss=RC, Rcssclass=c(RCC))
        if (!is.null(rhs.lab)) {
            Rcsstext(max(after[nowstim], before[nowstim]), nowy+pixel, rhs.lab[nowstim],
                     pos=4, Rcss=RC, Rcssclass=c(RCC, "plain"))
        }
        nowy = nowy+(1*pixel+pixel+spacer);
    }
    
    ## draw a manual label
    legend.y = ylim[2]*legend.y;
    Rcssrect(legend.x[1], legend.y, legend.x[2], legend.y+pixel,
             col="#999999",  Rcss=RC, Rcssclass=c(RCC, "highlight"))
    Rcssrect(legend.x[1], legend.y+pixel+spacer, legend.x[2], legend.y+spacer+(2*pixel),
             col="#999999", Rcss=RC, Rcssclass=c(RCC, "highlight", "before"))
    Rcsstext(legend.x[2], legend.y+(pixel/2), legend[2], pos=4, Rcss=RC, Rcssclass=RCC)
    Rcsstext(legend.x[2], legend.y+spacer+(3*pixel/2), legend[1], pos=4, Rcss=RC, Rcssclass=RCC)

    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    
    dev.off();
}






##' Draw a heatmap using Rcssplot.
##'
##' @param dat - matrix of colors.
##' @param legend - vector of colors. Legend shading. Warning: the user bares responsibility
##' to make sure the legend matches with the dat matrix. 
##' @param xlabels - logical. Toggle display of labels on x axis.
##' @param ylabels - logical. Toggle display of labels on y axis.
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title, above heatmap.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
E3PlotBasicHeat = function(dat, legend=NULL,
    ylabels=FALSE, xlabels=TRUE,
    xlab="", ylab="", main="", 
    Rcss="default", Rcssclass=c()) {
    
    RC=Rcss;
    RCC=Rcssclass;
    
    ## get dimensions of the heatmap
    xlim = c(0, ncol(dat))
    ylim = c(0, nrow(dat))

    ## find out how much below the x axis the labels should be
    xdown = RcssGetPropertyValueOrDefault(RC, "basicheat", "xdown", default=-0.1, Rcssclass=RCC)
    ##cat("xdown: ", xdown, "\n")
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
         xlab="", ylab="", frame=F, axes=F)
    
    ## for each segmented sample, for each chromosome, compute and draw smooth FC
    for (i in 1:ncol(dat)) {        
        for (j in 1:nrow(dat)) {
            nowcol = dat[j, i]
            Rcssrect((i-1), (j-1), i, j, col=nowcol, Rcss=RC, Rcssclass=RCC)
        }        
    }
    if (xlabels) {
        for (i in 1:ncol(dat)) {
            Rcsstext(i-0.5, ylim[2]*xdown, colnames(dat)[i], Rcss=RC, Rcssclass=c(RCC, "x"))  
        }
    }
    if (ylabels) {
        for (i in 1:nrow(dat)) {
            Rcsstext(0, i-0.5, rownames(dat)[i], pos=2, Rcss=RC, Rcssclass=c(RCC, "y"))
        }
    }
    
    ## draw divider lines
    for (i in 1:(ncol(dat)-1)) {
        Rcsslines(rep(i,2), ylim, Rcss=RC, Rcssclass=c(RCC, "divider"))
    }
    Rcsslines(c(xlim, rev(xlim), xlim[1]), c(rep(ylim, each=2), ylim[1]),
              Rcss=RC, Rcssclass=c(RCC, "box"))
    
    ## draw a legend in the corner
    E3DrawScaleLegend(legend=legend, xlim=xlim, ylim=ylim, RC=RC, RCC=RCC)
    
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    Rcssmtext(ylab, side=2, Rcss=RC, Rcssclass=c(RCC, "y"))
    Rcssmtext(xlab, side=1, Rcss=RC, Rcssclass=c(RCC, "x"))
    
}






##' Draw a set of violins (using Rcssplot, sm)
##' 
##' @param valset - named list of numeric vectors. These vectors are the datapoints
##' associated with each group and will be used to generate violins. The data points
##' should contain names themselves (to use highlight option below)
##' @param colors - vector of colors. One color for each violin
##' @param xlim - numeric vector of length two. Range of x values
##' @param xat - numeric vector. Determines where markers will be placed on the horizontal
##' axis
##' @param highlight - vector of item names to highlight
##' @param hidedots - logical. Set TRUE to avoid dots on the violins.
##' @param showconnect - logical. Set TRUE to connect highlighted items with lines
##' @param xlab - character string. Text to display beside x axis
##' @param ylab - character string. Text to display besdie y axis
##' @param main - character string. Text to display in the title of the plot
##' @param spacer - numeric. Determines space between adjacent violins
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotVioSet = function(valset, colors=NULL, xlim=c(-1,1), xat=c(-0.5,0,0.5),
    highlight=NULL, hidedots=FALSE, showconnect=FALSE, 
    xlab="", ylab="", main="", spacer=0.04, 
    Rcss="default", Rcssclass=c()) {
    
    RC = Rcss;
    RCC = Rcssclass;
    
    ylim = c(-2*spacer, length(valset)+2*spacer)
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)
    plot(0, 0, type="n", ylim=ylim, xlim=xlim,
         axes=FALSE, frame=FALSE, xaxs="i", yaxs="i", xlab="", ylab="");
    
    E3Axes(xlim, ylim, xat=xat, xlab=xlab, ylab=ylab, RC=RC, RCC=RCC, where="x")
    Rcssaxis(2, at=ylim, pos=xlim[1], labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "y"))

    highlight.x = list()
    if (!is.null(highlight)) {
        for (j in 1:length(highlight)) {
            highlight.x[[j]] = c(0)
        }
    }
    highlight.y = c()
    
    if (length(valset)>0) {
        for (i in 1:length(valset)) {
            nowname = names(valset)[i]
            nowvals = valset[[i]];
            nowy = (length(valset)-i+1);

            highlight.y = c(highlight.y, nowy-0.12, nowy-0.88)
            
            if (length(nowvals)>1) {
                E3Vioplot(nowvals, col=paste0(colors[nowname],"88"), at=nowy-0.5, add=T,
                        drawRect=F, lwd=0.3, wex=0.85, border=NA, horizontal=TRUE);
                
                if (!is.null(highlight)) {
                    plainvals = nowvals[!(names(nowvals)%in%highlight)]
                    highvals = nowvals[names(nowvals) %in% highlight]
                    if (!hidedots) {
                        Rcsspoints(plainvals, runif(length(plainvals), nowy-0.8, nowy-0.2),
                                   Rcss=RC, Rcssclass=RCC);
                    }
                    for (j in 1:length(highvals)) {
                        nowv = highvals[j];
                        Rcsslines(rep(nowv, 2), c(nowy-0.88, nowy-0.12),
                                  Rcss=RC, Rcssclass=c(RCC, "highlight", names(highvals[j])));
                        if (is.null(highlight.x[[j]])) {
                            highlight.x[[j]] = rep(nowv, 2)
                        } else {
                            highlight.x[[j]] = c(highlight.x[[j]], rep(nowv, 2))
                        }
                    }
                } else {
                    if (!hidedots) {
                        Rcsspoints(nowvals, runif(length(nowvals), nowy-0.8, nowy-0.2),
                                   Rcss=RC, Rcssclass=RCC);
                    }
                }
                
            }
            Rcsstext(xlim[1], nowy-0.5, nowname, pos=2, Rcss=RC, Rcssclass=c(RCC))
        }
    }

    if (!is.null(highlight) & showconnect) {
        for (j in 1:length(highlight)) {
            Rcsslines(highlight.x[[j]][-1], highlight.y, Rcss=RC,
                      Rcssclass=c(RCC, "highlight", "light", names(highvals[j])))
        }                
    }
    
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"));    
}








##' Draw a scatter plot using bubble of various sizes (using Rcssplot)
##'
##' @param x - numeric vector with names. X-coordinates
##' @param y - numeric vector with names. Y-coordinates
##' @param dsize - numeric vector with names. Determines size of the bubble
##' @param dp - numeric vector with names. Determines color shading of the bubble. These
##' values will be -log10 transformed before being converted into colors
##' (e.g. use p values associated with each item)
##' @param labels - character vector. Names of items which require text labels
##' @param xlab - character string. Text to display beside x axis
##' @param ylab - character string. Text to display besdie y axis
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotScatterBubbles = function(x, y, 
    dsize=NULL, dp=NULL, labels=NULL, xlab="", ylab="", 
    Rcss="default", Rcssclass=c()) {
    
    if (is.null(dsize)) {
        dsize = rep(1, length(x));
        names(dsize) = names(x)
    }
    if (is.null(dp)) {
        dp = rep(1, length(x));
        names(dp) = names(x)
    }

    RC = Rcss;
    RCC = Rcssclass;
    
    xnames = names(x)
    
    xlim = c(0, max(x)*1.05)
    ylim = c(0, max(y)*1.05)
    
    plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, axes=F, frame=F, 
         main="", xlab="", ylab="", xaxs="i", yaxs="i")
    
    KOgenes.fisher.bg = -log10(dp)
    KOgenes.fisher.bg[KOgenes.fisher.bg>5] = 5;
    KOgenes.fisher.bg = KOgenes.fisher.bg/5;
    KOgenes.fisher.bg = val2hex(KOgenes.fisher.bg)
    
    E3Axes(xlim, ylim, xlab=xlab, ylab=ylab, RC=RC, RCC=RCC, where=c("x","y"))
    Rcssaxis(2, at=ylim, labels=c("", ""), line=0, tck=0, Rcss=RC, Rcssclass=c(RCC, "y"))
    
    Rcsspoints(x, y[xnames], cex=dsize[xnames], pch=21, lwd=0.4,
               col="#000000", bg = paste0("#ff0000", KOgenes.fisher.bg), xpd=1,
               Rcss=RC, Rcssclass=RCC)
    if (!is.null(labels)) {
        Rcsstext(x[labels], y[labels], labels, cex=0.3, xpd=1,
                 Rcss=RC, Rcssclass=RCC)
    }
    
    
    Rcssmtext(xlab, side=1, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssmtext(ylab, side=1, Rcss=RC, Rcssclass=c(RCC, "y"))
    
}





##' Draw a contour plot with select markers
##'
##' @param x - numeric vector with names. Coordinates for x axis. 
##' @param y - numeric vector with names. Coordinates for y axis
##' @param highlight - character string. names of x,y to highlight.
##' @param highlight.label - character vector. Names of points to highlight on the plot
##' @param col - color. For highlighting points
##' @param label - logical. Set TRUE to draw text labels next to highlighted points
##' @param density.transformation - obsolete, do not change.
##' @param show.fit - logical. Set TRUE to display a fitted line through highlighted points
##' @param xylim - numeric vector of length two. Ranges for axes
##' @param xymajor - numeric vector. Location of major tick marks
##' @param xyminor - numeric vector. Location of minor tick marks
##' @param xlab - character string. Label for x axis
##' @param ylab - character string. Label for y axis
##' @param main - character strings. Title of the plot
##' @param breakrange - numeric vector of length 2. Lower and upper bounds for contour plot
##' @param numbreaks - integer. Number of partitions for contour matrix
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##'
##' @export
E3PlotContourScatter = function(x, y, highlight=NULL, col="#ff0000",
    highlight.label=FALSE,
    density.transformation = function(x) { x }, show.fit=FALSE, 
    xylim=c(0.1, 100),
    xymajor = 10^seq(-3, 6),
    xyminor = c(),
    xlab="", ylab="", main="", 
    breakrange=range(c(x, y)), numbreaks=20, 
    Rcss="default", Rcssclass="contour") {

    RC = Rcss;
    RCC = Rcssclass;
    
    ## transform the input data
    xlog = log10(x)
    ylog = log10(y)    
    xylimlog = log10(xylim)
    
    ## create a density from all the data
    breaks = seq(log10(breakrange[1]), log10(breakrange[2]), length=numbreaks+1)
    breakmids = (breaks[-1]+breaks[-length(breaks)])/2    
    dd = matrix(0, ncol=numbreaks, nrow=numbreaks)    
    for (i in 1:(numbreaks-1)) {
        inow = breaks[i]
        inext = breaks[i+1]
        for (j in 1:(numbreaks-1)) {
            jnow = breaks[j]
            jnext = breaks[j+1]            
            dd[i,j] = sum(xlog>=inow & xlog<inext & ylog>=jnow & ylog<jnext)            
        }
    }
    dd = density.transformation(dd)
    
    Rcsspar(Rcss=RC, Rcssclass=RCC)    
    plot(xylimlog, xylimlog, xlim=xylimlog, ylim=xylimlog, xlab="", ylab="", 
         xaxs="i", yaxs="i", type="n", frame=F, axes=F)

    ## draw the diagonal line
    Rcsslines(xylimlog, xylimlog, Rcss=RC, Rcssclass=c(RCC, "diagonal"))     
    
    ## draw the contour of the background distribution
    Rcsscontour(breakmids, breakmids, dd, add=TRUE, Rcss=RC, Rcssclass=RCC)
    ## draw the points to highlight
    if (!is.null(highlight)) {
        Rcsspoints(xlog[highlight], ylog[highlight], col=col, Rcss=RC, Rcssclass=RCC)
    }
    ## write text labels (can tune labels above/below the diagonal)
    if (!is.null(highlight.label)) {
        high.down = which(xlog[highlight.label]>=ylog[highlight.label])
        high.down = highlight.label[high.down]
        high.up = which(xlog[highlight.label]<ylog[highlight.label])
        high.up = highlight.label[high.up]
        if (length(high.down)>0) {
            Rcsstext(xlog[high.down], ylog[high.down], high.down, Rcss=RC, Rcssclass=c(RCC, "down"))
        }
        if (length(high.up)>0) {
            Rcsstext(xlog[high.up], ylog[high.up], high.up, Rcss=RC, Rcssclass=c(RCC, "up"))
        }        
    }
    ## draw a fitted line
    if (!is.null(highlight) & show.fit) {
        temp.model = lm(ylog[highlight]~xlog[highlight])$coefficients
        ymodel = temp.model[1] + (xylimlog*temp.model[2])
        ymodel = as.numeric(ymodel)
        Rcsslines(xylimlog, ymodel, xpd=0, Rcss=RC, Rcssclass=c(RCC, "fit"))
    }    
    
    ## draw the box around the data
    ##Rcssbox(Rcss=RC, Rcssclass=RCC)
    
    ## draw the axes and labels
    Rcssaxis(1, at=xylimlog, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(1, at=log10(xymajor), labels=rep("", length(xymajor)),
             Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssaxis(1, at=log10(xymajor), labels=xymajor, lwd=0, 
             Rcss=RC, Rcssclass=c(RCC, "x"))
    
    
    Rcssaxis(2, at=xylimlog, labels=c("", ""), tck=0, line=0, Rcss=RC, Rcssclass=c(RCC, "y"))
    Rcssaxis(2, at=log10(xymajor), labels=rep("", length(xymajor)), line=0, 
             Rcss=RC, Rcssclass=c(RCC, "y"))
    Rcssaxis(2, at=log10(xymajor), labels=xymajor, lwd=0, 
             Rcss=RC, Rcssclass=c(RCC, "y"))
    
    
    if (length(xyminor)>0) {
        Rcssaxis(1, at=log10(xyminor), labels=rep("", length(xyminor)),
                 Rcss=RC, Rcssclass=c(RCC, "x", "minor"))
        Rcssaxis(2, at=log10(xyminor), labels=rep("", length(xyminor)),
                 Rcss=RC, Rcssclass=c(RCC, "y", "minor"), line=0, labels=NA)                    
    }
    
    Rcssmtext(xlab, side=1, Rcss=RC, Rcssclass=c(RCC, "x"))
    Rcssmtext(ylab, side=2, Rcss=RC, Rcssclass=c(RCC, "y"))
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"));


    return(list(xlog=xlog, ylog=ylog))
    
}
