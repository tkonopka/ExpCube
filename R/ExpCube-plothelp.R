##
## ExpCube analysis
## Author: Tomasz Konopka
##
## Here are some helper functions used in plotting
##



##' Draw axes and labels on x/y, using Rcssplot
##'
##' @param xlim - two element vector. Range of x axis
##' @param ylim - two element vector. Range of y axis
##' @param xlab - character string. Text to write below x-axis
##' @param ylab - character string. Text to write below y-axis
##' @param RC - Rcss object. Style for formatting the axes using Rcssplot
##' @param RCC - character vector. Classes for tuning Rcssplot formating
##' @param where - vector containg "x" and/or "y". Determines where axes will be drawn
##' 
##' @export
E3Axes = function(xlim, ylim, xlab="", ylab="", xat=NULL,
    RC="default", RCC=c(), where=c("x", "y")) {

    ## display some axes
    if ("x" %in% where) {
        Rcssaxis(1, at=xlim, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "x"))
        if (is.null(xat)) {
            Rcssaxis(1, Rcss=RC, Rcssclass=c(RCC, "x"))
        } else {
            Rcssaxis(1, at=xat, labels=xat, Rcss=RC, Rcssclass=c(RCC, "x"))
        }
    }
    if ("y" %in% where) {
        Rcssaxis(2, at=ylim, labels=c("", ""), tck=0, Rcss=RC, Rcssclass=c(RCC, "y"), line=0)
        Rcssaxis(2, Rcss=RC, Rcssclass=c(RCC, "y"), lwd=0)
        Rcssaxis(2, Rcss=RC, Rcssclass=c(RCC, "y"), line=0, labels=NA)
    }
    
    ## display the axis labels
    if (xlab != "") {
        Rcssmtext(xlab, side=1, Rcss=RC, Rcssclass=c(RCC, "x"));
    }
    if (ylab != "") {
        Rcssmtext(ylab, side=2, Rcss=RC, Rcssclass=c(RCC, "y"));
    }
}




##' Convert a matrix of values (range [-Inf,Inf] to colors using transparency)
##'
##' See also E3Val2Col. Uses val2hex from package Rpipelines
##' 
##' @param x - numeric matrix.
##' @param col - vector of two colors in #XXXXXX format. First element determines
##' color associated with negative values. Second element determines color associated
##' with positive values
##' @param maxval - numeric. Value for which color reaches saturation
##'
##' @export
E3ValMat2ColMat = function(x, col=c("#0000ff", "#ff0000"), maxval=1) {            
    xpos = x>0;
    xneg = x<0;        
    ## force absolute values into range [0,1]
    y = abs(x);
    y[y>maxval]=maxval;
    yvals = as.vector(y)/maxval;
    ## make a matrix of colors
    temp = paste0(col[1], val2hex(yvals));
    if (sum(xpos)>0) {
        temp[xpos] = paste0(col[2], val2hex(yvals[xpos]));
    }
    ans = matrix(temp, ncol=ncol(x), nrow=nrow(x));
    rownames(ans) = rownames(x);
    colnames(ans) = colnames(x);
    ## replace full transparency with #ffffff
    ans[x==0] = "#ffffff";
    return(ans);
}




##' Convert a number into a color using  transparency
##'
##' See also E3ValMat2ColMat. Uses val2hex from package Rpipelines
##' 
##' @param x - numeric matrix.
##' @param col - vector of two colors in #XXXXXX format. First element determines
##' color associated with negative values. Second element determines color associated
##' with positive values
##' @param maxval - numeric. Value for which color reaches saturation
##'
##' @export
E3Val2Col = function(x, col=c("#0000ff", "#ff0000"), maxval=1) {        
    ## create a boolean vector of colors based on sign of x
    ans = rep(col[1], length(x))
    ans[x>0] = col[2]    
    ## make sure values are within [-1, 1] range
    x[x>maxval] = maxval;
    x[x<(-maxval)] = -maxval;
    x = x/maxval    
    ## append a transparency value and that's it
    ans = paste0(ans, val2hex(abs(x)))
    names(ans) = names(x)
    return(ans)
}




##' Draw a vertical recrangle with a color scale and labels
##'
##' This function uses Rcss selector "scalelegend" to determine legend position.
##' 
##' @param legend - named vector of colors. Colors defining a color scale. All
##' names associated with vector will be draw onto to the plot. (Set some names to ""
##' to avoid labeling every single color in a color gradient)
##' @param xlim - numeric vector of two elements. Gives range of x axis.
##' @param ylim - numeric vector of two elements. Gives range of y axis.
##' @param main - text two write above the legend
##' @param RC - Rcss object. Style to use for plotting, uses package Rcssplot
##' @param RCC - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
E3DrawScaleLegend = function(legend=NULL, xlim=c(0,1), ylim=c(0,1), main="",   
    RC="default", RCC=c()) {
    
    ## get features of the legend from the RC
    relposx =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relposx",
        default=-0.1, Rcssclass=RCC)
    relwidth =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relwidth",
        default=0.05, Rcssclass=RCC)
    relposy =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relposy",
        default=0.75, Rcssclass=RCC)    
    relheight =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relheight",
        default=0.15, Rcssclass=RCC)
    
    ## find distances on x/y axes
    xd = xlim[2]-xlim[1]
    yd = ylim[2]-ylim[1]
    ## define coordinates of legend corners
    bl = c(xlim[1]+(xd*relposx), ylim[1]+(yd*relposy))
    br = c(bl[1]+(xd*relwidth), bl[2])
    tl = c(bl[1], bl[2]+(yd*relheight))
    tr = c(br[1], tl[2])
    
    ysteps = seq(bl[2], tl[2], length=length(legend)+1)
    absv = yd*RcssGetPropertyValueOrDefault(RC, "scalelegend", "relv",
        default=0.05, Rcssclass=RCC)
    
    ## draw a legend in the corner
    if (!is.null(legend) & length(legend)>0) {        
        for (j in 1:length(legend)) {
            Rcssrect(bl[1], ysteps[j], tr[1], ysteps[j+1], col=legend[j], Rcss=RC,
                     Rcssclass=c(RCC,"scalelegend"))
        }
        Rcsslines(c(bl[1], br[1], tr[1], tl[1], bl[1]), c(bl[2], br[2], tr[2], tl[2], bl[2]), 
                  Rcss=RC, Rcssclass=c(RCC, "scalelegend"))
        Rcsstext(br[1], seq(bl[2], tl[2], length=length(legend)), names(legend),
                 Rcss=RC, Rcssclass=c(RCC, "scalelegend"))

        Rcsstext(tl[1], tl[2]+absv, main,
                 Rcss=RC, Rcssclass=c(RCC, "scalelegend", "main"))        
    }
    
    
}



##' Draw a vertical legend with markers (circles, boxes, and labels)
##'
##' This function uses Rcss object "markerlegend" to determine legend position.
##' 
##' @param legend - data frame with columns cex, pch, col, bg, and label. Order and contents
##' determines what markers are drawn
##' @param xlim - numeric vector of two elements. Range of x axis (used to place the legend)
##' @param ylim - numeric vector of two elements. Range of y axis (used to place the legend)
##' @param main - text to write above the legend
##' @param RC - Rcss object. Style to use for plotting, uses package Rcssplot
##' @param RCC - character vetor. Classes to tune Rcssplot formatting. 
##'
##' @export
E3DrawMarkerLegend = function(legend=NULL, xlim=c(0,1), ylim=c(0,1), main="",
    RC="default", RCC=c()) {
    
    if (is.null(legend)) {
        return()
    }
    if (sum(c("cex", "pch", "col", "bg", "label") %in% colnames(legend))!=5) {
        stop("legend object must contain columns cex, pch, col, bg, label\n")
    }
    
    ## get features of the legend from the RC
    relposx =  RcssGetPropertyValueOrDefault(RC, "markerlegend", "relposx",
        default=-0.1, Rcssclass=RCC)
    relwidth =  RcssGetPropertyValueOrDefault(RC, "markerlegend", "relwidth",
        default=0.05, Rcssclass=RCC)
    relposy =  RcssGetPropertyValueOrDefault(RC, "markerlegend", "relposy",
        default=0.75, Rcssclass=RCC)    
    
    ## find distances on x/y axes
    xd = xlim[2]-xlim[1]
    yd = ylim[2]-ylim[1]
    
    ## coordinates of top-left coorner
    tl = c(xlim[1]+(xd*relposx), ylim[1]+(yd*relposy))
    tr = c(tl[1]+(xd*relwidth), tl[2])
    tm = c(tl[1]+(xd*relwidth/2), tl[2])
    
    absv = yd*RcssGetPropertyValueOrDefault(RC, "scalelegend", "relv",
        default=0.05, Rcssclass=RCC)

    llen = nrow(legend)
    Rcsspoints(rep(tm[1], llen), tm[2]-(seq(0, llen-1) *absv),
               cex=legend[,"cex"], col=legend[,"col"], pch=legend[,"pch"], bg=legend[,"bg"],
               Rcss=RC, Rcssclass=c(RCC, "markerlegend"))
    Rcsstext(rep(tr[1], llen), tr[2]-(seq(0,llen-1)*absv),
             legend[,"label"], Rcss=RC, Rcssclass=c(RCC, "markerlegend"))
    
    
               
}




##' Draw violin onto a plot.
##'
##' This is a hacked version of vioplot from package vioplot. This version does not
##' "box()" every time it is called
##'
##' @param x - see vioplot from package vioplot
##' @param range
##' @param h
##' @param ylim
##' @param names
##' @param horizontal
##' @param col
##' @param border
##' @param lty
##' @param lwd
##' @param rectCol
##' @param colMed
##' @param pchMed
##' @param at
##' @param add
##' @param wex
##' @param drawRect
##' @param ... 
##' 
##' @export
E3Vioplot = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
  horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
  lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
  at, add = FALSE, wex = 1, drawRect = TRUE) {
  
  datas <- list(x, ...)
    n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) {
    args <- c(args, h = h)
  }
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  boxwidth <- 0.05 * wex

  if (!add) {
    plot.new()
  }
  
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    ##box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
                lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    ##box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
             boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}


