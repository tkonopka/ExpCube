##
## ExpCube 
## Author: Tomasz Konopka
##
## Functions used to modify gene expression intervals
##




##' Estimates an interval rescale factor for each feature/gene based on within-group variability.
##'
##' The goal is to find genes that are more variable in practice than the
##' uncertainty interval model may suggest.
##'
##' computes an empirical interval from values of a background sample set
##' computes a representative interval from values of the background sample set
##' computes the ratio of the empirical to representative - this is the output
##'
##' 
##' @param sampledata - a list with tables $expression, $expression.low, $expression.high
##' @param id.col - character strings. Name of column in sampledata that identifies
##' the gene/transcript species
##' @param bg.samples - list of sample names. Each vector should represent a group of related
##' samples (replicates). The function will compute variability within each group.
##' @param upper.quantile - numeric. The quantile used to equalize the input and output intervals.
##' @param min.factor - numeric. Minimum rescaling factor, i.e. when procedure gives
##' a suggestion to rescale by a factor lower than this, the output will in fact be min.factor.
##'
##' @export
E3CalcWithinGroupIntervalRescaleFactors = function(sampledata, id.col="Gene",
    bg.samples=NULL, upper.quantile=0.85, min.factor=1) {
    
    ## use tables from sampledata, but use shorthand code notation here 
    tab = sampledata$expression
    tab.high = sampledata$expression.high
    
    ## compute one set of renorm factors for each set of reference samples
    gene.res.vecs = lapply(bg.samples,
        function(x) {
            gene.rescaling = rep(1, nrow(tab))
            names(gene.rescaling) = as.character(tab[,id.col])
            
            ## compute the median and upper quantile expression for all genes
            bg.stats = t(apply(tab[,x], 1, quantile, p=c(0.5, upper.quantile)))
            ## compute the interval as the distance between the upper quantile and median,
            ## relative to median
            bg.interval = (bg.stats[,2]-bg.stats[,1])/bg.stats[,1]
            
            ## compute the actual median interval size
            now.interval = apply((tab.high[,x]-tab[,x]) / tab[,x], 1, median)
            
            gene.rescaling = bg.interval/now.interval
            gene.rescaling[!is.finite(gene.rescaling) | gene.rescaling==0] = 1
            gene.rescaling[gene.rescaling<min.factor] = min.factor
            gene.rescaling = signif(gene.rescaling, 4)
            names(gene.rescaling) = as.character(tab[,id.col])
            return(gene.rescaling)
        })
    
    gene.rescale.table = data.frame(Gene=as.character(tab[,id.col]), do.call(cbind, gene.res.vecs),
        stringsAsFactors=FALSE)
    colnames(gene.rescale.table)[1] = id.col
    
    if (length(bg.samples)==1) {
        gene.rescale.table[,"rescale.factor"] = gene.rescale.table[,names(bg.samples)]
    } else {
        qqq = function(x) { as.numeric(quantile(x, upper.quantile)) }
        gene.rescale.table[,"rescale.factor"] = 
            apply(gene.rescale.table[,names(bg.samples)], 1, qqq)
        rm(qqq)
    }
    gene.rescale.table[,"rescale.factor"] = signif(gene.rescale.table[,"rescale.factor"], 4)
        
    return(gene.rescale.table[,c(id.col, "rescale.factor")])
}




##' Estimate an interval rescale factor for each feature/gene based on across-group
##' variability.
##' 
##' This goal is to find genes that are prone to batch effects, i.e. that differ between
##' batches more than their intervals would suggest. 
##'
##' 
##' @param sampledata - a list with tables "expression", "expression.low", "expression.high"
##' @param id.col - character strings. Name of column in sampledata that identifies
##' the gene/transcript species
##' @param bg.samples - list of sample names. Each vector should represent a group of related
##' samples (replicates). The function will compute variability between across these groups.
##' @param upper.quantile - numeric. The quantile used to equalize the input and output intervals.
##' @param min.factor - numeric. Minimum rescaling factor, i.e. when procedure gives
##' a suggestion to rescale by a factor lower than this, the output will in fact be min.factor.
##'
##' @export
E3CalcAcrossGroupIntervalRescaleFactors = function(sampledata, id.col="Gene",
    bg.samples=NULL, upper.quantile=0.85, min.factor=1) {
    
    ## use tables from sampledata, but use shorthand code notation here 
    tab = sampledata$expression
    tab.high = sampledata$expression.high
    
    ## compute median values for genes within each bg.group
    gene.med.vecs = lapply(bg.samples,
        function(x) {
            gene.rescaling = apply(tab[,x], 1, median)
            names(gene.rescaling) = as.character(tab[,id.col])
            return(gene.rescaling)
        });
    gene.uq.vecs = lapply(bg.samples,
        function(x) {
            gene.uq = apply(tab.high[,x], 1, median)
            names(gene.uq) = as.character(tab[,id.col])
            return(gene.uq)
        });
    ## compute median across all groups
    gene.med.med = apply(do.call(cbind, gene.med.vecs), 1, median)
    
    ## compute how far off the within-group medians are from the global median
    ## in units of the within-group interval size
    gene.med.off = list()
    for (i in 1:length(bg.samples)) {
        x = gene.med.vecs[[i]]
        xuq = gene.uq.vecs[[i]]
        temp = abs(x-gene.med.med)/abs(xuq-x)
        temp[!is.finite(temp)] = 1
        temp[temp<min.factor] = 1
        gene.med.off[[i]] = temp
    }
    names(gene.med.off) = names(bg.samples);
    
    gene.rescale.table = data.frame(Gene=as.character(tab[,1]),
        do.call(cbind, gene.med.off),
        stringsAsFactors=FALSE)
    colnames(gene.rescale.table)[1] = id.col
    
    if (length(bg.samples)==1) {
        gene.rescale.table[, "rescale.factor"] = gene.rescale.table[,names(bg.samples)]
    } else {
        qqq = function(x) { as.numeric(quantile(x, upper.quantile, na.rm=TRUE)) }
        gene.rescale.table[, "rescale.factor"] =
            apply(gene.rescale.table[,names(bg.samples)], 1, qqq)
        rm(qqq)
    }
    gene.rescale.table[,"rescale.factor"] = signif(gene.rescale.table[,"rescale.factor"], 4)
    
    return(gene.rescale.table[,c(id.col, "rescale.factor")])
}






##' Apply interval rescale factors, i.e. change low and high expression estimates.
##'
##' When a gene expression estimate is, e.g.  (2, 4, 6) for (low, middle, high) and the rescale
##' factor is 1.5, the output expression estimate will be (1, 4, 7). The middle value stays
##' the same, the distances between the outer values and the middle is rescale by the factor.
##'
##' Function returns a similar list of table as that in sampledata. 
##' 
##' @param sampledata - a list with tables $expression, $expression.low, $expression.high
##' @param samplenames - character vector. Names of all samples to rescale
##' @param gene.rescaling - numeric vector. Rescales all intervals in both directions. 
##' The function assumes length of this vector is equal to rows of the expression data
##' 
##' @export
E3RescaleIntervals = function(sampledata, samplenames, gene.rescaling) {
    
    ## set up the basic structure of the answer, mimicking the input 
    tab = sampledata$expression
    tab.high = sampledata$expression.high
    
    ## get a new object with rescaled uncertainties
    ans = list()
    ans$expression = tab
    ans$expression.high = tab.high
    ans$expression.low = sampledata$expression.low
    
    ans$expression.high[,samplenames] =
        tab[,samplenames] + ((tab.high[,samplenames] - tab[,samplenames])*gene.rescaling)
    ans$expression.low[,samplenames] =
        tab[,samplenames] - ((tab.high[,samplenames] - tab[,samplenames])*gene.rescaling)
    
    ## make sure the low estimate are non-negative
    for (nowsample in samplenames) {
        temp = ans$expression.low[,nowsample]
        temp[temp<0] = 0
        ans$expression.low[,nowsample] = temp
    }
    
    return(ans)    
}


