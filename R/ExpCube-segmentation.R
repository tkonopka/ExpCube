##
## ExpCube - Segmentation analysis
## Author: Tomasz Konopka
##
##
## Functions relevant for segmentation analysis
##


##' Obtain genomic positions of genes from a refGene-style annotation table.
##'
##' @param annotab - data frame. A refGene style table with transcript definitions
##' @param chrnames - character vector. Use this to look at transcripts on just a few
##' chromosomes, e.g. to avoid genes on decoy contigs.
##' 
##' @export
E3GetGenePosition = function(annotab, chrnames=NULL) {
    
    if (is.null(chrnames)) {
        stop("must specify chromosome names")
    }
    ## select desired chromosomes
    annotab = annotab[annotab[,3] %in% chrnames, ];
    ## compute rough center of transcript based on txs start/end positions
    ans = data.frame(Gene = annotab[,13], chr=annotab[,3], pos=(annotab[,5]+annotab[,6])/2, 
        stringsAsFactors=F);

    ## convert between transcript level positions to gene level positions
    ans = split(ans, ans[,"Gene"]);
    ans = lapply(ans,
        function(x) {
            if (length(unique(x[,"chr"]))>1) {
                x[,"pos"] = NA;
                x[,"chr"] = NA;
                return(x[1, ,drop=FALSE]);
            } else {
                x[,"pos"] = round(median(x[,"pos"]));
                return(x[1,,drop=FALSE]);
            }
        })
    
    ans = data.frame(rbindlist(ans), stringsAsFactors=F)
    rownames(ans) = ans[,"Gene"];
    return(ans);    
}



##' Perform genomic segmentation.
##'
##' This function uses CBS (package DNAcopy) to segment genomic data. This wrapper
##' allows to start with a data matrix containing values associated with gene names. A
##' separated object mapping genes to genomic positions is used to arrange the data
##' before segmentation. 
##'
##' Values in fcdata will be log2 transformed before being sent to segmentation. Functions
##' returns a list with the resulting segmentation.
##' 
##' @param fcdata - numeric matrix with sample names in colnames and gene names in rownames.
##' @param gene.position - data frame with columns "chr", "pos", and rownames
##' containing gene names.
##' @param genenames - names of 
##' @param chrnames - character vector with names of chromosomes to include in segmentation.
##' Leave null, To use all the genes in the input matrix.
##' @param min.width - minimum number
##' @param verbose - logical. When TRUE, function prints dots indicating progress.
##' 
##' @export
E3MakeFCSegmentation = function(fcdata, gene.position=NULL, genenames=NULL,
    chrnames=NULL, min.width=5, verbose=TRUE) {
    
    if (is.null(chrnames)) {
        stop("must specify chromsome names")
    }
    if (!is.null(genenames)) {
        fcdata = fcdata[rownames(fcdata) %in% genenames,]
        gene.position = gene.position[rownames(gene.position) %in% genenames,]
    }
    
    gps = gene.position[gene.position[,"chr"] %in% chrnames,]
    fcdata = fcdata[rownames(fcdata) %in% rownames(gps),]
    
    segdata = list()
    for (nowname in colnames(fcdata)) {
        if (verbose) {
            cat(".")
        }

        ## create an object suitable for segmentation
        ## add one dummy point per chromosome to ensure that segments always start at
        ## position 1
        nowdata = rbind(data.frame(exp=log2(fcdata[,nowname]),
            chr=gps[rownames(fcdata), "chr"],
            pos=gps[rownames(fcdata), "pos"], stringsAsFactors=F),
            data.frame(exp=1, chr=chrnames, pos=1, stringsAsFactors=F))        
        ## remove up/down regulation that is so high it can skew the mean...
        nowdata[nowdata[,"exp"]>3, "exp"] = 3;
        nowdata[nowdata[,"exp"]<(-3), "exp"] = -3;
        ## segment the data using CBS
        nowdata = unique(nowdata);
        nowdata = nowdata[order(nowdata[,"chr"], nowdata[,"pos"]),]        
        nowseg = segment(
            CNA(nowdata[,"exp"], nowdata[,"chr"], nowdata[,"pos"], sampleid=nowname),
            verbose=0, min.width=min.width)
        segdata[[nowname]] = nowseg$output
    }
    if (verbose) {
        cat("\n")
    }
    return(segdata)
}

