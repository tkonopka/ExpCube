##
## ExpCube analysis
## Author: Tomasz Konopka
##
## Here are some frequently used "helper" functions
##



##' Create a list of vectors describing how samples are related into groups.
##' 
##' @param summarytable - data frame with columns Sample, Group, etc.
##' @param refgroups - character vector. Names of groups that contain "reference" samples.
##' These groups are relevant for splitting outlier samples into two sub-groups
##' @param outliers - names of samples to be considered as outliers. These samples will
##' be placed in separate groups
##' 
##' @export
E3MakeSampleGroups = function(summarytable, refgroups=c(), outliers=c()) {
    
    samplenames = summarytable[,"Sample"]

    ## check that outliers match the samplenames in the 
    outliers = outliers[outliers %in% samplenames];
    
    ## create a first attempt at a sampleGroups list (does not take outliers into account)
    ans = split(samplenames, summarytable[,"Group"])
    
    ## get a vector with names of reference samples (specified from GO.groups);
    ref.samples = unique(unlist(ans[unique(refgroups)]));
        
    ## find samples that are outliers AND are ref samples
    ## find samples that are outliers AND are NOT ref samples
    ref.outliers = intersect(outliers, ref.samples);
    other.outliers = outliers[!(outliers %in% ref.outliers)]
    
    ans = lapply(ans, function(x) {x[!(x%in%outliers)]} );
    ans = ans[sapply(ans, length)>0]
    if (length(ref.outliers)>0) {
        ans[["ZZZ_Outliers_Ref"]] = ref.outliers;
    }
    if (length(other.outliers)>0) {
        ans[["ZZZ_Outliers_Other"]] = other.outliers;
    }
    
    return(ans);
}




##' Determine genes that are expressed in a dataset. Expressed genes are those that
##' have expression in a range (min.expression, max.expression)
##' 
##' @param dlist - a list of dataframes, one of which must be $expression. Data frames
##' must have one identifier column and the remaining columns will be treated as samples.
##' @param id.col - character string. Name of column that contain gene identifier
##' @param min.expression - numeric. Minimal expression level
##' @param max.expression - numeric. Maximal expression level
##' @param min.count - integer. Number of samples that must lie in range
##' [min.expression, max.expression]. Use this to require a gene to be expressed
##' in several samples in a cohort.
##' 
##' @export
E3GetExpressedGenes = function(dlist, id.col="Gene", min.expression=1, max.expression=Inf,
    min.count=1) {

    ans = dlist$expression;
    samplenames = colnames(ans);
    samplenames = samplenames[samplenames != id.col]
    
    ## loop through all samples and determine if a gene is expressed above minimal required level
    for (nowsample in samplenames) {
        temp = dlist$expression[,nowsample]
        ans[,nowsample] = as.integer(temp>min.expression & temp<max.expression)
    }
    
    ## find genes that are expressed in at least 1 of the samples
    temp = apply(ans[,samplenames], 1, sum)>min.count;
    
    ## return a vector of gene names
    return(ans[temp, id.col]);  
}






##' Compute spearman pseudo-distances between all pairs of samples. Spearman
##' pseudo distance between two samples is defined as 1-spearman.rho
##' 
##'
##' @param dd - numeric matrix. Columns should contain sample names (next argument)
##' @param samplenames - character vector with column names
##'
##' @export
E3dist.spearman = function(dd, samplenames) {

  ## simplify the matrix
  dd = dd[,samplenames, drop=FALSE];

  ## compute the ranks once 
  ddranks = dd;
  for (i in 1:ncol(dd)) {
    ddranks[,i] = rank(dd[,i]);
  }
  
  ## make a matrix with all distances
  ans = matrix(0, ncol=length(samplenames), nrow=length(samplenames));
  colnames(ans) = samplenames;
  rownames(ans) = samplenames;
  
  for (i in 1:(length(samplenames)-1)) {
    for (j in (i+1):length(samplenames)) {      
      nowdist = 1-spearmanRhoFromRanks(dd[,i], dd[,j]);
      ans[i,j] = nowdist;
      ans[j,i] = nowdist;
    }
  }

  ## return an object with distances only
  return(as.dist(ans));  
}




##' Compute euclidean distances between all pairs of samples
##' 
##' @param dd - numeric matrix. Columns should contain sample names (next argument)
##' @param samplenames - character vector with column names in dd.
##' ##@param noise.sd - numeric. If set >0, the distances will contain a noise element.
##' ##If the computed distance is D, the output distance will be D + abs(rnorm(1, 0, noise.sd*D))
##' 
##' @export
E3dist.euclidean = function(dd, samplenames) {
    
    ## simplify the matrix
    dd = dd[,samplenames, drop=FALSE];
    
    ## make a matrix with all distances
    ans = matrix(0, ncol=length(samplenames), nrow=length(samplenames));
    colnames(ans) = samplenames;
    rownames(ans) = samplenames;
    
    for (i in 1:(length(samplenames)-1)) {
        for (j in (i+1):length(samplenames)) {      
            nowdist = sum((dd[,i]- dd[,j])*(dd[,i]-dd[,j]));
            ##if (noise.sd>0) {
            ##    nowdist = nowdist + abs(rnorm(1, 0, nowdist*noise.sd))
            ##}
            ans[i,j] = nowdist;
            ans[j,i] = nowdist;
        }
    }
    ans = sqrt(ans)
    
    ## return an object with distances only
    return(as.dist(ans));   
}




##' Compute a score summarizing a fold-change and a z-value into one number [-1, 1]
##'
##' This is a suggested scoring function that collapses two informative numbers, a fold
##' change conveying effect size and a z-score conveying significance, into a single number
##' in a fixed range [-1, 1]
##'
##' @param fc numeric. Values of fold change (e.g. 2 for 2-fold up-regualtion,
##' 0.25 for 4-fold down-regulation). 
##' @param z numeric. Values of z-score
##' @param fc.thresh numeric vector of length two. Should hold lower and upper thresholds for
##' interesting fold-changes. Values of fc below lower threshold will receive score 0. Values above
##' upper threshold will receive score 1. Intermediate values follow linear interpolation.
##' @param z.thresh numeric vector of length two. Analogous to fc.thresh, but here thresholds
##' applied to z score.
##'
##' @export
E3DEscore = function(fc, z, fc.thresh=c(1.25,1.75), z.thresh=c(1.25,1.75)) { 
    fc[fc==0]=1;
    fc = exp(abs(log(fc)));
    getScore = function(x, xlim) {
        x = (x - xlim[1])/(xlim[2]-xlim[1]);
        x[x<0]=0;
        x[x>1]=1;
        return(x);
    }
    return(sign(z)*getScore(fc, fc.thresh)*getScore(abs(z), z.thresh));
}
