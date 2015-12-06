##
## ExpCube 
## Author: Tomasz Konopka
##
## Functions for working with fold change matrices
##




## #######################################################################################
## Functions for analysis involving fold changes


## Helper function - obtain a mean fold change between a set of samples and 
## sampleexp - a data frame with colnames (SAMPLENAME) [can contain other columns, e.g. gene]
## samplist - a list with samplenames and then vectors of other samplenames to compare to
##          e.g. sample1 -> refsample1, refsample2
##               sample2 -> refsample2, refsample3
##              (note how reference samples can be the same or different each time)
##
## returns - a vector of average fold changes
##
E3GetTwoGroupFC = function(sampleexp, samplist, fun=mean) {
    
    ans = matrix(0, ncol=length(samplist), nrow=nrow(sampleexp))
    colnames(ans) = names(samplist);
    
    ## for each sample, compute reference samples,
    ## obtain "mean" FC wrt to each of the reference samples
    for (ns in names(samplist)) {
        nowref = samplist[[ns]]
        temp = sampleexp[,ns]/sampleexp[,nowref,drop=FALSE];
        ans[,ns] = apply(temp, 1, fun);
    }
    
    ## obtain mean FC over each of the nowsamples
    return(apply(ans, 1, fun))
}



## Helper function
## for a given group name, identify other groups that are worth comparing to
##
## returns:
## $group - the input groupname
## $refgroups - the direct comparison group as given by signature.ref
## $expectedgroups - the groups with similar expected gene expression
## $expectedrefgroups - the expected groups are relative to their own comparison groups
##
E3GetComparisonGroups = function(groupname, sampleGroups, signature.ref, signature.expected) {
    
    ans = list()
    ans$group = groupname;
    nowsamples = sampleGroups[[groupname]];
    ans$refgroups = unique(signature.ref[nowsamples]);
    
    ans$expectedgroups = unique(signature.expected[nowsamples])
    expected.samples = as.character(unlist(sampleGroups[ans$expectedgroups]))
    ans$expectedrefgroups = unique(signature.ref[expected.samples])
    
    return(ans);
}




##' Produce a score based on whether a genes is consistent with expected behavior (fold changes)
##'
##' Function returns a list with two tables. Values in one table are expected fold
##' changes. Values in the other tables are observed fold changes. See descriptions
##' of arguments to see how these are computed. Each table will have gene names in rownames
##' and group names in colnames.
##' 
##' @param sampleexp - data frame with columns (Gene, SAMPLENAMES)
##' @param sampleGroups - named list of character vectors (sample names). Names indicate
##' sample grouping, e.g. group WT_None3 containing samples c("WT_None3_R1", "WT_None2_R2")
##' @param signature.ref - named character vector. This should link samplenames to
##' a group with reference samples, e.g. sample KO_ACTA_R1 will link to a reference group
##' KO_None.
##' @param signature.expected - names character vector. This should link samplenames to
##' expected signature groups, e.g. sample KO_ACTA_R1 will be expected to be similar to
##' group WT_ACTA.
##' @param fun - function. Used to summarize data from multiple samples into one group value.
##' The mean summary function is usually appropriate and fast.
##' @param add.dark - numeric. a dark count for the expression values. This is useful to avoid
##' division by zeros.
##' @param id.col - name of column in sampleexp with gene identifier
##' 
##' @export
E3GetObservedExpectedGroupFC = function(sampleexp, sampleGroups,
    signature.ref, signature.expected, fun=mean, add.dark=1, id.col="Gene") {
    
    ## The answer object will contain two tables - expectedFC, observedFC
    ans = list()
    ## make a table scoring each gene for each group
    temp = matrix(0, ncol=length(sampleGroups), nrow=nrow(sampleexp))
    colnames(temp) = names(sampleGroups);
    rownames(temp) = sampleexp[,id.col];
    ans$observed = temp;
    ans$expected = temp;
    rm(temp)
    
    ## add a dark count to all the values in the sample expression matrix
    for (i in 1:ncol(sampleexp)) {
        if (class(sampleexp[,i])%in% c("integer", "numeric")) {
            sampleexp[,i] = add.dark + sampleexp[,i];
        }
    }
    rownames(sampleexp) = sampleexp[,id.col]
    
    ## shorthand for the number of rows in the answer matrix
    numgenes = nrow(sampleexp)

    ## make a cache of expected.FC
    expected.FC.cache = list();
    
    ## process each group in turn
    for (nowsamples.group in names(sampleGroups)) {
        cat(".")
        ## get names of groups that relevant for comparisons
        nowcompare = E3GetComparisonGroups(nowsamples.group,
            sampleGroups, signature.ref, signature.expected)
        
        ## get KO+stim samples
        nowsamples = sampleGroups[[nowsamples.group]]
        ## for each sample, get a list of their 
        nowsamples.refgroups = nowcompare$refgroups
        nowsamples.refsamples = lapply(as.list(nowsamples),
            function(x) {
                return(sampleGroups[[signature.ref[x]]])
            });
        names(nowsamples.refsamples) = nowsamples
        
        ## get WT+stim signature samples
        expected.groups = nowcompare$expectedgroups;
        expected.samples = as.character(unlist(sampleGroups[expected.groups]))        
        ## get WT+none samples
        expected.refgroups = nowcompare$expectedrefgroups;
        expected.refsamples = lapply(as.list(expected.samples),
            function(x) {
                return(sampleGroups[[signature.ref[x]]])
            })
        names(expected.refsamples) = expected.samples        
        
        nowsamples.refsamples = nowsamples.refsamples[!sapply(nowsamples.refsamples, is.null)]
        
        ## compute fold changes in this group and in the expected comparison
        nowsamples.FC = E3GetTwoGroupFC(sampleexp, nowsamples.refsamples, fun=fun)
        ## compute fold changes for the expected groups (use cache)
        cache.index = paste(c("CI",names(expected.refsamples)),collapse=":");
        if (cache.index %in% names(expected.FC.cache)) {
            expected.FC = expected.FC.cache[[cache.index]];
        } else {
            expected.FC = E3GetTwoGroupFC(sampleexp, expected.refsamples, fun=fun)
            expected.FC.cache[[cache.index]] = expected.FC;
        }
        
        ## record a measure of the discrepancy between the fold changes
        ans$observed[,nowsamples.group] = nowsamples.FC
        ans$expected[,nowsamples.group] = expected.FC
        
        rm(nowsamples.FC, expected.FC)
    }
    cat("\n")    
    
    return(ans);
}




