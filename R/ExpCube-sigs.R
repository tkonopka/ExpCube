
##
## ExpCube 
## Author: Tomasz Konopka
##
##
## Functions for manipulating and comparing signatures
##




##' Evaluate various metrics comparing sample and group signatures
##'
##' @param sampleGroups - named list of character vectors. Each vector should represent
##' a set of samples making up a group.
##' @param sample.sigs - named list of character vectors. Each vector should represent
##' a gene set (signature) associated with a sample.
##' @param group.sigs - named list of character vectors. Each vector should represent
##' a gene set (signature) associated with a group.
##' 
##' @export
E3MeasureGroupConsistency = function(sampleGroups, sample.sigs, group.sigs) {

    ans = data.frame(Group=names(sampleGroups), Replicates=0, Jaccard=0, Size=0, stringsAsFactors=F)
    rownames(ans) = names(sampleGroups)

    for (nowg in names(sampleGroups)) {

        nowgsig = group.sigs[[nowg]]
        
        ## for each sample in the group, compute the jaccard index with the
        ## group consensus signature
        nowsamples = sampleGroups[[nowg]]
        nowj = sapply(sample.sigs[nowsamples],
            function(x) {
                E3GetJaccard(x, nowgsig)
            })                
        ans[nowg, "Jaccard"] = mean(nowj)
        rm(nowj)
        
        ans[nowg, "Size"] = length(group.sigs[[nowg]])
        ans[nowg, "Replicates"] = length(nowsamples)
        
    }
    
    return(ans)    
}




##' Compute a distance matrix between gene sets based on the Jaccard index.
##'
##' The function will compute pairwise Jaccard Index (JI) between gene sets. The distance
##' between the gene sets will be 1-JI.
##' 
##' @param sigs - named list of character vectors. Each vector should represent a
##' gene set. The function will compute pairwise Jaccard Indexes between these vectors.
##' 
##' @export
E3GetJaccardDistance = function(sigs) {
    
    ans = matrix(0, ncol=length(sigs), nrow=length(sigs))
    
    for (i in 1:length(sigs)) {
        sigi = sigs[[i]]
        for (j in 1:length(sigs)) {
            sigj = sigs[[j]]
            ans[i,j] = 1-E3GetJaccard(sigi, sigj)
        }
    }
    rownames(ans) = names(sigs)
    colnames(ans) = names(sigs)

    ## eliminate NaN that can arise when signatures are empty
    ans[!is.finite(ans)] = 1
    
    return(ans)    
}




##' Collapse a per-sample table of covariates into a per-group table
##'
##' Function outputs a table with as many columns as covariates, as many rows as
##' groups defined in sampleGroups.
##' 
##' @param covariates - data frame. Table with samplenames in rownames, then columns
##' representing categorial and/or numeric covariates
##' 
##' @param sampleGroups - named list of character vectors. Each vector represents the
##' samples that make up one group
##' 
##' @export
E3MakeGroupCovariates = function(covariates=NULL, sampleGroups=NULL) {
    
    ans = list();
    for (nowgroup in names(sampleGroups)) {
        nowcov = covariates[sampleGroups[[nowgroup]],,drop=FALSE];

        for (j in colnames(nowcov)) {
            jclass = class(nowcov[,j]);
            if (jclass=="numeric" | jclass=="integer") {
                nowcov[,j] = median(nowcov[,j], na.rm=TRUE);                
            } else {
                if (length(unique(nowcov[,j]))>1) {
                    nowcov[,j] = "Complex";
                } 
            }            
        }
        nowcov = nowcov[1,];
        rownames(nowcov) = nowgroup;
        ans[[nowgroup]] = nowcov;
    }
    
    ans = data.frame(rbindlist(ans), stringsAsFactors=F);
    rownames(ans) = names(sampleGroups);
    return(ans);
}




##' Create GLM models linking Overlap with an expected signature to technical covariates.
##'
##' This function is rather particular for ExpCube. It creates a series of GLM models,
##' one for each "Expected" signature type. 
##' The function outputs a list of glm models and a new table analogous to sigmetrics.
##' 
##' @param sigmetrics - a data frame with columns "Group", "Expected", "Overlap"
##' @param group.covariates -
##' @param glm.covariate - character. Name of column in group.covariates to use
##' in the glm model
##' @param glm.readcounts - logical. Set TRUE to include a term in the glm model
##' for read depth.
##' @param depth.col - column identifier. Column in sigmetrics holding depth data.
##' 
##' @export
E3MakeGLMs = function(sigmetrics, group.covariates=NULL,
    glm.covariate=NULL, glm.readcounts=TRUE, depth.col="Depth_M_reads") {
    
    sigmetrics[,"GLM"] = 0;
    metrics.s = split(sigmetrics, sigmetrics[,"Expected"])
    
    ans = list(metrics=sigmetrics, models=list());
    
    ## now compute a score based on expected signature GLM
    for (nows in names(metrics.s)) {
        if (sum(is.finite(metrics.s[[nows]][,"Overlap"]))>2) {
            ## create a small data frame with the numbers for the glm
            nowxyz = metrics.s[[nows]][,c("Group", "Overlap")]
            if (glm.covariate %in% colnames(group.covariates)) {
                nowxyz[,"Covariate"] = group.covariates[nowxyz[,"Group"],glm.covariate]
            }
            nowxyz[,"Depth"] = group.covariates[nowxyz[,"Group"], depth.col]
            
            ## create the glm model
            if (glm.covariate %in% colnames(group.covariates)) {
                if (glm.readcounts) {
                    nowglm = glm(Overlap~Covariate+Depth, data=nowxyz)
                } else {
                    nowglm = glm(Overlap~Covariate, data=nowxyz)
                }                
            } else {
                if (glm.readcounts) {
                    nowglm = glm(Overlap~Depth, data=nowxyz)
                } else {
                    nowglm = list(residuals=rep(0, nrow(metrics.s[[nows]])))
                }                
            }
            
            ans$models[[nows]] = nowglm;
            metrics.s[[nows]][,"GLM"] = nowglm$residuals
            
        }
    }
    ans$metrics = data.frame(rbindlist(metrics.s), stringsAsFactors=F)
    rownames(ans$metrics) = ans$metrics[,"Group"]
    
    return(ans)
}





##' take two vectors A, B
##' returns a list with elements in A only, B only, and both A and B.
##' 
##' @param A - vector. First vectors to compare
##' @param B - vector. Second vector to compare
##' 
##' @export
E3CompareVectors = function(A=c(), B=c()) {
  ans = list();
  ans$A = A[!(A %in% B)];
  ans$B = B[!(B %in% A)];
  ans$both = A[A %in% B];
  ans$summary = c(A=length(ans$A), B=length(ans$B), both=length(ans$both));
  return(ans);  
}


##' Compute Jaccard index of two vectors (size of intersection divided by sum of union)
##'
##' @param x - vector
##' @param y -vector
##'
##' @export
E3GetJaccard = function(x, y) {
    temp = E3CompareVectors(x, y)$summary;
    return(as.numeric(temp["both"]/sum(temp)));    
}

##' Compute overlap (number of shared elements) between two vectors x, y
##' 
##' @param x - vector
##' @param y - vector
##' 
##' @export
E3GetOverlap = function(x, y) {
    temp = E3CompareVectors(x, y)$summary;
    return(as.numeric(temp["both"]));
}

##' Compute relative overlap (number of shared elements divided by size of first)
##' between two vectors x, y
##'
##' @param x - vector
##' @param y - vector
##' 
##' @export
E3GetRelativeOverlap = function(x, y) {
    temp = E3CompareVectors(x, y)$summary;
    return(as.numeric(temp["both"]/length(x)));    
}


##' Compute fisher p-value: selection of two vectors from a background set
##'
##' @param x - vector
##' @param y - vector
##' @param background - large vector from which x and y were picked
##' 
##' @export
E3GetFisherP = function(x,y,background) {
    temp1 = sum(x%in%y)
    temp2 = length(x)-temp1;
    temp3 = length(y)-temp1;
    temp4 = length(background[!(background %in% c(x, y))])
    return(fisher.test(matrix(c(temp1, temp2, temp3, temp4), nrow=2))$p.value)
}




##' Collect metrics (Jaccard, Overlap, Fisher) describing to what extent observed
##' gene signatures correspond to expected signatures.
##'
##' Metrics include overlap between observed and expected signatures, JI between observed
##' and expected signatures, and Fisher p value for selection of observed and expected
##' signatures from a background set.
##'
##' Warning: the Fisher p value is computed very crudely using a background set of fixed size.
##' This is in principle flawed, but the whole idea of using a Fisher test in such a scenario
##' is flawed because not all genes are equally powered to be part of signatures. Bottom line
##' is that the Fisher p values should be interpreted with care, they may be off even by
##' two/three orders of magnitude. 
##' 
##' @param sampleGroups - named list of character vectors. Each vector should represent
##' a set of samples making up a group.
##' @param group.signature - named list of character vectors. Each vector should represent
##' the gene set (signature) associated with a group
##' @param signature.ref - named character vector. Each element should link a sample name
##' to a reference group, e.g. KO_ACTA_R1 would be linked to group KO_None
##' @param signature.expected - named character vector. Each element should a sample
##' to an expected signature group, e.g. KO_ACTA_R1 would be linked to group WT_ACTA
##' @param numgenes - number of genes to use as background in calculation of Fisher metric.
##' 
##' 
##' @export
E3CompareSignatures = function(sampleGroups, group.signature=NULL, signature.ref=NULL,
    signature.expected=NULL, numgenes=5000) {
    
    ans = data.frame(Group=names(sampleGroups), Expected=NA, Jaccard=0,
        Overlap=0, Fisher=1, stringsAsFactors=F);
    rownames(ans) = names(sampleGroups);
    
    ## For each group, record Jaccard and Overlaps scores
    for (nowgroup in names(sampleGroups)) {
        nowcompare = E3GetComparisonGroups(nowgroup,
            sampleGroups, signature.ref, signature.expected)        
        nowgenes = group.signature[[nowgroup]];
        expgenes = unique(unlist(group.signature[nowcompare$expectedgroups]));
        
        ans[nowgroup, "Expected"] = paste(nowcompare$expectedgroups, collapse=", ")
        ans[nowgroup, c("Jaccard", "Overlap")] =
            c(E3GetJaccard(nowgenes, expgenes), E3GetRelativeOverlap(expgenes, nowgenes))
        
        nowcommon = sum(nowgenes %in% expgenes)
        nowfisher = matrix(
            c(nowcommon, length(nowgenes)-nowcommon,
              length(expgenes)-nowcommon,
              numgenes - length(unique(c(nowgenes, expgenes)))), nrow=2)        
        ans[nowgroup, "Fisher"] = fisher.test(nowfisher)$p.value
        
    }
    
    ## clean up the table a little and save it 
    for (zz in c("Jaccard", "Overlap", "Fisher")) {
        ans[,zz] = signif(ans[,zz], 3);        
    }    
    
    return(ans);
}




