##
## ExpCube 
## Author: Tomasz Konopka
##
## Functions used to create gene signatures
##


## look at pairs of samples and create a table with log FoldChange and Z
##
## return - two column data frame - logFC and Z (Gene names encoded as rownames)
##
E3MakeFCZ = function(sampledata, sampleA, sampleB, shiftval=0.01, id.col="Gene") {
    
    s.mid = sampledata$expression;
    s.high = sampledata$expression.high;
    s.low = sampledata$expression.low;
    
    sdnow = ((s.high[,sampleA]-s.mid[,sampleA])^2) + ((s.high[,sampleB]-s.mid[,sampleB])^2);
    sdnow = sqrt(sdnow);
    
    nowlogFC = log2((shiftval+s.mid[,sampleB])/(shiftval+s.mid[,sampleA]));
    nowlogFC[!is.finite(nowlogFC)] = 0.0;  
    nowZ = (s.mid[,sampleB]-s.mid[,sampleA])/sdnow;
    
    ans = data.frame(logFC = nowlogFC, z = nowZ, stringsAsFactors=F);
    rownames(ans) = s.mid[,id.col];
    
    return(ans);  
}




##' Compute DE scores for a single sample.
##'
##' Function returns a named vector with scores [-1,1]. Names are genes. Values are
##' scores that indicate whether a gene is DE based on FC and z criteria. 
##' 
##' @param sampledata - list with expression, expression.low, expression.high tables
##' @param sampleA - character string. The sample for which the signature is made
##' @param refsamples - character vector. Reference samples to which sampleA is compared.
##' @param scorefun - function. This function should accept two arguments (fc and z) and return
##' a numeric score.
##' @param rankrescale - logical. Set TRUE to penalize scores in an ad-hoc way for
##' genes that have a high score but lower fc or lower z than many other genes.
##' @param shiftval - numeric. Used when computing FC to avoid division by zero.
##' return - a vector with gene names and then a score [0,1]
##'          Score 1 means the gene is DE wrt all samples in groupB
##' 
##' @export
E3GetDEScores = function(sampledata, sampleA, refsamples,
    scorefun=E3DEscore, rankrescale=TRUE, shiftval=0.1) {
    
    ## make sure comparing a sample with other samples and not to itself
    refsamples = refsamples[!(refsamples %in% sampleA)]
    
    ## make tables holding scores for all comparisons
    numgenes = nrow(sampledata$expression)
    scores = matrix(0, ncol=length(refsamples), nrow=numgenes)
    colnames(scores) = refsamples
    rownames(scores) = sampledata$expression[,"Gene"]
    
    ## evaluate scores for each comparison
    for (nowB in refsamples) {        
        temp = E3MakeFCZ(sampledata, nowB, sampleA, shiftval=shiftval)
        temp[!is.finite(temp[,1]),1] = 0
        temp[!is.finite(temp[,2]),2] = 0
        
        if (rankrescale) {
            ## get rank-based dampening factor for the score
            ## highly ranked items will get a high score, others not
            numnonzero = sum(sampledata$expression[,sampleA] > 0)        
            rank.FC = rank(temp[,1])
            rank.Z = rank(temp[,2])
            rank.FC = 2*apply(cbind(rank.FC, numgenes-rank.FC+1), 1, min)
            rank.Z = 2*apply(cbind(rank.Z, numgenes-rank.Z+1), 1, min)
            rank.scale.FC = (1-(0.5*rank.FC/numnonzero))
            rank.scale.Z = (1-(0.5*rank.Z/numnonzero));       
            rank.scale.FC[rank.scale.FC < 0] = 0
            rank.scale.Z[rank.scale.Z < 0] = 0
            rank.scale = rank.scale.FC * rank.scale.Z
        } else {
            rank.scale = 1
        }
        
        scores[,nowB] = scorefun(2^temp[,1], temp[,2])*rank.scale       
    }
    
    ## average the scores and output
    return(apply(scores, 1, mean))  
}



##' Create a named list with gene signatures from a matrix with scors.
##' 
##' @param scores - a matrix with scores (gene names in rownames and
##' sample/group names in colnames required)
##' @param min.score - minimal score (absolute value). Genes with scores above the
##' threshold will be considered signature genes.
##' 
##' @export
E3MakeSignatureFromScores = function(scores, min.score=0.75) {
    ans = list();  
    allgenes = rownames(scores);
    samplenames = colnames(scores);
    for (i in 1:length(samplenames)) {
        ans[[i]] = allgenes[abs(scores[,i])>min.score];      
    }
    names(ans) = samplenames;  
    return(ans);
}



##' Average scores over several replicates to obtain one score per group
##'
##' @param scores - matrix with scores. Column names should be sample names.
##' @param sampleGroups - named list with character vectors. Each vector will
##' define a set of samples to average into a group
##' 
##' @export
E3MakeGroupSignatureScores = function(scores, sampleGroups) {
  ans = data.frame(matrix(0, ncol=length(sampleGroups), nrow=nrow(scores)));
  colnames(ans) = names(sampleGroups);
  rownames(ans) = rownames(scores);
  for (i in 1:length(sampleGroups)) {
    nowgroup = names(sampleGroups);
    nowsamples = sampleGroups[[i]];
    ans[,i] = apply(scores[,nowsamples, drop=FALSE], 1, mean);
  }
  return(ans);
}






##' Write signatures into gmt files (text format)
##'
##' @param sigs - named list with character vectors
##' @param expressed.genes - vector with all expressed genes (another signature)
##' @param sig.label - character string. Some text used to identify the gene signature within
##' the gmt file.
##' @param file - character string. output filename
##' 
##' @export
E3WriteSignatureGMT = function(sigs, expressed.genes=NULL, sig.label="ExpCube", file="out.gmt") {
    
    ans = rep("", length(sigs));
    for (i in 1:length(sigs)) {
        nowsample = names(sigs)[i];
        nowsig = sigs[[i]];
        ans[i] = paste0(nowsample,"\t",sig.label,":",nowsample,"\t",
               paste(nowsig, collapse="\t"));
    }
    if (!is.null(expressed.genes)) {
        ans[length(ans)+1] = paste0("ExpressedGenes\t",sig.label,":ExpressedGenes\t",
               paste(expressed.genes,collapse="\t"));
    }
    write(ans, file=file);
    
}



##' Read GMT signatures from a file and produce a list of vectors with gene names
##'
##' @param f - filename
##' 
##' @export
E3LoadGMT = function(f) {

  ## read the file from the disk
  fcon = file(f, open="r");
  flines = readLines(fcon);
  close(fcon);

  ## parse the files
  flines = strsplit(flines, "\t");
  ansempty = c("AA");
  ansempty = ansempty[ansempty=="b"];
  ans = lapply(flines,
    function(x) {
      if (length(x)>2) {
        return(x[3:length(x)]);
      } else {
        return(ansempty);
      }
    });  
  names(ans) = sapply(flines, function(x) {return(x[1])} );  
  return(ans);
}

