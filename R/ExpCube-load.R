##
## ExpCube analysis
## Author: Tomasz Konopka
##
## Here are functions to load data from disk into R objects
##



##' Load sample information from csv files.
##'
##' This function reads one or more csv files and concatenates them into a single data frame.
##' The function assumes the files contain certain columns, and ignore all others. 
##' 
##' @param filelist - character vector representing filenames. Function reads tables
##' from each file and merges them into one data frame.
##'
##' @import Rcssplot data.table sm DNAcopy
##' @export
E3LoadPlateConfig = function(filelist) {
    ans = list();
  for (i in 1:length(filelist)) {
      nowfile = filelist[i];
      cat(nowfile, "\n");
      temp = read.table(nowfile, header=T, stringsAsFactors=F);
      temp = temp[,c("Sample","Group","RefGroup","ExpectedSig","Stimulus","Lane",
          "Plate.96well","Row.96well","Column.96well",
          "RNA.concentration.ng.ul","X260.280","X260.230","Library.concentration.ng.ul","ERCC")];
      ans[[i]] = temp;
  }
  ans = data.frame(rbindlist(ans), stringsAsFactors=F);
  rownames(ans) = ans[,"Sample"];
  return(ans);
}





##' Load expression data for samples from Exp3p-like data files.
##'
##' @param samplenames - character vector with sample codes e.g. KO_None. The sample codes
##' must correspond to parts of the data path.
##' @param platenames - character vector with additional codes for each sample, e.g. a plate id.
##' These codes must correspond to parts of the data path.
##' @param path.template - character string. This is a generic path that can contain patterns
##' SAMPLE and PLATE. These substrings will be replaced by values from the previous two vectors.
##' @param column.template - character string marking the column in the data table containing
##' values to be interpreted as expression values. The pattern SAMPLE will be replaced by
##' values from the samplenames vector.
##' @param column.low.template - character string. Similar to column.template, but referring
##' to columns with a low estimate for expression values
##' @param column.high.template - character string. Similar to column.template, but referring
##' to columns with a high estimate for expression values
##' @param id.in - character string. Name of column in input data files containing the
##' gene or region identifier.
##' @param id.out - character string. Name of column in the output data frame containing
##' gene identifier
##' @param rescale.uncertainty - numeric. When different from one, the low/high estimates
##' will be modified by this factor. Use this value to read more conservative
##' expression intervals.
##' @param verbose - logical. Print out progress messages.
##'
##' @export
E3LoadExp3pGeneData = function(samplenames, platenames, path.template,
    column.template="SAMPLE.EPM",
    column.low.template="SAMPLE.EPM.low",
    column.high.template="SAMPLE.EPM.high",
    id.in="Region.id", id.out="Gene",
    rescale.uncertainty=1,
    verbose=TRUE) {
    
    ## if the processed data file does not exist, create it from scratch
    ## First, put simple tables into items in a list
    ans = list();
    
    previousfile = "";
    nowfile = "";
    nowdata = NULL;
    
    if (verbose) {
        cat("loading exp3p data\n");
    }
    
    for (i in 1:length(samplenames)) {
        nowsample = samplenames[i];
        nowplate = platenames[i];
        ## load the data for this sample from disk (or reuse previous table if loaded
        ## on previous step
        nowfile = gsub("PLATE", nowplate, gsub("SAMPLE",nowsample,path.template));
        if (verbose) {
            cat("loading from: ",nowfile,"\n");
        }
        if (!identical(nowfile, previousfile)) {
            nowdata = read.table(nowfile, header=T, sep="\t", stringsAsFactors=FALSE);
        }
        
        ## get names of columns
        nowcol = gsub("SAMPLE",nowsample, column.template);
        nowcol.low = gsub("SAMPLE",nowsample,column.low.template);
        nowcol.high = gsub("SAMPLE",nowsample,column.high.template);
        
        ## perhaps rename some columns with sample name starts with a number"
        nowcol.X = which(colnames(nowdata)==paste("X",nowcol,sep=""));
        nowcol.low.X = which(colnames(nowdata)==paste("X",nowcol.low,sep=""));
        nowcol.high.X = which(colnames(nowdata)==paste("X",nowcol.high, sep=""));    
        colnames(nowdata)[nowcol.X] = nowcol;
        colnames(nowdata)[nowcol.low.X] = nowcol.low;    
        colnames(nowdata)[nowcol.high.X] = nowcol.high;
        
        ## get the transcript and expression table
        temp = nowdata[,c(id.in, nowcol, nowcol.low, nowcol.high)];
        colnames(temp) = c(id.out, paste(nowsample,c("",".low",".high"),sep=""));
        temp[,1] = as.character(temp[,1]);
        for (k in 2:4) {
            temp[,k] = as.numeric(as.character(temp[,k]));
        }
        temp = temp[order(temp[,1]),];
        
        ans[[nowsample]] = temp;
        rm(temp);
        
        previousfile = nowfile;
    }
    
    ## make sure all the samples have the same number of columns and species names
    alltxs = ans[[1]][,id.out];
    for (nowsample in samplenames) {
        if (!identical(alltxs, ans[[nowsample]][,"Gene"])) {
            stop(paste("ids are not identical in samples (",nowsample,")\n",sep=""));
        }
    }
    
    ## make the master table
    sampledata = list();
    sampledata$expression = data.frame(matrix(0, ncol=1+length(samplenames), nrow=length(alltxs)));
    colnames(sampledata$expression) = c(id.out, samplenames);    
    sampledata$expression[,id.out] = as.character(alltxs);
    
    sampledata$expression.low = sampledata$expression;
    sampledata$expression.high = sampledata$expression;
    for (nowsample in samplenames) {
        temp = ans[[nowsample]];
        sampledata$expression[,nowsample] = temp[,nowsample];
        sampledata$expression.low[,nowsample] = temp[,paste0(nowsample,".low")];
        sampledata$expression.high[,nowsample] = temp[,paste0(nowsample,".high")];
    }
    
    ##     
    colnames(sampledata$expression.low) = c(id.out, paste0(samplenames));
    colnames(sampledata$expression.high) = c(id.out, paste0(samplenames));
    
    ## get rid of rows that have duplicate names
    badnames = table(sampledata$expression[,id.out])
    badnames = names(badnames[badnames>1])
    if (length(badnames)>0) {
        for (i in paste0("expression", c("", ".low", ".high"))) {
            sampledata[[i]] = sampledata[[i]][!(sampledata[[i]][,"Gene"] %in% badnames),]
        }
    }
    
    if (verbose) {
        cat("loading gene exp3p data - done\n");
    }
    
    ## finally, output the list with expression estimates
    return(sampledata);      
}



