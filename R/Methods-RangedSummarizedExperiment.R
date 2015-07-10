## reduce a RangedSummarizedExperiment by summing the counts of replicates runs of the same
## sample (e.g. on different lanes)
setMethod("reduce", "RangedSummarizedExperiment", function(x, by, collapse=";", FUN=rowSums){
    if(missing(by))
        stop("Argument by is required!")
    newX <- doReduce(x, by=by, collapse=collapse, FUN=FUN)
    ## a little heavy copying but doesn't work otherwise...
    newerX <- as(newX, "RangedSummarizedExperiment")
    rm(newX)
    rowRanges(newerX) <- rowRanges(x)
    return(newerX)
})

setMethod("reduce", "SummarizedExperiment0", function(x, by, collapse=";", FUN=rowSums){
    if(missing(by))
        stop("Argument by is required!")
    newX <- doReduce(x, by=by, collapse=collapse, FUN=FUN)
    return(newX)
})

doReduce <- function(x, by, collapse=";", FUN=rowSums){
    ## if length by matches length of colData then use that.
    cd <- colData(x)
    cn <- colnames(x)
    if(length(by)==1){
        ## if length by is 1 assume its a factor in colData
        if(any(colnames(cd)==by)){
            by <- as.character(cd[, by])
        }else{
            stop("If by is of length one it should correspond to a column name in colData(x) that can be used to merge/reduce replicates!")
        }
    }
    if(length(by)!=length(cn))
        stop("length of by does not match the number of columns!")
    ## OK now check the assays ...
    nassays <- length(assays(x))
    if(nassays==0)
        stop("No assays in the object! There is nothing to be done...")
    assaynames <- names(assays(x))
    newList <- vector("list", nassays)
    ## process the assays...
    if(class(FUN)!="list")
        FUN <- list(FUN)
    if(length(FUN)!=nassays){
        if(length(FUN)==1){
            ## that's trivial... just apply this one...
            FUN <- rep(FUN, nassays)
        }else{
            warning("length of FUN does not match the number of assays. Applying the first FUN to all assays!")
            FUN <- rep(FUN[1], nassays)
        }
    }
    for(i in seq_along(assaynames)){
        tmp <- do.call(cbind, lapply(unique(by), function(z){
            ## process the assay with the appropriate function...
            return(FUN[[i]](assays(x)[[i]][, by==z, drop=FALSE], na.rm=TRUE))
        }))
        colnames(tmp) <- unique(by)
        newList[[i]] <- tmp
    }
    names(newList) <- assaynames
    ## done with the assays, not process the colData...
    cd <- as.data.frame(cd)
    dataTypes <- unlist(lapply(cd, class))
    ## define the conversion functions...
    dataFuns <- paste0("as.", dataTypes)
    cd <- cbind(cd, collapsed=rownames(cd))
    res <- aggregate(cd, by=list(by), FUN=function(z){
        return(paste(unique(z), collapse=collapse))
    })
    rownames(res) <- res[, "Group.1"]
    res <- res[, -1]
    res <- res[unique(by), ]
    ## all those guys that do not contain a ; transform to what they are...
    for(i in 1:length(dataTypes)){
        if(length(grep(res[, i], pattern=collapse)) == 0){
            res[, i] <- do.call(dataFuns[i], list(res[, i]))
        }
    }
    res <- droplevels(res)
    cd <- as(res, "DataFrame")
    newX <- SummarizedExperiment(assays=SimpleList(newList),
                                 colData=cd)
    metadata(newX) <- metadata(x)
    return(newX)
}
