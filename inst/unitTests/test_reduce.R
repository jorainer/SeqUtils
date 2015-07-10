test_reduce <- function(){
    nrows <- 200
    ncols <- 6
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                         name=c("ChIP.a", "Input.a", "ChIP.b",
                                "Input.b", "ChIP.a", "Input.a"),
                         row.names=LETTERS[1:6])
    se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                colData=colData)
    ## reduce it...
    ser <- reduce(se0, by="name")
    unname <- unique(se0$name)
    for(i in length(unname)){
        tmp <- rowSums(assay(se0)[ , se0$name==unname[i], drop=FALSE])
        checkEquals(tmp, assay(ser)[, i])
    }

    ## now test an object with more assays...
    semulti <- SummarizedExperiment(assays=SimpleList(counts=counts,
                                                      other=counts,
                                                      evenanother=counts),
                                    colData=colData)
    semultir <- reduce(semulti, by="name")
    for(i in length(unname)){
        tmp <- rowSums(assays(semulti)[[1]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[1]][, i])
        ## second
        tmp <- rowSums(assays(semulti)[[2]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[2]][, i])
        ## third
        tmp <- rowSums(assays(semulti)[[3]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[3]][, i])
    }
    ## next testing what happens if we use rowSums, rowMeans and rowMedians
    semultir <- reduce(semulti, by="name", FUN=c(rowSums, rowMeans, rowMedians))
    for(i in length(unname)){
        tmp <- rowSums(assays(semulti)[[1]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[1]][, i])
        ## second
        tmp <- rowMeans(assays(semulti)[[2]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[2]][, i])
        ## third
        tmp <- rowMedians(assays(semulti)[[3]][, semulti$name==unname[i], drop=FALSE])
        checkEquals(tmp, assays(semultir)[[3]][, i])
    }
}

