if(!isGeneric("overview",))
    setGeneric("overview", function(x, ...)
               standardGeneric("overview"))
setMethod("overview", "QA", function(x, ...){
    df <- as(values(x@src), "data.frame")
    return(df)
})



#########################################
##
##   Distribution of average read quality.
##
setMethod("plot", "QAReadQuality", function(x, y, pal=brewer.pal(8, "Set1"), subset=NULL, ...){
    df <- as(values(x), "data.frame")
    if(!is.null(subset)){
        if(!is.numeric(subset))
            stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
        if(any(!(subset %in% as.numeric(unique(df$Id)))))
            stop("The indices submitted with the argument subset are outside of the range of allowed values!")
        ## doing it with split, so the index can be in any ordering...
        df.splitted <- split(df, f=df$Id)
        df.splitted <- df.splitted[as.character(subset)]
        names(df.splitted) <- NULL
        df <- do.call(rbind, df.splitted)
        df <- droplevels(df)
    }
    lvl <- levels(df$Id)
    flag <- lvl[x@flag]
    df$Id <- factor(df$Id, levels=c(lvl[!lvl %in% flag], flag))
    xmin <- min(df$Score)
    ymax <- max(df$Density)
    col <- c(rep("gray", length(lvl) - length(flag)),
             pal[1 + (seq_along(flag) - 1) %% 8])
    plt <-
        xyplot(Density ~ Score, group=Id, df,
               type = "l", xlab = "Average (calibrated) base quality",
               nylab = "Proportion of reads", col = col, strip=FALSE,
               key = list(space = "top",
                   lines = list(col=tail(col, length(flag)), size=3L, lwd=2),
                   text = list(lab=tail(lvl, length(flag))),
                   columns=min(length(col), 10L), cex=.6))
    return(plt)
})



###############################################
##
## Nucleotide count by Cycle
##
## y is not used.
setMethod("plot", "QANucleotideByCycle", function(x, y,
                                                   dnaCol=c("#1F78B4", "#33A02C", "#B2DF8A", "#A6CEE3", "#E31A1C"),
                                                   subset=NULL, ...){
              df <- as(values(x), "data.frame")
              if(!is.null(subset)){
                  if(!is.numeric(subset))
                      stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
                  if(any(!(subset %in% as.numeric(unique(df$Id)))))
                      stop("The indices submitted with the argument subset are outside of the range of allowed values!")
                  ## doing it with split, so the index can be in any ordering...
                  df.splitted <- split(df, f=df$Id)
                  df.splitted <- df.splitted[as.character(subset)]
                  names(df.splitted) <- NULL
                  df <- do.call(rbind, df.splitted)
                  df <- droplevels(df)
              }
              df <- df[ df$Base!="N" & df$Base!="-", ]
              df$Base <- factor(df$Base)
              return(plotNucleotideCountByCycle(x=df, dnaCol=dnaCol, ...))
          })
## plots the nucleotide count by cycle
plotNucleotideCountByCycle <- function(x, dnaCol=c("#1F78B4", "#33A02C", "#B2DF8A", "#A6CEE3", "#E31A1C"), ...){
    if(class(x)!="data.frame")
        stop("x has to be a data.frame!")
    ## check for required colnames:
    if(sum(c("Id", "Cycle", "Base", "Count") %in% colnames(x))!=4)
        stop("One or more required columns: Id, Cycle, Byse, Count are missing in x!")

    xmax <- max(x$Cycle)
    ymax <- log10(max(x$Count))
    plt <-
        xyplot(log10(Count) ~ as.integer(Cycle) | Id,
               group = factor(Base),
               x[order(x$Id, x$Base, x$Cycle),],
               panel = function(..., subscripts) {
                   lbl <- as.character(unique(x$Id[subscripts]))
                   ltext(xmax, ymax, lbl, adj = c(1, 1))
                   panel.xyplot(..., subscripts = subscripts)
               },
               type = "l", col = dnaCol[1:4],
               key = list(
                   space = "top", lines = list(col = dnaCol[1:4]),
                   text = list(lab = levels(x$Base)),
                   columns = length(levels(x$Base))),
               xlab = "Cycle", aspect = 2, strip = FALSE)
        return(plt)
}


######################################################
## base call Quality by Cycle
##
setMethod("plot", "QAQualityByCycle", function(x, y, type="levelplot", subset=NULL, ...){
    return(plotQualByCycle(x=x, type=type, subset=subset, ...))
})

plotQualByCycle <- function(x, type="levelplot", subset=NULL){
    type <- match.arg(type, c("levelplot", "boxplot"))
    calc_means <- function(x, y, z)
        rowsum(y * z, x)/rowsum(z, x)
    calc_quantile <- function(x, y, z, q = c(0.25, 0.5, 0.75)){
        by(list(y, z), x, function(x) {
            scoreRle <- Rle(x[[1]], x[[2]])
            quantile(scoreRle, q)
        })
    }
    calc_boxplot <- function(x, y, z, ...){
        by(list(y, z), x, function(x){
            Vals <- rleboxplot.stats(Rle(x[[1]], x[[2]]), ...)
            return(Vals$stats)
        })
    }
    df <- as(values(x), "data.frame")
    if(!is.null(subset)){
        if(!is.numeric(subset))
            stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
        if(any(!(subset %in% as.numeric(unique(df$Id)))))
            stop("The indices submitted with the argument subset are outside of the range of allowed values!")
        ## doing it with split, so the index can be in any ordering...
        df.splitted <- split(df, f=df$Id)
        df.splitted <- df.splitted[as.character(subset)]
        names(df.splitted) <- NULL
        df <- do.call(rbind, df.splitted)
        df <- droplevels(df)
    }
    Id <- df$Id
    pal <- c("#66C2A5", "#FC8D62")
    lvlPal <- c("#F5F5F5", "black")
    rng <- range(df$Count)
    at <- seq(rng[1], rng[2], length.out = 512)
    np <- length(unique(Id))
    nrow <- ceiling(np/4)
    layout <- c(ceiling(np/nrow), nrow)
    ymin <- min(df$Score)
    if(type=="levelplot"){
        plt <- xyplot(Score ~ Cycle | Id, df,
                      panel = function(x, y, ..., subscripts) {
                          z <- df$Count[subscripts]
                          mean <- calc_means(x, y, z)
                          qtiles <- calc_quantile(x, y, z)
                          sxi <- sort(unique(x))
                          panel.levelplot(x, y, z, subscripts = TRUE, at = at,
                                          col.regions = colorRampPalette(lvlPal))
                          llines(sxi, mean, type = "l", col = pal[[1]], lwd = 1)
                          llines(sxi, sapply(qtiles, "[[", 1), type = "l",
                                 col = pal[[2]], lwd = 1, lty = 3)
                          llines(sxi, sapply(qtiles, "[[", 2), type = "l",
                                 col = pal[[2]], lwd = 1)
                          llines(sxi, sapply(qtiles, "[[", 3), type = "l",
                                 col = pal[[2]], lwd = 1, lty = 3)
                          lbl <- as.character(unique(df$Id[subscripts]))
                          ltext(1, ymin, lbl, adj = c(0, 0))
                      }, ylab = "Quality Score", layout = layout, strip = FALSE)
    }
    if(type=="boxplot"){
        plt <- xyplot(Score ~ Cycle | Id, df,
                      panel = function(x, y, ..., subscripts) {
                          ## x: Cycle
                          ## y: Score
                          ## z: Count
                          ## -> from y and z i can generate a Rle.
                          z <- df$Count[subscripts]
                          mean <- calc_means(x, y, z)
                          qtiles <- calc_quantile(x, y, z)
                          boxplot.vals <- calc_boxplot(x, y, z, do.conf=FALSE)
                          sxi <- sort(unique(x))
                          ## values I get are: 1, 2, 3, 4, 5: lowe wisker, box...
                          ## draw arrows... bummer, use a loop.
                          for(i in 1:length(sxi)){
                              grid.lines(x=unit(rep(sxi[ i ], 2), "native"),
                                         y=unit(c(boxplot.vals[[ i ]][ 1 ], boxplot.vals[[ i ]][ 5 ]), "native"),
                                         arrow=arrow(angle=90, length=unit(0.3, "native"), ends="both"),
                                         gp=gpar(col="grey")
                                        )
                          }
                          ## grid.polyline(x=unit(rep(sxi, each=2), "native"),
                          ##             y=unit(rep(c(10, 20), length(sxi)), "native")
                          ##            )
                          grid.rect(x=unit(sxi, "native"),
                                    y=unit(sapply(boxplot.vals, "[[", 2), "native"),
                                    width=unit(rep(0.6, length(sxi)), "native"),
                                    height=unit(sapply(boxplot.vals, "[[", 4) - sapply(boxplot.vals, "[[", 2), "native"),
                                    just="bottom",
                                    gp=gpar(fill="white", col="grey")
                                   )
                          ## at last plot the median...
                          grid.points(x=unit(sxi, "native"),
                                      y=unit(sapply(qtiles, "[[", 2), "native"),
                                      size=unit(0.6, "native"),
                                      pch=16
                                     )
                          llines(sxi, mean, type = "l", col = pal[[1]], lwd = 2)
                          lbl <- as.character(unique(df$Id[subscripts]))
                          ltext(1, ymin, lbl, adj = c(0, 0))
                      }, ylab = "Quality Score", layout = layout, strip = FALSE)
    }
    return(plt)
}

###################################
## QANucleotideUse
setMethod("plot", "QANucleotideUse", function(x, y, dnaCol=c("#1F78B4", "#33A02C", "#B2DF8A", "#A6CEE3", "#E31A1C"), subset=NULL, ...){
    return(plotNucleotideUse(x=x, dnaCol=dnaCol, subset=subset, ...))
})
plotNucleotideUse <- function(x, dnaCol=c("#1F78B4", "#33A02C", "#B2DF8A", "#A6CEE3", "#E31A1C"), subset=NULL, ...){
    df <- as(values(x), "data.frame")
    if(!is.null(subset)){
        if(!is.numeric(subset))
            stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
        if(any(!(subset %in% as.numeric(unique(df$Id)))))
            stop("The indices submitted with the argument subset are outside of the range of allowed values!")
        ## doing it with split, so the index can be in any ordering...
        df.splitted <- split(df, f=df$Id)
        df.splitted <- df.splitted[as.character(subset)]
        names(df.splitted) <- NULL
        df <- do.call(rbind, df.splitted)
        df <- droplevels(df)
    }
    df$Id <- as.integer(as.character(df$Id))
    plt <-
        dotplot(Id ~ Count|factor(ifelse(df$Nucleotide == "N", "N", "O")),
                group=Nucleotide, df,
                base=df$Nucleotide,
                type = "b", pch = 20, col = dnaCol,
                key = list(space = "top", lines = list(col = dnaCol),
                    text = list(lab = levels(values(x)[["Nucleotide"]])),
                    columns = 5L),
                strip=FALSE, scale=list(relation="free"),
                par.settings=list(layout.widths = list(panel = c(1, 2))))
    return(plt)
}
setMethod("table", "QANucleotideUse", function(x, ...){
    df <- as.data.frame(values(x))
    Res <- split(df, f=df$Id)
    T <- lapply(Res, function(x){
        tmp <- c(x$Count, x$Count * 100 / sum(x$Count))
        names(tmp) <- c(as.character(x$Nucleotide), paste("%", x$Nucleotide))
        return(tmp)
    })
    return(do.call(rbind, T))
})


###################################
## QAQualityUse
plotQualityUse <- function(x, subset=NULL, ...){
    df <- as(values(x), "data.frame")
    if(!is.null(subset)){
        if(!is.numeric(subset))
            stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
        if(any(!(subset %in% as.numeric(unique(df$Id)))))
            stop("The indices submitted with the argument subset are outside of the range of allowed values!")
        ## doing it with split, so the index can be in any ordering...
        df.splitted <- split(df, f=df$Id)
        df.splitted <- df.splitted[as.character(subset)]
        names(df.splitted) <- NULL
        df <- do.call(rbind, df.splitted)
        df <- droplevels(df)
    }
    id <- df[["Id"]]
    q <- df[["Quality"]]
    q <- factor(q, levels=levels(q)[min(as.integer(q)):max(as.integer(q))])
    df[["Quality"]] <- q
    df <- df[order(df$Id, df$Quality),]
    df[["Proportion"]] <-
        with(df, unlist(Map("/",
                            lapply(split(Count, Id), cumsum),
                            lapply(split(Count, Id), sum)),
                        use.names=FALSE))
    pal <-       # brewer.pal(9, "RdYlBu")
        c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
          "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")
    col <- colorRampPalette(pal)(length(levels(q)))
    plt <- dotplot(Id ~ Proportion, group=Quality, df,
                   type = "b", pch = 20, col = col,
                   xlab="Cummulative Proportion",
                   key = list(space = "top",
                       lines = list(col = col, size=3L),
                       text = list(lab = levels(df[["Quality"]])),
                       columns = min(length(col), 10L), cex=.6))
    return(plt)
}
setMethod("plot", "QAQualityUse", function(x, y, subset=NULL, ...){
    return(plotQualityUse(x=x, subset=subset, ...))
})
plotQualityUseDensity <- function(x, xlab="Used base-call quality", ylab="density", col,
                                  pal=brewer.pal(8, "Set1"), lty=1, subset=NULL, ...){
    df <- as(values(x), "data.frame")
    if(!is.null(subset)){
        if(!is.numeric(subset))
            stop("Argument subset has to be a numeric vector representing the indices of the files to plot!")
        if(any(!(subset %in% as.numeric(unique(df$Id)))))
            stop("The indices submitted with the argument subset are outside of the range of allowed values!")
        ## doing it with split, so the index can be in any ordering...
        df.splitted <- split(df, f=df$Id)
        df.splitted <- df.splitted[as.character(subset)]
        names(df.splitted) <- NULL
        df <- do.call(rbind, df.splitted)
        df <- droplevels(df)
    }
    lvl <- levels(df$Id)
    if(missing(col)){
        col <- pal[1 + (seq_along(lvl) - 1) %% 8]
    }
    densities <- lapply(split(df, df$Id), function(x){
        return(density(as.numeric(Rle(x$Score, x$Count))))
    })
    all.x <- do.call(cbind, lapply(densities, function(x){
        x$x
    }))
    all.y <- do.call(cbind, lapply(densities, function(x){
        x$y
    }))
    matplot(all.x, all.y, type="l", xlab=xlab, ylab=ylab, col=col, lty=1, ...)
}
if(!isGeneric("plotDensity",))
    setGeneric("plotDensity", function(x, ...)
               standardGeneric("plotDensity"))

setMethod("plotDensity", "QAQualityUse", function(x, xlab="Used base-call quality", ylab="density",
                                                   col, pal=brewer.pal(8, "Set1"), lty=1, subset=NULL, ...){
    plotQualityUseDensity(x, xlab=xlab, ylab=ylab, col=col, pal=pal, lty=lty, subset=subset, ...)
})


###################################
## QASequenceUse
plotSequenceUse <- function(x, ...){
    tmp <- as(values(x), "data.frame")
    df <- with(tmp, {
        nOccur <- tapply(Occurrences, Id, c)
        cumulative <- tapply(Occurrences * Reads, Id, function(elt) {
            cs <- cumsum(elt)
            (cs - cs[1] + 1)/(diff(range(cs)) + 1L)
        })
        id <- tapply(Id, Id, c)
        data.frame(Occurrences = unlist(nOccur),
                   Cumulative = unlist(cumulative),
                   Id = unlist(id), row.names = NULL)
    })
    ## df <- with(values(x), {
    ##     nOccur <- tapply(Occurrences, Id, c)
    ##     cumulative <- tapply(Occurrences * Reads, Id, function(elt) {
    ##         cs <- cumsum(elt)
    ##         (cs - cs[1] + 1)/(diff(range(cs)) + 1L)
    ##     })
    ##     id <- tapply(Id, Id, c)
    ##     data.frame(Occurrences = unlist(nOccur),
    ##                Cumulative = unlist(cumulative),
    ##                Id = unlist(id), row.names = NULL)
    ## })
    xmax <- log10(max(df$Occurrences))
    plt <- xyplot(Cumulative ~ log10(Occurrences) | factor(Id), df,
                  xlab = expression(paste("Number of occurrences of each sequence (",
                      log[10], ")", sep = "")),
                  ylab = "Cumulative proportion of reads",
                  aspect = 2, panel = function(x, y, ..., subscripts, type) {
                      lbl <- unique(df$Id[subscripts])
                      ltext(xmax, 0.05, lbl, adj = c(1, 0))
                      type <- if (1L == length(x)) "p" else "l"
                      panel.xyplot(x, y, ..., type = type)
                  }, strip = FALSE)
    return(plt)
}
setMethod("plot", "QASequenceUse", function(x, y, ...){
    return(plotSequenceUse(x=x, ...))
})


###################################
## QAFrequentSequence
## return a data.frame with:
## Id: the sample id
## Records: number of reads
## Seq_count: number of times the sequence was present
## Sequence: the sequence
tableFrequentSequence <- function(x){
    ## well, have to guess what's in that object.
    df <- values(x)
    maxNoSeq <- x@n
    ## TopCount is a IntegerList, sequences are names.
    Res <- lapply(split(df, df$Id), function(z){
        Seqnames <- unname(names(z$TopCount[[1]]))
        Seqcount <- unname(z$TopCount[[1]])
        return(data.frame(Id=rep(z$Id, length(Seqnames)),
                           Records=rep(z$Records, length(Seqnames)),
                           Seq_count=Seqcount,
                           Sequence=Seqnames))
    })
    return(do.call(rbind, Res))
}
setMethod("table", "QAFrequentSequence", function(x, ...){
    return(tableFrequentSequence(x))
})


###################################
## QAAdapterContamination
tableAdapterContamination <- function(x){
    df <- values(x)
    return(df)
}
setMethod("table", "QAAdapterContamination", function(x, ...){
    return(tableAdapterContamination)
})


