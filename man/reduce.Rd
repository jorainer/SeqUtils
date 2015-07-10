\name{reduce}
\Rdversion{1.1}
\docType{methods}
\alias{reduce}
\alias{reduce,RangedSummarizedExperiment-method}
\alias{reduce,SummarizedExperiment0-method}

\title{Reducing a SummarizedExperiment by combinding replicates}
\description{
  This function combines values in the \code{assays} slot for replicates
  e.g. of counts of the same biological sample sequenced on different
  lanes.
}
\usage{

\S4method{reduce}{RangedSummarizedExperiment}(x, by, collapse=";", FUN=rowSums)

\S4method{reduce}{SummarizedExperiment0}(x, by, collapse=";", FUN=rowSums)

}
\arguments{

  \item{x}{
    The \code{RangedSummarizedExperiment} or
    \code{SummarizedExperiment0} (defined in the
    \code{SummarizedExperiment} package) that should be reduced.
  }

  \item{by}{
    Either a character of length 1 specifying the column in
    \code{colData(x)} or a vector of the same length than
    \code{colnames(x)} that defines the replicates. Replicates that
    should be merged are assumed to be identified by e.g. the same value
    in \code{by}.
  }

  \item{collapse}{
    A character of length 1 to be used as a separator for values in
    columns of \code{colData(x)} that are not unique for the same
    replicate.
  }

  \item{seed}{
    The seed for the random number generator.
  }

  \item{x}{
    The BAM file name or a \code{BamFile} object.
  }

}
\section{Methods and Functions}{
  \describe{
    \item{sampleBam}{
      the function samples the specified fraction of reads from the Bam
      file. The returned object is a list of lists. For more details
      please see the help of the \code{\link[Rsamtools]{scanBam}} function.
    }
  }
}
\details{
  This function combines values from replicates in a
  SummarizedExperiment identified by the argument \code{by}. By default,
  for each row, the values of replicates are summed (assuming count
  data). The argument \code{FUN} however allows to apply a different
  aggregation function. Also, by default, the same aggregation function
  is applied to all assays in \code{x}. Different aggregation functions
  can be applied by specifying several functions, where the number of
  specified functions has to match \code{length(assays(x))}; functions
  are applied in sequential order, i.e. the first function is applied to
  the first assay (\code{assays(x)[1]}), the second to the second and
  so on.

  In addition to the combination of assay data, also the \code{colData}
  is reduced to fit the number of unique replicates
  (\code{length(unique(by))}). Values in columns with non-unique values for the
  same replicate are collapsed using the characted defined by
  \code{collapse} as separator between the unique values. The original
  colnames of the object are collapsed in a column added to the
  \code{colData} named \code{"collapsed"}.
}
\value{
  A \code{SummarizedExperiment0} or \code{RangedSummarizedExperiment}.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link[SummarizedExperiment]{SummarizedExperiment}}
}
\examples{

nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     name=c("ChIP.a", "Input.a", "ChIP.b",
                            "Input.b", "ChIP.a", "Input.a"),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            colData=colData)

## sum assay values of replicates identified by column "name"
## in the colData.
ser <- reduce(se0, by="name")
colData(ser)
colnames(ser)

## alternatively use a vector of length=length(colnames(se0))
set <- reduce(se0, by=c(1, 2, 1, 2, 3, 4))

}
\keyword{classes}

