\name{summarizeAlignment}
\alias{summarizeAlignment}
\title{
  Summarize alignment results
}
\description{
  This function summarized the results from an alignment by reading
  through the corresponding BAM file and counting the number of uniquely
  aligned or multi-mapped reads.
}
\usage{
summarizeAlignment( x, index=character(), sbp=ScanBamParam( what=c( "flag", "mapq" ), tag=c( "NH", "HI" ) ), yieldSize=1e+6, v=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{
    The file name of the BAM file to summarize.
  }

  \item{index}{
    The file name of the index file for the BAM (usually the BAM
    filename with suffix \code{".bai"}).
  }

  \item{sbp}{
    The \code{ScanBamParam} specifying which fields to read.
  }

  \item{yieldSize}{
    The number of entries that should be sequentially processed.
  }

  \item{v}{
    verbosity.
  }

}
\details{
  This function scans throught the full BAM file and records the number
  of uniquely or multi-aligned reads along with their mapping quality.
}
\value{
  A list with 3 elements:
  \itemize{
    \item{\code{flag}}{
      a list of 3 matrices \code{no_algn}, \code{unique_algn} and
      \code{multi_algn} with rows corresponding to the BAM flags and
      (one) column with the number of reads for a certain flag.
    }
    \item{mapq}{
      a list of 3 matrices \code{no_algn}, \code{unique_algn} and
      \code{multi_algn} with 256 rows corresponding to all possible
      mapping qualities (Sanger Phred encoded; from 0 to 255) and
      (one) column with the number of reads with a certain map quality.
    }
    \item{counts}{
      a list of 3 numeric values \code{no_algn}, \code{unique_algn} and
      \code{multi_algn} representing the total number of non-,
      unique- or mult-mapping reads, respectively.
    }
  }
}
\author{
  Johannes Rainer
}

\examples{

dirpath <- system.file( package="SeqUtils", "extdata" )
bamfile <- dir( dirpath, "example.bam", full=TRUE )

Res <- summarizeAlignment( bamfile )
Res$counts

}
\keyword{ data }
