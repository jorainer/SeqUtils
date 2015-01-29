## what do we do:
## sum of reads with NH==0 -> not aligned.
## sum of reads with NH==1 -> unique alignment.
## sum of reads with NH>1 -> multiple alignments; count only those with HI==1 (first alignment)
## mapq of the 3 categories above.
## flag summaries of the 3 categories above.
summarizeAlignment <- function( x, index=character(), sbp=ScanBamParam( what=c( "flag", "mapq" ), tag=c( "NH", "HI" ) ), yieldSize=1e+6, v=FALSE ){
    ## open the bamfile.

    stream <- open( BamFile( x, index=index, obeyQname=TRUE, yieldSize=yieldSize ) )
    on.exit( close( stream ) )
    ## define result matrices: na (no alignment), ua (unique alignment), ma (multiple alignments)
    FC.na <- matrix( ncol=1, nrow=length( FLAG_BITNAMES ), 0 )  ## flag count
    rownames( FC.na ) <- FLAG_BITNAMES
    SPL <- unlist( strsplit( x, split=.Platform$file.sep, fixed=TRUE ) )
    colnames( FC.na ) <- SPL[ length( SPL ) ]
    FC.ua <- FC.na
    FC.ma <- FC.ua
    ## map quality matrix. expect values from 0 to 256 (or from 1 to 8???)
    MAPQ.na <- matrix( ncol=1, nrow=256, 0 )
    colnames( MAPQ.na ) <- colnames( FC.na )
    rownames( MAPQ.na ) <- 0:255
    MAPQ.ua <- MAPQ.na
    MAPQ.ma <- MAPQ.na
    ## counter for the alignments...
    counts.na <- 0
    counts.ua <- 0
    counts.ma <- 0
    counta <- 0
    repeat{
        res <- scanBam( stream, param=sbp )
        if( length( res[[1]]$flag ) == 0 )
            break
        ## unique aligning reads:
        bm <- res[[1]]$tag$NH==1
        bm[ is.na( bm ) ] <- FALSE
        if( any( bm ) ){
            flagcounts <- colSums( bamFlagAsBitMatrix( res[[1]]$flag[ bm ] ) )
            FC.ua <- FC.ua + flagcounts[ rownames( FC.ua ) ]
            ## process mapping quality:
            MAPQ.bam <- table( res[[1]]$mapq[ bm ] )
            MAPQ.ua[ as.character( names( MAPQ.bam ) ), 1 ] <- MAPQ.ua[ as.character( names( MAPQ.bam ) ), 1 ] + as.numeric( MAPQ.bam )
            counts.ua <- counts.ua + sum( bm )
        }
        ## multi aligning reads; use only first aligned read...
        bm <- res[[1]]$tag$NH>1 & res[[1]]$tag$HI==1
        bm[ is.na( bm ) ] <- FALSE
        if( any( bm ) ){
            ## Note: in this case we're only checking the first reported alignment!
            flagcounts <- colSums( bamFlagAsBitMatrix( res[[1]]$flag[ bm ] ) )
            FC.ma <- FC.ma + flagcounts[ rownames( FC.ma ) ]
            ## process mapping quality:
            MAPQ.bam <- table( res[[1]]$mapq[ bm ] )
            MAPQ.ma[ as.character( names( MAPQ.bam ) ), 1 ] <- MAPQ.ma[ as.character( names( MAPQ.bam ) ), 1 ] + as.numeric( MAPQ.bam )
            counts.ma <- counts.ma + sum( bm )
        }
        ## do we have no-aligning reads???
        bm <- res[[1]]$tag$NH==0
        bm[ is.na( bm ) ] <- FALSE
        if( any( bm ) ){
            flagcounts <- colSums( bamFlagAsBitMatrix( res[[1]]$flag[ bm ] ) )
            FC.na <- FC.na + flagcounts[ rownames( FC.na ) ]
            ## process mapping quality:
            MAPQ.bam <- table( res[[1]]$mapq[ bm ] )
            MAPQ.na[ as.character( names( MAPQ.bam ) ), 1 ] <- MAPQ.na[ as.character( names( MAPQ.bam ) ), 1 ] + as.numeric( MAPQ.bam )
            counts.na <- counts.na + sum( bm )
        }
        counta <- counta + length( res[[1]]$flag )
        if( v )
            cat( paste0( "processed ", counta, " entries (", counts.na, ", ", counts.ua, ", ", counts.ma ,") number of (no-, unique, multi) aligned reads)\n" ) )
    }
    return( list( flag=list( no_algn=FC.na, unique_algn=FC.ua, multi_algn=FC.ma ),
                 mapq=list( no_algn=MAPQ.na, unique_algn=MAPQ.ua, multi_algn=MAPQ.ma ),
                 counts=list( no_algn=counts.na, unique_algn=counts.ua, multi_algn=counts.ma ))
           )
}

