## function that uses system.call to sample from a BAM file using the samtools
## samtools view -sb
## paired end reads: according to https://groups.google.com/forum/#!topic/bedtools-discuss/gf0KeAJN2Cw
## Samtools's view -s option selects reads according to scores based on (a hash calculated from) their read names.  Hence the in or out decision comes out the same way for all reads with the same name, so read pairs are indeed kept intact.
## maintains read pairs intact.
## option -s: -s FLOAT  Integer part is used to seed the random number generator [0]. Part after the decimal point sets the fraction of templates/pairs to subsample [no  subsampling].
if( !isGeneric( "sampleBam", ) )
    setGeneric( "sampleBam", function( x, ... )
               standardGeneric( "sampleBam" ))
## x: file name (including path?) of the BAM file
## seed: random seed to be used; to allow reproducible sampling.
## fraction: fraction of the BAM file to sample
## ... additional parameters to be submitted to the scanBam function
setMethod( "sampleBam", "character", function( x, seed=sample( 1:1e+6, 1 ), fraction=0.001, outfile=tempfile(), ... ){
    ResFile <- .doSampleBam( x, seed=seed, fraction=fraction, outfile=outfile )
    ## read the file....
    SB <- scanBam( ResFile, ... )
    ## remove all those NA entries:
    Names <- names( SB[[1]] )
    idx <- 1:length( Names )
    if( any( Names=="tag" ) ){
        idx <- idx[ -which( Names=="tag" ) ]
    }
    IsNa <- is.na( SB[[1]][[ idx[ 1 ] ]] )
    if( any( IsNa ) ){
        for( i in idx ){
            SB[[1]][[ i ]] <- SB[[1]][[ i ]][ !IsNa ]
        }
        ## also remove entries from tag... if tag is present:
        if( any( Names=="tag" ) ){
            for( i in 1:length( SB[[1]]$tag ) ){
                SB[[1]]$tag[[ i ]] <- SB[[1]]$tag[[ i ]][ !IsNa ]
            }
        }
    }
    ## remove the file
    file.remove( ResFile )
    return( SB )
} )

setMethod( "sampleBam", "BamFile", function( x, seed=sample( 1:1e+6, 1 ), fraction=0.001, outfile=tempfile(), ... ){
    return( sampleBam( path( x ), seed=seed, fraction=fraction, outfile=outfile, ... ) )
} )

.doSampleBam <- function( bamfile, fraction=0.001, seed=sample( 1:1e+6, 1 ), outfile=tempfile() ){
    if( fraction > 0.5 )
        warning( "depending on its version, samtools might have a problem when fraction is larger 0.5!" )
    if( fraction > 1 )
        stop( "fraction has to be smaller than 1!" )
    fract <- unlist( strsplit( as.character( fraction ), split=".", fixed=TRUE ) )[ 2 ]
    samcall <- paste0( "samtools view -s ", seed, ".", fract, " -b ", bamfile, " > ", outfile )
    Res <- system( samcall )
    if( Res!=0 )
        stop( "Error running samtools!" )
    return( outfile )
}

