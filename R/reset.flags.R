#' A utility script to reset to FALSE (or TRUE) the locus metric flags after
#' some individuals or populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,
#' PICSnp, and AvgPIC included, but the allelic composition will change when
#' some individuals are removed from the dataset and so the initial statistics
#' will no longer apply. This applies also to some variable calculated by dartR
#' (e.g. maf). This script resets the locus metrics flags to FALSE to indicate
#' that these statistics in the genlight object are no longer current. The
#' verbosity default is also set, and in the case of SilcoDArT, the flags PIC
#' and OneRatio are also set.
#'
#' If the locus metrics do not exist then they are added to the genlight object
#'  but not populated. If the locus metrics flags do not exist, then they are
#'  added to the genlight object and set to FALSE (or TRUE).
#'
#' @param x Name of the genlight object containing the SNP data or
#' tag presence/absence data (SilicoDArT) [required].
#' @param set Set the flags to TRUE or FALSE [default FALSE].
#' @noRd

reset.flags <- function(x,
                        set = FALSE) {

    # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics)) {
      x@other$loc.metrics <- as.data.frame(array(NA, nLoc(x)))
    }
    # Check if the x@other$loc.metrics.flags slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics.flags)) {
      x@other$loc.metrics.flags <- as.data.frame(array(NA, 1))
    }
    # loc.metric should be a dataframe
    x@other$loc.metrics <- as.data.frame(x@other$loc.metrics)
    # AvgPIC
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$AvgPIC <- set
    # OneRatioRef
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$OneRatioRef <- set
    # OneRatioSnp
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$OneRatioSnp <- set
    # PICRef
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$PICRef <- set
    # PICSnp
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$PICSnp <- set
    # CallRate
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$CallRate <- set
    # FreqHomRef
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$FreqHomRef <- set
    # FreqHomSnp
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$FreqHomSnp <- set
    # FreqHets
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$FreqHets <- set
    # monomorphs
    if (is.null(x@other$loc.metrics$monomorphs)) {
      x@other$loc.metrics$monomorphs <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$monomorphs <- set
    # maf
    if (is.null(x@other$loc.metrics$maf)) {
      x@other$loc.metrics$maf <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$maf <- set
    # OneRatio
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$OneRatio <- FALSE
    # PIC
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA, nLoc(x))
    }
    x@other$loc.metrics.flags$PIC <- FALSE
    # monomorphs
    x@other$loc.metrics.flags$monomorphs <- set
    # allna
    x@other$loc.metrics.flags$allna <- set

  return(x)

}
