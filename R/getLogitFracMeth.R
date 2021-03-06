#' Helper function for compartment inference
#'
#' Want an object with nominally Gaussian error for compartment inference, so
#' this function uses 'suitable' (defaults to to 3 or more reads in 2 or more 
#' samples) measurements. Using Dirichlet smoothing (adding 'k' reads to M
#' and U), these measurements are then turned into lightly moderated,
#' logit-transformed methylated-fraction estimates for compartment calling.
#'
#' @param x        A bsseq object with methylated and total reads
#' @param minCov   Minimum read coverage for landmarking samples (DEFAULT: 3)
#' @param minSamp  Minimum landmark samples with at least minCov (DEFAULT: 2)
#' @param k        Pseudoreads for smoothing (DEFAULT: 0.1)
#' @param r        Regions to collapse over - if NULL, do it by CpG
#'                   (DEFAULT: NULL)
#'
#' @return         Smoothed logit(M / Cov) GRanges with coordinates as row names
#'
#' @import gtools
#' @import bsseq
#' @importFrom methods is
#'
#' @aliases getMvals
#'
#' @examples
#'
#'   orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
#'                           package="biscuiteer")
#'   orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
#'                           package="biscuiteer")
#'   bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
#'                       merged = FALSE)
#'
#'   reg <- GRanges(seqnames = rep("chr11",5),
#'                  strand = rep("*",5),
#'                  ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
#'                                   end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
#'                  )
#'
#'   frac <- getLogitFracMeth(bisc, minSamp = 1, r = reg)
#'
#' @export
#'
getLogitFracMeth <- function(x,
                             minCov = 3,
                             minSamp = 2,
                             k = 0.1,
                             r = NULL) {

  # do any loci/regions have enough read coverage in enough samples? 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    covgs <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal")
  } else { 
    covgs <- getCoverage(x, type="Cov", what="perBase")
  } 

  usable <- DelayedMatrixStats::rowSums2(covgs >= minCov) >= minSamp
  if (!any(usable)) stop("No usable loci/regions ( >= minCov in >= minSamp )!")
    
  # construct a subset of the overall BSseq object with smoothed mvalues 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    out <- subset(sort(r), usable)
    smoothed <- getSmoothedLogitFrac(x, k=k, minCov=minCov, r=out)
  } else { 
    out <- subset(x, usable) # Leave as BSseq object
    smoothed <- getSmoothedLogitFrac(out, k=k, minCov=minCov)
    out <- granges(out) # Pull out GRanges portion for output
  } 

  mcols(out) <- as.data.frame(smoothed)
  return(out)

}

# Helper function to find smoothed logit fraction
# x: bsseq object
# k: pseudoreads for smoothing (DEFAULT: 0.1)
# minCov: minimum coverage (DEFAULT: 3)
# maxFrac: maximum fraction of NAs allowed (DEFAULT 0.5)
# r: regions to collapse over
getSmoothedLogitFrac <- function(x,
                                 k = 0.1,
                                 minCov = 3,
                                 maxFrac = 0.5,
                                 r = NULL) {

  if (!is.null(r) && is(r, "GenomicRanges")) {
    M <- getCoverage(x, sort(r), type="M", what="perRegionTotal")
    U <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal") - M 
    rnames <- as.character(sort(r))
  } else { 
    M <- getCoverage(x, type="M", what="perBase")
    U <- getCoverage(x, type="Cov", what="perBase") - M 
    rnames <- as.character(granges(x))
  } 

  res <- logit((M + k) / ((M + k) + (U + k))) 
  rownames(res) <- rnames 

  makeNA <- ((M + U) < minCov)
  maxPct <- paste0(100 * maxFrac, "%")
  tooManyNAs <- (DelayedMatrixStats::colSums2(makeNA)/nrow(x)) > maxFrac
  if (any(tooManyNAs)) {
    message(paste(colnames(x)[tooManyNAs],collapse=", ")," are >",maxPct," NA!")
  }
  res[ makeNA ] <- NA
  return(res)

}

#' @describeIn getLogitFracMeth Alias for getLogitFracMeth
#'
getMvals <- getLogitFracMeth
