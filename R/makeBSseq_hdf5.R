#' make a HDF5-backed BSseq object from a biscuit BED
#'
#' @param tbl       a tibble (from read_tsv) or a data.table (from fread)
#' @param params    parameters from checkBiscuitBED
#' @param grid      RegularArrayGrid for writing HDF5 file
#' @param M_sink    HDF5RealizationSink for M matrix
#' @param Cov_sink  HDF5RealizationSink for Cov matrix
#' @param sink_lock IPC lock to avoid race condition when writing to HDF5 file
#' @param simplify  simplify sample names by dropping .foo.bar.hg19 and similar (FALSE)
#' @param verbose   be verbose about what is happening (FALSE)
#'
#' @return          a HDF5-backed BSseq object
#'
#' @import GenomicRanges
#' @import bsseq
#'
#' @seealso makeBSseq
#'
#' @export
#'
makeBSseq_hdf5 <- function(tbl,
                           params,
                           grid,
                           M_sink,
                           Cov_sink,
                           sink_lock,
                           simplify = FALSE,
                           verbose = FALSE) {

    # helper function for coercing vector to matrix
    # x      : vector to coerce
    # verbose: be verbose about what is happening
    # returns: a matrix
    matMe <- function(x, verbose) {
        if (!is(x, "matrix")) {
            if (verbose) message("Turning a vector in a matrix...")
            x <- as.matrix(x)
        }
        return(x)
    }

    # helper function to fix row and column names
    # x      : data.frame-like object that needs fixing
    # gr     : GenomicRanges object for rownames
    # what   : what matrix needs the names fixed?
    # verbose: be verbose about what is happening
    # returns: data.frame-like object with fixed names
    fixNames <- function(x, gr, what=c("M","Cov"), verbose) {
        if (is.null(rownames(x))) {
            if (verbose) message("Adding rownames...")
            rownames(x) <- as.character(gr)
        }
        colnames(x) <- base::sub("beta", match.arg(what), colnames(x))
        return(x)
    }

    gr <- resize(makeGRangesFromDataFrame(tbl[, 1:3]), 1)

    # Create M_tmp and Cov_tmp matrices from tbl
    if (params$how == "data.table") {
        # deal with data.table weirdness
        betas <- match(params$betaCols, names(tbl))
        covgs <- match(params$covgCols, names(tbl))
        if (verbose) message("Generating M matrix")
        M_tmp <- matMe(x=fixNAs(round(tbl[,..betas]*tbl[,..covgs]), y=0, params$sparse),
                   verbose=verbose)
        if (verbose) message("Generating Cov matrix")
        Cov_tmp <- matMe(x=fixNAs(tbl[, ..covgs], y=0, params$sparse),
                     verbose=verbose)
    } else { 
        M_tmp <- with(params, 
                  matMe(x=fixNAs(round(tbl[,betaCols]*tbl[,covgCols]), y=0, sparse),
                        verbose=verbose))
        Cov_tmp <- with(params, 
                    matMe(x=fixNAs(tbl[, covgCols], y=0, sparse), 
                          verbose=verbose))
    }

    # Write HDF5 to disk
    for (i in 1:params$nSamples) {
        if (verbose) message("Writing sample, ", i, " now...")
        viewport <- grid[[i]] #TODO: figure out how to properly loop over grid
        ipclock(sink_lock)
        write_block(x = M_sink, viewport = viewport, block = M_tmp)
        write_block(x = Cov_sink, viewport = viewport, block = Cov_tmp)
        ipcunlock(sink_lock)
    }

    # Create M and Cov matrices from M_sink and Cov_sink
    if (verbose) message("Make M a DelayedArray")
    M <- as(M_sink, "DelayedArray")
    if (verbose) message("Make Cov a DelayedArray")
    Cov <- as(Cov_sink, "DelayedArray")

    if (verbose) message("Make the SummarizedExperiment now...")
    sum.expt <- SummarizedExperiment(assays = list(M = M, Cov = Cov),
                                     rowRanges = gr)
                                     #colData = params$colData)
    #if (verbose) message("Fixing names now...")
    #M_tmp <- fixNames(M_tmp, gr, what="M", verbose=verbose)
    #Cov_tmp <- fixNames(Cov_tmp, gr, what="Cov", verbose=verbose)
    #colnames(Cov_tmp) <- colnames(M_tmp) <- params$pData$sampleNames

    if (verbose) message("Create the BSseq object!")
    bsseq <- new2("BSseq", sum.expt, check = FALSE)

    x <- bsseq
    x@assays <- HDF5Array:::.shorten_assay2h5_links(x@assays)
    saveRDS(x, file = params$rds_path)
    bsseq
}
