#' a bsseq loader for Biscuit output (BED-like format, 2 or 3 cols/sample)
#' 
#' e.g. P01-028-T06.joint.ch.hg19.bed.gz has 3 samples, and 9 columns total,
#' while P01-028-T06.joint.cg.merged.hg19.bed.gz has 3 samples and 12 columns.
#' Note: the defaults assume alignment against hg19 (use genome=xyz to override)
#' Note 2: if a BED has no header, a VCF header can be used to autodetect names.
#'
#' @param BEDfile     the file (compressed or not, doesn't matter) to load
#' @param VCFfile     the file (compressed and tabixed, with header) to load
#' @param merged      are CpG sites merged?
#' @param sampleNames if NULL, create; if VCF, read; if data.frame, make pData
#' @param simplify    simplify sample names by dropping .foo.bar.hg19 or similar
#' @param genome      what genome assembly were the runs aligned against? (hg19)
#' @param how         how to load the data? "data.table" (default) or "readr"
#' @param hdf5        make the object HDF5-backed? (FALSE; use in-core storage) 
#' @param hdf5dir     if hdf5 is TRUE, where should HDF5 files be stored? (NULL)
#' @param replace     replace hdf5dir if is already exists (FALSE)
#' @param chunkdim    chunk dimensions for writing HDF5 file to disk (NULL)
#' @param level       compression level for writing HDF5 to disk (NULL)
#' @param sparse      are there a lot of zero-coverage sites? (default is FALSE)
#' @param clumpSize   number of rows before readr reading becomes clumped (1e6)
#' @param chr         load a specific chromosome (to rbind() later)? (NULL)
#' @param which       a GRanges of regions to load (default NULL, load them all)
#' @param verbose     be verbose? (FALSE) 
#' 
#' @return            a bsseq::BSseq object, possibly Matrix- or HDF5-backed
#'
#' @import data.table
#' @import readr
#' @import bsseq
#' @import rlang
#'
#' @seealso BSseq
#' @aliases load.biscuit
#' @seealso checkBiscuitBED
#'
#' @export
read.biscuit <- function(BEDfile, 
                         VCFfile, 
                         merged, 
                         sampleNames=NULL, 
                         simplify=FALSE, 
                         genome="hg19",
                         how=c("data.table","readr"),
                         hdf5=FALSE, 
                         hdf5dir=NULL,
                         replace=FALSE,
                         chunkdim=NULL,
                         level=NULL,
                         sparse=FALSE,
                         clumpSize=1e6, 
                         chr=NULL,
                         which=NULL,
                         verbose=FALSE) { 

  # Check if required inputs are missing
  # Print more useful messages if they are
  if (is_missing(BEDfile)) stop("Tabix'ed BED file from biscuit is required.\n")
  if (is_missing(VCFfile)) {
    err_message <- paste("Tabix'ed VCF file from biscuit is required.",
                         "Header information is used to set up column names.\n")
    stop(err_message)
  }
  if (is_missing(merged)) {
    err_message <- paste("merged flag is required.",
                         "merged = TRUE if 'biscuit mergecg' was run after 'biscuit vcf2bed'.",
                         "Otherwise use merged = FALSE.\n")
    stop(err_message)
  }

  how <- match.arg(how)
  params <- checkBiscuitBED(BEDfile=BEDfile, VCFfile=VCFfile, how=how, chr=chr,
                            sampleNames=sampleNames, clumpSize=clumpSize, hdf5=hdf5,
                            hdf5dir=hdf5dir, replace=replace,
                            sparse=sparse, merged=merged)
  message("Reading ", ifelse(params$merged, "merged", "unmerged"), 
          " input from ", params$tbx$path, "...")

  if (params$how == "data.table") {
    # {{{
    select <- grep("\\.context", params$colNames, invert=TRUE)
    if (is.null(which)) {
      cmd <- paste("gunzip -c", params$tbx$path) # for mac compatibility
    } else { 
      tmpBed <- tempfile(fileext=".bed")
      export(which, tmpBed)
      cmd <- paste("tabix -R", tmpBed, params$tbx$path)
    }
    tbl <- fread(cmd=cmd, sep="\t", sep2=",", fill=TRUE, na.string=".", 
                 select=select)
    if (params$hasHeader == FALSE) names(tbl) <- params$colNames[select]
    names(tbl) <- sub("^#", "", names(tbl))
    # }}}
  } else if (params$how == "readr") {
    # {{{
    if (params$passes > 1) { 
      f <- function(x, pos) {
        message("Reading line ", pos, "...")
        return(x)
      }
      message("Making ",params$passes," passes of ",clumpSize," loci each...")
      tbl <- with(params,
                  read_tsv_chunked(tbx$path, DataFrameCallback$new(f), na=".",
                                   skip=as.numeric(params$hasHeader), 
                                   col_names=colNames, col_types=colSpec, 
                                   chunk_size=clumpSize))
    } else { 
      message("If the following is slow, you may need to decrease clumpSize")
      message("from ",clumpSize," to something smaller & do multiple passes.")
      tbl <- with(params,
                  read_tsv(tbx$path, na=".", comment="#",
                           skip=as.numeric(params$hasHeader), 
                           col_names=colNames, col_types=colSpec))
    }
    # }}}
  }

  # shift from 0-based to 1-based coordinates  
  tbl[, 2] <- tbl[, 2] + 1 # FIXME: can this be done automagically?

  # Remove CpG sites with zero-coverage
  if(!params$sparse) {
    message("Excluding CpG sites with uniformly zero coverage...")
    tbl <- tbl[rowSums(is.na(tbl)) == 0, ]
  }

  # Make HDF5-backed BSseq object
  message("Loaded ", params$tbx$path, ". Creating bsseq object...")
  if (hdf5) {
    # Set up RegularArrayGrid for reading in data
    # Variable Description:
    #   - rag_nrow: number of rows in data
    #   - rag_ncol: number of samples in the data
    #   - rag_dim : matrix dimensions for the data (nrow x ncol)
    #   - grid    : refdim - matrix dimensions
    #               spacings - how to read fill/read the data (currently reads single-column-wise)
    # TODO: Investigate best number of columns to load (change 1L to nSamples????)
    # TODO: Investigate other grid shapes (ArbitraryArrayGrid, BlockGrid)
    rag_nrow <- as.integer(nrow(tbl))
    rag_ncol <- as.integer(params$nSamples)
    rag_dim <- c(rag_nrow, rag_ncol)
    grid <- RegularArrayGrid(refdim = rag_dim, spacings = c(rag_nrow, 1L))

    # Initialize HDF5RealizationSink
    hdf5_path <- file.path(hdf5dir, "assays.h5")
    M_sink <- HDF5RealizationSink(dim = rag_dim,
                                  dimnames = NULL,
                                  type = "integer",
                                  filepath = hdf5_path,
                                  name = "M",
                                  chunkdim = chunkdim,
                                  level = level)
    on.exit(close(M_sink), add = TRUE)
    Cov_sink <- HDF5RealizationSink(dim = rag_dim,
                                    dimnames = NULL,
                                    type = "integer",
                                    filepath = hdf5_path,
                                    name = "Cov",
                                    chunkdim = chunkdim,
                                    level = level)
    on.exit(close(Cov_sink), add = TRUE)
    sink_lock <- ipcid()
    on.exit(ipcremove(sink_lock), add = TRUE)
    res <- makeBSseq_hdf5(tbl, params, grid = grid,
                          M_sink = M_sink, Cov_sink = Cov_sink,
                          sink_lock = sink_lock, simplify = simplify,
                          verbose = verbose)
  } else {
    # Write into memory rather than to disk
    res <- makeBSseq(tbl, params, simplify=simplify, verbose=verbose)
  }
  metadata(res)$vcfHeader <- params$vcfHeader
  genome(rowRanges(res)) <- genome

  return(res)

}


#' @export
load.biscuit <- read.biscuit
