#' Function to load HDF5-backed SummarizedExperiment from read.biscuit
#'
#' Local wrapper for HDF5Array::loadHDF5SummarizedExperiment
#' 
#' @param   dir     path to the HDF5-backed file
#' @param   prefix  optional prefix to add to filenames inside ‘dir’; default is ""
#' 
#' @return        a SummarizedExperiment-derived object (such as a BSseq object)
#' 
#' @import HDF5Array
#' 
#' @export
#'
loadHDF5 <- function(dir, prefix = "") {
    # Check directory exists before running load function
    if (!dir.exists(dir)) {
        stop(dir, " does not exist. Enter a valid directory")
    }
    loadHDF5SummarizedExperiment(dir = dir, prefix = prefix)
}
