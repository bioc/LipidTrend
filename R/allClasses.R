setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))

#' @title Class for storing analyzeLipidRegion analysis results
#' @description This class extends SummarizedExperiment to store
#' analyzeLipidRegion analysis results
#' @slot split_chain Logical. Indicates whether chains were split by
#' parity (even vs. odd).
#'   When TRUE, results are stored in even_chain_result and
#'   odd_chain_result slots.
#'   When FALSE, results are stored in the result slot.
#' @slot result Data frame of analysis results (for non-split data)
#' @slot even_chain_result Data frame of results for even chains
#' (when split_chain is TRUE)
#' @slot odd_chain_result Data frame of results for odd chains
#' (when split_chain is TRUE)
#' @slot abund_weight Logical. Logical. Set to TRUE to incorporate average
#' lipid abundance as a weight in the test statistic.
#' @importFrom methods setClass setClassUnion setValidity
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @exportClass LipidTrendSE
setClass(
    "LipidTrendSE", contains="SummarizedExperiment",
    slots=c(
        split_chain="logical", result="data.frameOrNULL",
        even_chain_result="data.frameOrNULL",
        odd_chain_result="data.frameOrNULL", abund_weight="logical"),
    prototype=list(
        split_chain=FALSE, result=NULL,
        even_chain_result=NULL, odd_chain_result=NULL, abund_weight=TRUE)
)

#' @title Validate LipidTrendSE object
#' @name LipidTrendSE-validity
#' @param object A LipidTrendSE object to validate
#' @return TRUE if valid, otherwise an error message
setValidity("LipidTrendSE", function(object) {
    if (isTRUE(object@split_chain)) {
        if (is.null(object@even_chain_result) && is.null(
            object@odd_chain_result)) {
            return("Split chain is TRUE but no chain results provided")
        }
    } else {
        if (is.null(object@result)) {
            return("Split chain is FALSE but no result provided")
        }
    }
    return(TRUE)
})
