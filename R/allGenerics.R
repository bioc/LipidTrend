#' @title Get Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing analysis results. The result table includes
#' the following columns:
#' \enumerate{
#'   \item{Feature columns: Lipid feature values from the input \code{rowData},
#'   such as chain length or double bond count. Column names vary depending on
#'   input.}
#'   \item{avg.abund: Mean abundance of each lipid across all samples. For
#'   one-dimensional analysis, this may also include \code{avg.abund.ctrl} and
#'   \code{avg.abund.case} for group-wise means.}
#'   \item{direction: Sign of the smoothed statistic:
#'     \itemize{
#'       \item{+ : Trend increases in the case group.}
#'       \item{- : Trend decreases in the case group.}
#'     }
#'   }
#'   \item{smoothing.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   region-based permutation test.}
#'   \item{marginal.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   marginal test (per lipid).}
#'   \item{log2.FC: Log2 fold-change in abundance between case and control
#'   groups.}
#'   \item{significance: Overall significance label based on smoothed test and
#'   FC direction:
#'     \itemize{
#'       \item{Increase: Significant positive trend in case group.}
#'       \item{Decrease: Significant negative trend in case group.}
#'       \item{NS: Not significant.}
#'     }
#'   }
#'}
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result
#' results <- result(res_se)
#' @export
setGeneric("result", function(object) standardGeneric("result"))

#' @rdname result
#' @exportMethod result
setMethod("result", "LipidTrendSE", function(object) {
    object@result
})

#' @title Get Even Chain Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing even chain result. The result table includes
#' the following columns:
#' \enumerate{
#'   \item{Feature columns: Lipid feature values from the input \code{rowData},
#'   such as chain length or double bond count. Column names vary depending on
#'   input.}
#'   \item{avg.abund: Mean abundance of each lipid across all samples. For
#'   one-dimensional analysis, this may also include \code{avg.abund.ctrl} and
#'   \code{avg.abund.case} for group-wise means.}
#'   \item{direction: Sign of the smoothed statistic:
#'     \itemize{
#'       \item{+ : Trend increases in the case group.}
#'       \item{- : Trend decreases in the case group.}
#'     }
#'   }
#'   \item{smoothing.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   region-based permutation test.}
#'   \item{marginal.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   marginal test (per lipid).}
#'   \item{log2.FC: Log2 fold-change in abundance between case and control
#'   groups.}
#'   \item{significance: Overall significance label based on smoothed test and
#'   FC direction:
#'     \itemize{
#'       \item{Increase: Significant positive trend in case group.}
#'       \item{Decrease: Significant negative trend in case group.}
#'       \item{NS: Not significant.}
#'     }
#'   }
#' }
#' @examples
#' data("lipid_se_CL")
#' sub <- lipid_se_CL[seq_len(10), ]
#' res_se <- analyzeLipidRegion(
#'     lipid_se=sub, ref_group="sgCtrl", split_chain=TRUE,
#'     chain_col="chain", radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result summary
#' results <- even_chain_result(res_se)
#' @export
setGeneric(
    "even_chain_result",
    function(object) standardGeneric("even_chain_result"))

#' @rdname even_chain_result
#' @aliases even_chain_result,LipidTrendSE-method
#' @exportMethod even_chain_result
setMethod("even_chain_result", "LipidTrendSE", function(object) {
    object@even_chain_result
})

#' @title Get Odd Chain Result from LipidTrendSE
#' @param object A LipidTrendSE object
#' @return A data frame containing odd chain result. The result table includes
#' the following columns:
#' \enumerate{
#'   \item{Feature columns: Lipid feature values from the input \code{rowData},
#'   such as chain length or double bond count. Column names vary depending on
#'   input.}
#'   \item{avg.abund: Mean abundance of each lipid across all samples. For
#'   one-dimensional analysis, this may also include \code{avg.abund.ctrl} and
#'   \code{avg.abund.case} for group-wise means.}
#'   \item{direction: Sign of the smoothed statistic:
#'     \itemize{
#'       \item{+ : Trend increases in the case group.}
#'       \item{- : Trend decreases in the case group.}
#'     }
#'   }
#'   \item{smoothing.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   region-based permutation test.}
#'   \item{marginal.pval.BH: Benjamini–Hochberg adjusted p-value from the
#'   marginal test (per lipid).}
#'   \item{log2.FC: Log2 fold-change in abundance between case and control
#'   groups.}
#'   \item{significance: Overall significance label based on smoothed test and
#'   FC direction:
#'     \itemize{
#'       \item{Increase: Significant positive trend in case group.}
#'       \item{Decrease: Significant negative trend in case group.}
#'       \item{NS: Not significant.}
#'     }
#'   }
#' }
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=TRUE,
#'     chain_col="chain", radius=3, own_contri=0.5, permute_time=100)
#' # Get complete result summary
#' results <- odd_chain_result(res_se)
#' @export
setGeneric(
    "odd_chain_result",
    function(object) standardGeneric("odd_chain_result"))

#' @rdname odd_chain_result
#' @aliases odd_chain_result,LipidTrendSE-method
#' @exportMethod odd_chain_result
setMethod("odd_chain_result", "LipidTrendSE", function(object) {
    object@odd_chain_result
})


setGeneric(
    ".split_chain",
    function(object) standardGeneric(".split_chain"))
setMethod(".split_chain", "LipidTrendSE", function(object) {
    object@split_chain
})


setGeneric(
    ".abund_weight",
    function(object) standardGeneric(".abund_weight"))
setMethod(".abund_weight", "LipidTrendSE", function(object) {
    object@abund_weight
})


#' @title Show method for LipidTrendSE objects
#' @param object A LipidTrendSE object
#' @return LipidTrendSE object information
#' @importMethodsFrom SummarizedExperiment show
#' @aliases show,LipidTrendSE-method
#' @exportMethod show
setMethod("show", "LipidTrendSE", function(object) {
    callNextMethod()
    cat("\nLipidTrend Results:\n")
    cat("------------------------\n")
    cat("Split chain analysis:",
        if(.split_chain(object)) "Yes" else "No", "\n")

    if (.split_chain(object)) {
        even_results <- even_chain_result(object)
        if (!is.null(even_results)) {
            cat(
                "Even chain result: ", nrow(even_results),
                " features\n", sep="")
        }
        odd_results <- odd_chain_result(object)
        if (!is.null(odd_results)) {
            cat(
                "Odd chain result: ",
                nrow(odd_results), " features\n", sep="")
        }
    } else {
        result <- result(object)
        if (!is.null(result)) {
            cat("Result: ", nrow(result), " features\n", sep="")
        }
    }
})
