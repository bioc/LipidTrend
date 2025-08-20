#' @title Conduct statistical to analyze lipid features tendencies
#' @description
#' Performs region-based statistical analysis to identify lipidomic trends
#' between groups. This function applies a two-stage procedure:
#' \enumerate{
#'      \item \bold{Marginal Test}:
#'  Each lipid feature is first tested individually using
#'  either a t-test (with glog10 transformation) or a Wilcoxon test to obtain
#'  marginal statistics and p-values.
#'      \item \bold{Region-Based Permutation Test with Smoothing}:
#'  Marginal statistics are smoothed using a Gaussian kernel that incorporates
#'  neighborhood information (e.g., proximity in chain length or double bond
#'  count). Statistical significance is assessed by comparing the smoothed
#'  statistic against a null distribution generated through permutation.
#' }
#'
#' This approach enhances statistical robustness, especially for small-sample
#' datasets. It supports both one-dimensional and two-dimensional analyses
#' depending on the number of features provided in the input.
#'
#' The input must be a \code{SummarizedExperiment} object, and the output is a
#' \code{LipidTrendSE} object, which can be used for result visualization or
#' further downstream analysis.
#' @param lipid_se A SummarizedExperiment object. Must contain following data:
#' \enumerate{
#'      \item{Assay: A matrix containing lipid abundance data, where rows
#'      represent lipids and columns represent samples.}
#'      \item{RowData: A data frame of lipid features (e.g., double bond count,
#'      chain length), with rows as lipids and columns as lipid features (
#'      limited to 1 or 2 columns). The order of lipids must match the abundance
#'      data. If RowData contains one column, a one-dimensional analysis will
#'      be performed. If \code{rowData} includes two columns, a two-dimensional
#'      analysis will be conducted.}
#'      \item{ColData: A data frame containing group information, where rows
#'      represent sample names and columns must include sample name, label name,
#'      and group, arranged accordingly.}
#' }
#' @param ref_group Character. Group name of the reference group. It must be
#' one of the group names in the \code{colData$group} column.
#' @param split_chain Logical. If \code{TRUE} the results will split to shown by
#'  odd and even chain. Default is \code{FALSE}.
#' @param chain_col Character. The column name in \code{rowData} that specifies
#' chain length. Must be provided if \code{split_chain=TRUE}, otherwise should
#' be set to \code{NULL}. Default is \code{NULL}.
#' @param radius Numeric. Distance of neighboring features to be included in
#' the smoothing kernel. Default is \code{3}.
#' @param own_contri Numeric. Proportion of self-contribution when smoothing.
#' Default is \code{0.5}. Recommended range: 0.5–1 to avoid over-emphasizing
#' neighbors.
#' @param test Character. Type of statistical test: either \code{"t.test"}
#' or \code{"Wilcoxon"}. Default is \code{"t.test"}.
#' @param abund_weight Logical. Whether to use average abundance as a weight
#' in calculating the region-based smoothed statistic. When set to TRUE, lipid
#' species with higher mean abundance contribute more to the smoothed trend.
#' Default is \code{TRUE}.
#' @param permute_time Integer. Number of permutations used to calculate
#' empirical p-values in the region-based permutation test. Default is
#' \code{100000}. For the Wilcoxon test (i.e., \code{test="Wilcoxon"}),
#' we recommend setting \code{permute_time} to fewer than 10,000 to ensure
#' reasonable runtime.
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom stats pt p.adjust
#' @importFrom dplyr mutate
#' @importFrom matrixTests col_wilcoxon_twosample
#' @importFrom MKmisc glog10
#' @return A LipidTrendSE object containing lipidomic feature testing result.
#' @examples
#' data("lipid_se_CL")
#' res_se <- analyzeLipidRegion(
#'     lipid_se=lipid_se_CL, ref_group="sgCtrl", split_chain=FALSE,
#'     chain_col=NULL, radius=3, own_contri=0.5, test="t.test",
#'     abund_weight=TRUE, permute_time=100)
#' @export
#' @seealso
#' \code{\link{plotRegion1D}} for one-dimensional visualization
#' \code{\link{plotRegion2D}} for two-dimensional visualization
analyzeLipidRegion <- function(
        lipid_se, ref_group, split_chain=FALSE, chain_col=NULL, radius=3,
        own_contri=0.5, test="t.test", abund_weight=TRUE, permute_time=100000){
    if (isFALSE(inherits(lipid_se, "SummarizedExperiment")) ) {
        stop("Invalid lipid_se input.")
    }
    if (!is.numeric(radius) || radius <= 0) {
        stop("radius must be a positive numeric value.")
    }
    if (!is.numeric(own_contri) || own_contri < 0 || own_contri > 1) {
        stop("own_contri must be between 0 and 1.")
    }
    X <- assay(lipid_se) %>% t()
    if (ncol(X) < 2) {
        stop(
            "Lipid characteristic (e.g., chain length, double bond) must ",
            "contain at least two features.")
    }
    X.info <- rowData(lipid_se) %>% as.data.frame()
    group_info <- colData(lipid_se) %>% as.data.frame()
    if (is.null(ref_group) || !(ref_group %in% group_info$group)) {
        stop("ref_group must be one of the groups in colData group column.")
    }
    group_info %<>% mutate(
        lipidTrend_group=ifelse(group_info$group==ref_group, 0, 1))
    group <- group_info$lipidTrend_group
    if (isTRUE(split_chain)) {
        if (!is.null(chain_col) && chain_col %in% colnames(X.info)) {
            output.df <- .analyzeSplitChain(
                X, X.info, group, chain_col, radius, own_contri, test,
                permute_time, abund_weight)
            return(new(
                "LipidTrendSE", lipid_se, split_chain=TRUE,
                even_chain_result=output.df$even,
                odd_chain_result=output.df$odd, abund_weight=abund_weight))
        } else {
            stop("chain_col is not present in the column name of rowData.")
        }
    } else if (isFALSE(split_chain)) {
        output.df <- .stat_lipidTrend(
            X, X.info, group, radius, own_contri, test, permute_time,
            abund_weight)
        return(new(
            "LipidTrendSE", lipid_se, split_chain=FALSE, result=output.df,
            abund_weight=abund_weight))
    } else {
        stop("split_chain must be logical value.")
    }
}

.analyzeSplitChain <- function(
        X, X.info, group, chain_col, radius, own_contri, test,
        permute_time, abund_weight) {
    chain_parity <- X.info[[chain_col]] %% 2
    chain_groups <- list(even=chain_parity == 0, odd=chain_parity == 1)
    output.df <- lapply(names(chain_groups), function(chain_type) {
        sub_feature <- chain_groups[[chain_type]]
        if (sum(sub_feature) < 2) {
            warning(
                sprintf(
                    "Fail to proceed with the %s-chain analysis. ", chain_type),
                sprintf("Detect less than two %s-chain lipids.", chain_type)
            )
            return (NULL)
        }
        .stat_lipidTrend(
            X=X[, sub_feature], X.info=subset(X.info, sub_feature),
            group, radius, own_contri, test, permute_time, abund_weight)
    })
    names(output.df) <- names(chain_groups)
    return(output.df)
}

.regionStat <- function(X, Y, test="t.test", ...){
    n1 <- sum(Y[,1])
    n2 <- sum(1-Y[,1])
    G1 <- Y / n1
    G2 <- (1-Y) / n2
    mean.group1 <- t(X) %*% G1
    mean.group2 <- t(X) %*% G2
    if (test == "t.test") {
        var.group1 <- (t(X^2) %*% Y  - n1*mean.group1^2)/(n1 - 1)
        var.group2 <- (t(X^2) %*% (1-Y)  - n2*mean.group2^2)/(n2 - 1)
        pooled.var <- ((n1-1)*var.group1 + (n2-1)*var.group2)/(n1+n2-2)
        t.stat <- (mean.group1-mean.group2) / sqrt((1/n1+1/n2)*pooled.var)
        df <- n1 + n2 - 2
        t.test.pval <- 2*pt(abs(t.stat), df = df, lower.tail = FALSE)
        res <- -log10(t.test.pval) * sign(t.stat)
    }
    if (test == "Wilcoxon") {
        wilcox.pval <- apply(Y, 2, function(y) {
            group1 <- X[y == 1, , drop=FALSE]
            group2 <- X[y == 0, , drop=FALSE]
            col_wilcoxon_twosample(
                x=group1, y=group2, alternative="two.sided")$pvalue
        })
        res <- -log10(wilcox.pval) * sign(mean.group1-mean.group2)
    }
    return(res)
}

.stat_lipidTrend <- function(
        X, X.info, group, radius, own_contri, test, permute_time, abund_weight){
    if (test == 't.test') {
        glog.X <- glog10(X)
        region.stat.obs <- .regionStat(X=glog.X, Y=as.matrix(group), test=test)
    } else {
        region.stat.obs <- .regionStat(X=X, Y=as.matrix(group), test=test)
    }
    dis_res <- .countDistance(X.info)
    dist_input <- dis_res$normalized_matrix
    x.distance <- dis_res$x_distance
    y.distance <- dis_res$y_distance
    dim <- dis_res$dimensions
    dist.mat <- as.matrix(dist(dist_input, method = "euclidean"))
    dist.mat[dist.mat>=radius] <- Inf
    stat_res <- .smooth_permutation(
        X, group, dist.mat, own_contri, test, abund_weight, permute_time,
        region.stat.obs)
    smooth.stat <- stat_res$smooth.stat
    smooth.stat.permute <- stat_res$smooth.stat.permute
    output.df <- .resultSummary(
        X, X.info, group, smooth.stat, smooth.stat.permute, region.stat.obs,
        dim, permute_time)
    if(all(is.na(output.df$smoothing.pval.BH))){
        stop(
            "Gaussian kernel smoothing may fail due to some lipids
            having zero ", "variance in abundance between the two groups or
            having identical ", "sample values. Please remove those lipids and
            re-run the function.")
    }
    return(output.df)
}

.smooth_permutation <- function(
        X, group, dist.mat, own_contri, test, abund_weight, permute_time,
        region.stat.obs) {
    dist.idx <- which.max(colSums(exp(-dist.mat^2)))
    radius.seq <- seq(0.01, 5, 0.01)
    own_contri.seq <- vapply(
        radius.seq, function(x) 1/sum(exp(-dist.mat[,dist.idx]^2*x)),
        FUN.VALUE=numeric(1))
    radius.idx <- sum(own_contri.seq <= own_contri)
    if (radius.idx + 1 > length(radius.seq)) {
        stop(sprintf(
            "Based on the data structure and 'radius' setting, \n
        the 'own_contri' should smaller than %0.4f", max(own_contri.seq)))
    }
    kernel.radius <- radius.seq[radius.idx + 1]
    exp.dist <- exp(-dist.mat^2*kernel.radius)
    rescaled.X <- X-min(X)
    rescaled.X <- 4*rescaled.X/max(rescaled.X) + 1
    avg.log.abund <- colMeans(rescaled.X)
    if (abund_weight) {
        smooth.stat <- exp.dist %*% (region.stat.obs * avg.log.abund)
    } else {
        smooth.stat <- exp.dist %*% region.stat.obs
    }
    Y.permute <- vapply(
        seq_len(permute_time), function(x) sample(group, replace=FALSE),
        FUN.VALUE=numeric(length(group)))
    if (test == 't.test') {
        glog.X <- glog10(X)
        region.stat.permute <- .regionStat(X=glog.X, Y=Y.permute, test=test)
    } else {
        region.stat.permute <- .regionStat(X=X, Y=Y.permute, test=test)
    }
    if (abund_weight) {
        smooth.stat.permute <- exp.dist %*% (region.stat.permute*avg.log.abund)
    } else {
        smooth.stat.permute <- exp.dist %*% region.stat.permute
    }
    return(list(
        smooth.stat=smooth.stat, smooth.stat.permute=smooth.stat.permute))
}

.resultSummary <- function(
        X, X.info, group, smooth.stat, smooth.stat.permute,
        region.stat.obs, dimension, permute_time, ...){
    smooth.stat.all <- cbind(smooth.stat, smooth.stat.permute)
    smooth.stat.pval <- apply(
        smooth.stat.all, 1, function(x)
            mean(abs(round(x[1], digits=10)) < abs(round(x[-1], digits=10))))
    smooth.stat.pval[smooth.stat.pval == 0] <- 1/permute_time
    smooth.stat.pval.BH <- p.adjust(smooth.stat.pval, method="BH")
    direction <- ifelse(sign(smooth.stat) > 0, "+", "-")
    marginal.pval <- 10^-abs(region.stat.obs)
    marginal.pval.BH <- p.adjust(marginal.pval, method="BH")
    log2.FC <- apply(
        X, 2, function(x) log2(mean(x[group == 1])/mean(x[group == 0])))
    FC.direction <- ifelse(log2.FC > 0, "+", "-")
    smooth.stat.pval.BH[FC.direction != direction] <- 1
    result <- data.frame(
        X.info, direction=direction, smoothing.pval.BH=smooth.stat.pval.BH,
        marginal.pval.BH=marginal.pval.BH, log2.FC=log2.FC) %>%
        mutate(significance=ifelse(
            smoothing.pval.BH < 0.05 & log2.FC > 0, 'Increase',
            ifelse(smoothing.pval.BH < 0.05 & log2.FC < 0, 'Decrease', 'NS')))
    if (dimension == 2) {
        result %<>%
            mutate(avg.abund=colMeans(X), .before="direction")
    }
    if (dimension == 1) {
        result %<>%
            mutate(
                avg.abund.ctrl=colMeans(X[group == 0,]),
                avg.abund.case=colMeans(X[group == 1,]),
                .before="direction")
    }
    return(result)
}


