library(testthat)
library(SummarizedExperiment)

if (!exists(".countDistance", envir=asNamespace("LipidTrend"))) {
    .countDistance <- function(X.info) {
        if ("y" %in% colnames(X.info)) {
            normalized_matrix <- as.matrix(X.info[, c("x", "y")])
            list(
                normalized_matrix=normalized_matrix,
                x_distance=min(diff(sort(unique(X.info$x)))),
                y_distance=min(diff(sort(unique(X.info$y)))), dimensions=2)
        } else {
            # 1D case
            normalized_matrix <- as.matrix(X.info[, "x", drop=FALSE])
            list(
                normalized_matrix=normalized_matrix,
                x_distance=min(diff(sort(unique(X.info$x)))),
                y_distance=NULL, dimensions=1)
        }
    }
    assignInNamespace(".countDistance", .countDistance, "LipidTrend")
}

# Helper function to create mock SummarizedExperiment object
create_mock_se <- function(n_samples=6, n_features=10, split_chain=FALSE) {
    if (split_chain) {
        feature_names <- as.character(seq(33, 33 + n_features - 1, by = 1))
        assay_data <- matrix(
            rnorm(n_samples * n_features, mean=5, sd=1),
            nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
        row_data <- data.frame(
            chain=as.numeric(feature_names), x=as.numeric(feature_names),
            row.names=feature_names)
    } else {
        base_chains <- seq(33, 33 + n_features - 1)
        double_bonds <- rep(1:2, length.out = n_features)
        feature_names <- paste("TG", paste0(base_chains, ":", double_bonds))
        assay_data <- matrix(
            rnorm(n_samples * n_features, mean=5, sd=1),
            nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
        row_data <- data.frame(
            x=base_chains, y=double_bonds, Total.C=base_chains,
            Total.DB=double_bonds, row.names=feature_names)
    }
    col_data <- data.frame(
        sample_name=paste0("sample", seq_len(n_samples)),
        label_name=paste0("label", seq_len(n_samples)),
        group=rep(c("control", "case"), length.out=n_samples),
        row.names=paste0("sample", seq_len(n_samples)))
    SummarizedExperiment(
        assays=list(abundance=assay_data), rowData=row_data, colData=col_data)
}

set.seed(1234)

test_that("analyzeLipidRegion handles basic input correctly", {
    se <- create_mock_se()
    expect_no_error(
        result <- analyzeLipidRegion(
            lipid_se=se, ref_group="control", split_chain=FALSE, radius=3,
            permute_time=100)
    )
    expect_true(inherits(result, "LipidTrendSE"))
    expect_false(.split_chain(result))
    res_df <- result(result)
    expect_true(is.data.frame(res_df))
    expect_true(
        all(c(
            "direction", "smoothing.pval.BH", "log2.FC") %in% colnames(res_df)))
})

test_that("analyzeLipidRegion handles split chain analysis", {
    se <- create_mock_se(split_chain=TRUE)
    expect_no_error(
        result <- analyzeLipidRegion(
            lipid_se=se, ref_group="control", split_chain=TRUE,
            chain_col="chain", radius=3,
            permute_time=100)
    )
    expect_true(.split_chain(result))
    expect_true(!is.null(
        even_chain_result(result)) || !is.null(odd_chain_result(result)))
})

test_that("analyzeLipidRegion handles empty chain groups correctly", {
    n_samples <- 6
    n_features <- 5
    feature_names <- as.character(seq(33, 33 + 2*(n_features-1), by = 2))
    assay_data <- matrix(
        rnorm(n_samples * n_features, mean=5, sd=1),
        nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
    row_data <- data.frame(
        chain=as.numeric(feature_names), row.names=feature_names)
    col_data <- data.frame(
        sample_name=paste0("sample", seq_len(n_samples)),
        label_name=paste0("label", seq_len(n_samples)),
        group=rep(c("control", "case"), length.out = n_samples),
        row.names=paste0("sample", seq_len(n_samples)))
    se_odd_only <- SummarizedExperiment(
        assays=list(abundance=assay_data), rowData=row_data, colData=col_data)
    expect_warning(
        analyzeLipidRegion(
            lipid_se=se_odd_only, ref_group="control", split_chain=TRUE,
            chain_col="chain", radius=3, permute_time=100)
    )
    # Create the opposite case - only even chains
    feature_names_even <- as.character(seq(34, 34 + 2*(n_features-1), by = 2))
    assay_data_even <- matrix(
        rnorm(n_samples * n_features, mean=5, sd=1), nrow=n_features,
        ncol=n_samples, dimnames=list(feature_names_even, NULL))
    row_data_even <- data.frame(
        chain=as.numeric(feature_names_even), row.names=feature_names_even)
    se_even_only <- SummarizedExperiment(
        assays=list(abundance=assay_data_even), rowData=row_data_even,
        colData=col_data)
    warning_message <-
    expect_warning(
        analyzeLipidRegion(
            lipid_se=se_even_only, ref_group="control", split_chain=TRUE,
            chain_col="chain", radius=3, permute_time=100)
    )
})

test_that("analyzeLipidRegion validates input parameters", {
    se <- create_mock_se()
    expect_error(
        analyzeLipidRegion(se, ref_group="nonexistent"),
        "ref_group must be one of the groups in colData group column"
    )
    expect_error(
        analyzeLipidRegion(se, ref_group="control", split_chain="invalid"),
        "split_chain must be logical value"
    )
    expect_error(
        analyzeLipidRegion(se, ref_group="control", split_chain=TRUE),
        "chain_col is not present in the column name of rowData"
    )
})

test_that(".regionStat calculates statistics correctly", {
    X <- matrix(rnorm(20), nrow = 4, ncol = 5)
    Y <- matrix(c(rep(1, 2), rep(0, 2)), ncol = 1)
    t_result <- .regionStat(X, Y, test="t.test")
    expect_true(is.numeric(t_result))
    expect_equal(length(t_result), ncol(X))
    wilcox_result <- .regionStat(X, Y, test="Wilcoxon")
    expect_true(is.numeric(wilcox_result))
    expect_equal(length(wilcox_result), ncol(X))
})

test_that(".smooth_permutation handles edge cases", {
    X <- matrix(rnorm(20), nrow=4, ncol=5)
    group <- c(1, 1, 0, 0)
    dist.mat <- as.matrix(dist(t(X)))
    region.stat.obs <- rnorm(5)
    expect_no_error(
        .smooth_permutation(
            X, group, dist.mat, test="t.test", abund_weight=TRUE,
            permute_time=100, region.stat.obs=region.stat.obs)
    )
    expect_no_error(
        .smooth_permutation(
            X, group, dist.mat, test="t.test",
            abund_weight=FALSE, permute_time=100,
            region.stat.obs=region.stat.obs)
    )
    result <- .smooth_permutation(
        X, group, dist.mat, test="t.test", abund_weight=TRUE,
        permute_time=100, region.stat.obs=region.stat.obs)
    expect_true(all(is.finite(result$smooth.stat)))
    expect_true(all(is.finite(result$smooth.stat.permute)))
    result <- .smooth_permutation(
        X, group, dist.mat, test="t.test", abund_weight=FALSE,
        permute_time=100, region.stat.obs=region.stat.obs)
    expect_true(all(is.finite(result$smooth.stat)))
    expect_true(all(is.finite(result$smooth.stat.permute)))
})

test_that(".resultSummary formats output correctly", {
    X <- matrix(abs(rnorm(20, mean=10, sd=2)), nrow=4, ncol=5)
    X.info <- data.frame(feature=paste0("feature", seq_len(5)), x=seq_len(5))
    group <- c(1, 1, 0, 0)
    smooth.stat <- abs(rnorm(5))
    smooth.stat.permute <- matrix(abs(rnorm(500)), nrow=5)
    region.stat.obs <- abs(rnorm(5))
    dimension <- 1
    result <- .resultSummary(
        X, X.info, group, smooth.stat, smooth.stat.permute,
        region.stat.obs, dimension, permute_time=100)
    expect_true(is.data.frame(result))
    expect_true(all(c(
        "direction", "smoothing.pval.BH", "log2.FC") %in% colnames(result)))
    expect_equal(nrow(result), ncol(X))
    expect_true(all(is.finite(result$log2.FC)))
    expect_true(all(
        result$smoothing.pval.BH >= 0 & result$smoothing.pval.BH <= 1))
    expect_true(all(
        result$marginal.pval.BH >= 0 & result$marginal.pval.BH <= 1))
})

test_that("analyzeLipidRegion handles invalid input formats", {
    expect_error(
        analyzeLipidRegion(lipid_se=matrix(1:10), ref_group="control"),
        "Invalid lipid_se input"
    )
    se <- create_mock_se()
    expect_error(
        analyzeLipidRegion(
            lipid_se=se, ref_group="ctrl", split_chain=FALSE, chain_col=NULL,
            radius=3, permute_time=100),
        "ref_group must be one of the groups in colData group column"
    )
    expect_error(
        analyzeLipidRegion(
            lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
            radius=-1, permute_time=100),
        "radius must be a positive numeric value"
    )
})

test_that(".smooth_permutation handles small radius_values", {
    X <- matrix(abs(rnorm(20)), nrow=4, ncol=5)
    group <- c(1, 1, 0, 0)
    dist.mat <- as.matrix(dist(t(X)))
    region.stat.obs <- rnorm(5)
    radius_values <- 0.1
    expect_error(
        .smooth_permutation(
            X, group, dist.mat * radius_values, test="t.test",
            abund_weight=TRUE, permute_time=100,
            region.stat.obs=region.stat.obs),
        "cannot be achieved with the current"
    )
})

test_that("analyzeLipidRegion requires exactly 2 groups", {
    se <- create_mock_se(n_samples=9)
    colData(se)$group <- rep(c("a", "b", "c"), length.out=9)
    expect_error(
        analyzeLipidRegion(
            se, ref_group="a", split_chain=FALSE, chain_col=NULL,
            radius=3, permute_time=100),
        "analyzeLipidRegion\\(\\) requires exactly 2 groups"
    )
})
