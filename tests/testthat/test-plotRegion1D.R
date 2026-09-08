library(testthat)
library(SummarizedExperiment)

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

test_that("plotRegion1D creates valid plot", {
    set.seed(1234)
    se <- create_mock_se()
    se_split <- create_mock_se(split_chain=TRUE)
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3, permute_time=100)
    expect_no_error(plot <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot, "ggplot"))
    # positive regions
    result_df <- result(result)
    result_df$direction <- "+"
    result_df$smoothing.pval.BH <- 0.01
    result_df$avg.expr.ctrl <- 1
    result_df$avg.expr.case <- 2
    attr(result, "result") <- result_df
    expect_no_error(plot_pos <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot_pos, "ggplot"))
    # negative regions
    result_df$direction <- "-"
    result_df$avg.expr.ctrl <- 2
    result_df$avg.expr.case <- 1
    attr(result, "result") <- result_df
    expect_no_error(plot_neg <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot_neg, "ggplot"))
    # mixed positive and negative regions
    result_df$direction <- rep(c("+", "-"), length.out=nrow(result_df))
    result_df$smoothing.pval.BH <- rep(c(0.01, 0.1), length.out=nrow(result_df))
    attr(result, "result") <- result_df
    expect_no_error(plot_mixed <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot_mixed, "ggplot"))
})

test_that("plotRegion1D handles split chain analysis correctly", {
    set.seed(1234)
    se_split <- create_mock_se(split_chain=TRUE)
    result <- analyzeLipidRegion(
        lipid_se=se_split, ref_group="control", split_chain=TRUE,
        chain_col="chain", radius=3, permute_time=100)
    # both even and odd chain
    expect_no_error(plot <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot$even_result, "ggplot"))
    expect_true(inherits(plot$odd_result, "ggplot"))
    # one chain type has no significant regions
    even_results <- even_chain_result(result)
    even_results$smoothing.pval.BH <- 1
    attr(result, "even_chain_results") <- even_results
    expect_no_error(plot_one_sig <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
})

test_that("plotRegion1D handles invalid input and edge cases", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3, permute_time=100)
    # invalid p_cutoff values
    expect_error(
        plotRegion1D(result, p_cutoff=-1, y_scale='identity'),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    expect_error(
        plotRegion1D(result, p_cutoff=2, y_scale='identity'),
        "p_cutoff must be a numeric value between 0 and 1"
    )
    # invalid regional_p_cutoff values
    expect_error(
        plotRegion1D(result, regional_p_cutoff=-1),
        "regional_p_cutoff must be a numeric value between 0 and 1"
    )
    expect_error(
        plotRegion1D(result, regional_p_cutoff=2),
        "regional_p_cutoff must be a numeric value between 0 and 1"
    )
    # no significant regions
    result_df <- result(result)
    result_df$smoothing.pval.BH <- 1
    attr(result, "result") <- result_df
    expect_no_error(plot_no_sig <- plotRegion1D(
        result, p_cutoff=0.05, y_scale='identity'))
    expect_true(inherits(plot_no_sig, "ggplot"))
    # missing results
    bad_result <- result
    attr(bad_result, "result") <- NULL
    expect_error(plotRegion1D(bad_result, p_cutoff=0.05, y_scale='identity'))
    # invalid split chain status
    bad_split_result <- result
    attr(bad_split_result, "split_chain") <- "invalid"
    expect_error(plotRegion1D(
        bad_split_result, p_cutoff=0.05, y_scale='identity'))
})

test_that("plotRegion1D handles various p-value cutoffs", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, permute_time=100)
    p_cutoffs <- c(0.01, 0.05, 0.1)
    for(p_cut in p_cutoffs) {
        expect_no_error(plot <- plotRegion1D(
            result, p_cutoff=p_cut, y_scale='identity'))
        expect_true(inherits(plot, "ggplot"))
    }
})

# Builds a dataset with a strong, engineered group difference so that both
# the Increase and Decrease regions reliably end up with more than 5
# significant features -- this is what triggers plotRegion1D()'s new
# paired-test-based coloring branch (see regionalTestFC()). Plain
# create_mock_se() noise does not reliably land on either side of that
# n.features > 5 boundary, so it cannot be used to test this branch.
create_strong_signal_se <- function(n_features=16, n_samples=8) {
    base_chains <- seq(33, 33 + n_features - 1)
    feature_names <- paste0("TG_", base_chains)
    assay_data <- matrix(
        rnorm(n_samples * n_features, mean=5, sd=0.3),
        nrow=n_features, ncol=n_samples, dimnames=list(feature_names, NULL))
    group <- rep(c("control", "case"), each=n_samples / 2)
    case_idx <- which(group == "case")
    assay_data[1:7, case_idx] <- assay_data[1:7, case_idx] + 3
    assay_data[8:14, case_idx] <- assay_data[8:14, case_idx] - 3
    row_data <- data.frame(x=base_chains, row.names=feature_names)
    col_data <- data.frame(
        sample_name=paste0("s", seq_len(n_samples)),
        label_name=paste0("l", seq_len(n_samples)), group=group,
        row.names=paste0("s", seq_len(n_samples)))
    SummarizedExperiment(
        assays=list(abundance=assay_data), rowData=row_data, colData=col_data)
}

test_that("plotRegion1D switches to paired-test coloring when both directions have >5 significant features", {
    set.seed(42)
    se <- create_strong_signal_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3, permute_time=500)
    tab <- regionalTestFC(result)
    # sanity check: the engineered signal really does clear the n>5 gate
    expect_true(all(tab$n.features > 5))
    decision <- .regionColorDecision(tab, 0.05)
    expect_equal(unname(decision), unname(
        tab$regional.test.pval[match(names(decision), tab$direction)] < 0.05))
    expect_no_error(plot <- plotRegion1D(result, p_cutoff=0.05))
    expect_true(inherits(plot, "ggplot"))
})

test_that("plotRegion1D's regional_p_cutoff independently controls the paired-test validation", {
    set.seed(42)
    se <- create_strong_signal_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, chain_col=NULL,
        radius=3, permute_time=500)
    tab <- regionalTestFC(result)
    expect_true(all(tab$n.features > 5))
    # derive cutoffs from the actual p-values so this doesn't depend on how
    # extreme the engineered signal's p-values happen to be
    strict_cutoff <- min(tab$regional.test.pval) / 100
    loose_cutoff <- min(1, max(tab$regional.test.pval) * 10)
    expect_true(all(!.regionColorDecision(tab, strict_cutoff)))
    expect_true(all(.regionColorDecision(tab, loose_cutoff)))
    expect_no_error(
        plotRegion1D(result, p_cutoff=0.05, regional_p_cutoff=strict_cutoff))
    expect_no_error(
        plotRegion1D(result, p_cutoff=0.05, regional_p_cutoff=loose_cutoff))
})

test_that("plotRegion1D handles various y_scale", {
    se <- create_mock_se()
    result <- analyzeLipidRegion(
        lipid_se=se, ref_group="control", split_chain=FALSE, permute_time=100)
    y_scales <- c('identity', 'log2', 'log10', 'sqrt')
    for(y_scale in y_scales) {
        expect_no_error(plot <- plotRegion1D(
            result, p_cutoff=0.05, y_scale=y_scale))
        expect_true(inherits(plot, "ggplot"))
    }
    expect_error(
        plotRegion1D(result, p_cutoff=0.05, y_scale=NULL),
        "y_scale must be one of identity, log2, log10, or sqrt"
    )
})
