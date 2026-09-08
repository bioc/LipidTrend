library(testthat)
library(SummarizedExperiment)

create_minimal_se <- function(n_features=4, n_samples=4) {
    assay_data <- matrix(
        rnorm(n_features * n_samples), nrow=n_features,
        dimnames=list(paste0("f", seq_len(n_features)), NULL))
    row_data <- data.frame(x=seq_len(n_features), row.names=rownames(assay_data))
    col_data <- data.frame(
        sample_name=paste0("s", seq_len(n_samples)),
        label_name=paste0("l", seq_len(n_samples)),
        group=rep(c("control", "case"), length.out=n_samples),
        row.names=paste0("s", seq_len(n_samples)))
    SummarizedExperiment(
        assays=list(abundance=assay_data), rowData=row_data, colData=col_data)
}

test_that("result()/even_chain_result()/odd_chain_result() return stored data", {
    se <- create_minimal_se()
    res_df <- data.frame(x=seq_len(4), direction=c("+", "+", "-", "-"))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=res_df,
        abund_weight=TRUE)
    expect_identical(result(lt_se), res_df)
    expect_null(even_chain_result(lt_se))
    expect_null(odd_chain_result(lt_se))

    lt_se_split <- new(
        "LipidTrendSE", se, split_chain=TRUE, even_chain_result=res_df,
        odd_chain_result=NULL, abund_weight=FALSE)
    expect_identical(even_chain_result(lt_se_split), res_df)
    expect_null(odd_chain_result(lt_se_split))
})

test_that(".split_chain()/.abund_weight() return the stored slot values", {
    se <- create_minimal_se()
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=data.frame(x=1),
        abund_weight=TRUE)
    expect_false(.split_chain(lt_se))
    expect_true(.abund_weight(lt_se))
})

test_that("show() prints result summary when split_chain=FALSE", {
    se <- create_minimal_se()
    res_df <- data.frame(x=seq_len(4), direction=c("+", "+", "-", "-"))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=res_df,
        abund_weight=TRUE)
    output <- capture.output(show(lt_se))
    expect_true(any(grepl("LipidTrend Results:", output)))
    expect_true(any(grepl("Split chain analysis: No", output)))
    expect_true(any(grepl("Result: 4 features", output)))
})

test_that("show() prints even/odd chain summaries when both present", {
    se <- create_minimal_se()
    res_df <- data.frame(x=seq_len(4), direction=c("+", "+", "-", "-"))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=TRUE, even_chain_result=res_df,
        odd_chain_result=res_df[1:2, ], abund_weight=TRUE)
    output <- capture.output(show(lt_se))
    expect_true(any(grepl("Split chain analysis: Yes", output)))
    expect_true(any(grepl("Even chain result: 4 features", output)))
    expect_true(any(grepl("Odd chain result: 2 features", output)))
})

test_that("show() skips the odd chain line when odd_chain_result is NULL", {
    se <- create_minimal_se()
    res_df <- data.frame(x=seq_len(4), direction=c("+", "+", "-", "-"))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=TRUE, even_chain_result=res_df,
        odd_chain_result=NULL, abund_weight=TRUE)
    output <- capture.output(show(lt_se))
    expect_true(any(grepl("Even chain result: 4 features", output)))
    expect_false(any(grepl("Odd chain result:", output)))
})
