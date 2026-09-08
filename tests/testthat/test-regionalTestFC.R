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

# Builds a minimal chain_result data.frame directly (bypassing
# analyzeLipidRegion's permutation test) so n.features / direction can be
# set exactly, instead of hoping random data happens to land on a boundary.
craft_chain_result <- function(
        se, incr_idx=integer(0), decr_idx=integer(0), sig_pval=0.001) {
    n <- nrow(se)
    direction <- rep("+", n)
    direction[decr_idx] <- "-"
    smoothing.pval.BH <- rep(1, n)
    smoothing.pval.BH[c(incr_idx, decr_idx)] <- sig_pval
    log2.FC <- rep(0.05, n)
    log2.FC[incr_idx] <- 2
    log2.FC[decr_idx] <- -2
    data.frame(
        row.names=rownames(se), direction=direction,
        smoothing.pval.BH=smoothing.pval.BH, log2.FC=log2.FC)
}

set.seed(1234)

test_that("regionalTestFC returns the expected columns and counts", {
    se <- create_mock_se(n_features=12)
    chain_result <- craft_chain_result(se, incr_idx=1:3, decr_idx=4:5)
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    tab <- regionalTestFC(lt_se)
    expect_true(is.data.frame(tab))
    expect_equal(
        colnames(tab),
        c("direction", "n.features", "regional.test.pval", "regional.FC",
          "log2.regional.FC"))
    expect_equal(tab$direction, c("Increase", "Decrease"))
    expect_equal(tab[tab$direction == "Increase", "n.features"], 3L)
    expect_equal(tab[tab$direction == "Decrease", "n.features"], 2L)
    # the two scales must agree: regional.FC == 2^log2.regional.FC
    expect_equal(tab$regional.FC, round(2^tab$log2.regional.FC, 2))
})

test_that("regionalTestFC returns NA regional.test.pval when n.features is 0 or 1", {
    se <- create_mock_se(n_features=12)
    # 0 increase, 1 decrease
    chain_result <- craft_chain_result(se, incr_idx=integer(0), decr_idx=4)
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    tab <- regionalTestFC(lt_se)
    incr_row <- tab[tab$direction == "Increase", ]
    decr_row <- tab[tab$direction == "Decrease", ]
    expect_equal(incr_row$n.features, 0L)
    expect_true(is.na(incr_row$regional.test.pval))
    expect_true(is.na(incr_row$regional.FC))
    expect_true(is.na(incr_row$log2.regional.FC))
    expect_equal(decr_row$n.features, 1L)
    expect_true(is.na(decr_row$regional.test.pval))
    expect_false(is.na(decr_row$regional.FC))
})

test_that("regionalTestFC returns NA regional.test.pval when 2 <= n.features <= 5", {
    se <- create_mock_se(n_features=12)
    # n.features = 5: technically enough for a paired t-test, but still
    # underpowered -- should stay NA, unlike the old n < 2 threshold.
    chain_result <- craft_chain_result(se, incr_idx=1:5, decr_idx=integer(0))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    tab <- regionalTestFC(lt_se)
    incr_row <- tab[tab$direction == "Increase", ]
    expect_equal(incr_row$n.features, 5L)
    expect_true(is.na(incr_row$regional.test.pval))
    expect_false(is.na(incr_row$regional.FC))
})

test_that("regionalTestFC computes a paired t-test when n.features > 5", {
    se <- create_mock_se(n_features=12)
    chain_result <- craft_chain_result(se, incr_idx=1:6, decr_idx=7:12)
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    tab <- regionalTestFC(lt_se)
    incr_row <- tab[tab$direction == "Increase", ]
    decr_row <- tab[tab$direction == "Decrease", ]
    expect_false(is.na(incr_row$regional.test.pval))
    expect_true(incr_row$regional.test.pval >= 0 && incr_row$regional.test.pval <= 1)
    expect_false(is.na(decr_row$regional.test.pval))
})

test_that("regionalTestFC respects p_cutoff", {
    se <- create_mock_se(n_features=12)
    chain_result <- craft_chain_result(se, incr_idx=1:3, decr_idx=integer(0))
    chain_result$smoothing.pval.BH[1:3] <- 0.03
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    tab_default <- regionalTestFC(lt_se, p_cutoff=0.05)
    tab_strict <- regionalTestFC(lt_se, p_cutoff=0.01)
    expect_equal(tab_default[tab_default$direction == "Increase", "n.features"], 3L)
    expect_equal(tab_strict[tab_strict$direction == "Increase", "n.features"], 0L)
})

test_that("regionalTestFC validates p_cutoff", {
    se <- create_mock_se(n_features=12)
    chain_result <- craft_chain_result(se, incr_idx=1:3, decr_idx=4:5)
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    expect_error(
        regionalTestFC(lt_se, p_cutoff=-1),
        "p_cutoff must be a numeric value between 0 and 1")
    expect_error(
        regionalTestFC(lt_se, p_cutoff=2),
        "p_cutoff must be a numeric value between 0 and 1")
})

test_that("regionalTestFC requires exactly 2 groups", {
    se <- create_mock_se(n_features=12, n_samples=9)
    colData(se)$group <- rep(c("a", "b", "c"), length.out=9)
    chain_result <- craft_chain_result(se, incr_idx=1:3, decr_idx=4:5)
    lt_se <- new(
        "LipidTrendSE", se, split_chain=FALSE, result=chain_result,
        abund_weight=TRUE)
    expect_error(regionalTestFC(lt_se), "requires exactly 2 groups")
})

test_that("regionalTestFC handles split_chain=TRUE", {
    se <- create_mock_se(n_features=12, split_chain=TRUE)
    even_result <- craft_chain_result(se, incr_idx=1:3, decr_idx=4:5)
    odd_result <- craft_chain_result(se, incr_idx=integer(0), decr_idx=integer(0))
    lt_se <- new(
        "LipidTrendSE", se, split_chain=TRUE, even_chain_result=even_result,
        odd_chain_result=odd_result, abund_weight=TRUE)
    tab <- regionalTestFC(lt_se)
    expect_named(tab, c("even_result", "odd_result"))
    expect_equal(
        tab$even_result[tab$even_result$direction == "Increase", "n.features"],
        3L)
    expect_equal(tab$odd_result$n.features, c(0L, 0L))
})

test_that(".regionColorDecision keeps original coloring unless both directions have a usable p-value", {
    # both_usable mirrors what .pairedFCTest() actually returns when
    # n.features > 5 for both directions: a real (non-NA) regional.test.pval.
    both_usable <- data.frame(
        direction=c("Increase", "Decrease"), n.features=c(6, 6),
        regional.test.pval=c(0.2, 0.01), regional.FC=c(2, 0.5),
        log2.regional.FC=c(1, -1))
    decision <- .regionColorDecision(both_usable, 0.05)
    expect_equal(unname(decision["Increase"]), FALSE)
    expect_equal(unname(decision["Decrease"]), TRUE)

    # one_unusable mirrors n.features <= 5 for Decrease: .pairedFCTest()
    # would report NA there, which must fall back the whole figure even
    # though Increase on its own has a usable, non-significant p-value.
    one_unusable <- data.frame(
        direction=c("Increase", "Decrease"), n.features=c(6, 3),
        regional.test.pval=c(0.2, NA_real_), regional.FC=c(2, 0.5),
        log2.regional.FC=c(1, -1))
    decision2 <- .regionColorDecision(one_unusable, 0.05)
    expect_true(all(decision2))
})

test_that(".regionColorDecision respects a custom regional_p_cutoff", {
    both_usable <- data.frame(
        direction=c("Increase", "Decrease"), n.features=c(6, 6),
        regional.test.pval=c(0.03, 0.03), regional.FC=c(2, 0.5),
        log2.regional.FC=c(1, -1))
    # default-like 0.05: both pass
    expect_true(all(.regionColorDecision(both_usable, 0.05)))
    # stricter than 0.03: both fail
    expect_false(any(.regionColorDecision(both_usable, 0.01)))
    # looser than 0.03: both still pass
    expect_true(all(.regionColorDecision(both_usable, 0.1)))
})
