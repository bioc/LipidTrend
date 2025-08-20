#' Example Dataset for one-dimensional data
#'
#' @docType data
#' @usage data(lipid_se_CL)
#' @return A `SummarizedExperiment` object with the following slots:
#' \describe{
#'   \item{\code{colData}}{ A data frame with 6 observations on the following
#'    3 variables, containing sample name, label name, and group name}
#'   \item{\code{assay}}{ A 29*6 matrix containing lipid abundance data}
#'   \item{\code{rowData}}{ A data frame with 29 observations on the following
#'   1 variables containing chain characteristic information}
#' }
#'
#' @examples
#' data(lipid_se_CL)
"lipid_se_CL"

#' Example Dataset for two-dimensional data
#'
#' @docType data
#' @usage data(lipid_se_2D)
#' @return A `SummarizedExperiment` object with the following slots:
#' \describe{
#'   \item{\code{colData}}{ A data frame with 6 observations on the following
#'    3 variables, containing sample name, label name, and group name}
#'   \item{\code{assay}}{ A 137*6 matrix containing lipid abundance data}
#'   \item{\code{rowData}}{ A data frame with 137 observations on the following
#'   2 variables containing total chain length and total double bond
#'   characteristic information}
#' }
#'
#' @examples
#' data(lipid_se_2D)
"lipid_se_2D"

#' Example lipid abundance data for one-dimensional LipidTrend analysis
#'
#' @docType data
#' @usage data(abundance_CL)
#' @return A `matrix` object of lipid abundance with 29 lipids over 6 samples
#'
#' @source Tomoyuki Shiota et al. ,Hepatoviruses promote very-long-chain fatty
#' acid and sphingolipid synthesis for viral RNA replication and quasi-enveloped
#' virus release. Sci Adv. 9(42):eadj4198
#'  \url{https://www.science.org/doi/10.1126/sciadv.adj4198}.
#' @examples
#' data(abundance_CL)
"abundance_CL"

#' Example lipid abundance data for two-dimensional LipidTrend analysis
#'
#' @docType data
#' @usage data(abundance_2D)
#' @return A `matrix` object of lipid abundance with 137 lipids over 6 samples
#'
#' @source Tomoyuki Shiota et al. ,Hepatoviruses promote very-long-chain fatty
#' acid and sphingolipid synthesis for viral RNA replication and quasi-enveloped
#' virus release. Sci Adv. 9(42):eadj4198
#'  \url{https://www.science.org/doi/10.1126/sciadv.adj4198}.
#' @examples
#' data(abundance_2D)
"abundance_2D"

#' Example lipid characteristics table for one-dimensional LipidTrend analysis
#'
#' @docType data
#' @usage data(char_table_CL)
#' @return A `data.frame` object of chain characteristics over 29 lipids
#'
#' @source Tomoyuki Shiota et al. ,Hepatoviruses promote very-long-chain fatty
#' acid and sphingolipid synthesis for viral RNA replication and quasi-enveloped
#' virus release. Sci Adv. 9(42):eadj4198
#'  \url{https://www.science.org/doi/10.1126/sciadv.adj4198}.
#' @examples
#' data(char_table_CL)
"char_table_CL"

#' Example lipid characteristics table for two-dimensional LipidTrend analysis
#'
#' @docType data
#' @usage data(char_table_2D)
#' @return A `data.frame` object of total chain length and total double bond
#' characteristics over 137 lipids
#'
#' @source Tomoyuki Shiota et al. ,Hepatoviruses promote very-long-chain fatty
#' acid and sphingolipid synthesis for viral RNA replication and quasi-enveloped
#' virus release. Sci Adv. 9(42):eadj4198
#'  \url{https://www.science.org/doi/10.1126/sciadv.adj4198}.
#' @examples
#' data(char_table_2D)
"char_table_2D"

#' Example group information table for LipidTrend analysis
#'
#' @docType data
#' @usage data(group_info)
#' @return A `data.frame` object of sample name, lable name, and group name
#' over 6 samples
#'
#' @source Tomoyuki Shiota et al. ,Hepatoviruses promote very-long-chain fatty
#' acid and sphingolipid synthesis for viral RNA replication and quasi-enveloped
#' virus release. Sci Adv. 9(42):eadj4198
#'  \url{https://www.science.org/doi/10.1126/sciadv.adj4198}.
#' @examples
#' data(group_info)
"group_info"
