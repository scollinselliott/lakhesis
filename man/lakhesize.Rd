% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lakhesize.R
\name{lakhesize}
\alias{lakhesize}
\alias{lakhesize.strands}
\alias{lakhesize.default}
\title{Lakhesize}
\usage{
lakhesize(strands, crit = "cor_sq", pbar = TRUE)

\method{lakhesize}{strands}(strands, crit = "cor_sq", pbar = TRUE)

\method{lakhesize}{default}(strands, crit = "cor_sq", pbar = TRUE)
}
\arguments{
\item{strands}{A \code{list} of \code{strands} class (see \code{\link[lakhesis]{strand_add}}).}

\item{crit}{The criterion used to assess the seration resulting from two strands \eqn{x} and \eqn{y}:
\itemize{
\item \code{"cor_sq"} Computes a squared correlation coefficent for the incidences of 1s as points \eqn{(i,j)}. See \code{\link[lakhesis]{cor_sq}}. Higher values are more optimal. This is the default option.
\item \code{"conc_wrc"} Computes weighted row-column concentration as the optimality criterion. See \code{\link[lakhesis]{conc_wrc}}. Lower values are more optimal.
}}

\item{pbar}{Displaying a progress bar. Default is \code{TRUE}.}
}
\value{
A \code{list} of class \code{lakhesis} containing the following:
\itemize{
\item \code{row} A seriated vector of row elements.
\item \code{col} A seriated vector of column elements
\item \code{coef}  A \code{data frame} containing the following columns:
\itemize{
\item \code{Strand} The number of the strand.
\item \code{Agreement} The measure of agreement, i.e., how well each strand accords with the consensus seriation. Using the square of Spearman's rank correlation coefficient, \eqn{\rho^2}, between each strand and the consensus ranking, agreement is computed as the product of \eqn{\rho^2} for their row and column rankings, \eqn{\rho_r^2}\eqn{\rho_c^2}.
\item \code{Criterion} Criterion of the optimality of each strand (per the \code{"crit"} option above).
}
\item \code{im_seriated} The seriated incidence matrix, of class \code{incidence_matrix}.
}
}
\description{
This function returns the row and column consensus seriation for a \code{list} object of the \code{strands} class, containing their rankings, coefficients of association, and criterion. Consensus seriation is achieved by iterative simple linear regression to handle \code{NA} vales in each strand. To initialize, a regression is performed pairwise, with every strand as the dependent \eqn{y} variate and every other strand as the independent \eqn{x} variate. The independent variate's rankings are then regressed onto \eqn{f(x) = \hat{\beta}_1 x + \hat{\beta}_0}. If \eqn{y \neq f(x)}, the mean of \eqn{y} and \eqn{f(x)} is used. Then, the values of the dependent variate and those of the regressed independent varaite are re-ranked together, to form a combined ranking, which serves as the dependent variate on the next iteration. The pair of strands is chosen which minimizes a specified optimality criterion. The process is repeated until all strands have been regressed and re-ranked into a single consensus seriation.
}
\examples{
data("qf_strands")
L <- lakhesize(qf_strands, pbar = FALSE)
# summary(L) 

}
