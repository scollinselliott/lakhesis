% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criterion.R
\name{conc_wrc}
\alias{conc_wrc}
\alias{conc_wrc.matrix}
\alias{conc_wrc.incidence_matrix}
\title{Optimality Criterion: Weighted Row-Column Concentration}
\usage{
conc_wrc(obj)

\method{conc_wrc}{matrix}(obj)

\method{conc_wrc}{incidence_matrix}(obj)
}
\arguments{
\item{obj}{A seriated binary matrix.}
}
\value{
The weighted row-column coefficient of concentration.
}
\description{
Extends the Kendall-Doran (column) measure of concentration (see \code{\link[lakhesis]{conc_c}}) to include rows and then weights the total measure by the total sum of values in the matrix.
}
\examples{
data("quattrofontanili")
conc_wrc(quattrofontanili)

}
