% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strand.R
\name{strand_suppress}
\alias{strand_suppress}
\alias{strand_suppress.strands}
\alias{strand_suppress.default}
\title{Suppress Element from Strands}
\usage{
strand_suppress(strands, ...)

\method{strand_suppress}{strands}(strands, ...)

\method{strand_suppress}{default}(strands, elements)
}
\arguments{
\item{strands}{A \code{list} of class \code{strands}.}

\item{elements}{A vector of one or more row or column ids to suppress.}
}
\value{
A list of the strands.
}
\description{
Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca_procrustes_ser}}) to generate a new list of strands that exclude the specified elements. If the resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
}
\examples{
data("quattrofontanili")
data("qfStrands")
strand_suppress(qfStrands, "QF II 15-16")

strand_suppress(qfStrands, c("QF II 15-16", "I", "XIV"))

}
