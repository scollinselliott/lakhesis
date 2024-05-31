#' Quattro Fontanili
#'
#' The seriation of tombs from necropoleis at Veii, primarily Quattro Fontanili, but also Valle la Fata, Vaccareccia, and Picazzano, in southern Etruria, established by \insertCite{close-brooks_veii_1979;textual}{lakhesis}.
#'
#' @format
#' A seriated incidence matrix of 81 rows (tombs) and 82 columns (types).
#' \describe{
#' Data entered from \insertCite{close-brooks_veii_1979;textual}{lakhesis}, an English translation of the authors' original publication in \emph{Notizie degli Scavi} (1963). Descriptions of types may be found in that paper.
#' }
#' @usage data("quattrofontanili")
#' @examples
#' data("quattrofontanili")
#' print(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
"quattrofontanili"


#' Quattro Fontanili - Lakhesis Results
#'
#' The results of \code{\link[lakhesis]{lakhesize}} applied to \code{\link[lakhesis]{quattrofontanili}} data (exported and reimported into the Calculator). Seven strands were selected by the package author as an example for the documentation of functions. Consult the \code{\link[lakhesis]{lakhesize}} function for the structure of the object.
#'
#' @format
#' A list containing the following objects output by \code{\link[lakhesis]{lakhesize}}.
#' }
#' @usage data("qfLakhesis")
#' @examples
#' data("qfLakhesis")
#' print(qfLakhesis)
#' 
"qfLakhesis"



