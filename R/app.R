#' Lakhesis Calculator
#'
#' Launch Lakhesis Calculator, a graphical interface to explore binary matrices via correspondence analysis, select potentially well-seriated sequences, and perform consensus seriation. Interface is made with \code{ggplot2}, \code{shiny}, \code{shinydashboard}, and \code{bslib} \insertCite{wickham_ggplot2:_2016,chang_shiny_2024,chang_shinydashboard_2021,sievert_bslib_2024}{lakhesis}. 
#' 
#' Input is done in the calculator, via a "long" format a two-column \code{.csv} file giving pairs of row and column incidences. See \code{\link[lakhesis]{im.csv.read}} for details. Conversion of a pre-existing incidence matrix to long format can be performed with \code{\link[lakhesis]{im.long}}.
#' 
#' Results can be downloaded from the calculator as an \code{.rds} file containing a `list` of the following:
#' 
#' * `results` The consensus seriation, PCA, and coefficients of agreement and concentration 
#' (\code{\link[lakhesis]{lakhesize}}).
#' * `strands` The strands selected by the investigator.
#' * `im.seriated` The incidence matrix of the consensus seriation.
#' 
#' @returns Opens the Lakhesis Calculator. 

#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
LC <- function() {
    shiny::runApp(appDir = system.file('app', package = 'lakhesis') )
}

