#' Lakhesis Calculator
#'
#' Launch Lakhesis Calculator, a graphical interface to explore binary matrices via correspondence analysis, 
#' select potentially well-seriated sequences, and perform consensus seriation. 
#' Interface is made with \code{ggplot2}, \code{shiny}, \code{shinydashboard}, \code{bslib} 
#' \insertCite{wickham_ggplot2:_2016,chang_shiny_2024,chang_shinydashboard_2021,sievert_bslib_2024}{lakhesis}. 
#' 
#' @returns Opens the Lakhesis Calculator. CSV files are imported (see \link[lakhesis]{im.csv.read} for formatting). 
#' Results are download as an \code{/rds} file containing a list object `results` of  the following:
#' * `results` The consensus seriation achieved, including coefficients of agreement and concentration 
#' (\code{\link[lakhesis]{lakhesize}}).
#' * `strands` The strands selected by the investigator.
#' * `im.seriated` The incidence matrix of the consensus seriation.
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
LC <- function() {
    shiny::runApp(appDir = system.file("app", package = "lakhesis"))
}

