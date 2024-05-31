#' Spearman Correlation Squared
#' 
#' The square of Spearman's rank correlation coefficient applied to two rankings \insertCite{spearman_proof_1904}{lakhesis}. Rows with `NA` values are automatically removed.
#' 
#' @param r1,r2 Two vectors of paired ranks.
#' @returns The square of Spearman's rank correlation coefficient with NA values removed.
#' 
#' @examples
#' # e.g., for two partial seriations:
#' x <- c(1, 2, 3, 4, NA, 5, 6, NA, 7.5, 7.5, 9)
#' y <- c(23, 17, 19, NA, 21, 22, 25, 26, 27, 36, 32)
#' spearman.sq(x, y)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
spearman.sq <- function(r1, r2) {
    dat <-  stats::na.omit( data.frame(r1,r2) )
    r <- (stats::cor(dat)[1,2])^2
    return(r)
    }


#' Kendall-Doran Concentration
#'
#' The Kendall-Doran measure of concentration \insertCite{kendall_statistical_1963,doran_computer_1971}{lakhesis}. In a seriated matrix, this function computes the total number cells between the first and last non-zero value, column by column.
#' 
#' @param obj A seriated binary matrix.
#' @returns The measure of concentration.
#' 
#' @examples 
#' data("quattrofontanili")
#' concentration.col(quattrofontanili)
#' 
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
concentration.col <- function(obj) {
conc <- numeric(ncol(obj))
    for (j in 1:ncol(obj)) {
        conc[j] <- max(which(obj[,j] != 0)) - min(which(obj[,j] != 0)) + 1
    }
    return(sum(conc))
}



#' Kappa Concentration
#'
#' The concentration coefficient \eqn{\kappa}, which extends the Kendall-Doran measure of concentration to include rows and then weights the total measure by the total sum of values in the matrix. See \code{\link[lakhesis]{concentration.col}}.
#' 
#' @param obj A seriated binary matrix.
#' @returns The \eqn{\kappa} coefficient of concentration.
#' 
#' @examples 
#' data("quattrofontanili")
#' kappa.coef(quattrofontanili)
#' 
#' @export kappa.coef
kappa.coef <- function(obj) {
     nu <- sum(obj)
     k <- ( concentration.col(obj) + concentration.col(t(obj)) ) / (2 * nu)
     return(k)
}



#' Evaluating Element Fit
#'
#' Performs a goodness-of-fit test on individual row and column elements using deviance, using a quadratic-logistic model to fit row and column occurrences. In the case of perfect separation of 0/1 values, an `NA` value is assigned. Results are reported as \eqn{p} values for each row and column.
#' 
#' @param obj A seriated binary matrix.
#' @returns A `list` containing results in frames for row and column elements:
#' 
#' * `RowFit` a data frame containing
#'   * `id` Row element
#'   * `p.val` \eqn{p} values of the row elements
#' * `ColFit` a data frame containing
#'   * `id` Column element
#'   * `p.val` \eqn{p} values of the column elements
#' 
#' @examples 
#' data("quattrofontanili")
#' element.eval(quattrofontanili)
#' 
#' @export
element.eval <- function(obj) {
    dev.c <- numeric(ncol(obj))
    dev.c[] <- NA
    for (j in 1:ncol(obj)) {
        conc_j <- max(which(obj[,j] != 0)) - min(which(obj[,j] != 0)) + 1
        conc_js <- conc_j / sum(obj[,j])
        if (conc_js != 1) {
            x <- 1:length(obj[,j])
            y <- obj[,j]
            suppressWarnings({
            fit <- stats::glm(y ~ x + I(x^2), family = stats::binomial(link = "logit"))
            })
            dev.res <- fit$deviance
            dev.nul <- fit$null.deviance
            df.res <- fit$df.residual
            df.nul <- fit$df.null

            pval <- 1 - stats::pchisq(dev.nul - dev.res, df.nul - df.res)
            dev.c[j] <- pval
        } 
    }
    ColFit = data.frame(id = colnames(obj), p.val = dev.c)
    ColFit <- ColFit[order(ColFit$p.val, decreasing = TRUE) , ]

    dev.r <- numeric(nrow(obj))
    dev.r[] <- NA
    for (i in 1:nrow(obj)) {
        conc_i <- max(which(obj[i,] != 0)) - min(which(obj[i,] != 0)) + 1
        conc_is <- conc_i / sum(obj[i,])
        if (conc_is != 1) {
            x <- 1:length(obj[i,])
            y <- obj[i,]
            suppressWarnings({
            fit <- stats::glm(y ~ x + I(x^2), family = stats::binomial(link = "logit"))
            })
            dev.res <- fit$deviance
            dev.nul <- fit$null.deviance
            df.res <- fit$df.residual
            df.nul <- fit$df.null

            pval <- 1-stats::pchisq(dev.nul - dev.res, df.nul - df.res)
            dev.r[i] <- pval
        } 
    }
    RowFit <- data.frame(id = rownames(obj), p.val = dev.r)
    RowFit <- RowFit[order(RowFit$p.val, decreasing = TRUE) , ]

    results <- list(RowFit = RowFit, ColFit = ColFit)
    return(results)
}






#' Strand Extract
#' 
#' From the results of \code{\link[lakhesis]{ca.procrustes.curve}}, extrect two matrices containing the ranks of the rows and columns. The row/column elements are contained in the rows, and the strands are contained in the columns. NA values are entered where a given row/column element is missing from that strand.
#' 
#' @param strands A list of `strands`, which are data frames returned by \code{\link[lakhesis]{ca.procrustes.curve}}.
#' @param obj The intial incidence matrix.
#' 
#' @return A list of.
#' * `Row` A matrix of the ranks of the row elements.
#' * `Col` A matrix of teh ranks of the column elements.
#' 
#' @export
strand.extract <- function(strands, obj) {

    rowranks <- matrix(NA, nrow = nrow(obj),  ncol = length(strands)) 
    rownames(rowranks) <- rownames(obj)

    colranks <- matrix(NA, nrow = ncol(obj),  ncol = length(strands)) 
    rownames(colranks) <- colnames(obj)

    for (i in 1:length(strands)) {
        strand <- strands[[i]]

        strand.r <- strand[strand$Type == "row", ]
        strand.c <- strand[strand$Type == "col",]

        rowranks[ match(rownames(strand.r), rownames(rowranks)) , i] <- strand.r$Rank
        colranks[ match(rownames(strand.c), rownames(colranks)) , i] <- strand.c$Rank
    }

    rowranks <- rowranks[rowSums(is.na(rowranks)) != ncol(rowranks), ]
    colranks <- colranks[rowSums(is.na(colranks)) != ncol(colranks), ]

    results <- list()
    results[["Row"]] <- rowranks
    results[["Col"]] <- colranks
    return(results)
}



#' Suppress Element from Strands
#' 
#' Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca.procrustes.curve}}) to generate a new list of strands. The row or column element id is stored in an output list in case it should be added back into the strand. Column and row names are required to be unique. If the resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
#' 
#' @param strands A list of strands, which are data frames returned by \code{\link[lakhesis]{ca.procrustes.curve}}.
#' @param obj The intial incidence matrix.
#' @param elements A vector of one or more row or column ids to suppress.
#' 
#' @return A list of the strands.
#' 
#' @export
strand.suppress <- function(strands, obj, elements) {
    new <- list()
    for (i in 1:length(strands)) {
        rows <- rownames(strands[[i]])[ !(rownames(strands[[i]]) %in% elements) & (strands[[i]]$Type == "row") ]
        cols <- rownames(strands[[i]])[ !(rownames(strands[[i]]) %in% elements) & (strands[[i]]$Type == "col") ]
        obj.copy <- obj
        obj.copy <- obj.copy[rows, ]
        obj.copy <- obj.copy[ ,cols]
        obj.copy <- obj.copy[rowSums(obj.copy) !=0 , ]
        obj.copy <- obj.copy[ , colSums(obj.copy) !=0]
        newstrand <- NULL
        suppressWarnings({
            newstrand <- ca.procrustes.curve(obj.copy)
        })
        if (!is.null(newstrand)) { 
           new[[i]] <- newstrand
        }
    }
    return(new)
}


