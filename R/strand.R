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
    dat <-  na.omit( data.frame(r1,r2) )
    r <- (cor(dat)[1,2])^2
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
#' # example using Hodson's data from Münsingen-Rain
#' library(seriation)
#' data(Munsingen)
#' munsingen <- Munsingen
#' 
#' # row and column names should be unique
#' rownames(munsingen) <- paste("Context", rownames(munsingen))
#' colnames(munsingen) <- paste("Find", colnames(munsingen))
#' 
#' concentration.col(munsingen)
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
#' The concentration coefficient \eqn{\kappa}, which extends the Kendall-Doran measure of concentration to include rows, and then weights the total measure by the total sum of the matrix. See \code{\link[lakhesis]{concentration.col}}.
#' 
#' @param obj A seriated binary matrix.
#' @returns The \eqn{\kappa} coefficient of concentration.
#' 
#' @examples 
#' # using the Münsingen data from the Münsingen-Rain
#' library(seriation)
#' munsingen <- data(Munsingen)
#' 
#' # row and column names should be unique
#' rownames(munsingen) <- paste("Context", rownames(munsingen))
#' colnames(munsingen) <- paste("Find", colnames(munsingen))
#' 
#' kappa.coef(munsingen)
#' 
#' @export kappa.coef
kappa.coef <- function(obj) {
     nu <- sum(obj)
     k <- ( concentration.col(obj) + concentration.col(t(obj)) ) / (2 * nu)
     return(k)
}




#' Evaluating Element Fit
#'
#' Performs a deviance test using a quadratic logistic regression on the rows and columns of a seriated incidence matrix, returning the \eqn{p} value for each row or column. Rows or columns which exhibit perfect separation in 0/1 values are assigned an NA value. 
#' 
#' @param obj A seriated binary matrix.
#' @returns A list of the \eqn{p} values for each row and column 
#' 
#' @examples 
#' # using the Münsingen data from the Münsingen-Rain
#' library(seriation)
#' munsingen <- data(Munsingen)
#' 
#' # row and column names should be unique
#' rownames(munsingen) <- paste("Context", rownames(munsingen))
#' colnames(munsingen) <- paste("Find", colnames(munsingen))
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
            fit <- glm(y ~ x + I(x^2), family = binomial(link = "logit"))
            })
            dev.res <- fit$deviance
            dev.nul <- fit$null.deviance
            df.res <- fit$df.residual
            df.nul <- fit$df.null

            pval <- 1-pchisq(dev.nul - dev.res, df.nul - df.res)
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
            fit <- glm(y ~ x + I(x^2), family = binomial(link = "logit"))
            })
            dev.res <- fit$deviance
            dev.nul <- fit$null.deviance
            df.res <- fit$df.residual
            df.nul <- fit$df.null

            pval <- 1-pchisq(dev.nul - dev.res, df.nul - df.res)
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
#' From the results of `ca.procrustes.curve()`, extrect two matrices containing the ranks of the rows and columns. The row/column elements are contained in the rows, and the strands are contained in the columns. NA values are entered where a given row/column element is missing from that strand.
#' 
#' @param strands A list of `strands`, which are data frames returned by `ca.procrustes.curve()`.
#' @param obj The intial incidence matrix.
#' 
#' @return A list of the following:.
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
#' Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca.procrustes.curve}}) to generate a new list of strands. The row or column element id is stored in an output list in case it should be added back into the strand. Column and row names are required to be unique. If the resulting strand lacks sufficient points, the entire strand is deleted in the output.
#' 
#' @param strands A list of `strands`, which are data frames returned by `ca.procrustes.curve()`.
#' @param obj The intial incidence matrix.
#' @param elements The name of one or more row or column ids to suppress.
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


