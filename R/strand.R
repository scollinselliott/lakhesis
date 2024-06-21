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
#' spearman_sq(x, y)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
spearman_sq <- function(r1, r2) {
    UseMethod("spearman_sq")
}

#' @export 
spearman_sq.numeric <- function(r1, r2) {
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
#' conc_col(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
conc_col <- function(obj) {
    UseMethod("conc_col")
}

#' @export
conc_col.matrix <- function(obj) {
    conc <- numeric(ncol(obj))
    for (j in 1:ncol(obj)) {
        conc[j] <- max(which(obj[,j] != 0)) - min(which(obj[,j] != 0)) + 1
    }
    return(sum(conc))
}

#' @export
conc_col.incidence_matrix <- function(obj) {
    conc_col.matrix(obj)
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
#' conc_kappa(quattrofontanili)
#' 
#' @export
conc_kappa <- function(obj) {
    UseMethod("conc_kappa")
}

#' @export
conc_kappa.matrix <- function(obj) {
    nu <- sum(obj)
    k <- ( conc_col(obj) + conc_col(t(obj)) ) / (2 * nu)
    return(k)
}

#' @export
conc_kappa.incidence_matrix <- function(obj) {
    conc_kappa.incidence_matrix(obj)
}


#' Evaluating Element Fit
#'
#' Performs a goodness-of-fit test on individual row and column elements using deviance, using a quadratic-logistic model to fit row and column occurrences. In the case of perfect separation of 0/1 values, an `NA` value is assigned. Results are reported as \eqn{p} values for each row and column.
#' 
#' @param obj A seriated binary matrix.
#' @returns A \code{list} containing results in data frames for row and column elements:
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
#' element_eval(quattrofontanili)
#' 
#' @export
element_eval <- function(obj) {
    UseMethod("element_eval")
}


#' @export
element_eval.incidence_matrix <- function(obj) {
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
    class(results) <- c("list")
    return(results)
}







#' Strand Extract
#' 
#' From a `list` of strands produced by \code{\link[lakhesis]{ca.procrustes.curve}}, extract two matrices containing the ranks of the rows and columns. The row/column elements are contained in the rows, and the strands are contained in the columns. NA values are entered where a given row/column element is missing from that strand.
#' 
#' @param strands A \code{list} of class \code{strands}.
#' 
#' @return A list of two matrices: 
#' * `Row` A matrix of the ranks of the row elements.
#' * `Col` A matrix of the ranks of the column elements.
#' 
#' @examples
#' data("quattrofontanili")
#' data("qfStrands")
#' strand_extract(qfStrands)
#' 
#' @export
strand_extract <- function(strands, ...) {
    UseMethod("strand_extract")
}  

#' @export
strand_extract.strands <- function(strands) {
    obj <- strands[[1]]$im_seriated
    for (i in 2:length(strands)) {
        obj <- im_merge(obj, strands[[i]]$im_seriated)
    }

    rowranks <- matrix(NA, nrow = nrow(obj),  ncol = length(strands)) 
    rownames(rowranks) <- rownames(obj)

    colranks <- matrix(NA, nrow = ncol(obj),  ncol = length(strands)) 
    rownames(colranks) <- colnames(obj)

    for (i in 1:length(strands)) {
        strand <- strands[[i]]$dat

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
#' Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca.procrustes.curve}}) to generate a new list of strands that exclude the specified elements. If the resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
#' 
#' @param strands A \code{list} of class \code{strands}.
#' @param elements A vector of one or more row or column ids to suppress.
#' 
#' @return A list of the strands.
#' 
#' @examples
#' data("quattrofontanili")
#' data("qfStrands")
#' strand_suppress(qfStrands, "QF II 15-16")
#' 
#' strand_suppress(qfStrands, c("QF II 15-16", "I", "XIV"))
#' 
#' @export
strand_suppress <- function(strands, ...) {
    UseMethod("strand_suppress")
}


#' @export
strand_suppress.strands <- function(strands, ...) {
    strand_suppress.default(strands, ...)
}

#' @rdname strand_suppress
#' @method strand_suppress default
#' @export
strand_suppress.default <- function(strands, elements) {
    new <- list()
        obj <- strands[[1]]$im_seriated
        for (i in 2:length(strands)) {
            obj <- im_merge(obj, strands[[i]]$im_seriated)
        }
    for (i in 1:length(strands)) {
        strand <- strands[[i]]$dat
        rows <- rownames(strand)[ !(rownames(strand) %in% elements) & (strand$Type == "row") ]
        cols <- rownames(strand)[ !(rownames(strand) %in% elements) & (strand$Type == "col") ]
        obj.copy <- obj
        obj.copy <- obj.copy[rows, ]
        obj.copy <- obj.copy[ ,cols]
        obj.copy <- obj.copy[rowSums(obj.copy) !=0 , ]
        obj.copy <- obj.copy[ , colSums(obj.copy) !=0]
        newstrand <- NULL
        suppressWarnings({
            newstrand <- ca_procrustes_ser(obj.copy)
        })
        if (!is.null(newstrand)) { 
           new[[i]] <- newstrand
        }
    }
    class(new) <- c("strands", "list")
    return(new)
}




#' Add Strand to List of Strands
#' 
#' Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca.procrustes.curve}}) to generate a new list of strands that exclude the specified elements. If the resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
#' 
#' @param strand An object of class \code{strand} returned by \code{\link[lakhesis]{ca.procrustes.curve}}.
#' @param strands A \code{list} of strands.
#'
#' @return A \code{list} of class \code{strands}.
#' 
#' @export 
strand_add <- function(strand, ...) {
    UseMethod("strand_add")
}


#' @export
strand_add.strand <- function(strand, strands) {
    strands[[length(strands) + 1]] <- strand
    class(strands) <- c("strands, list")
    return(strands)
}


