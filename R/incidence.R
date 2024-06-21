#' Read csv File to Incidence Matrix
#'
#' Wrapper around the \code{\link[readr]{read_csv}} function from the \code{readr} package \insertCite{wickham_readr_2024}{lakhesis}. Read a \code{.csv} file in which the first column represents row elements and the second column represents column elements, and convert it into an incidence matrix.
#'  
#' @param filename The filename to uploaded (must be in \code{.csv} format).
#' @param header If the \code{.csv} file contains a header. Default is `FALSE`.
#' @param characterencoding File encoding as used by \code{\link[readr]{locale}}. Default is `"iso-8859-1"` to handle special characters. 
#' @param remove.hapax Remove any row or column which has a sum of 1 (i.e., is only attested once), since they do not directly contribute to the result of the seriation. Default is `FALSE`.
#' @returns A matrix of binary values (0 = row/column occurrence is absence; 1 = row/column occurrence is present). 
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
im_read_csv <- function(filename, header = FALSE, characterencoding = "iso-8859-1", remove.hapax = FALSE) {
    dat <- readr::read_csv(filename, col_names = header, show_col_types = FALSE, locale = readr::locale(encoding = characterencoding))
    dat <- data.frame(dat)

    finds <- levels(factor(dat[,2]))
    contexts <- levels(factor(dat[,1]))
    obj <- matrix(0, nrow = length(contexts), ncol = length(finds))
    rownames(obj) <- contexts
    colnames(obj) <- finds

    for (k in 1:nrow(dat)) {
        j <- which(finds == dat[k, 2]) 
        i <- which(contexts == dat[k, 1])
        obj[i, j] <- 1
    }

    if (remove.hapax == TRUE) {
        obj <- obj[, colSums(obj) > 1]
        obj <- obj[rowSums(obj) > 1, ]
    }

    class(obj) <- c("incidence_matrix", "matrix", "array")
    return(obj)
}



#' Create Reference Matrix
#'
#' Create an ideal reference matrix of well-seriated values of the same size as the input matrix. 
#'
#' @param obj A matrix of size \eqn{n \times k}.
#' @returns A matrix of size \eqn{n \times k} with 1s along the diagonal. If \eqn{n > k}, 1s are placed from cell \eqn{(i,i)} to \eqn{(i,i+k-n)}, with 0 in all other cells.
#' @examples
#' im_ref(matrix(NA, 5, 5))
#' im_ref(matrix(1, 7, 12))
#' 
#' @export
im_ref <- function(obj) {
    UseMethod("im_ref")
}


#' @export
im_ref.matrix <- function(obj) { 
        ref <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))
        if (nrow(ref) > ncol(ref)) {
            r <- nrow(ref) - ncol(ref) 
            for (j in 1:ncol(ref)) {
                end <- j+r
                ref[j:end,j] <- 1
            } 
        } else if (nrow(ref) == ncol(ref)) {
            diag(ref) <- 1
            r <- nrow(ref) - ncol(ref) 
            for (j in 1:(ncol(ref)-1)) {
                end <- j+1
                ref[j:end,j] <- 1
            } 
        } else {
            r <- ncol(ref) - nrow(ref) 
            for (i in 1:nrow(ref)) {
                end <- i+r
                ref[i, i:end] <- 1
            }
        }
        if (is.null(colnames(ref))) {
            rownames(ref) <- paste("R", 1:nrow(ref), sep = '')
            colnames(ref) <- paste("C", 1:ncol(ref), sep = '')
        }
    class(ref) <- c("incidence_matrix", "matrix", "array")
    return(ref)
}


#' @export
im_ref.incidence_matrix <- function(obj) { 
    im_ref.matrix(obj)
}


#' Convert Incidence Matrix to Pairs (Long Format)
#' 
#' Take an incidence matrix and convert it to a data frame of two columns, where the first column represents the row elements of the incidence matrix and the second column represents the column elements of the incidence matrix. Each row pair represents the incidence (or occurrence) of that row and column element together.
#'
#' @param obj An incidence matrix. 
#' @returns A data frame of two columns (row and column of the incidence matrix), in which row of the data frame represents a pair of an 
#' @examples 
#' data(quattrofontanili)
#' qf <- im_long(quattrofontanili)
#' 
#' # to export for uploading into the Lakhesis Calculator, use write.table() to 
#' # remove both row and column names:
#' 
#' # write.table(qf, file = 'qf.csv', row.names = FALSE, col.names = FALSE, sep = ",")
#' 
#' @export
im_long <- function(obj) {
    UseMethod("im_long")
}


#' @export
im_long.matrix <- function(obj) { 
    dat <- data.frame(Row = c(), Col = c())
    for (i in 1:nrow(obj)) {
        row <- obj[i, ]
        rows <- rownames(obj)[i]
        cols <- names(row[row == 1])
        tmp <- matrix(c(rep(rows, length(cols)) , cols), nrow = length(cols), ncol = 2)
        dat <- rbind(dat, tmp)
    }
    colnames(dat) <- c("Row", "Col")
    class(dat) <- "data.frame"
    return(dat)
}


#' @export
im_long.incidence_matrix <- function(obj) { 
    im_long.matrix(obj)
}

#' Merge Two Incidence Matrices
#'
#' From two incidience matrices, create a single incidence matrix. Matrices may contain same row or column elements. 
#'
#' @param obj1,obj2 Two incidence matrices of any size.
#' @returns A single incidence matrix.
#' @examples 
#' data(quattrofontanili)
#' qf1 <- quattrofontanili[1:20, 1:40]
#' qf1 <- qf1[rowSums(qf1) != 0, colSums(qf1) != 0]
#' 
#' qf2 <- quattrofontanili[30:50, 20:60]
#' qf2 <- qf2[rowSums(qf2) != 0, colSums(qf2) != 0]
#' 
#' im_merge(qf1, qf2)
#' 
#' @export
im_merge <- function(obj1, obj2) {
    UseMethod("im_merge")
}

#' @export
im_merge.matrix <- function(obj1, obj2) {
    dat <- rbind(im_long(obj1), im_long(obj2))

    finds <- levels(factor(dat[,2]))
    contexts <- levels(factor(dat[,1]))
    obj <- matrix(0, nrow = length(contexts), ncol = length(finds))
    rownames(obj) <- contexts
    colnames(obj) <- finds

    for (k in 1:nrow(dat)) {
        j <- which(finds == dat[k, 2]) 
        i <- which(contexts == dat[k, 1])
        obj[i, j] <- 1
    }
    class(obj) <- c("incidence_matrix", "matrix", "array")
    return(obj)
}

#' @export
im_merge.incidence_matrix <- function(obj1, obj2) {
    im_merge.matrix(obj1, obj2)
}

