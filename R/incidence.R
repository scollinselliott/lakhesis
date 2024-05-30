#' Read csv File to Incidence Matrix
#'
#' Wrapper around the \code{\link[readr]{read_csv}} function from the \code{\link[readr]{readr}} package \insertCite{wickham_readr_2024}{lakhesis}. Read a CSV file in which the first column represents row elements and the second column represents column elements, and convert it into an incidence matrix.
#'  
#' @param filename The filename to uploaded (must be in CSV format).
#' @param header If the CSV file contains a header. Default is `FALSE`.
#' @param characterencoding File encoding as used by \code{\link[readr::locale]{locale}}. Default is `"iso-8859-1"` to handle special characters. 
#' @param remove.hapax Remove any row or column which has a sum of 1 (i.e., is only attested once), since they do not directly contribute to the result of the seriation. Default is `TRUE`.
#' @returns An incidence matrix of binary values (0 = column element is lacking for row element; 1 = column element is present for row element).
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
im.csv.read <- function(filename, header = FALSE, characterencoding = "iso-8859-1", remove.hapax = TRUE) {
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
    return(obj)
}



#' Create Reference Matrix
#'
#' Create an ideal reference matrix of well-seriated values of the same size as the input matrix. 
#'
#' @param obj A matrix of size n x k.
#' @returns A matrix with 1s along the diagonal. If n > k, 1s are placed from (i,i) to (i,i+k-n), with 0 in all other cells.
#' @examples
#' im.ref(matrix(NA, 5, 5))
#' im.ref(matrix(1, 7, 12))
#' 
#' @export
im.ref <- function(obj) {
    if (is.matrix(obj)) {
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
    return(ref)
    } else {
        print("Input into im.ref() is not a matrix.")
    }
}



#' Convert Incidence Matrix to Pairs (Long Format)
#' 
#' Take an incidence matrix and convert it to a data frame of two columns, where the first column represents the row elements of the incidence matrix and the second column represents the column elements of the incidence matrix. Each row pair represents the incidence (or occurrence) of that row and column element together.
#'
#' @param obj An incidence matrix. 
#' @returns A data frame of two columns (row and column of the incidence matrix), in which row of the data frame represents a pair of an incidence.
#' 
#' @export
im.long <- function(obj) {
    dat <- data.frame(Row = c(), Col = c())
    for (i in 1:nrow(obj)) {
        row <- obj[i, ]
        rows <- rownames(obj)[i]
        cols <- names(row[row == 1])
        tmp <- matrix(c(rep(rows, length(cols)) , cols), nrow = length(cols), ncol = 2)
        dat <- rbind(dat, tmp)
    }
    colnames(dat) <- c("Row", "Col")
    return(dat)
}



#' Merge Two Incidence Matrices
#'
#' Create an ideal reference matrix of well-seriated values of the same size as the input matrix. 
#'
#' @param obj1,obj2 Two incidence matrices, which may contain overlapping row or column values.
#' @returns A single incidence matrix.
#' 
#' @export
im.merge <- function(obj1, obj2) {
    dat <- rbind(im.long(obj1), im.long(obj2))

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
    return(obj)
}


