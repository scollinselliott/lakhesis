#' Strand Extract
#' 
#' From a \code{list} of strands produced by \code{\link[lakhesis]{ca_procrustes_ser}}, extract two matrices containing the ranks of the rows and columns. The row/column elements are contained in the rows, and the strands are contained in the columns. NA values are entered where a given row/column element is missing from that strand.
#' 
#' @param strands A \code{list} of class \code{strands}.
#' 
#' @return A list of two matrices: 
#' * \code{Row} A matrix of the ranks of the row elements.
#' * \code{Col} A matrix of the ranks of the column elements.
#' 
#' @examples
#' data("quattrofontanili")
#' data("qf_strands")
#' strand_extract(qf_strands)
#' 
#' @export
strand_extract <- function(strands) {
    UseMethod("strand_extract")
}  

#' @export
strand_extract.strands <- function(strands) {
    obj <- strands[[1]]$im_seriated
    for (i in 2:length(strands)) {
        obj <- im_merge(obj, strands[[i]]$im_seriated)
    }
    if (sum(is.na(match(rownames(obj), colnames(obj)))) != nrow(obj))  {
        stop("Strands improperly produced. Row names and column names of input matrix be unique.")
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
#' Given a \code{list} of \code{strands} produced by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca_procrustes_ser}}), remove one or more row or column elements, re-seriating each strand. This generates a new list of strands that exclude the specified elements. If a resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
#' 
#' @param strands A \code{list} of class \code{strands}.
#' @param elements A vector of one or more row or column ids to suppress.
#' 
#' @return A list of the strands.
#' 
#' @examples
#' data("qf_strands")
#' strand_suppress(qf_strands, "QF II 15-16")
#' 
#' strand_suppress(qf_strands, c("QF II 15-16", "I", "XIV"))
#' 
#' @export
strand_suppress <- function(strands, elements) {
    UseMethod("strand_suppress")
}

#' @rdname strand_suppress
#' @export
strand_suppress.strands <- function(strands, elements) {
    strand_suppress.default(strands, elements)
}

#' @rdname strand_suppress
#' @export
strand_suppress.default <- function(strands, elements) {
    new <- list()
    obj <- strands[[1]]$im_seriated
    for (s in 2:length(strands)) {
        obj <- im_merge(obj, strands[[s]]$im_seriated)
    }
    j <- 1
    for (s in 1:length(strands)) {
        strand <- strands[[s]]$dat
        rows <- rownames(strand)[ !(rownames(strand) %in% elements) & (strand$Type == "row") ]
        cols <- rownames(strand)[ !(rownames(strand) %in% elements) & (strand$Type == "col") ]
        obj_copy <- obj
        obj_copy <- obj_copy[rows, ]
        obj_copy <- obj_copy[ ,cols]
        obj_copy <- obj_copy[rowSums(obj_copy) != 0 , ]
        obj_copy <- obj_copy[ , colSums(obj_copy) != 0]
        newstrand <- NULL
        suppressWarnings({
            newstrand <- ca_procrustes_ser(obj_copy)
        })
        if (!is.null(newstrand)) { 
           new[[j]] <- newstrand
           j <- j + 1
        }
    }
    class(new) <- c("strands", "list")
    return(new)
}




#' Add Strand to List of Strands
#' 
#' Given a list of strands, remove a row or column element and re-run seriation by correspondence analysis with Procrustes fitting (\code{\link[lakhesis]{ca_procrustes_ser}}) to generate a new list of strands that exclude the specified elements. If the resulting strand lacks sufficient points to perform correspondence analysis, that strand is deleted in the output.
#' 
#' @param strand An object of class \code{strand} returned by \code{\link[lakhesis]{ca_procrustes_ser}}.
#' @param strands A \code{list} of strands.
#'
#' @return A \code{list} of class \code{strands}.
#' 
#' @export 
strand_add <- function(strand, strands) {
    UseMethod("strand_add")
}

#' @rdname strand_add
#' @export
strand_add.strand <- function(strand, strands) {
    strands[[length(strands) + 1]] <- strand
    class(strands) <- c("strands", "list")
    return(strands)
}


#' Create List of Strands
#' 
#' Given one or more individual strand objects, create a single \code{list} of class \code{strands}.
#' 
#' @param strands A \code{list} of strands.
#'
#' @return A \code{list} of class \code{strands}.
#' 
#' @export 
strands_create <- function(strands) {
    UseMethod("strands_create")
}

#' @rdname strands_create
#' @export
strands_create.list <- function(strands) {
    class(strands) <- c("strands", "list")
    return(strands)
}



#' Create Strand Object from Seriated Incidence Matrix
#' 
#' Given a seriated incidence matrix with unique row and column names, create a \code{strand} object.
#' 
#' @param obj A \code{list} of strands.
#' @param method The method used to create the strand (optional).
#'
#' @return A \code{list} of class \code{strands}.
#' 
#' @export 
strand_create <- function(obj, method = NULL) {
    UseMethod("strand_create")
}

#' @rdname strand_create
#' @export
strand_create.matrix <- function(obj, method = NULL) {

    ord_row <- data.frame(Rank = 1:nrow(obj))
    ord_row$Type <- factor("row")
    rownames(ord_row) <- rownames(obj)
    ord_col <- data.frame(Rank = 1:ncol(obj))
    ord_col$Type <- factor("col")
    rownames(ord_col) <- colnames(obj)

    ranking <- rbind(ord_row,ord_col)
    ranking$sel <- FALSE

    class(obj) <- c("incidence_matrix", "matrix", "array")

    res <- list(dat = ranking, im_seriated = obj, method = method)
    class(res) <- c("strand", "list")
    return(res)
}

#' @rdname strand_create
#' @export
strand_create.incidence_matrix <- function(obj, method = NULL) {
    strand_create.matrix(obj, method)
}
