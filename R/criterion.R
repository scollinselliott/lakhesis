#' Optimality Criterion: Squared Correlation
#' 
#' Treating each incidence of 1 in an element \eqn{(i,j)} of a seriated matrix as an \eqn{(x,y)} point, computes the squared correlation coefficient \insertCite{@see @mccormick_identification_1969, 147-148}{lakhesis}.
#' 
#' @param obj A seriated binary matrix.
#' @returns Spearman's rank correlation coefficient.
#' 
#' @examples
#' data("quattrofontanili")
#' cor_sq(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
cor_sq <- function(obj) {
    UseMethod("cor_sq")
}

#' @rdname cor_sq
#' @export 
cor_sq.matrix <- function(obj) {
    x <- c()
    y <- c()
    for (i in 1:nrow(obj)) {
        yi <- which(obj[i, ] == 1)
        xi <- rep(i, length(yi))
        x <- c(x, xi)
        y <- c(y, yi)
    }
    r_ <- ( stats::cor(x,y, method = "pearson")  )^2
    return(r_)
}

#' @rdname cor_sq
#' @export
cor_sq.incidence_matrix <- function(obj) {
    cor_sq.matrix(obj)
}



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

#' @rdname spearman_sq
#' @export 
spearman_sq.numeric <- function(r1, r2) {
    dat <-  stats::na.omit( data.frame(r1,r2) )
    r <- (stats::cor(dat)[1,2])^2 
    return(r)
}



#' Optimality Criterion: Kendall-Doran (Column) Concentration
#'
#' The Kendall-Doran measure of concentration \insertCite{kendall_statistical_1963,doran_computer_1971}{lakhesis}. In a seriated matrix, this function computes the total number cells between the first and last non-zero value, column by column.
#' 
#' @param obj A seriated binary matrix.
#' @returns The measure of concentration.
#' 
#' @examples 
#' data("quattrofontanili")
#' conc_c(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
conc_c <- function(obj) {
    UseMethod("conc_c")
}

#' @rdname conc_c
#' @export
conc_c.matrix <- function(obj) {
    conc <- numeric(ncol(obj))
    for (j in 1:ncol(obj)) {
        conc[j] <- max(which(obj[,j] != 0)) - min(which(obj[,j] != 0)) + 1
    }
    return(sum(conc))
}

#' @rdname conc_c
#' @export
conc_c.incidence_matrix <- function(obj) {
    conc_c.matrix(obj)
}


#' Optimality Criterion: Weighted Row-Column Concentration
#'
#' Extends the Kendall-Doran (column) measure of concentration (see \code{\link[lakhesis]{conc_c}}) to include rows and then weights the total measure by the total sum of values in the matrix.
#' 
#' @param obj A seriated binary matrix.
#' @returns The weighted row-column coefficient of concentration.
#' 
#' @examples 
#' data("quattrofontanili")
#' conc_wrc(quattrofontanili)
#' 
#' @export
conc_wrc <- function(obj) {
    UseMethod("conc_wrc")
}

#' @rdname conc_wrc
#' @export
conc_wrc.matrix <- function(obj) {
    nu <- sum(obj)
    k <- ( conc_c(obj) + conc_c(t(obj)) ) / (2 * nu)
    return(k)
}

#' @rdname conc_wrc
#' @export
conc_wrc.incidence_matrix <- function(obj) {
    conc_wrc.matrix(obj)
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


#' @rdname element_eval
#' @export
element_eval.matrix <- function(obj) {
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

#' @rdname element_eval
#' @export
element_eval.incidence_matrix <- function(obj) {
    element_eval.matrix(obj)
}



