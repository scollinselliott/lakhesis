#' Correspondence Analysis with Procrustes Fitting
#'
#' Fit scores of correspondence analysi on an incidence matrix to those produced by reference matrix which contain an ideal seriation using a Procrustes method (on the reference matrix, see \code{\link[lakhesis]{im_ref}}).  Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @param symmetric Whether to use standard scores for both rows and columns. Default is \code{TRUE}. Setting this to \code{FALSE} will result in a asymmetric map of standard scores for rows and principal scores for columns.
#' @returns A \code{list} object of class \code{strand} containing the following:
#' * `ref` The Procrustes-fit coordinates of the scores of the reference seriation.
#' * `x` The coordinates of the row standard scores of the data.
#' * `y` The coordinates of the column principal scores of the data.
#' * `x_pr` The Procrustes-fit coordinates of the row standard scores of the data.
#' * `y_pr` The Procrustes-fit coordinates of the column column scores of the data.
#' 
#' @examples 
#' data("quattrofontanili")
#' s <- ca_procrustes(quattrofontanili)
#' # print(s)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
ca_procrustes <- function(obj, symmetric = TRUE) {
    UseMethod("ca_procrustes")
}
#'
#' @export
ca_procrustes.default <- function(obj, symmetric = TRUE) {
    stop(paste('ca_procrustes cannot take input of class "', class(obj), '" ', sep =''))
}

#' @rdname ca_procrustes
#' @export
ca_procrustes.matrix <- function(obj, symmetric = TRUE) {
    if (sum(is.na(match(rownames(obj), colnames(obj)))) != nrow(obj))  {
        stop("Row names and column names of input matrix be unique.")
    }
    results <- list()

    obj <- obj[sort(rownames(obj)), sort(colnames(obj))]
    transposed <- FALSE
    if (nrow(obj) > ncol(obj)) {
        obj <- t(obj)
        transposed <- TRUE
    }

    if (symmetric == TRUE) {
        ref <- im_ref(obj)

        G.r <- ca::ca(obj)$rowcoord[,1:2]
        G.c <- ca::ca(obj)$colcoord[,1:2]
        R.r <- ca::ca(ref)$rowcoord[,1:2]
        R.c <- ca::ca(ref)$colcoord[,1:2]
    } else {
        ref <- im_ref(obj)
        S <- diag(ca::ca(obj)$sv)
        Sref <- diag(ca::ca(ref)$sv)

        G.r <- ca::ca(obj)$rowcoord[,1:2]
        G.c <- ca::ca(obj)$colcoord %*% S
        G.c <- G.c[,1:2]
        R.r <- ca::ca(ref)$rowcoord[,1:2]
        R.c <- ca::ca(ref)$colcoord %*% Sref
        R.c <- R.c[,1:2]
    }

    #data
    G.rc <- G.r - t(replicate( nrow(G.r) , apply(G.r, 2, mean))) # mean center
    G.rc.rad.max <- max( sqrt(G.rc[,1]^2 + G.rc[,2]^2) ) # rescale
    x.r <- G.rc / G.rc.rad.max

    G.cc <- G.c - t(replicate( nrow(G.c) , apply(G.r, 2, mean))) # mean center with row data
    x.c <- G.cc / G.rc.rad.max  # rescale with row data

    #reference
    R.rc <- R.r - t(replicate( nrow(R.r) , apply(R.r, 2, mean))) # mean center
    R.rc.rad.max <- max( sqrt(R.rc[,1]^2 + R.rc[,2]^2) ) # rescale
    ref_r <- R.rc / R.rc.rad.max

    R.cc <- R.c - t(replicate( nrow(R.c) , apply(R.r, 2, mean))) # mean center with row data
    ref.c <- R.cc / R.rc.rad.max  # rescale with row data


    # procrustes rotation without landmark points

    # minimize rss of one-to-many distance from data to reference curve
    mid <- round(nrow(ref_r)/2)
    ref_mid_rad <- atan2(ref_r[mid,2], ref_r[mid,1])
    rss <- rss_rotation(x.r, ref_r, ref_mid_rad)

    # rotate data using row points
    idx <- which.min(rss)

    theta <- atan2(x.r[idx,2], x.r[idx,1]) - ref_mid_rad
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
    x.r.rot <- x.r %*% rotate
    x.c.rot <- x.c %*% rotate

    results[['ref']] <- ref_r
    if (transposed == FALSE) {
        results[['x']] <- G.r
        results[['y']] <- G.c
        results[['x_pr']] <- x.r.rot
        results[['y_pr']] <- x.c.rot
    } else {
        results[['y']] <- G.r
        results[['x']] <- G.c
        results[['y_pr']] <- x.r.rot
        results[['x_pr']] <- x.c.rot
    }

    class(results) <- c("procrustean", "list")

    return(results)
}

#' @rdname ca_procrustes
#' @export
ca_procrustes.incidence_matrix <- function(obj, symmetric = TRUE) {
    ca_procrustes.matrix(obj, symmetric)
}



#' Seriate Procrustes-Fit CA Scores
#' 
#' Obtain a ranking of row and column scores projected onto a reference curve of an ideal seriation (row and column scores are ranked separately). Scores of correspondence analysis have been fit to those produced by reference matrix contain an ideal seriation using a Procrustes method, projecting them. Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @param projection Which projection to use:
#'   * \code{ca1} The first axis for CA scores.
#'   * \code{ca2} The second axis for CA scores.
#'   * \code{procrustes1} The first axis of Procrustes-fit CA scores.
#'   * \code{procrustes2} The second axis of Procrustes-fit CA scores.
#'   * \code{curve} The reference curve of an ideal seriation, using Procrustes fitting (this is the default).
#' @param samples Number of samples to use for plotting points along polynomial curve. Default is \code{10^5}.
#' @return A \code{list} of class \code{strand} containing the following:.
#' * \code{$dat} A data frame with the following columns:
#'   * \code{Procrustes1, Procrustes2} The location of the point on the biplot after fitting.
#'   * \code{CurveIndex} The orthogonal projection of the point onto the reference curve, given as the index of the point sampled along \eqn{y = \beta_2 x^2 + \beta_0}. 
#'   * \code{Distance} The squared Euclidean distance of the point to the nearest point on the reference curve.
#'   * \code{Rank} The ranking of the row or column, a range of `1:nrow`` and `1:ncol``.
#'   * \code{Type} Either `row` or `col`.
#'   * \code{sel} Data frame column used in \code{shiny} app to indicate whether point is selected in biplot/curve projection.
#' * \code{$im_seriated} The seriated incidence matrix, of class \code{incidence_matrix}.
#' 
#' @examples
#' data("quattrofontanili")
#' s <- ca_procrustes_ser(quattrofontanili)
#' # print(s)
#' # summary(s)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
ca_procrustes_ser <- function(obj, projection = "curve", samples = 10^5, symmetric = TRUE) {
    UseMethod("ca_procrustes_ser")
}

#' @rdname ca_procrustes_ser
#' @export
ca_procrustes_ser.incidence_matrix <- function(obj, projection = "curve", samples = 10^5, symmetric = TRUE) {
    ca_procrustes_ser.matrix(obj, projection, samples, symmetric)
}

#' @rdname ca_procrustes_ser
#' @export
ca_procrustes_ser.matrix <- function(obj, projection = "curve", samples = 10^5, symmetric = TRUE) {
    if ((nrow(obj) > 2 ) & (ncol(obj) > 2) ) {
        if (!(projection %in% c("ca1", "ca2", "procrustes1", "procrustes2", "curve"))) {
            message("Invalid projection method.")
        } else {
            obj <- obj[sort(rownames(obj)), sort(colnames(obj))]
            obj_procca <- ca_procrustes(obj, symmetric = symmetric)
            ref_r <- obj_procca$ref
            x_r_rot <- obj_procca$x_pr
            x_c_rot <- obj_procca$y_pr
            x <- ref_r[,1]
            y <- ref_r[,2]

            # polynomial fitting
            xx <- seq(-1, 1, length.out = samples)
            fit <- stats::lm(y ~ x + I(x^2))
            yy <- stats::predict(fit, data.frame(x = xx))
            ref2_r <- as.matrix(cbind(xx,yy))

            x_ <- orth_proj_fit(x_r_rot, ref2_r)
            x_dat <- data.frame(index = x_[,2], x_fit = x_[,3], y_fit = x_[,4], dist = x_[,1])
            rownames(x_dat) <- rownames(obj)

            # col points
            y_ <- orth_proj_fit(x_c_rot, ref2_r)
            y_dat <- data.frame(index = y_[,2], x_fit = y_[,3], y_fit = y_[,4], dist = y_[,1])
            rownames(y_dat) <- colnames(obj)

            if (projection == "curve") {
                r_idx <- x_dat$index
                c_idx <- y_dat$index
            } else if (projection == "procrustes1") {
                r_idx <- x_r_rot[,1]
                c_idx <- x_c_rot[,1]
            } else if (projection == "procrustes2") {
                r_idx <- x_r_rot[,2]
                c_idx <- x_c_rot[,2]
            } else if (projection == "ca1") {
                r_idx <- obj_procca$x[,1]
                c_idx <- obj_procca$y[,1]
            } else if (projection == "ca2") {
                r_idx <- obj_procca$x[,2]
                c_idx <- obj_procca$y[,2]
            }

            r_rank <- rank(r_idx)
            c_rank <- rank(c_idx)

            ord_row1 <- data.frame(CA1 = obj_procca$x[,1], CA2 = obj_procca$x[,2],Procrustes1 = obj_procca$x_pr[,1], Procrustes2 = obj_procca$x_pr[,2], CurveIndex = r_idx, Distance = x_dat$dist, Rank = r_rank)
            rownames(ord_row1) <- rownames(x_dat)
            ord_row1$Type <- factor("row")
            ord_col1 <- data.frame(CA1 = obj_procca$y[,1], CA2 = obj_procca$y[,2],Procrustes1 = obj_procca$y_pr[,1], Procrustes2 = obj_procca$y_pr[,2], CurveIndex = c_idx, Distance = y_dat$dist, Rank = c_rank)
            rownames(ord_col1) <- rownames(y_dat)
            ord_col1$Type <- factor("col")
            ranking <- rbind(ord_row1,ord_col1)
            ranking$sel <- FALSE

            R <- rownames(x_dat)[order(r_rank)]
            C <- rownames(y_dat)[order(c_rank)]

            im_seriated <- obj[R, C]
            class(im_seriated) <- c("incidence_matrix", "matrix", "array")
            strand <- list(dat = ranking, im_seriated = im_seriated, method = projection)
            class(strand) <- c("strand", "list")
            return(strand)
        }
    } else {
        message("Insufficient rows and columns for correspondence analysis.")
    }
}



