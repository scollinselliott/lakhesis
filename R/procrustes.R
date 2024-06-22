#' Correspondence Analysis with Procrustes Fitting
#'
#' Fit scores of correspondence analysis on an incidence matrix to those produced by reference matrix which contain an ideal seriation using a Procrustes method (on the reference matrix, see \code{\link[lakhesis]{im_ref}}). Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @returns A \code{list} object of class \code{strand} containing the following:
#' * `ref` The Procrustes-fit coordinates of the scores of the reference seriation.
#' * `x` The Procrustes-fit coordinates of the row scores of the data.
#' * `y` The Procrustes-fit coordinates of the column scores of the data.
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
ca_procrustes <- function(obj) {
    UseMethod("ca_procrustes")
}
#'
#' @export
ca_procrustes.default <- function(obj) {
    stop(paste('ca_procrustes cannot take input of class "', class(obj), '" ', sep =''))
}

#' @rdname ca_procrustes
#' @export
ca_procrustes.matrix <- function(obj) {
    transposed <- FALSE
    if (nrow(obj) > ncol(obj)) {
        obj <- t(obj)
        transposed <- TRUE
    }

    ref <- im_ref(obj)

    G.r <- ca::ca(obj)$rowcoord[,1:2]
    G.c <- ca::ca(obj)$colcoord[,1:2]
    R.r <- ca::ca(ref)$rowcoord[,1:2]
    R.c <- ca::ca(ref)$colcoord[,1:2]

    #data
    G.rc <- G.r - t(replicate( nrow(G.r) , apply(G.r, 2, mean))) # mean center
    G.rc.rad.max <- max( sqrt(G.rc[,1]^2 + G.rc[,2]^2) ) # rescale
    x.r <- G.rc / G.rc.rad.max

    G.cc <- G.c - t(replicate( nrow(G.c) , apply(G.r, 2, mean))) # mean center with row data
    x.c <- G.cc / G.rc.rad.max  # rescale with row data

    #reference
    R.rc <- R.r - t(replicate( nrow(R.r) , apply(R.r, 2, mean))) # mean center
    R.rc.rad.max <- max( sqrt(R.rc[,1]^2 + R.rc[,2]^2) ) # rescale
    ref.r <- R.rc / R.rc.rad.max

    R.cc <- R.c - t(replicate( nrow(R.c) , apply(R.r, 2, mean))) # mean center with row data
    ref.c <- R.cc / R.rc.rad.max  # rescale with row data

    # procrustes rotation without landmark points

    # minimize rss of one-to-many distance from data to reference curve

    mid <- round(nrow(ref.r)/2)
    ref.mid.rad <- atan2(ref.r[mid,2], ref.r[mid,1])

    rss <- c()
    for (j in 1:nrow(x.r)) {

        theta <- atan2(x.r[j,2], x.r[j,1]) - ref.mid.rad
        rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
        x.r.rot <- x.r %*% rotate
        #x.c.rot <- x.c %*% rotate

        rss1 <- 0
        for (k in 1:nrow(x.r)) {
        rss1 <- rss1 + min( rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(x.r.rot)))  - ref.r)^2) )
        }

        rss <- c(rss, rss1)
    }






    # rotate data using row points
    idx <- which.min(rss)

    theta <- atan2(x.r[idx,2], x.r[idx,1]) - ref.mid.rad
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
    x.r.rot <- x.r %*% rotate
    x.c.rot <- x.c %*% rotate

    results <- list()

    results[['ref']] <- ref.r
    if (transposed == FALSE) {
        results[['x']] <- x.r.rot
        results[['y']] <- x.c.rot
    }
    if (transposed == TRUE) {
        results[['y']] <- x.r.rot
        results[['x']] <- x.c.rot
    }

    class(results) <- c("procrustean", "list")

    return(results)
}

#' @rdname ca_procrustes
#' @export
ca_procrustes.incidence_matrix <- function(obj) {
    ca_procrustes.matrix(obj)
}



#' Seriate Using Reference Curve for Procrustes-Fit CA Scores
#' 
#' Obtain a ranking of row and column scores projected onto a reference curve of an ideal seriation (row and column scores are ranked separately). Scores of correspondence analysis have been fit to those produced by reference matrix contain an ideal seriation using a Procrustes method, projecting them. Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
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
ca_procrustes_ser <- function(obj, samples = 10^5) {
    UseMethod("ca_procrustes_ser")
}

#' @rdname ca_procrustes_ser
#' @export
ca_procrustes_ser.incidence_matrix <- function(obj, samples = 10^5) {
    ca_procrustes_ser.matrix(obj, samples = 10^5)
}

#' @rdname ca_procrustes_ser
#' @export
ca_procrustes_ser.matrix <- function(obj, samples = 10^5) {
    if ((nrow(obj) > 2 ) & (ncol(obj) > 2) ) {

        obj.procca <- ca_procrustes(obj)
        ref.r <- obj.procca$ref
        x.r.rot <- obj.procca$x
        x.c.rot <- obj.procca$y
        x <- ref.r[,1]
        y <- ref.r[,2]

        # polynomial fitting
        xx <- seq(-1, 1, length.out = samples)
        fit <- stats::lm(y ~ x + I(x^2))
        yy <- stats::predict(fit, data.frame(x = xx))
        ref2.r <- as.matrix(cbind(xx,yy))

        # row points
        index <- c()
        rss <- c()
        x.prc <- c()
        y.prc <- c()
        for (k in 1:nrow(x.r.rot)) {
            rss1 <- min( rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(ref2.r)))  - ref2.r)^2) )
            idx <- which.min(rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(ref2.r)))  - ref2.r)^2))
            index <- c(index, idx)
            x.prc <- c(x.prc, ref2.r[idx,1])
            y.prc <- c(y.prc, ref2.r[idx,2])
            rss <- c(rss, rss1)
        }
        x.dat <- data.frame(index = index, x.prc = x.prc, y.prc = y.prc, dist = rss)
        rownames(x.dat) <- rownames(obj)

        # col points
        index <- c()
        rss <- c()
        x.prc <- c()
        y.prc <- c()
        for (k in 1:nrow(x.c.rot)) {
            rss1 <- min( rowSums(( t(matrix(x.c.rot[k,],  nrow = 2, ncol = nrow(ref2.r)))  - ref2.r)^2) )
            idx <- which.min(rowSums(( t(matrix(x.c.rot[k,],  nrow = 2, ncol = nrow(ref2.r)))  - ref2.r)^2))
            index <- c(index, idx)
            x.prc <- c(x.prc, ref2.r[idx,1])
            y.prc <- c(y.prc, ref2.r[idx,2])
            rss <- c(rss, rss1)
        }
        y.dat <- data.frame(index = index, x.prc = x.prc, y.prc = y.prc, dist = rss)
        rownames(y.dat) <- colnames(obj)

        r.idx <- x.dat$index
        c.idx <- y.dat$index

        r.rank <- rank(r.idx)
        c.rank <- rank(c.idx)
        ord.row1 <- data.frame(Procrustes1 = obj.procca$x[,1], Procrustes2 = obj.procca$x[,2], CurveIndex = r.idx, Distance = x.dat$dist, Rank = r.rank)
        rownames(ord.row1) <- rownames(x.dat)
        ord.row1$Type <- factor("row")
        ord.col1 <- data.frame(Procrustes1 = obj.procca$y[,1], Procrustes2 = obj.procca$y[,2], CurveIndex = c.idx, Distance = y.dat$dist, Rank = c.rank)
        rownames(ord.col1) <- rownames(y.dat)
        ord.col1$Type <- factor("col")
        proc.ranking <- rbind(ord.row1,ord.col1)
        proc.ranking$sel <- FALSE
        class(proc.ranking) <- c("data.frame")
        R <- rownames(x.dat)[order(r.rank)]
        C <- rownames(y.dat)[order(c.rank)]
        im_seriated <- obj[R, C]
        class(im_seriated) <- c("incidence_matrix", "matrix", "array")
        strand <- list(dat = proc.ranking, im_seriated = im_seriated)
        class(strand) <- c("strand", "list")
        return(strand)
    } else {
        message("Insufficient rows and columns for correspondence analysis.")
    }
}



