#' Correspondence Analysis with Procrustes Fitting
#'
#' Fit scores of correspondence analysis on an incidence matrix to those produced by reference matrix which contain an ideal seriation using a Procrustes method (on the reference matrix, see \code{\link[lakhesis]{im.ref}}). Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @returns A list of the following:.
#' 
#' * `ref` The Procrustes-fit coordinates of the scores of the reference seriation.
#' * `x` The Procrustes-fit coordinates of the row scores of the data.
#' * `x.dat` A data frame containing the following information related to the fit of the row score after Procrustes fitting: `index`, the row name, `match`, the reference point nearest to the row score, and `dist`, the Euclidean distance between the row score and reference score point.
#' * `y` The Procrustes-fit coordinates of the column scores of the data.
#' * `y.dat` A data frame containing the same information as `x.dat`, but with respect to the column scores.
#' 
#' @examples 
#' # Quattro Fontanili
#' data(quattrofontanili)
#' ca.procrustes(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
ca.procrustes <- function(obj) {
    transposed = FALSE
    if (nrow(obj) > ncol(obj)) {
        obj <- t(obj)
        transposed = TRUE
    }

    ref <- im.ref(obj)

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

    argmin <- c()
    rss <- c()
    for (k in 1:nrow(x.r)) {
        rss1 <- min( rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(ref.r)))  - ref.r)^2) )
        idx <- which.min(rowSums(( t(matrix(x.r.rot[k,],  nrow = 2, ncol = nrow(ref.r)))  - ref.r)^2))

        argmin <- c(argmin, idx)
        rss <- c(rss, rss1)
    }

    dat.x <- data.frame(index = rownames(x.r.rot), match = argmin, dist = rss)

    argmin.c <- c()
    rss.c <- c()
    for (k in 1:nrow(x.c)) {
        rss1 <- min( rowSums(( t(matrix(x.c.rot[k,],  nrow = 2, ncol = nrow(ref.r)))  - ref.r)^2) )
        idx <- which.min(rowSums(( t(matrix(x.c.rot[k,],  nrow = 2, ncol = nrow(ref.r)))  - ref.r)^2))

        argmin.c <- c(argmin.c, idx)
        rss.c <- c(rss.c, rss1)
    }

    dat.y <- data.frame(index = rownames(x.c.rot), match = argmin.c, dist = rss.c)

    results <- list()

    results[['ref']] <- ref.r
    if (transposed == FALSE) {
    results[['x']] <- x.r.rot
    results[['y']] <- x.c.rot
    results[['x.dat']] <- dat.x
    results[['y.dat']] <- dat.y
    }
    if (transposed == TRUE) {
    results[['y']] <- x.r.rot
    results[['x']] <- x.c.rot
    results[['y.dat']] <- dat.x
    results[['x.dat']] <- dat.y
    }

    return(results)
}



#' Projection onto Reference Curve
#'
#' Performs a polynomial regression on the row reference scores and orthogonally projects data points on to the reference curve. Sampling can be increased to refine ranking and avoid ties, but default is largely sufficient. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @param resolution Number of samples to use for plotting points along polynomial curve (default is 10000).
#' @returns A list of the following:.
#' * `ref`` The Procrustes-fit coordinates of the scores of the reference seriation.
#' * `x` The Procrustes-fit coordinates of the row scores of the data.
#' * `x.dat` A data frame containing the following information related to the fit of the row score after Procrustes fitting: `index`, the row name, `match`, the reference point nearest to the row score, and `dist`, the Euclidean distance between the row score and reference score point.
#' * `y` The Procrustes-fit coordinates of the column scores of the data.
#' * `y.dat` A data frame containing the same information as `x.dat`, but with respect to the column scores.
#' 
#' @examples
#' # Quattro Fontanili
#' data(quattrofontanili)
#' ca.procrustes.poly(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
ca.procrustes.poly <- function(obj, resolution = 10000) {
    obj.procca <- ca.procrustes(obj)    
    ref.r <- obj.procca$ref
    x.r.rot <- obj.procca$x
    x.c.rot <- obj.procca$y
    x <- ref.r[,1]
    y <- ref.r[,2]

    # polynomial fitting
    xx <- seq(-1, 1, length.out = resolution)
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

    results <- list()
    results$proc <- obj.procca
    results$fit <- fit
    results$rank.x <- x.dat
    results$rank.y <- y.dat
    return(results)
}



#' Seriate Using Reference Curve
#' 
#' Obtain a ranking of row and column scores projected onto a reference curve of an ideal seriation (row and column scores are ranked separately). Scores of correspondence analysis have been fit to those produced by reference matrix contain an ideal seriation using a Procrustes method, projecting them. Rotation is determined by minimizing Euclidean distance from each row score to the nearest reference row score. Correspondence analysis is performed using the \code{\link[ca]{ca}} package \insertCite{nenadic_correspondence_2007}{lakhesis}.
#'
#' @param obj An incidence matrix of size n x k.
#' @param resolution Number of samples to use for plotting points along polynomial curve (default is 10000).
#' @return A data frame of the following:.
#' * `Procrustes1,Procrustes2` The location of the point on the biplot after fitting.
#' * `CurveIndex` The orthogonal projection of the point onto the reference curve, given as the index of the point sampled along \eqn{y = \beta_2 x^2 + \beta_0}. 
#' * `Distance` The squared Euclidean distance of the point to the nearest point on the reference curve.
#' * `Rank` The ranking of the row or column, a range of `1:nrow`` and `1:ncol``.
#' * `Type` Either `row` or `col`.
#' * `sel` Data frame column used in `shiny` app to indicate whether point is selected in biplot/curve projection.
#'
#' @examples
#' # Quattro Fontanili
#' data(quattrofontanili)
#' ca.procrustes.curve(quattrofontanili)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
ca.procrustes.curve <- function(obj, resolution = 10000) {
    if ((nrow(obj) > 2 ) & (ncol(obj) > 2) ) {
        mat.proc.poly <- ca.procrustes.poly(obj, resolution)
        mat.procca <- mat.proc.poly$proc
        r.idx <- mat.proc.poly$rank.x$index
        c.idx <- mat.proc.poly$rank.y$index
        r.rank <- rank(r.idx)
        c.rank <- rank(c.idx)
        ord.row1 <- data.frame(Procrustes1 = mat.procca$x[,1], Procrustes2 = mat.procca$x[,2], CurveIndex = r.idx, Distance = mat.proc.poly$rank.x$dist, Rank = r.rank)
        rownames(ord.row1) <- rownames(mat.proc.poly$rank.x)
        ord.row1$Type <- factor("row")
        ord.col1 <- data.frame(Procrustes1 = mat.procca$y[,1], Procrustes2 = mat.procca$y[,2], CurveIndex = c.idx, Distance = mat.proc.poly$rank.y$dist, Rank = c.rank)
        rownames(ord.col1) <- rownames(mat.proc.poly$rank.y)
        ord.col1$Type <- factor("col")
        proc.ranking <- rbind(ord.row1,ord.col1)
        proc.ranking$sel <- FALSE
        return(proc.ranking)
    } else {
        print("Insufficient points for correspondence analysis.")
    }
}
