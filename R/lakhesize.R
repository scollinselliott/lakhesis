#' Lakhesize
#'
#' This function returns the exploratory row and column consensus seriation for a `list` object of the `strands` class, containing their rankings, the results of their PCA, and coefficients of association and concentration.
#' 
#' Consensus seriation is achieved by iterative, multi-step linear regression using random simulation. On one iteration, strands are chosen at random, omitting incomplete or missing pairs, using PCA to determine the best-fitting line for their rankings. Both strands' rankings are then regressed onto that line to determine missing values, and then re-ranked, repeating until all strands have been regressed. PCA of the simulated rankings is then used to determine the final sequence of the row and column elements.
#' 
#' Two \code{method} options are available:
#' 
#' * \code{"exploratory"} is used in the Lakhesis calculator (\code{\link[lakhesis]{LC}}), with 100 iterations (\code{sim = 100, opt = 1}).
#' * \code{"optimize"} is recommended for post-exploratory analysis, to determine if a more optimal seriations can be found for the initial "exploratory" result. The number of iterations (\code{iter}) and simulation runs of iterations (\code{sim}) can be increased. The seriation obtained which yields the lowest (most optimal) measure of concentration \eqn{\kappa} (see \code{\link[lakhesis]{conc_kappa}}) chosen as the output. 
#'
#' @param strands A `list` of the `strands` class, comprised of `strand` objects returned by \code{\link[lakhesis]{ca.procrustes.curve}}.
#' @param obj The intial incidence matrix.
#' @param method Either \code{"exploratory"} or \code{"optimized"}. The \code{"exploratory"} method is teh default, used in the Lakhesis Calculator (\code{\link[lakhesis]{LC}}).
#' @param iter The number of iterations to use in lakhesizing the \code{strands}. The default is set at 100, also used for \code{"exploratory"}.
#' @param sim The number of simulation runs of iterations. The default is set at 1, also used for \code{"exploratory"}. Increasing the number of simulations will better ensure an optin 
#' @return A `list` of class `lakhesis` containing the following:
#' 
#' * \code{row} A seriated vector of row elements. 
#' * \code{col} A seriated vector of column elements
#' * \code{rowPCA} The results of \code{\link[stats]{prcomp}} performed on the row elements of \code{strands}.
#' * \code{colPCA} The results of \code{\link[stats]{prcomp}} performed on the column elements of \code{strands}.
#' * \code{coef}  A \code{data frame} containing the following columns:
#'   * \code{Strand} The number of the strand.
#'   * \code{Agreement} The measure of agreement, i.e., how well each strand accords with the consensus seriation. Using the square of Spearman's rank correlation coefficient, \eqn{\rho^2}, between each strand and the consensus ranking, agreement is computed as the product of \eqn{\rho^2} for their row and column rankings, \eqn{\rho_r^2}\eqn{\rho_c^2}. 
#'   * \code{Concentration} the concentration coefficient \eqn{\kappa}, which provides a measure of the optimality of each strand (see \code{\link[lakhesis]{kappa.coef}}).
#' * \code{im_seriated} The seriated incidence matrix, of class \code{incidence_matrix}.
#' 
#' @examples
#' data("qfStrands")
#' x <- lakhesize(qfStrands, method = "exploratory")
#' # summary(x) 
#'
#' @export
lakhesize <- function(strands, ...) {
    UseMethod("lakhesize")
}

#' @export 
lakhesize.strands <- function(strands, method = "exploratory", iter = 100, sim = 1) {
    lakhesize.default(strands, method, iter, sim)
}

#' @export
lakhesize.default <- function(strands, method = "exploratory", iter = 100, sim = 1) {    
    if (!(method %in% c("exploratory", "optimize"))) {
        stop('Method needs to be either "exploratory" or "optimize" ')
    } else {
        obj <- strands[[1]]$im_seriated
        for (i in 2:length(strands)) {
            obj <- im_merge(obj, strands[[i]]$im_seriated)
        }
        if (method == "exploratory") {
            iter <- 100
            sim <- 1
        }
        if (length(strands) > 1) {
            strand.mat <- strand_extract(strands)

            rowranks <- strand.mat[["Row"]]
            colranks <- strand.mat[["Col"]]

            if (method == "optimize") {
                pb <- utils::txtProgressBar(min = 0, max = sim, style = 3)
            }

            k <- 99999

            for (mOpt in 1:sim) {

                R <- matrix(NA, nrow = nrow(rowranks), ncol = 0)
                C <- matrix(NA, nrow = nrow(colranks), ncol = 0)
                
                for (m in 1:iter) { 
                    remaining <- 1:ncol(rowranks)
                    start <- sample(remaining, 1)
                    regressed <- data.frame(Row = rowranks[,start])
                    remaining <- remaining[-start]
                    
                    while (length(remaining) != 0) {
                        m2 <- sample(remaining, 1)
                        dat <- data.frame(regressed$Row, rowranks[,m2])
                        if (sum(!is.na(rowSums(dat))) > 3)  {
                            dat2 <- dat[!is.na(rowSums(dat)),]
                            y <- stats::prcomp(dat2)$x[,1]
                            fit1 <- stats::lm(y ~ dat2[,1])
                            fit2 <- stats::lm(y ~ dat2[,2])
                            regr1 <- dat[,1] * fit1$coef[2] + fit1$coef[1]
                            regr2 <- dat[,2] * fit2$coef[2] + fit2$coef[1]
                            regr <- matrix(c(regr1,regr2), ncol = 2)
                            merged <- rowMeans(regr, na.rm = TRUE)
                            merged <- rank(merged, na.last = "keep")
                            regressed$Row <- merged
                            remaining <- remaining[!(remaining %in% m2)]
                        }     
                    }
                    R <- cbind(R, regressed)
                }
                R <- R[!is.na(rowSums(R)),]

                for (m in 1:iter) { 
                    remaining <- 1:ncol(colranks)
                    start <- sample(remaining, 1)
                    regressed <- data.frame(Col = colranks[,start])
                    remaining <- remaining[-start]

                    while (length(remaining) != 0) {
                        m2 <- sample(remaining, 1)
                        dat <- data.frame(regressed$Col, colranks[,m2])
                        if (sum(!is.na(rowSums(dat))) > 3)  {
                            dat2 <- dat[!is.na(rowSums(dat)),]
                            y <- stats::prcomp(dat2, scale. = FALSE)$x[,1]
                            fit1 <- stats::lm(y ~ dat2[,1])
                            fit2 <- stats::lm(y ~ dat2[,2])
                            regr1 <- dat[,1] * fit1$coef[2] + fit1$coef[1]
                            regr2 <- dat[,2] * fit2$coef[2] + fit2$coef[1]
                            regr <- matrix(c(regr1,regr2), ncol = 2)
                            merged <- rowMeans(regr, na.rm = TRUE)
                            merged <- rank(merged, na.last = "keep")
                            regressed$Col <- merged
                            remaining <- remaining[!(remaining %in% m2)]
                        }     
                    }
                    C <- cbind(C, regressed)
                }
                C <- C[!is.na(rowSums(C)),]        

                consensus.row.pca <- stats::prcomp(R, scale. = FALSE)
                consensus.col.pca <- stats::prcomp(C, scale. = FALSE)
                consensus.r0 <- rank(consensus.row.pca$x[,1])
                consensus.c0 <- rank(consensus.col.pca$x[,1])

                RL0 <- names(consensus.r0)[order(consensus.r0)]
                CL0 <- names(consensus.c0)[order(consensus.c0)] 

                k0 <- conc_kappa(obj[RL0, CL0])

                if (k0 < k) {
                    k <- k0
                    Rpca <- consensus.row.pca
                    Cpca <- consensus.col.pca
                    RL <- RL0
                    CL <- CL0
                    consensus.r <- consensus.r0
                    consensus.c <- consensus.c0
                }
                if (method == "optimize") {
                    utils::setTxtProgressBar(pb, mOpt)
                }
            }

            if (method == "optimize") {
                utils::close(pb)
            }

            # concentration of each strand
            strand.k.c <- c()
            for (i in 1:length(strands)) {
                ctx <- stats::na.omit(rowranks[,i])
                fnd <- stats::na.omit(colranks[,i])
                roworder <- names(ctx[order(ctx)])
                colorder <- names(fnd[order(fnd)])
                strand.im <- obj[roworder,colorder]
                k.c <- conc_kappa(strand.im)
                strand.k.c <- c(strand.k.c, k.c)
            }

            # agreement with conensus seration
            assoc <- c()
            for (i in 1:length(strands)) {
                sp.f <- spearman_sq( rowranks[,i], consensus.r )
                sp.c <- spearman_sq( colranks[,i], consensus.c )
                assoc <- c(assoc,  sp.f * sp.c)
            }

            coefs <- data.frame( Strand = 1:length(strands), Agreement = assoc, Concentration = strand.k.c)

            results <- list()

            results[["row"]] <- RL #consensus.row.dat
            results[["col"]] <- CL #consensus.col.dat
            results[["rowPCA"]] <- Rpca
            results[["colPCA"]] <- Cpca
            results[["coef"]] <- coefs
            im <- obj[RL,CL]
            class(im) <- c("incidence_matrix", "matrix")
            results[["im_seriated"]] <- im
            class(results) <- c("lakhesis", "list")
            return(results)
        } else {
            warning("Need more than 1 strand to create consensus seration.")
        }
    }
}



