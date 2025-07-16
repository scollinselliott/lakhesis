#' Lakhesize
#'
#' This function returns the row and column consensus seriation for a \code{list} object of the \code{strands} class, containing their rankings, coefficients of association, and criterion. Consensus seriation is achieved by iterative simple linear regression to handle \code{NA} vales in each strand. To initialize, a regression is performed pairwise, with every strand as the dependent \eqn{y} variate and every other strand as the independent \eqn{x} variate. The independent variate's rankings are then regressed onto \eqn{f(x) = \hat{\beta}_1 x + \hat{\beta}_0}. If \eqn{y \neq f(x)}, the mean of \eqn{y} and \eqn{f(x)} is used. Then, the values of the dependent variate and those of the regressed independent varaite are re-ranked together, to form a combined ranking, which serves as the dependent variate on the next iteration. The pair of strands is chosen which minimizes a specified optimality criterion. The process is repeated until all strands have been regressed and re-ranked into a single consensus seriation.
#'
#' @param strands A \code{list} of \code{strands} class (see \code{\link[lakhesis]{strand_add}}).
#' @param crit The criterion used to assess the seration resulting from two strands \eqn{x} and \eqn{y}: 
#'   * \code{"cor_sp"} Computes Spearman's correlation coefficent for the incidences of 1s as points \eqn{(i,j)}. See \code{\link[lakhesis]{cor_sp}}.  Higher values are more optimal. This is the default option.
#'   * \code{"conc_wrc"} Computes weighted row-column concentration as the optimality criterion. See \code{\link[lakhesis]{conc_wrc}}. Lower values are more optimal.
#' 
#' @param pbar Displaying a progress bar. Default is \code{TRUE}.
#' @return A \code{list} of class \code{lakhesis} containing the following:
#' 
#' * \code{row} A seriated vector of row elements. 
#' * \code{col} A seriated vector of column elements
#' * \code{coef}  A \code{data frame} containing the following columns:
#'   * \code{Strand} The number of the strand.
#'   * \code{Agreement} The measure of agreement, i.e., how well each strand accords with the consensus seriation. Using the square of Spearman's rank correlation coefficient, \eqn{\rho^2}, between each strand and the consensus ranking, agreement is computed as the product of \eqn{\rho^2} for their row and column rankings, \eqn{\rho_r^2}\eqn{\rho_c^2}. 
#'   * \code{Criterion} Criterion of the optimality of each strand (per the \code{"crit"} option above).
#' * \code{im_seriated} The seriated incidence matrix, of class \code{incidence_matrix}.
#' 
#' @examples
#' data("qf_strands")
#' L <- lakhesize(qf_strands, pbar = FALSE)
#' # summary(L) 
#'
#' @export
lakhesize <- function(strands, crit = "cor_sp", pbar = TRUE) {
    UseMethod("lakhesize")
}

#' @rdname lakhesize
#' @export 
lakhesize.strands <- function(strands, crit = "cor_sp", pbar = TRUE) {
    lakhesize.default(strands, crit, pbar)
}

#' @rdname lakhesize
#' @export
lakhesize.default <- function(strands, crit = "cor_sp", pbar = TRUE) {    
    if (!(crit %in% c("cor_sp", "conc_wrc"))) {
        stop('Criterion must be "cor_sp" or "conc_wrc".')
    }
    obj <- strands[[1]]$im_seriated
    for (i in 2:length(strands)) {
        obj <- im_merge(obj, strands[[i]]$im_seriated)
    }
    ns <- length(strands)

    if (ns > 1) {
        strand.mat <- strand_extract(strands)

        rowranks <- strand.mat[["Row"]]
        colranks <- strand.mat[["Col"]]

        check_matrix <- matrix(NA, ns, ns)
        for (i in 1:ns) {
            for (j in i:ns) {
                if (i != j) {
                    dat <- rowranks[,c(i,j)]
                    check_matrix[i,j] <- nrow( dat[!is.na(rowSums(dat)),] ) 
                }
            }
        }

        check_matrix <- check_matrix > 3

        if (ns > 2) {
            check <- colSums(check_matrix[ ,2:ncol(check_matrix)], na.rm = TRUE) + rowSums(check_matrix[1:(ncol(check_matrix)-1) ,], na.rm = TRUE) 
        }

        if (ns == 2) {
            check <- check_matrix[1,2]
        }


        if (0 %in% check ) {
            warning("Each strand must share at least four joint elements with another strand.")
            return(NULL)
        } else {
            if (pbar == TRUE) {
                pb <- utils::txtProgressBar(min = 0, max = ns, style = 3)
            }   
            regressed <- c()
            kappa_matrix <- matrix(NA, ns, ns)

            for (i in 1:ns) {
                for (j in i:ns) {
                    if (!(i == j)) {
                        rdat <- data.frame(rowranks[,i], rowranks[,j])
                        rdat1 <- rdat[!is.na(rowSums(rdat)),]
                        cdat <- data.frame(colranks[,i], colranks[,j])
                        cdat1 <- cdat[!is.na(rowSums(cdat)),]
                        if ((nrow(rdat1) > 3) & (nrow(cdat1) > 3) ) {
                            ry <- rdat1[,2]
                            rx <- rdat1[,1]
                            rfit <- stats::lm(ry ~ rx)

                            rdat[,1] <- rdat[,1] * rfit$coef[2] + rfit$coef[1]

                            rdat <- rowMeans(rdat, na.rm = TRUE)
                            rdat <- rank(rdat, na.last = "keep")

                            cy <- cdat1[,2]
                            cx <- cdat1[,1]
                            cfit <- stats::lm(cy ~ cx)

                            cdat[,1] <- cdat[,1] * cfit$coef[2] + cfit$coef[1]

                            cdat <- rowMeans(cdat, na.rm = TRUE)
                            cdat <- rank(cdat, na.last = "keep")

                            R <- names(rdat)[order(rdat)]
                            C <- names(cdat)[order(cdat)]
                            if (crit == "conc_wrc") {
                                K <- conc_wrc(obj[R, C])
                            } else if (crit == "cor_sp") {
                                K <- cor_sp(obj[R, C])
                            }
                            kappa_matrix[i,j] <- K
                        }
                    }
                }
            }

            regressed <- rev( arrayInd(which.min(kappa_matrix), dim(kappa_matrix)))
            remaining <- 1:ns
            remaining <- remaining[-regressed]

            rdat <- data.frame(rowranks[,regressed[2]], rowranks[,regressed[1]])
            rdat1 <- rdat[!is.na(rowSums(rdat)),]
            cdat <- data.frame(colranks[,regressed[2]], colranks[,regressed[1]])
            cdat1 <- cdat[!is.na(rowSums(cdat)),]

            ry <- rdat1[,2]
            rx <- rdat1[,1]
            rfit <- stats::lm(ry ~ rx)

            rdat[,1] <- rdat[,1] * rfit$coef[2] + rfit$coef[1]

            rdat <- rowMeans(rdat, na.rm = TRUE)
            rdat_y <- rank(rdat, na.last = "keep")

            cy <- cdat1[,2]
            cx <- cdat1[,1]
            cfit <- stats::lm(cy ~ cx)

            cdat[,1] <- cdat[,1] * cfit$coef[2] + cfit$coef[1]

            cdat <- rowMeans(cdat, na.rm = TRUE)
            cdat_y <- rank(cdat, na.last = "keep")

            while (length(remaining) > 0) {
                kappa_check = numeric(ns)
                kappa_check[] <- NA
                for (i in remaining) {
                    rdat <- data.frame(rowranks[,i], rdat_y)
                    rdat1 <- rdat[!is.na(rowSums(rdat)),]
                    cdat <- data.frame(colranks[,i], cdat_y)
                    cdat1 <- cdat[!is.na(rowSums(cdat)),]
                    if ((nrow(rdat1) > 3) & (nrow(cdat1) > 3) ) {
                        ry <- rdat1[,2]
                        rx <- rdat1[,1]
                        rfit <- stats::lm(ry ~ rx)

                        rdat[,1] <- rdat[,1] * rfit$coef[2] + rfit$coef[1]

                        rdat <- rowMeans(rdat, na.rm = TRUE)
                        rdat <- rank(rdat, na.last = "keep")

                        cy <- cdat1[,2]
                        cx <- cdat1[,1]
                        cfit <- stats::lm(cy ~ cx)

                        cdat[,1] <- cdat[,1] * cfit$coef[2] + cfit$coef[1]

                        cdat <- rowMeans(cdat, na.rm = TRUE)
                        cdat <- rank(cdat, na.last = "keep")

                        R <- names(rdat)[order(rdat)]
                        C <- names(cdat)[order(cdat)]

                        if (crit == "conc_wrc") {
                            K <- conc_wrc(obj[R, C])
                        } else if (crit == "cor_sp") {
                            K <- cor_sp(obj[R, C])
                        }

                        kappa_check[i] <- K
                    }
                }
                idx_next <- which.min(kappa_check)
                regressed <- c(regressed, idx_next)
                remaining <- remaining[-which(remaining == idx_next)]

                rdat <- data.frame(rowranks[,idx_next], rdat_y)
                rdat1 <- rdat[!is.na(rowSums(rdat)),]
                cdat <- data.frame(colranks[,idx_next], cdat_y)
                cdat1 <- cdat[!is.na(rowSums(cdat)),]

                ry <- rdat1[,2]
                rx <- rdat1[,1]
                rfit <- stats::lm(ry ~ rx)

                rdat[,1] <- rdat[,1] * rfit$coef[2] + rfit$coef[1]

                rdat <- rowMeans(rdat, na.rm = TRUE)
                rdat_y <- rank(rdat, na.last = "keep")

                cy <- cdat1[,2]
                cx <- cdat1[,1]
                cfit <- stats::lm(cy ~ cx)

                cdat[,1] <- cdat[,1] * cfit$coef[2] + cfit$coef[1]

                cdat <- rowMeans(cdat, na.rm = TRUE)
                cdat_y <- rank(cdat, na.last = "keep")

                if (pbar == TRUE) {
                    utils::setTxtProgressBar(pb, length(regressed))
                }   
            }

            RL <- names(rdat_y)[order(rdat_y)]
            CL <- names(cdat_y)[order(cdat_y)]

            if (pbar == TRUE) {
                close(pb)
            }

            # criterion of each strand
            strand.k.c <- c()
            if (crit == "conc_wrc") {
                for (i in 1:length(strands)) {
                    ctx <- stats::na.omit(rowranks[,i])
                    fnd <- stats::na.omit(colranks[,i])
                    roworder <- names(ctx)[order(ctx)]
                    colorder <- names(fnd)[order(fnd)]
                    strand.im <- obj[roworder,colorder]
                    k.c <- conc_wrc(strand.im)
                    strand.k.c <- c(strand.k.c, k.c)
                }
            } else if (crit == "cor_sp") {
                for (i in 1:length(strands)) {
                    ctx <- stats::na.omit(rowranks[,i])
                    fnd <- stats::na.omit(colranks[,i])
                    roworder <- names(ctx)[order(ctx)]
                    colorder <- names(fnd)[order(fnd)]
                    strand.im <- obj[roworder,colorder]
                    k.c <- cor_sp(strand.im)
                    strand.k.c <- c(strand.k.c, k.c)
                }
            }

            # agreement with conensus seration
            assoc <- c()
            for (i in 1:length(strands)) {
                sp.f <- spearman_sq( rowranks[,i], rdat_y )
                sp.c <- spearman_sq( colranks[,i], cdat_y )
                assoc <- c(assoc,  sp.f * sp.c)
            }

            coefs <- data.frame( Strand = 1:length(strands), Agreement = assoc, Criterion = strand.k.c)

            results <- list()

            results[["row"]] <- RL #consensus.row.dat
            results[["col"]] <- CL #consensus.col.dat
            results[["coef"]] <- coefs
            im <- obj[RL,CL]
            class(im) <- c("incidence_matrix", "matrix")
            results[["im_seriated"]] <- im
            class(results) <- c("lakhesis", "list")
            return(results)
        } 
    } else {
        warning("Need more than 1 strand to create consensus seration.")
    }
}



