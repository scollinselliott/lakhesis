#' Lakhesize
#'
#' This function returns the row and column consensus serations for a data frame of strands, containing their rankings,
#' the results of their regression, and the Lakhesis coefficient of each strand. See documentation 
#' on \code{\link[lakhesis]{kappa.coef}} and  \code{\link[lakhesis]{spearman.sq}}.
#'
#' @param strands A list of `strands`, which are data frames returned by \code{\link[lakhesis]{ca.procrustes.curve}}.
#' @param obj The intial incidence matrix.
#' @return A list of the following:
#' * `RowConsensus` Data frame of the consensus seriation of the row elements in the order of their projection on the first principal axis.
#' * `ColConsensus` Data frhe consensus seriation of the column elements in the order of their project onto the first principal axis.
#' * `RowPCA` The results of \code{prcomp} performed on the regressed row elements of strands.
#' * `ColPCA` The results of \code{prcomp} performed on the regressed column elements of strands.
#' * `Coef`  Coefficient of agreement with the consensus seriation (the product of the square of Spearman's rank correlation coefficient, \eqn{\rho}, of the \eqn{i}th row and column strand with those of consensus, 
#' and the concentration coefficient \eqn{\kappa} (\code{\link[lakhesis]{kappa.coef}}).
#' 
#' @export
lakhesize <- function(strands, obj) {
    if (length(strands) > 1) {
        strand.mat <- strand.extract(strands, obj)

        rowranks <- strand.mat[["Row"]]
        colranks <- strand.mat[["Col"]]

        R <- matrix(NA, nrow = nrow(rowranks), ncol = 0)
        C <- matrix(NA, nrow = nrow(colranks), ncol = 0)

        for (m in 1:100) { 
            remaining <- 1:ncol(rowranks)
            start <- sample(remaining, 1)
            regressed <- data.frame(Row = rowranks[,start])
            remaining <- remaining[-start]
            
            while (length(remaining) != 0) {
                m2 <- sample(remaining, 1)
                dat <- data.frame(regressed$Row, rowranks[,m2])
                if (sum(!is.na(rowSums(dat))) > 3)  {
                    dat2 <- dat[!is.na(rowSums(dat)),]
                    y <- prcomp(dat2)$x[,1]
                    fit1 <- lm(y ~ dat2[,1])
                    fit2 <- lm(y ~ dat2[,2])
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

        for (m in 1:100) { 
            remaining <- 1:ncol(colranks)
            start <- sample(remaining, 1)
            regressed <- data.frame(Col = colranks[,start])
            remaining <- remaining[-start]

            while (length(remaining) != 0) {
                m2 <- sample(remaining, 1)
                dat <- data.frame(regressed$Col, colranks[,m2])
                if (sum(!is.na(rowSums(dat))) > 3)  {
                    dat2 <- dat[!is.na(rowSums(dat)),]
                    y <- prcomp(dat2)$x[,1]
                    fit1 <- lm(y ~ dat2[,1])
                    fit2 <- lm(y ~ dat2[,2])
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


        # for(i in 1:ncol(rowranks.imputed)){
        #     rowranks.imputed[is.na(rowranks.imputed[,i]), i] <- mean(rowranks.imputed[,i], na.rm = TRUE)
        # }
        # for(i in 1:ncol(colranks.imputed)){
        #     colranks.imputed[is.na(colranks.imputed[,i]), i] <- mean(colranks.imputed[,i], na.rm = TRUE)
        # }
        consensus.row.pca <- prcomp(R, scale. = FALSE)
        consensus.col.pca <- prcomp(C, scale. = FALSE)
        consensus.row.scores <- consensus.row.pca$x[,1]
        consensus.col.scores <- consensus.col.pca$x[,1]
        consensus.row.ranks <- rank(consensus.row.scores)
        consensus.col.ranks <- rank(consensus.col.scores)

        consensus.row.dat <- data.frame(Row = names(consensus.row.scores)[order(consensus.row.scores)] )
        consensus.col.dat <- data.frame(Col = names(consensus.col.scores)[order(consensus.col.scores)] )

        # concentration of each strand
        strand.k.c <- c()
        for (i in 1:length(strands)) {
            ctx <- na.omit(rowranks[,i])
            fnd <- na.omit(colranks[,i])
            roworder <- names(ctx[order(ctx)])
            colorder <- names(fnd[order(fnd)])
            strand.im <- obj[roworder,colorder]
            k.c <- kappa.coef(strand.im)
            strand.k.c <- c(strand.k.c, k.c)
        }

        # agreement with conensus seration
        assoc <- c()
        for (i in 1:length(strands)) {
            sp.f <- spearman.sq( rowranks[,i], consensus.row.ranks )
            sp.c <- spearman.sq( colranks[,i], consensus.col.ranks )
            assoc <- c(assoc,  sp.f * sp.c)
        }

        #lakh.K <- assoc / strand.k.c
        coefs <- data.frame( Strand = 1:length(strands), Concentration.Kappa = strand.k.c, Consensus.Spearman.Sq = assoc) #, Lakhesis.Coef = lakh.K )
        #coefs <- coefs[order(coefs$Lakhesis.Coef) , ]

        results = list()

        results[["RowConsensus"]] <- consensus.row.dat
        results[["ColConsensus"]] <- consensus.col.dat
        results[["RowPCA"]] <- consensus.row.pca
        results[["ColPCA"]] <- consensus.col.pca
        results[["Coef"]] <- coefs
        return(results)
    } else {
        print("Need more than 1 strand to create consensus seration.")
    }
}


