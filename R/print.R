
#' @export 
print.procrustean <- function(result) {
    cat("Procrustes-fit CA scores of class", class(result), "\n")
    cat("$ref : matrix of reference matrix scores\n")
    cat("$x, $y : principal scores\n")
}




#' @rdname ca_procrustes
#' @export 
plot.procrustean <- function(result) {
    ref <- data.frame(Procrustes1 = result$ref[,1], Procrustes2 = result$ref[,2])
    r.dat <- data.frame(Procrustes1 = result$x[,1], Procrustes2 = result$x[,2])
    c.dat <- data.frame(Procrustes1 = result$y[,1], Procrustes2 = result$y[,2])
    Type <- c(rep("row", nrow(result$x)), rep("col", nrow(result$y)) )
    dat <- cbind( rbind(r.dat, c.dat), Type)
    #ref <- im_ref( matrix(NA, nrow = sum(strand$Type == 'row'),  ncol = sum(strand$Type == 'col')  ))
    #refcurve <- ca_procrustes_curve(ref)
    curve.plot <- ggplot2::ggplot() +
        ggplot2::geom_point(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = Type), size = 2) +
        ggplot2::geom_text(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = Type ), label = rownames(dat), size = 3, hjust = 0.025, nudge_x = 0.015, check_overlap = TRUE) + 
        ggplot2::geom_line(data = ref,  ggplot2::aes(x = Procrustes1, y = Procrustes2), linewidth=1, alpha=0.4, linetype=1) +
        ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1, legend.position="none") 
    curve.plot
}




#' @export 
print.strand <- function(strand) {
    cat("Procrustes-fit CA of class:", class(strand), "\nUse summary() for more information. \n")
}




#' @export 
summary.strand <- function(strand) {
    s <- strand$dat
    cat("Procrustes-fit CA of class:", class(strand), "\n")
    cat("\n")
    cat("Ranking of row elements: \n")
    rows <- s[s$Type == "row", ]
    cat(rownames(rows)[order(rows$Rank)], fill = TRUE)
    cat("\n")
    cat("Ranking of column elements: \n")
    cols <- s[s$Type == "col", ]
    cat(rownames(cols)[order(cols$Rank)], fill = TRUE)
    cat("\n")
    cat("$dat contains the following columns:\n")
    cat("   $Procrustes1, $Procrustes2 : x,y coords for the Procrustes-fit CA principal scores\n")
    cat("   $CurveIndex: index of nearest point on refrence curve to the Procrustes-fit CA score point\n")
    cat("   $Distance: distance of principal score point to nearest point on reference curve\n")
    cat("   $Rank: ranking of score points projected onto the reference curve\n")
    cat('   $Type: "row" or "col"\n')
    cat("$im_seriated: the seriated incidence matrix")
    cat('\n')
}



#' @rdname ca_procrustes_ser
#' @export 
plot.strand <- function(strand, display = "ca") {
    #ref <- im_ref( matrix(NA, nrow = sum(strand$Type == 'row'),  ncol = sum(strand$Type == 'col')  ))
    #refcurve <- ca_procrustes_curve(ref)
    if (display == "ca") {
        dat <- strand$dat
        curve.plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = Type), size = 2) +
            ggplot2::geom_text(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = Type ), label = rownames(dat), size = 3, hjust = 0.025, nudge_x = 0.015, check_overlap = TRUE) + 
            ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1, legend.position="none") 
            #ggplot2::geom_line(data = refcurve,  ggplot2::aes(x = Procrustes1, y = Procrustes2), linewidth=1, alpha=0.4, linetype=1) 
        print(curve.plot)
    } else if (display == "ref") {
        dat <- strand$dat
        ord.plot <- ggplot2::ggplot() + 
            ggplot2::geom_point(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = Type ), size = 2) + 
            ggplot2::geom_text(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = Type ), label = rownames(dat), hjust = "left", nudge_y = 0.0005, size = 3,  angle = 90,  check_overlap = TRUE) + 
            ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = .3, legend.position="none")
        print(ord.plot)
    } else if (display == "im_seriated") {
        plot(strand$im_seriated)
        #strand2 <- strand$im_seriated
        #im.Image <- t(strand2[nrow(strand2):1 , ])
        #image(im_Image, col=c('white','black'), xaxt='n', yaxt='n')
    } else {
        stop('Display option must be "ca", "ref", "im_serated"')
    }
}




#' @export 
print.strands <- function(strands) {
    cat('List of', length(strands), 'strands\n')
}





#' @export 
print.lakhesis <- function(result) {
    cat("Lakhesis analysis of class:", class(result), "\nUse summary() for more information. \n")
}




#' @export 
summary.lakhesis <- function(result) {
    cat("Lakhesis analysis of class:", class(result), "\n")
    cat("\n")
    cat("Ranking of row elements: \n")
    cat(result$row, fill = TRUE)
    cat("\n")
    cat("Ranking of column elements: \n")
    cat(result$col, fill = TRUE)
    cat("\n")
    cat("PCA results contained in $rowPCA, $colPCA:\n")
    cat("Seriated incidience matrix in $im_seriated\n")
    cat("Coefficients in $coef\n")
    cat("kappa = ", conc_kappa(result$im_seriated))
}



#' @rdname lakhesize
#' @export 
plot.lakhesis <- function(result, display = "im_seriated") {
    lakhcoef <- result$coef
    if (display == "im_seriated") {
        im_seriated <- result$im_seriated
        k <- conc_kappa(im_seriated)
        ttl <- paste(format(nrow(im_seriated))," x ",format(ncol(im_seriated)),"; kappa = ",format(k), sep = "")
        im_Image <- t(im_seriated[nrow(im_seriated):1 , ])
        graphics::image(im_Image, col=c('white','black'), xaxt='n', yaxt='n', main = ttl)
    } else if (display == "rowPCA") {
        biplot(result$rowPCA)
    } else if (display == "colPCA") {
        biplot(result$colPCA)
    } else if (display == "agreement") {
        lakhcoef$Strand <- factor(lakhcoef$Strand, levels = lakhcoef$Strand[order(lakhcoef$Agreement, decreasing = FALSE)]) 
        plot_agreement <- ggplot2::ggplot(lakhcoef, ggplot2::aes(x=Strand, y=Agreement)) + 
            ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()
        print( plot_agreement )
    } else if (display == "concentration") {
        lakhcoef$Strand <- factor(lakhcoef$Strand, levels = lakhcoef$Strand[order(lakhcoef$Concentration, decreasing = TRUE)]) 
        plot_concentration <- ggplot2::ggplot(lakhcoef, ggplot2::aes(x=Strand, y=Concentration)) + 
            ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()
        print( plot_concentration )
    } else {
        stop("Choose a valid display option: im_seriated, rowPCA, colPCA, agreement, concentration")
    }
}



#' @rdname im_read_csv
#' @export 
plot.incidence_matrix <- function(im_seriated) {
    k <- conc_kappa(im_seriated)
    ttl <- paste(format(nrow(im_seriated))," x ",format(ncol(im_seriated)),"; kappa = ",format(k), sep = "")
    im_Image <- t(im_seriated[nrow(im_seriated):1 , ])
    graphics::image(im_Image, col=c('white','black'), xaxt='n', yaxt='n', main = ttl)
}

