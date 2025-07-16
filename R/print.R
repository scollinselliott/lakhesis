#' @export 
print.procrustean <- function(x, ...) {
    cat("Procrustes-fit CA scores of class", class(x), "\n")
    cat("  $ref : principal scores: reference points\n")
    cat("  $x : principal scores: row points\n")
    cat("  $y : principal scores: column points\n")
}



#' @export 
plot.procrustean <- function(x, ...) {
    Procrustes1 <- Procrustes2 <- NULL
    ref <- data.frame(Procrustes1 = x$ref[,1], Procrustes2 = x$ref[,2])
    r.dat <- data.frame(Procrustes1 = x$x[,1], Procrustes2 = x$x[,2])
    c.dat <- data.frame(Procrustes1 = x$y[,1], Procrustes2 = x$y[,2])
    Type <- c(rep("row", nrow(x$x)), rep("col", nrow(x$y)) )
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
print.strand <- function(x, ...) {
    cat("Procrustes-fit CA of class:", class(x), "\nUse summary() for more information. \n")
}



#' @export 
summary.strand <- function(object, ...) {
    s <- object$dat
    cat("Procrustes-fit CA of class:", class(object), "\n")
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



#' @export 
plot.strand <- function(x, display = "ca", ...) {
    CurveIndex <- Distance <- Type <- NULL
    #ref <- im_ref( matrix(NA, nrow = sum(strand$Type == 'row'),  ncol = sum(strand$Type == 'col')  ))
    #refcurve <- ca_procrustes_curve(ref)
    if (display == "ca") {
        plot(ca_procrustes(x$im_seriated))
    } else if (display == "ref") {
        dat <- x$dat
        ord.plot <- ggplot2::ggplot() + 
            ggplot2::geom_point(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = Type ), size = 2) + 
            ggplot2::geom_text(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = Type ), label = rownames(dat), hjust = "left", nudge_y = 0.0005, size = 3,  angle = 90,  check_overlap = TRUE) + 
            ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = .3, legend.position="none")
        print(ord.plot)
    } else if (display == "im_seriated") {
        plot(x$im_seriated)
        #strand2 <- strand$im_seriated
        #im.Image <- t(strand2[nrow(strand2):1 , ])
        #image(im_Image, col=c('white','black'), xaxt='n', yaxt='n')
    } else {
        stop('Display option must be "ca", "ref", "im_serated"')
    }
}




#' @export 
print.strands <- function(x, ...) {
    cat('List of', length(x), 'strands\n')
}





#' @export 
print.lakhesis <- function(x, ...) {
    cat("Lakhesis analysis of class:", class(x), "\nUse summary() for more information. \n")
}




#' @export 
summary.lakhesis <- function(object, ...) {
    cat("Lakhesis analysis of class:", class(object), "\n")
    cat("\n")
    cat("Ranking of row elements: \n")
    cat(object$row, sep = ", ", fill = TRUE)
    cat("\n")
    cat("Ranking of column elements: \n")
    cat(object$col, sep = ", ", fill = TRUE)
    cat("\n")
    cat("Seriated incidience matrix in $im_seriated\n")
    cat("Coefficients in $coef\n")
    cat("cor_sp = ", cor_sp(object$im_seriated), "\n")
    cat("conc_wrc = ", conc_wrc(object$im_seriated), "\n")
}



#' @export 
plot.lakhesis <- function(x, display = "im_seriated", ...) {
    Strand <- Agreement <- Criterion <- NULL
    lakhcoef <- x$coef
    if (display == "im_seriated") {
        im_seriated <- x$im_seriated
        k_sp <- cor_sp(im_seriated)
        k_wrc <- conc_wrc(im_seriated)
        ttl <- paste(format(nrow(im_seriated))," x ",format(ncol(im_seriated)),"; cor_sp = ",format(round(k_sp,3)), "; conc_wrc = ",format(round(k_wrc,3)), sep = "")
        im_Image <- t(im_seriated[nrow(im_seriated):1 , ])
        graphics::image(im_Image, col=c('white','black'), xaxt='n', yaxt='n', main = ttl)
    } else if (display == "agreement") {
        lakhcoef$Strand <- factor(lakhcoef$Strand, levels = lakhcoef$Strand[order(lakhcoef$Agreement, decreasing = FALSE)]) 
        plot_agreement <- ggplot2::ggplot(lakhcoef, ggplot2::aes(x=Strand, y=Agreement)) + 
            ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()
        print( plot_agreement )
    } else if (display == "criterion") {
        lakhcoef$Strand <- factor(lakhcoef$Strand, levels = lakhcoef$Strand[order(lakhcoef$Criterion, decreasing = TRUE)]) 
        plot_criterion <- ggplot2::ggplot(lakhcoef, ggplot2::aes(x=Strand, y=Criterion)) + 
            ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()
        print( plot_criterion )
    } else {
        stop("Choose a valid display option: im_seriated, agreement, criterion")
    }
}



#' @export 
plot.incidence_matrix <- function(x, ...) {
    k_wrc <- conc_wrc(x)
    k_sp <- cor_sp(x)
    ttl <- paste(format(nrow(x))," x ",format(ncol(x)),"; cor_sp = ", format(round(k_sp,3)), "; conc_rc = ",format(round(k_wrc,3)), sep = "")
    im_Image <- t(x[nrow(x):1 , ])
    graphics::image(im_Image, col=c('white','black'), xaxt='n', yaxt='n', main = ttl)
}

