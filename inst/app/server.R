options(shiny.maxRequestSize=50*1024^2)

server <- function(input, output, session) { 

    # creating reactive values
    mats <- reactiveValues(mat = NULL, mat_initial = NULL, caproc = NULL, caproc_ref = NULL, nr = NULL, devr = NULL, devc = NULL)

    sup <- reactiveValues(rowsel = NULL, colsel = NULL)

    mats$mat <- im_ref(matrix(NA, 20, 20))
    mats$mat_initial <- im_ref(matrix(NA, 20, 20))
    mats$caproc <- ca_procrustes_ser(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
    mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
    mats$nr <- nrow(isolate(mats$caproc)$dat)
    mats$devr <- data.frame(Message = c("Lakhesize and run deviance test."))
    mats$devc <- data.frame(Message = c("Lakhesize and run deviance test."))
    sup$rowsel <- rownames(isolate(mats$mat_initial))
    sup$colsel <- colnames(isolate(mats$mat_initial))

    observeEvent(input$datafile, {
        mats$mat <- im_read_csv(input$datafile$datapath, char = input$char, remove.hapax = as.logical(input$hapax))
        mats$mat_initial <- mats$mat
        mats$caproc <- ca_procrustes_ser(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
        mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
        mats$nr <- nrow(isolate(mats$caproc)$dat)
        selections$caplot <- rep(FALSE, isolate(mats$nr))
        selections$biplot <- rep(FALSE, isolate(mats$nr))
        selections$curve <- rep(FALSE, isolate(mats$nr))

        sup$rowsel <- rownames(isolate(mats$mat_initial))
        sup$colsel <- colnames(isolate(mats$mat_initial))

        updateSelectInput(session,
            "suppressrows",
            choices = isolate(sup$rowsel)
        )
    
        updateSelectInput(session,
            "suppresscols",
            choices = isolate(sup$colsel)
        )
    })

    updateSelectInput(session,
        "suppressrows",
        choices = isolate(sup$rowsel)
    )
    
    updateSelectInput(session,
        "suppresscols",
        choices = isolate(sup$colsel)
    )

    # interactive plots
    selections <- reactiveValues(caplot = NULL, biplot = NULL, curve = NULL )
    selections$caplot <- rep(FALSE, isolate(mats$nr))
    selections$biplot <- rep(FALSE, isolate(mats$nr))
    selections$curve <- rep(FALSE, isolate(mats$nr))

    observeEvent(input$plot_brush_caplot, {
        brushed_caplot <- brushedPoints(mats$caproc$dat, input$plot_brush_caplot, allRows = TRUE)$selected_
        selections$caplot <- (brushed_caplot | selections$caplot)
    })

    observeEvent(input$plot_brush_biplot, {
        brushed_biplot <- brushedPoints(mats$caproc$dat, input$plot_brush_biplot, allRows = TRUE)$selected_
        selections$biplot <- (brushed_biplot | selections$biplot)
    })

    observeEvent(input$plot_brush_curve, {
        brushed_curve <- brushedPoints(mats$caproc$dat, input$plot_brush_curve, allRows = TRUE)$selected_
        selections$curve <- (brushed_curve | selections$curve)        
    })

    observeEvent(input$plot_reset, {
        nr <- mats$nr
        selections$caplot <- rep(FALSE, nr)
        selections$biplot <- rep(FALSE, nr)
        selections$curve <- rep(FALSE, nr)
    })

    output$caplot <- renderPlot({
        dat <- mats$caproc$dat
        dat$sel <- selections$caplot

        ca_plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = dat, ggplot2::aes(x = CA1, y = CA2, color = interaction(Type, sel, ":") ), size = 2) +
            ggplot2::geom_text(data = dat, ggplot2::aes(x = CA1, y = CA2, color = interaction(Type, sel, ":") ), label = rownames(dat), size = 3, hjust = 0.025, nudge_x = 0.015, check_overlap = TRUE) + 
            ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1, legend.position="none")  
        ca_plot
    },  res = 96)

    output$biplot <- renderPlot({
        dat <- mats$caproc$dat
        dat$sel <- selections$biplot

        refcurve <- isolate(mats$caproc_ref)
        refcurve <- refcurve$ref
        refcurve <- data.frame(Procrustes1 = refcurve[,1], Procrustes2 = refcurve[,2])

        curve.plot <- ggplot2::ggplot() +
            ggplot2::geom_point(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = interaction(Type, sel, ":") ), size = 2) +
                ggplot2::geom_text(data = dat, ggplot2::aes(x = Procrustes1, y = Procrustes2, color = interaction(Type, sel, ":") ), label = rownames(dat), size = 3, hjust = 0.025, nudge_x = 0.015, check_overlap = TRUE) + 
                ggplot2::geom_line(data = refcurve,  ggplot2::aes(x = Procrustes1, y = Procrustes2), linewidth=1, alpha=0.4, linetype=1) + 
                ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1, legend.position="none")  
        curve.plot
    },  res = 96)

    output$procrustesplot <- renderPlot({
        dat <- mats$caproc$dat
        dat$sel <- selections$curve #selected_curve()

        ord.plot <- ggplot2::ggplot() + 
            ggplot2::geom_point(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = interaction(Type, sel, ":") ), size = 2) + 
                ggplot2::geom_text(data = dat, ggplot2::aes(x = CurveIndex, y = Distance, color = interaction(Type, sel, ":") ), label = rownames(dat), hjust = "left", nudge_y = 0.0005, size = 3,  angle = 90,  check_overlap = TRUE) + 
                ggplot2::theme_bw() +  ggplot2::theme(aspect.ratio = .3, legend.position="none")
        ord.plot
    }, res = 96)

    # re-run CA-Procrustes on selected points
    observeEvent(input$rerun, {
        m <- mats$mat
        keepRow <- !(rownames(m) %in% input$removeRow)
        keepCol <- !(colnames(m) %in% input$removeCol)
        keep.rc <- c(keepRow, keepCol)
        keep <- ((selections$caplot | selections$biplot)  | selections$curve ) & keep.rc
        dat <- mats$caproc$dat

        if ((length(rownames(dat)[keep & (dat$Type == "row") ]) > 2) & (length(rownames(dat)[keep & (dat$Type == "col") ]) > 2)) {

            tmp <- mats$mat[rownames(dat)[keep & (dat$Type == "row") ], rownames(dat)[keep & (dat$Type == "col") ] ]  
            tmp <- tmp[rowSums(tmp) != 0 , ]
            tmp <- tmp[,  colSums(tmp) != 0]

            if ( all( is.finite( tmp ) ) ) {

            mats$mat <- tmp
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)

            selections$caplot <- rep(FALSE, isolate(mats$nr))
            selections$biplot <- rep(FALSE, isolate(mats$nr))
            selections$curve <- rep(FALSE, isolate(mats$nr))
            } else {
                message("Infinite values in selection.")
            }
        } else {
            message("Insufficient number of points selected.")
        }
    })

    # re-run CA on selected points
    observeEvent(input$rerun, {
        m <- mats$mat
        keepRow <- !(rownames(m) %in% input$removeRow)
        keepCol <- !(colnames(m) %in% input$removeCol)
        keep.rc <- c(keepRow, keepCol)
        keep <- ((selections$caplot | selections$biplot)  | selections$curve ) & keep.rc
        dat <- mats$caproc$dat

        if ((length(rownames(dat)[keep & (dat$Type == "row") ]) > 2) & (length(rownames(dat)[keep & (dat$Type == "col") ]) > 2)) {

            tmp <- mats$mat[rownames(dat)[keep & (dat$Type == "row") ], rownames(dat)[keep & (dat$Type == "col") ] ]  
            tmp <- tmp[rowSums(tmp) != 0 , ]
            tmp <- tmp[,  colSums(tmp) != 0]

            if ( all( is.finite( tmp ) ) ) {

            mats$mat <- tmp
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)

            selections$caplot <- rep(FALSE, isolate(mats$nr))
            selections$biplot <- rep(FALSE, isolate(mats$nr))
            selections$curve <- rep(FALSE, isolate(mats$nr))
            } else {
                message("Infinite values in selection.")
            }
        } else {
            message("Insufficient number of points selected.")
        }
    })

    # restart the Lakhesis Calculator with the starting matrix
    observeEvent(input$reinitialize, {
        selections$caplot[] <- FALSE
        selections$biplot[] <- FALSE
        selections$curve[] <- FALSE

        mats$mat <- mats$mat_initial
        mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
        mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
        mats$nr <- nrow(mats$caproc$dat)
        selections$caplot <- rep(FALSE, isolate(mats$nr))
        selections$biplot <- rep(FALSE, isolate(mats$nr))
        selections$curve <- rep(FALSE, isolate(mats$nr))
    })

    # processing results

    # reactives
    results <- reactiveValues(strands = list(), strand_hold = list,  strand_backup = NULL, lakhesis_results = NULL, lakhesized = FALSE)

    # record the observed seriation as a strand
    observeEvent(input$log_ca1, {
        mats$caproc <- ca_procrustes_ser(mats$mat, projection = "ca1", symmetric = isolate(as.logical(input$sym)))
        results$strands <- strand_add(isolate(mats$caproc), results$strands)
        results$strand_backup <- isolate(results$strands)
        #print(isolate(results$strands))
    })
    observeEvent(input$log_ca2, {
        mats$caproc <- ca_procrustes_ser(mats$mat, projection = "ca2", symmetric = isolate(as.logical(input$sym)))
        results$strands <- strand_add(isolate(mats$caproc), results$strands)
        results$strand_backup <- isolate(results$strands)
        #print(isolate(results$strands))
    })
    observeEvent(input$log_proc1, {
        mats$caproc <- ca_procrustes_ser(mats$mat, projection = "procrustes1", symmetric = isolate(as.logical(input$sym)))
        results$strands <- strand_add(isolate(mats$caproc), results$strands)
        results$strand_backup <- isolate(results$strands)
        #print(isolate(results$strands))
    })
    observeEvent(input$log_proc2, {
        mats$caproc <- ca_procrustes_ser(mats$mat, projection = "procrustes2", symmetric = isolate(as.logical(input$sym)))
        results$strands <- strand_add(isolate(mats$caproc), results$strands)
        results$strand_backup <- isolate(results$strands)
        #print(isolate(results$strands))
    })
    observeEvent(input$log_curve, {
        mats$caproc <- ca_procrustes_ser(mats$mat, projection = "curve", symmetric = isolate(as.logical(input$sym)))
        results$strands <- strand_add(isolate(mats$caproc), results$strands)
        results$strand_backup <- isolate(results$strands)
        #print(isolate(results$strands))
    })

    # lachesize strands
    observeEvent(input$lakhesize, {    
        if (length(isolate(results$strands)) > 1) {
            s <- isolate(results$strands)
            results$strand_backup <- isolate(results$strands)

            results$lakhesis_results <- lakhesize(s, crit = input$crits, pbar = FALSE) 

            if ("lakhesis" %in% class(isolate(results$lakhesis_results))) {
                results$lakhesized <- TRUE
            }

            selections$caplot[] <- FALSE
            selections$biplot[] <- FALSE
            selections$curve[] <- FALSE

            mats$mat <- mats$mat_initial
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)

            selections$caplot <- rep(FALSE, mats$nr)
            selections$biplot <- rep(FALSE, mats$nr)
            selections$curve <- rep(FALSE, mats$nr)
            
        }
    })

    # # plots output after lakhesize() is run
    # output$consensusrowplot <- renderPlot({
    #     if (results$lakhesized == TRUE) {
    #         L <- results$lakhesis_results
    #         plot(L, display = "rowPCA")
    #     } else {
    #         plot(0,0, pch = " ", xlab = "Lakhesize to produce consensus PCA biplot.", ylab = " " )
    #     }
    # })

    # output$consensuscolplot <- renderPlot({
    #     if (results$lakhesized == TRUE) {
    #         L <- results$lakhesis_results
    #         plot(L, display = "colPCA")
    #     } else {
    #         plot(0,0, pch = " ", xlab = "Lakhesize to produce consensus PCA biplot.", ylab = " " )
    #     }
    # })

    output$consensusmatrixplot <- renderPlot({
        if (results$lakhesized == TRUE) {
            L <- results$lakhesis_results
            plot(L, display = "im_seriated")

        } else {
            plot(0,0, pch = " ", xlab = "Lakhesize to produce consensus matrix plot.", ylab = " " )
        }
    })

    # coefficient plots
    output$consensusplot <- renderPlot({
        if (results$lakhesized == TRUE) {
            L <- results$lakhesis_results
            #lakhcoef$Strand <- factor(lakhcoef$Strand, levels = lakhcoef$Strand[order(lakhcoef$Agreement, decreasing = FALSE)]) 
            suppressWarnings({
                plot(L, display = "agreement")
            })
            } else {
                plot(0,0, pch = " ", xlab = "Lakhesize to render plot of coefficents for strands.", ylab = " " )
        } 
    })

    output$kappaplot <- renderPlot({
        if (results$lakhesized == TRUE) {
            L <- results$lakhesis_results
            suppressWarnings({
                plot(L, display = "criterion")
            })
            } else {
                plot(0,0, pch = " ", xlab = "Lakhesize to render plot of coefficents for strands.", ylab = " " )
        }   
    })

    # delete a strand
    observeEvent(input$submitStrandDelete, {
        if (results$lakhesized == TRUE) {

            delValue <- as.integer(input$strandDelete)
            
            results$strand_hold <- isolate(results$strands[delValue])

            results$strands[delValue] <- NULL

            if (length(results$strands) > 1) {

            mats$mat <- mats$mat_initial
            suppressWarnings({
               results$lakhesis_results <- lakhesize(isolate(results$strands), crit = input$crits, pbar = FALSE) 
            })
            selections$caplot[] <- FALSE
            selections$biplot[] <- FALSE
            selections$curve[] <- FALSE

            mats$mat <- mats$mat_initial
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)
            } else {
                mats$mat <- mats$mat_initial
                mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
                mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
                mats$nr <- nrow(mats$caproc$dat)
            }
            selections$caplot[] <- FALSE
            selections$biplot <- rep(FALSE, mats$nr)
            selections$curve <- rep(FALSE, mats$nr)
        }
    })

    observeEvent(input$submitStrandDeleteUndo, {
        if (results$lakhesized == TRUE) {

            results$strands[[length(results$strands) + 1]] <- isolate(results$strand_hold[[1]])

            results$strand_hold <- NULL

            mats$mat <- mats$mat_initial
            suppressWarnings({
            results$lakhesis_results <- lakhesize(isolate(results$strands), crit = input$crits, pbar = FALSE) 
            })
            selections$caplot[] <- FALSE
            selections$biplot[] <- FALSE
            selections$curve[] <- FALSE

            mats$mat <- mats$mat_initial
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)

            selections$caplot[] <- FALSE
            selections$biplot <- rep(FALSE, mats$nr)
            selections$curve <- rep(FALSE, mats$nr)
        }
    })

    # deviance test results
    observeEvent(input$deviancetest, {
        if (results$lakhesized == TRUE) {
            L <- isolate(results$lakhesis_results)

            dev <- element_eval(L$im_seriated)
            mats$devc <- head(dev$Col, n = 10L)
            mats$devr <- head(dev$Row, n = 10L)
        }
    })

    output$deviancerows <- renderTable(mats$devr)
    output$deviancecols <- renderTable(mats$devc)

    # export data to a .rds object
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("Lakhesis", format(Sys.Date(), '%Y%m%d'), format(Sys.time(), '%H%M%S'), ".rds", sep="")
        },
        content = function(file) {
            # reset plot when exporting
            results$lakhesized <- TRUE

            selections$caplot[] <- FALSE
            selections$biplot[] <- FALSE
            selections$curve[] <- FALSE

            # lakhesize to ensure that all row/col elements from strands are in results if user has not performed this action
            s <- isolate(results$strands)
            suppressWarnings({
            results$lakhesis_results <- lakhesize(s, crit = input$crits, pbar = FALSE) 
            })
            lr <- isolate(results$lakhesis_results)
            results <- list(consensus = lr, strands = s)
            saveRDS(results, file = file)
          
            mats$mat <- mats$mat_initial
            mats$caproc <- ca_procrustes_ser(mats$mat, symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(mats$caproc$dat)

            selections$caplot[] <- FALSE
            selections$biplot <- rep(FALSE, mats$nr)
            selections$curve <- rep(FALSE, mats$nr)
        }
    )

    observeEvent(input$suppressrows, {
        if (!is.null(input$suppressrows)) {
            mi <- isolate(mats$mat_initial)
            sr <- isolate(input$suppressrows)
            mi <- mi[(!(rownames(mi) %in% sr)), ]
            mi <- mi[ , (colSums(mi) != 0)]
            mats$mat <- mi
            mats$caproc <- ca_procrustes_ser(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(isolate(mats$caproc$dat))
            selections$caplot <- rep(FALSE, isolate(mats$nr))
            selections$biplot <- rep(FALSE, isolate(mats$nr))
            selections$curve <- rep(FALSE, isolate(mats$nr))
        }
    })

    observeEvent(input$suppresscols, {
        if (!is.null(input$suppresscols)) {
            mi <- isolate(mats$mat_initial)
            sc <- isolate(input$suppresscols)
            mi <- mi[, (!(colnames(mi) %in% sc)) ]
            mi <- mi[(rowSums(mi) != 0), ]
            mats$mat <- mi
            mats$caproc <- ca_procrustes_ser(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$caproc_ref <- ca_procrustes(isolate(mats$mat), symmetric = isolate(as.logical(input$sym)))
            mats$nr <- nrow(isolate(mats$caproc$dat))
            selections$caplot <- rep(FALSE, isolate(mats$nr))
            selections$biplot <- rep(FALSE, isolate(mats$nr))
            selections$curve <- rep(FALSE, isolate(mats$nr))
        }
    })

    observeEvent(input$end, {
            session$close()
            stopApp()
    })

}



