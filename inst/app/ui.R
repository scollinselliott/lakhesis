ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Lakhesis Calculator"),  
  shinydashboard::dashboardSidebar(
        shiny::column(12,align = "center",offset = 0,
        shiny::radioButtons("hapax", label = "Selection:", choices = c("Use All Data" = FALSE, "Remove Hapax" = TRUE)),
        shiny::radioButtons("char", label = "Encoding:", choices = c("ISO-8859-1" = "iso-8859-1", "UTF-8" = "utf-8")),
        shiny::fileInput('datafile', 'Choose CSV File:',
                accept=c('csv', 'comma-separated-values','.csv')),
        shiny::radioButtons("sym", label = "CA Scores:", choices = c("Symmetric" = TRUE, "Asymmetric" = FALSE)),
        shiny::actionButton("reinitialize", "Reinitialize", width = "90%"),
        shiny::actionButton("rerun", "Replot with Selection", width = "90%"),
        # shiny::actionButton("log", "Save Seriation as Strand", width = "90%"),

        shiny::br(),
        shiny::downloadButton("downloadData", "Export Data", width = "90%"),
        shiny::br(),shiny::br(),
        shiny::tags$img(src = file.path('img', 'logo.png'), height = "50%", width = "50%", align="center"),
        shiny::br(),shiny::br(),
        shiny::actionButton("end", "End Calculator", width = "90%")
        )
        ),
    shinydashboard::dashboardBody(  
        shiny::fluidRow(
            shinydashboard::tabBox(title = "Seriation Explorer",
                id = "pointselector", height = "500px", width = 5, 
                    shiny::tabPanel("CA Plot",
                    shiny::plotOutput(outputId="caplot", brush = shiny::brushOpts(id="plot_brush_caplot", resetOnNew = TRUE), dblclick = "plot_reset"),
                    #column(12,align = "center",offset = 0,
                    #shiny::actionButton("save_curve", "Gather Plot Points", width = "90%"))
                ),
                shiny::tabPanel("CA-Procrustes Plot",
                    shiny::plotOutput(outputId="biplot", brush = shiny::brushOpts(id="plot_brush_biplot", resetOnNew = TRUE), dblclick = "plot_reset"),
                    #column(12,align = "center",offset = 0,
                    #shiny::actionButton("save_biplot", "Gather Plot Points", width = "90%"))
                ),
                shiny::tabPanel("Curve Plot",
                    shiny::plotOutput(outputId="procrustesplot", brush = shiny::brushOpts(id="plot_brush_curve", resetOnNew = TRUE), dblclick = "plot_reset"),
                    #column(12,align = "center",offset = 0,
                    #shiny::actionButton("save_curve", "Gather Plot Points", width = "90%"))
                )
            ),
            shinydashboard::tabBox(title = "Save Strand",
                id = "pcaplot", height = "500px", width = 2, 
                shiny::tabPanel("Projection",
                shiny::column(12,align = "center",offset = 0,
                p("Select a strand by choosing which axis to use for its seriation."),
                shiny::actionButton("log_ca1", "CA1", width = "90%"),
                shiny::actionButton("log_ca2", "CA2", width = "90%"),
                shiny::actionButton("log_proc1", "Procrustes1", width = "90%"),
                shiny::actionButton("log_proc2", "Procrustes2", width = "90%"),
                shiny::actionButton("log_curve", "Curve", width = "90%"),
                shiny::br(),
                shiny::radioButtons("crits", label = "Criterion:", choices = c("cor_sq" = "cor_sq", "conc_wrc" = "conc_wrc")),
                shiny::actionButton("lakhesize", "Lakhesize (Consensus)", width = "90%"),
                shiny::actionButton("deviancetest", "Run Deviance Test", width = "90%")

                )
                )
            ),
            shinydashboard::tabBox(title = "Consensus Seriation",
                id = "pcaplot", height = "500px", width = 5, 
                shiny::tabPanel("Matrix",
                    shiny::plotOutput(outputId="consensusmatrixplot")
                )
            )
        ),
        shiny::fluidRow(
            shinydashboard::tabBox(title = "Diagnostics",
                id = "results",
                # shiny::tabPanel("Lakhesis Coefficient",   
                #     shiny::plotOutput(outputId="lakhesiscoefplot")
                # ),
                shiny::tabPanel("Agreement",   
                    shiny::plotOutput(outputId="consensusplot")
                ),
                shiny::tabPanel("Criterion",   
                    shiny::plotOutput(outputId="kappaplot")
                ),
                shiny::tabPanel("Deviance (Rows)",   
                    shiny::tableOutput("deviancerows")
                ),
                shiny::tabPanel("Deviance (Columns)",
                    shiny::tableOutput("deviancecols")
                )
            ),
            shinydashboard::box(title = "Modify", 
                selectInput("suppressrows", label = "Select Rows to Suppress", choices = "rowsel", multiple = TRUE),
                selectInput("suppresscols", label = "Select Columns to Suppress", choices = "colsel", multiple = TRUE),
                shiny::p("Currently suppressed row or column values will not be included when saving strands, and are not removed from previously saved strands."),
                shiny::br(),
                shiny::textInput("strandDelete", "Delete Strand (Indices Will Resort Automatically, Undo Will Re-Add as Last Strand)"),
                shiny::actionButton("submitStrandDelete", "Delete"),
                shiny::actionButton("submitStrandDeleteUndo", "Undo")
                )
        )
    )
)

