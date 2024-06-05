ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Lakhesis Calculator"),  
  shinydashboard::dashboardSidebar(
        shiny::column(12,align = "center",offset = 0,
        shiny::fileInput('datafile', 'Choose CSV file',
                accept=c('csv', 'comma-separated-values','.csv')),
        shiny::br(),
        shiny::actionButton("reinitialize", "Reinitialize", width = "90%"),
        shiny::actionButton("rerun", "Recompute with Selection", width = "90%"),
        shiny::actionButton("log", "Save Seriation as Strand", width = "90%"),
        shiny::actionButton("lakhesize", "Lakhesize Strands", width = "90%"),
        shiny::actionButton("deviancetest", "Run Deviance Test", width = "90%"),
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
                id = "pointselector", height = "600px",
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
            shinydashboard::tabBox(title = "Consensus Seriation",
                id = "pcaplot", height = "600px",
                shiny::tabPanel("Row PCA",
                    shiny::plotOutput(outputId="consensusrowplot")
                ),
                shiny::tabPanel("Column PCA",
                    shiny::plotOutput(outputId="consensuscolplot")
                ),
                shiny::tabPanel("Matrix",
                    shiny::plotOutput(outputId="consensusmatrixplot")
                )
            )
        ),
        shiny::fluidRow(
            shinydashboard::tabBox(title = "Criteria",
                id = "results",
                # shiny::tabPanel("Lakhesis Coefficient",   
                #     shiny::plotOutput(outputId="lakhesiscoefplot")
                # ),
                shiny::tabPanel("Agreement",   
                    shiny::plotOutput(outputId="consensusplot")
                ),
                shiny::tabPanel("Concentration",   
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
                shiny::p("Removing row or column values may also entail the automatic suppression of any resulting columns and rows containing all zeros. Saving a new strand or lakhesizing will permanently eliminate suppressed row and column values for all past selection of strands."),
                shiny::br(),
                shiny::textInput("strandDelete", "Delete Strand (Indices Will Resort Automatically, Undo Will Re-Add as Last Strand)"),
                shiny::actionButton("submitStrandDelete", "Delete"),
                shiny::actionButton("submitStrandDeleteUndo", "Undo")
                )
        )
    )
)

