 install.packages("shinyjs")

 source('global.R')

ui <- fluidPage(theme = shinytheme("darkly"),
  navbarPage("Single Cell RNA Seq Analyzer",
    tabPanel("Home Page",
    includeMarkdown("./markdown/about.md")),
    tabPanel("scRNAseq Analyzer",
      fileInput("file", "Upload File (.rds)", multiple = FALSE, accept = c('.rds')),
      actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
      actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
      tabsetPanel(
        tabPanel("Instructions", includeMarkdown("./markdown/app-instructions.md"))
      )
    )
  ),
  dashboardHeader(title = " "),
  dashboardSidebar(),
  dashboardBody()
)

server <- function(input, output, session) {
    options(shiny.maxRequestSize=300*1024^2)

values <- reactiveValues()

    shinyjs::disable("run")

    observe({
    if(is.null(input$file) != TRUE) {
        shinyjs::enable("run")
    } else {
        shinyjs::disable("run")
    }
    })
 observeEvent(input$reset, {
        shinyjs::disable("file")
        shinyjs::disable("run")
 })

observeEvent(input$run, {
        shinyjs::disable("run")

show_modal_progress_line(text = "Preparing...") # show the modal window
update_modal_progress(0.2) # update progress bar value

 obj <- load_seurat_obj(input$file$datapath)
        if (is.vector(obj)){
            showModal(modalDialog(
                title = "Error with file",
                HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                    paste(unlist(obj), collapse = "<br><br>"))
            ))
             shinyjs::enable("run")

        } else { 

            output$umap <- renderPlot({
                if (!is.null(input$metadata_col)) {
                    create_metadata_UMAP(obj, input$metadata_col)
                }

            })

                     output$tSNE <- renderPlot({
                if (!is.null(input$gene)) {
                    create_tSNE_plot(obj, input$metadata_col)
                }
            })

            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "UMAP",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'umap'),
                        downloadButton("download_umap", "Download UMAP")
                    ),
                    column(
                        width = 4,
                        selectizeInput("metadata_col", 
                            "Metadata Column", 
                            colnames(obj@meta.data)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                ),
                select = TRUE

                )
                        
                  insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "tSNE",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'tSNE'),
                        
                    ),
                    column(
                        width = 4,
                        selectizeInput("metadata_col", 
                            "Metadata Column", 
                            colnames(obj@meta.data)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                ),
                )
                
                remove_modal_progress_line()
                shinyjs::enable("run")
        }

})
observeEvent(input$reset, {
        shinyjs::reset("file")
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "tSNE")
        shinyjs::disable("run")
    })

}


shinyApp(ui, server)
