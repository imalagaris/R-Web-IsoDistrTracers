ui <- fluidPage(
    with(tags, head(link(rel = "stylesheet", type = "text/css", href = "style.css"))),
    sidebarLayout(
        sidebarPanel(width = 3, style = "background-color: lightgray;",
            helpText("Select type and number of unlabeled atoms"),
            with(tags,
                table(class = "aek",
                    tr(th("Atoms"),
                        th("# unlabeled"),
                        th("# labeled"),
                        th("% isotopic purity")),
                    tr(td("Carbon"),
                        td(numInput("carbon", 17, 1e4)),
                        td(numInput("carbon1", 0, 1e4)),
                        td(numInput("carbon2", 100, 100, step = 0.1))),
                    tr(td("Hydrogen"),
                        td(numInput("hydrogen", 30, 1e4)),
                        td(numInput("hydrogen1", 0, 1e4)),
                        td(numInput("hydrogen2", 100, 100, step = 0.1))),
                    tr(td("Nitrogen"),
                        td(numInput("nitrogen", 1, 1e4)),
                        td(numInput("nitrogen1", 0, 1e4)),
                        td(numInput("nitrogen2", 100, 100, step = 0.1))),
                    tr(td("Oxygen"),
                        td(numInput("oxygen", 2, 1e4)),
                        td(""),
                        td("")),
                    tr(td("Silicon"),
                        td(numInput("silicon", 2, 1e4)),
                        td(""),
                        td("")),
                    tr(td("Iron"),
                        td(numInput("iron", 0, 1e4)),
                        td(""),
                        td("")),
                    tr(td("Sulfur"),
                        td(numInput("sulfur", 0, 1e4)),
                        td(""),
                        td(""))                
                ) # end table
            ), # end with
            hr(),
            br(),
            helpText("Options"),
            fluidRow(
                column(
                    width = 5,
                    numericInput("digits", "Display Digits", 4, 1, 7, NA, "100%")
                ), #end column
                column(
                    width = 5, 
                    offset = 1, 
                    numericInput("detect_lim", "Detection Limit (%)", 0.1 , 0, 1, 0.0001, "100%")
                ) # end column
            ), # end fluidRow
            br(),
            actionButton(inputId = "go", label = "update")
        ), # end sidebarPanel
        
        mainPanel(
            fluidRow(
                column(
                    width = 4, 
                    style = "padding-left: 5%; padding-right: 5%; background-color: lightgray;",
                    conditionalPanel(
                        condition = "output.a",
                        helpText("Chemical Formula"),
                        hr(),
                        uiOutput("chemical_formula"),
                        hr(),
                        helpText("Mass Distribution"),
                        tableOutput("stats"),
                        hr(),
                        downloadButton("download", "Download")
                    ) # end conditionalPanel
                ), # end column
                column(
                    width = 8, 
                    imageOutput("plot1")
                ) # end column
            ) # end fluidRow
        ) # end main Panel
    ) # sidebarLayout
) # end fluidPage