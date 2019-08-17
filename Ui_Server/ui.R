# Elements ####
Main_Table <- with(
    tags,
    div(class = "MainTable",
        table(
            tr(th(""),
               td("unlabeled"),
               td("labeled"),
               td("% purity")),
            tr(th("Carbon"),
               td(numInput("carbon", 17, 400)),
               td(numInput("carbon1", 0, 10)),
               td(numInput("carbon2", 100, 100, step = 0.1))),
            tr(th("Hydrogen"),
               td(numInput("hydrogen", 30, 400)),
               td(numInput("hydrogen1", 0, 10)),
               td(numInput("hydrogen2", 100, 100, step = 0.1))),
            tr(th("Nitrogen"),
               td(numInput("nitrogen", 1, 400)),
               td(numInput("nitrogen1", 0, 10)),
               td(numInput("nitrogen2", 100, 100, step = 0.1))),
            tr(th("Oxygen"),
               td(numInput("oxygen", 2, 400)),
               td(""),
               td("")),
            tr(th("Silicon"),
               td(numInput("silicon", 2, 5)),
               td(""),
               td("")),
            tr(th("Iron"),
               td(numInput("iron", 0, 4)),
               td(""),
               td("")),
            tr(th("Sulfur"),
               td(numInput("sulfur", 0, 8)),
               td(""),
               td(""))
        )
    )
)
options <- with(
    tags,
    table(
        class = "OptionsTable",
        tr(th("Display Digits"),
           td(numInput("digits", 3, 6, 1, NA))),
        tr(th("Detection Limit (%)"),
           td(numInput("detect_lim", 0.1 , 1, 0, 0.0001)))
    )
)

MainPanel.conditionalFormula <- conditionalPanel(
    condition = "output.Condition",
    tags$div(
        class = "OutPutClass",
        h3("Chemical Formula"),
        uiOutput("chemical_formula")
    )
)

MainPanel.conditionalPanel <- conditionalPanel(
    condition = "output.Condition",
    tags$div(
        h4("Mass Distribution"),
        tableOutput("stats"),
        downloadButton("download", "Download")
    )
)
MainPanel.Plot <- tags$div(
    class = "OutPutClass",
    imageOutput("DistributionPlot")
)
MainText <- with(
    tags,
    div(
        class = "MainText",
        includeMarkdown("../www/newMark.md")
    )
)
# Combined elements ####
MiddleContent <- with(
    tags,
    div(
        class ="MiddleContent",
        h4("Select type and number of atoms"),
        Main_Table,
        hr(),
        h5("Options"),
        options,
        actionButton(inputId = "go", label = "update")
    )
)
tabs <- tabsetPanel(
    type = "tabs",
    tabPanel("Table", MainPanel.conditionalPanel),
    tabPanel("Plot", MainPanel.Plot)
)
Footer <- with(
    tags,
    p(br(), br(), br(), br(), br(), br(),
      br(), br(), br(), br(), br(), br())
)
# rows ####
row1 <- fluidRow(column(width = 12,  MainText))
row2 <- fluidRow(
    column(width = 2),
    column(width = 4, with(tags, div(br(), h3("Input"), br())), fluidRow(column(12, MiddleContent))),
    column(width = 4, offset = 1, MainPanel.conditionalFormula, fluidRow(column(12, tabs))),
    column(width = 2)
)
row3 <- fluidRow(
    column(width = 12, Footer)
)
# User Interface
ui <- fluidPage(
    HeadCss("../www/mystyle.css"),
    ShinyHeader("Isotopic Distribution of Tracers"),
    row1,
    row2,
    row3
)
