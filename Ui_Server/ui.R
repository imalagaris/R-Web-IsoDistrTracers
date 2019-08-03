docHead <- tags$head(includeCSS("../www/mystyle.css"))

Main_Table <- with(tags,
    table(class = "MainTable",
        tr(th(""),
            td("unlabeled"),
            td("labeled"),
            td("% purity")),
        tr(th("Carbon"),
            td(numInput("carbon", 17, 1e4)),
            td(numInput("carbon1", 0, 1e4)),
            td(numInput("carbon2", 100, 100, step = 0.1))),
        tr(th("Hydrogen"),
            td(numInput("hydrogen", 30, 1e4)),
            td(numInput("hydrogen1", 0, 1e4)),
            td(numInput("hydrogen2", 100, 100, step = 0.1))),
        tr(th("Nitrogen"),
            td(numInput("nitrogen", 1, 1e4)),
            td(numInput("nitrogen1", 0, 1e4)),
            td(numInput("nitrogen2", 100, 100, step = 0.1))),
        tr(th("Oxygen"),
            td(numInput("oxygen", 2, 1e4)),
            td(""),
            td("")),
        tr(th("Silicon"),
            td(numInput("silicon", 2, 1e4)),
            td(""),
            td("")),
        tr(th("Iron"),
            td(numInput("iron", 0, 1e4)),
            td(""),
            td("")),
        tr(th("Sulfur"),
            td(numInput("sulfur", 0, 1e4)),
            td(""),
            td(""))
    )
)

options <- with(
    tags,
    table(
        class = "OptionsTable",
        tr(
            th("Display Digits"),
            td(numInput("digits", 4, 7, 1, NA))),
        tr(
            th("Detection Limit (%)"),
            td(numInput("detect_lim", 0.1 , 1, 0, 0.0001)))
    )
)

MainPanel.conditionalPanel <- conditionalPanel(
    condition = "output.Condition",
    tags$div(
        class = "OutPutClass",
        h3("Chemical Formula"),
        uiOutput("chemical_formula"),
        hr(),
        h4("Mass Distribution"),
        tableOutput("stats"),
        downloadButton("download", "Download")
    )
)

MainPanel.Plot <- tags$div(
    class = "OutPutClass",
    imageOutput("DistributionPlot")
)

ui <- fluidPage(
    docHead,
    titlePanel("Isotopic distribution of tracer molecules"),
    sidebarLayout(

        sidebarPanel(
            h3("Select type and number of atoms"),
            Main_Table,
            hr(),
            h4("Options"),
            options,
            actionButton(inputId = "go", label = "update")
        ), # end sidebarPanel

        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Table", MainPanel.conditionalPanel),
                tabPanel("Plot", MainPanel.Plot )
            ) # end fluidRow
        ) # end fluidRow
    ) # sidebarLayout
) # end fluidPage
