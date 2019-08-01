library(shiny)
source("~/Google Drive/R/shiny/thesis/helpers.R")
isotope <- readRDS(file = "./thesis/data/isotope.RDS")
dat <- readRDS(file = "./thesis/data/dat.RDS")

ui <- fluidPage(style = "padding: 5%;",
 with(tags,
     head(style(HTML('
table {
    font-family: arial;
    border-style: 1px solid black;
    vertical-align: middle;
}
th {
    padding: 0px;
    border-bottom: 5px double;
    text-align: center;
    font-size: 80%;
    width: auto;
    vertical-align: middle;
}
 
th:nth-child(1) {
    text-align: left;
}
tr {
    padding: 0px;
    height: 30px;
    vertical-align: middle;
}
td {
    padding: 0px;
    margin: 0px;
    vertical-align: middle;
}

hr {
    margin-top: 0.5em;
    margin-bottom: 0.5em;
    margin-left: auto;
    margin-right: auto;
    border-style: inset;
    border-width: 3px;
}

.aek div.form-group.shiny-input-container {
    padding: 1%;
    margin: 0 auto;
    width: 75%;
    vertical-align: middle;
}
.aek .form-control {
    padding: 0px;
    text-align: center;
    font-weight: bold;
    width: 100%;
    vertical-align: middle;
}

div.shiny-html-output {
width: 20%;
font-size: 80%;
}

#digits.shiny-input-container  {
margin: 0px;
}

#digits.form-control  {
text-align: center;
font-weight: bold;
}

#detect_lim.form-control  {
text-align: center;
font-weight: bold;

}

label[for = "digits"] {
font-size: 80%;
text-align: center;
}

label[for = "detect_lim"] {
font-size: 80%;
}


div.shiny-html-output  tr{
height: auto;
}
')
))),
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
                    td(numericInput(inputId = "carbon",    label = NULL, 
                        value = 17,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "carbon1",   label = NULL, 
                        value = 0,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "carbon2",   label = NULL, 
                        value = 100, min = 0, max = 100, step = 0.1))),
                tr(td("Hydrogen"),
                    td(numericInput(inputId = "hydrogen",  label = NULL, 
                        value = 30,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "hydrogen1", label = NULL, 
                        value = 0,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "hydrogen2", label = NULL, 
                        value = 100, min = 0, max = 100, step = 0.1))),
                tr(td("Nitrogen"),
                    td(numericInput(inputId = "nitrogen",  label = NULL, 
                        value = 1,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "nitrogen1", label = NULL, 
                        value = 0,   min = 0, max = 1e4)),
                    td(numericInput(inputId = "nitrogen2", label = NULL, 
                        value = 100, min = 0, max = 100, step = 0.1))),
                tr(td("Oxygen"),
                    td(numericInput(inputId = "oxygen",    label = NULL, 
                        value = 2,   min = 0, max = 1e4)),
                    td(""),
                    td("")),
                tr(td("Silicon"),
                    td(numericInput(inputId = "silicon",   label = NULL, 
                        value = 2,   min = 0, max = 1e4)),
                    td(""),
                    td("")),
                tr(td("Iron"),
                    td(numericInput(inputId = "iron",   label = NULL, 
                        value = 0,   min = 0, max = 1e4)),
                    td(""),
                    td("")),
                tr(td("Sulfur"),
                    td(numericInput(inputId = "sulfur",   label = NULL, 
                        value = 0,   min = 0, max = 1e4)),
                    td(""),
                    td(""))                
                )
                ),
            hr(),
            br(),
            helpText("Options"),
            fluidRow(
                column(5,
            numericInput(inputId = "digits", label = "Display Digits", width = "100%", value = 4, min = 1, max = 7)),
                column(5, offset = 1, 
            numericInput(inputId = "detect_lim", label = "Detection Limit (%)", width = "100%", value = 0.1, min = 0, max = 1, step = 0.0001))),
            br(),
            actionButton(inputId = "go", label = "update")),

    
    mainPanel(
        fluidRow(
            column(4, style = "padding-left: 5%; padding-right: 5%; background-color: lightgray;",
                conditionalPanel(condition = "output.a",
                helpText("Chemical Formula"),
                hr(),
                uiOutput("chemical_formula"),
                hr(),
                helpText("Mass Distribution"),
                tableOutput("stats"),
                hr(),
                downloadButton("download", "Download"))),
            column(8, 
                imageOutput("plot1"))))))
    


server <- function(input, output) {
    arguments <- eventReactive(input$go, {
        temp <-  list(
            "C"  = input$carbon, 
            "H"  = input$hydrogen, 
            "N"  = input$nitrogen,
            "O"  = input$oxygen,
            "Si" = input$silicon,
            "Fe" = input$iron,
            "S"  = input$sulfur)
        sub_iso_list(temp)
        })
    
    tracers <- eventReactive(input$go, {
        temp <- list(
            "C" = c(input$carbon1,   input$carbon2   / 100),
            "H" = c(input$hydrogen1, input$hydrogen2 / 100),
            "N" = c(input$nitrogen1, input$nitrogen2 / 100))
        sub_iso_list(temp)
        })
    
    output$chemical_formula <- renderUI(
        {
            nn_arg      <- names(arguments())
            nn_tr       <- names(tracers())
            n_atoms_arg <- unlist(lapply(arguments(), '[', 1), use.names = FALSE)
            n_atoms_tr  <- unlist(lapply(tracers(),   '[', 1), use.names = FALSE)
            temp        <- data.frame(
                    "atoms" = c("C", "C", "H", "H", "N", "N", "O", "Si", "Fe", "S"), 
                    "supscr" = c("{^{13}","", "{^{2}", "", "{^{15}", "", "", "", "", ""),
                    "type" = c(1, 1, 2, 2, 3, 3, 4, 5, 6, 7),
                    "ind_tr" = c(1, 2, 1, 2, 1, 2, 2, 2, 2, 2))
            index_tr    <- with(temp, which((ind_tr == 1) & (atoms %in% nn_tr)))
            index_arg   <- with(temp, which((ind_tr == 2) & (atoms %in% nn_arg)))
            temp_tr <- temp[index_tr, ]
            temp_tr[["n_atoms"]] <- n_atoms_tr
            temp_arg <- temp[index_arg, ]
            temp_arg[["n_atoms"]] <- n_atoms_arg
            
            temp     <- rbind(temp_tr, temp_arg)
            temp <- temp[order(temp$type, temp$ind_tr), ]
            
            temp[["chemm"]] <- NA
            for (i in 1:nrow(temp)) {
                if (temp$ind_tr[i] == 1) {
                    temp[i, "chemm"] <- paste0("(", temp$supscr[i], temp$atoms[i], "_", "{",temp$n_atoms[i], "}", ")", "}", collapse = '')
                } else {
                    nnn <- ""
                    if (temp$n_atoms[i] > 1) nnn <- temp$n_atoms[i]
                    temp[i, "chemm"] <- paste0(temp$atoms[i], "_", "{",nnn, "}", collapse = '')
                }
            }
            chem <- paste0("$$", paste0(temp[["chemm"]], collapse = ''), "$$", collapse = '')
            withMathJax(h5(chem))
        })
            
    
    decimals <- eventReactive(input$go, {input$digits})
    detect   <- eventReactive(input$go, {input$detect_lim * 10^-2})
    
    isotope_table <- reactive({myfft1(arguments(), tracers(), detect())})
    output$stats <- renderTable({isotope_table()},
        digits =  decimals, rownames = TRUE, spacing = "xs", display = c("s", "d", "f", "f"))
    
    output$download <- downloadHandler(
    filename = "isotope_table.txt", 
    content = function(fname) { 
        write.table(isotope_table(), fname, sep = "\t", col.names = NA)
    }
  )
    output$a <- reactive({is.matrix(isotope_table())})
    outputOptions(output, "a", suspendWhenHidden = FALSE)
    aek <- eventReactive(input$go, {nrow(isotope_table())})
    # output$plot1 <- renderPlot({
    #     par(mar = c(5, 8, 2, 6))
    #     plot(isotope_table()[, 2], yaxt = "n", ylab = "", bty = "n",
    #         type = "h",
    #         xaxt = "n", xlab = "", lwd = 4)
    #     axis(1, at = 1:nrow(isotope_table()), labels = row.names(isotope_table()))
    #     axis(2, las = 1,  cex.axis = 1.5)
    #     mtext(text = "Probability", side = 2, line = 5, cex = 2, font = 2)
    #     par(new = TRUE)
    #     plot(isotope_table()[, 3], las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", dpi = 300,
    #         pch = "", bty = "n")
    #     axis(4, las = 1, cex.axis = 1.5)
    #     mtext(text = "Normalized", side = 4, line = 5, cex = 2, font = 2)
    #     }) 
    output$plot1 <- renderImage({
        if (aek() < 5) {
            hh <- 400;
            ww <- 400;
        } else if (aek() < 15) {
            hh <- 500;
            ww <- 500;
        } else if (aek() < 30) {
            ww <- 600;
            hh <- 600;
        } else if (aek() < 50) {
            ww <- 700;
            hh <- 700;
        } else if (aek() < 60) {
            ww <- 800;
            hh <- 800;
        } else {
            ww <- 1000
            hh <- 800
        }
        
        outfile <- tempfile(fileext = '.png')
        
        png(outfile, width = ww, height = hh)
        par(mar = c(5, 8, 2, 6))
        plot(isotope_table()[, 2], yaxt = "n", ylab = "", bty = "n",
            type = "h",
            xaxt = "n", xlab = "", lwd = 2)
        axis(1, at = 1:nrow(isotope_table()), labels = row.names(isotope_table()))
        axis(2, las = 1,  cex.axis = 1.5)
        mtext(text = "Probability", side = 2, line = 5, cex = 2, font = 2)
        par(new = TRUE)
        plot(isotope_table()[, 3], las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
            pch = "", bty = "n")
        axis(4, las = 1, cex.axis = 1.5)
        mtext(text = "Normalized", side = 4, line = 5, cex = 2, font = 2)
        dev.off()
        
         list(src = outfile,
         contentType = 'image/png',
         width = ww,
         height = hh,
         alt = "This is alternate text")
        
        }, deleteFile = TRUE) 
    }
    
shinyApp(ui = ui, server = server)


















