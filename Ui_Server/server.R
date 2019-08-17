server <- function(input, output) {
    arguments <- eventReactive(input$go, { extract_atoms(input) })
    tracers <- eventReactive(input$go, { extract_atoms(input, "labeled") })
    output$chemical_formula <- renderUI({
        temp <-  prepareDataChemFormula(arguments(), tracers())
        chem <- createChemString(temp)
        withMathJax(h2(chem))
    })
    decimals <- eventReactive(input$go, {input$digits})
    detect   <- eventReactive(input$go, {input$detect_lim * 10^-2})
    isotope_table <- reactive(wrapperClass$new(arguments(), tracers())$result)
    output$stats <- renderTable({isotope_table()},
        digits =  decimals,
        rownames = TRUE,
        spacing = "xs",
        display = c("s", "f", "f", "f"))
    output$download <- downloadHandler(
        filename = "isotope_table.txt",
        content = function(fname) {
            write.table(isotope_table(), fname, sep = ",", col.names = NA)
       })
    output$Condition <- reactive({is.matrix(isotope_table())})
    outputOptions(output, "Condition", suspendWhenHidden = FALSE)
    NumberOfPeaks <- eventReactive(input$go, {nrow(isotope_table())})
    output$DistributionPlot <- renderImage({
        size <- image_size(NumberOfPeaks())
        ww <- size$ww
        hh <- size$hh
        outfile <- tempfile(fileext = '.png')
        png(outfile, width = ww, height = hh)
        par(mar = c(5, 8, 2, 6))
        plot(isotope_table()[,1], isotope_table()[, 2], yaxt = "n", ylab = "", bty = "n",
            type = "h",
            xlab = "", lwd = 2)
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
