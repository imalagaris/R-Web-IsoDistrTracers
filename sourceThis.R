library(shiny)
ScriptPath = "./Rscripts/"
Scripts = list.files(ScriptPath)
for (file in Scripts) {
    source(paste0(ScriptPath, file))
}
runApp("./Ui_Server")

