library(shiny)
library(R6)
library(partitions)
ScriptPath = "./Rscripts/"
Scripts = list.files(ScriptPath)
for (file in Scripts) {
    source(paste0(ScriptPath, file))
}
runApp("./Ui_Server")

