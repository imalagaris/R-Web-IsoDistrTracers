loadPkg <- function(package) {
    packageString <- deparse(substitute(package))
    if (!require(packageString, character.only = TRUE, quietly = TRUE)) {
        install.packages(packageString, dependencies = TRUE)
    }
    library(packageString, character.only = TRUE, quietly = TRUE)
}
loadPkg(shiny)
loadPkg(R6)
library(shiny)
library(R6)
ScriptPath = "./Rscripts/"
Scripts = list.files(ScriptPath)
for (file in Scripts) {
    source(paste0(ScriptPath, file))
}
runApp("./Ui_Server")
