mjType <- "text/javascript"
mjScript <- "http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
chemDF <- data.frame(
    "atoms" = c("C", "C", "H", "H", "N", "N", "O", "Si", "Fe", "S"),
    "supscr" = c("{^{13}", "", "{^{2}", "", "{^{15}", "", "", "", "", ""),
    "type" = c(1, 1, 2, 2, 3, 3, 4, 5, 6, 7),
    "ind_tr" = c(1, 2, 1, 2, 1, 2, 2, 2, 2, 2)
)
atoms <- list()
atoms[["unlabeled"]] <- c(
    "C" = "carbon",
    "H" = "hydrogen",
    "N" = "nitrogen",
    "O" = "oxygen",
    "Si" = "silicon",
    "Fe" = "iron",
    "S" = "sulfur"
)
atoms[["labeled"]] <- c(
    "C1" = "carbon1",
    "C2" = "carbon2",
    "H1" = "hydrogen1",
    "H2" = "hydrogen2",
    "N1" = "nitrogen1",
    "N2" = "nitrogen2"
)

allowedRange <- list()
allowedRange[["unlabeled"]] <- list(
    "C" = c("low"= 0, "high" = 400),
    "H" = c("low"= 0, "high" = 400),
    "N" = c("low"= 0, "high" = 400),
    "O" = c("low"= 0, "high" = 400),
    "Si" = c("low"= 0, "high" = 5),
    "Fe" = c("low"= 0, "high" = 4),
    "S" = c("low"= 0, "high" = 8)
)
allowedRange[["labeled"]] <- list(
    "C1" = c("low"= 0, "high" = 10),
    "C2" = c("low"= 0, "high" = 100),
    "H1" = c("low"= 0, "high" = 10),
    "H2" = c("low"= 0, "high" = 100),
    "N1" = c("low"= 0, "high" = 10),
    "N2" = c("low"= 0, "high" = 100)
)
