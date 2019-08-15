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
