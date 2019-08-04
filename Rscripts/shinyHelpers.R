HeadCss <- function(path) {
    tags$head(includeCSS(path), tags$script(type = mjType, src=mjScript))
}

ShinyHeader <- function(Title) {
    fluidRow(
        column(12, tags$div(id = "MainTitle", h1(Title)))
    )
}

numInput = function(inputId, value, max, min = 0, step = NA, label = NULL, width = NULL) {
    allArgs = formals(numInput)
    thisCall = as.list(match.call())[-1]
    allArgs[names(thisCall)] <-  thisCall
    do.call(numericInput, allArgs)
}

sub_iso_list <- function(alist) {
    len <- length(alist)
    index <- rep(FALSE, len)
    for (i in 1:len) {
        if (alist[[i]][1] != 0 & !is.na(alist[[i]][1])) index[i] <- TRUE
    }
    out <- alist[index]
    if (length(out) == 0) {
        return(NULL)
    } else {
        return(out)
    }
}

extract_atoms <- function(input, type = "unlabeled") {
    out = list();
    tmp = atoms[[type]]
    for (nam in names(tmp)) {
        out[[nam]] = input[[tmp[nam]]]
    }
    if (type == "labeled") {
        out[["C"]] = unlist(out[c("C1", "C2")], use.names = F) * c(1, 0.01)
        out[["H"]] = unlist(out[c("H1", "H2")], use.names = F) * c(1, 0.01)
        out[["N"]] = unlist(out[c("N1", "N2")], use.names = F) * c(1, 0.01)
        out = out[c("C", "H", "N")]
    }
    return(sub_iso_list(out))
}

prepareDataChemFormula <- function(arguments, tracers) {
    nn_arg      <- names(arguments)
    nn_tr       <- names(tracers)
    n_atoms_arg <- unlist(lapply(arguments, '[', 1), use.names = FALSE)
    n_atoms_tr  <- unlist(lapply(tracers,   '[', 1), use.names = FALSE)
    index_tr    <- with(chemDF, which((ind_tr == 1) & (atoms %in% nn_tr)))
    index_arg   <- with(chemDF, which((ind_tr == 2) & (atoms %in% nn_arg)))
    temp_tr <- chemDF[index_tr, ]
    temp_tr[["n_atoms"]] <- n_atoms_tr
    temp_arg <- chemDF[index_arg, ]
    temp_arg[["n_atoms"]] <- n_atoms_arg
    temp     <- rbind(temp_tr, temp_arg)
    temp <- temp[order(temp$type, temp$ind_tr), ]
    temp[["chemm"]] <- NA
    return(temp)
}

createChemString <- function(data) {
    for (i in 1:nrow(data)) {
        if (data$ind_tr[i] == 1) {
            data[i, "chemm"] <- paste0("(", data$supscr[i], data$atoms[i], "_", "{",data$n_atoms[i], "}", ")", "}", collapse = '')
        } else {
            nnn <- ""
            if (data$n_atoms[i] > 1) nnn <- data$n_atoms[i]
            data[i, "chemm"] <- paste0(data$atoms[i], "_", "{",nnn, "}", collapse = '')
        }
    }
    paste0("$$", paste0(data[["chemm"]], collapse = ''), "$$", collapse = '')
}

image_size <- function(x) {
    if (x < 5) {
        hh <- 400;
        ww <- 400;
    } else if (x < 15) {
        hh <- 500;
        ww <- 500;
    } else if (x < 30) {
        ww <- 600;
        hh <- 600;
    } else if (x < 50) {
        ww <- 700;
        hh <- 700;
    } else if (x < 60) {
        ww <- 800;
        hh <- 800;
    } else {
        ww <- 1000
        hh <- 800
    }
    list("hh" = hh, "ww" = ww)
}
