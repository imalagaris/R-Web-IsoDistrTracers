library(shiny)
isotope <- readRDS(file = "./data/isotope.RDS")
dat <- readRDS(file = "./data/dat.RDS")
dat_wt <- dat[dat$atoms != "", c("atoms", "weight")]

chemDF <- data.frame(
    "atoms" = c("C", "C", "H", "H", "N", "N", "O", "Si", "Fe", "S"),
    "supscr" = c("{^{13}","", "{^{2}", "", "{^{15}", "", "", "", "", ""),
    "type" = c(1, 1, 2, 2, 3, 3, 4, 5, 6, 7),
    "ind_tr" = c(1, 2, 1, 2, 1, 2, 2, 2, 2, 2))

atoms = list();
atoms[["unlabeled"]] = c(
    "C"  = "carbon",
    "H"  = "hydrogen",
    "N"  = "nitrogen",
    "O"  = "oxygen",
    "Si" = "silicon",
    "Fe" = "iron",
    "S"  = "sulfur"
)
atoms[["labeled"]] = c(
    "C1" = "carbon1",
    "C2" = "carbon2",
    "H1" = "hydrogen1",
    "H2" = "hydrogen2",
    "N1" = "nitrogen1",
    "N2" = "nitrogen2"
)

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

numInput = function(inputId, value, max, min = 0, step = NA, label = NULL, width = NULL) {
    allArgs = formals(numInput)
    thisCall = as.list(match.call())[-1]
    allArgs[names(thisCall)] <-  thisCall
    do.call(numericInput, allArgs)
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

fft_unit <- function(x, n_atoms) {
    if (n_atoms == 1) return(x)
    find_len <- function(n, k) {n * (k - 1L) + 1}
    n_species <- max(which(x > 0))
    x1 <- x[1:n_species]
    len <- find_len(n_atoms, n_species)
    x <- numeric(len)
    x[1:n_species] <- x1
    x <- fft(x)
    x <- x^n_atoms
    final <- Re(fft(x / len, inverse = TRUE))
    return(final)
}

myfft1 <- function(arguments, tracers = NULL, detect_limit = 0.1 * 10^-2) {

    tracers_exist <- !is.null(tracers)
    if (tracers_exist) {
        nn_tracer <- names(tracers)
        n_tracers <- length(tracers)
        prob_tracer <- list()
        n_atoms_tr <- numeric(n_tracers)
        for (i in 1L:n_tracers) {
            n_atoms_tr[i] = k <- tracers[[i]][1]
            p <- tracers[[i]][2]
            prob_tracer[[i]] <- dbinom(0:k, k, p)
        }
        ind_wt <- dat_wt[["atoms"]] %in% nn_tracer
        base_wt_tracer <- sum(dat_wt[["weight"]][ind_wt] * n_atoms_tr)
    } else {
        prob_tracer <- NULL
        base_wt_tracer <- NULL
    }

    arguments_exist <- !is.null(arguments)
    if (arguments_exist) {
        nn_unlab <- names(arguments)
        n_arg <- length(arguments)
        probs_unlab <- lapply(nn_unlab, function(x) eval(isotope[[x]]))
        index <- lapply(seq_len(n_arg), function(i) max(which(probs_unlab[[i]] > 0)))
        n_species <- lapply(index, sum)
        n_atoms_unlab <- numeric(n_arg)
        for (i in seq_len(n_arg)) {
            n_atoms_unlab[i] = k <- arguments[[i]][1]
            if (n_species[[i]] == 2) {
                probs_unlab[[i]] <- dbinom(0L:k, k, probs_unlab[[i]][2])
            } else if (n_species[[i]] > 2) {
                probs_unlab[[i]] <- fft_unit(probs_unlab[[i]], k)
            }
        }
        ind_wt <- dat_wt[["atoms"]] %in% nn_unlab
        base_wt <- sum(dat_wt[["weight"]][ind_wt] * n_atoms_unlab)
    } else {
        probs_unlab <- NULL
        base_wt <- NULL
    }

    probs <- c(probs_unlab, prob_tracer)
    base_wt <- sum(c(base_wt, base_wt_tracer))

    list_len <- length(probs)
    # if (length(probs) == 1) return(probs[[1]])
    list_len <- length(probs)
    if (list_len == 0) {
        stupid <- data.frame("Mass" = NA, "Probability" = NA, "Normalized" = NA)
        rownames(stupid) <- "M0"
        return(stupid)
    } else if (list_len == 1) {
        probs <- unlist(probs, use.names = FALSE)
    } else if (list_len > 1) {
        slots <- sum(unlist(lapply(probs, function(x) max(which(x > 1e-16)))))
        current_len <- lapply(probs, length)
        for (i in 1:length(probs)) {
            if (current_len[[i]] > slots) {
                probs[[i]] <- probs[[i]][seq_len(slots)]
            } else {
                probs[[i]] <- c(probs[[i]], rep(0, slots - current_len[[i]]))
            }
        }
        probs <- lapply(probs, fft)
        for (i in 2:length(probs)) {
            probs[[1]] <- probs[[1]] * probs[[i]]
        }
        probs <- Re(fft(probs[[1]] / slots, inverse = TRUE))
    }
    mass_number <- base_wt + 0:(length(probs) - 1)
    final <- matrix(c(mass_number, probs), ncol = 2,
        dimnames = list(c(paste0("M",0:(length(probs) - 1))), c("Mass", "Probability")))
    final[final[,2] < 3e-16, 2] <- 0
    highest_peak <- max(final[, 2])
    limit <- highest_peak * detect_limit
    index_min <- min(which(final[, 2] > limit))
    index_max <- max(which(final[, 2] > limit))
    index <- index_min:index_max
    final <- final[index, , drop = FALSE]
    Normalized <- final[, 2] / max(final[, 2]) * 100
    final <- cbind(final, Normalized)
    return(final)
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
