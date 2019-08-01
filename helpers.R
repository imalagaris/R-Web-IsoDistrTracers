isotope <- readRDS(file = "./thesis/data/isotope.RDS")
dat <- readRDS(file = "./thesis/data/dat.RDS")
dat_wt <- dat[dat$atoms != "", c("atoms", "weight")]

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