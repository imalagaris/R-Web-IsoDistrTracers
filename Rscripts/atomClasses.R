# class containing isotope info table
IsotopeData <- R6Class(
    classname = "IsotopeData",
    public = list(
        List = NULL,
        filepath = "../data/new_iso3.RDS",
        keepThem = c("type", "int_mass", "mass", "neuEx", "massEx", "prob"),
        initialize = function() {
            dat <- readRDS(self$filepath)
            self[["List"]] <- split(dat[self$keepThem], dat[["E"]])
        }
    )
)
# class returning atom specific object
AtomData <- R6Class(
    classname = "AtomData",
    public = list(
        type = NULL,
        int_mass = NULL,
        mass = NULL,
        massEx = NULL,
        neuEx = NULL,
        peak = NULL,
        prob = NULL,
        n = NULL,
        p = NULL,
        len = NULL,
        wt = NULL,
        low = NULL,
        high = NULL,
        print = function(...) {
            out <- list()
            for (nn in names(self)) {
                tmp <- class(self[[nn]])
                cond <- tmp != "function" && tmp != "environment"
                if (cond) out[[nn]] <- self[[nn]]
            }
            print(out)
        },
        resetLen = function() {
            self$len <- length(self$p)
        },
        reduceProbVec = function() {
            self$p <- self$p[which(self$p > 1e-7)]
            self$low <- as.integer(names(self$p)[1])
            self$high <- as.integer(names(self$p)[length(self$p)])
            self$resetLen()
        },
        findIndex = function(id, allow = 1e-7) {
            qbinom(1 - allow, self$n, self$prob[id])
        },
        getWt = function() {
            nams <- names(self$p)
            n0 <- as.numeric(nams[1L])
            baseWt <- self$n * self$mass[1L] + self$massEx[2L] * n0
            wt <- cumsum(c(baseWt, rep(self$massEx[2L], length(self$p) - 1L)))
            self$wt <- setNames(wt, nams)
        },
        initialize = function(List, n, purity = NULL) {
            for (nn in names(List)) {
                self[[nn]] <- List[[nn]]
            }
            if (!is.null(purity)) {
                self[["prob"]] <- c(1 - purity, purity)
            }
            self[["n"]] <- n
            self[["peak"]] <- ceiling(sum(self[["prob"]] * self[["neuEx"]] * n))
            private$fft_unit()
        }
    ),
    private = list(
        setNms = function(vec, len = NULL) {
            if (is.null(len)) len <- length(vec)
            nam <- as.character(0L:(len - 1L))
            setNames(vec, nam)
        },
        find_index = function(vec, peak) {
            ProbableVec$new(vec, peak)$one2HighVec()
        },
        atomGroupDist = function() {
            k <- self$findIndex(2L)
            out <- dbinom(0L:k, self$n, self$prob[2L])
            self$p <- private$setNms(out)
            self[["len"]] <- k + 1L
        },
        fft_unit = function() {
            if (self$n == 1L) {
                self$p <- private$setNms(self$prob)
                self[["len"]] <- self[["type"]]
            } else {
                private$atomGroupDist()
            }
        }
    )
)
AtomDataType3 <- R6Class(
    classname = "AtomDataType3",
    inherit = AtomData,
    public = list(
        findIndex = function(id, allow = 1e-7) {
            low <- qbinom(allow, self$n, self$prob[id])
            high <- qbinom(1 - allow, self$n, self$prob[id])
            c(low, high)
        },
        getWt = function() {
            ranges <- matrix(0L, nrow = 2L, ncol = 3L)
            for (i in seq(3L)) {
                ranges[, i] <- self$findIndex(i)
            }
            low <- ranges[1L, ]
            high <- ranges[2L, ]
            blocks <- blockparts(high - low, self$n - sum(low))
            blocks <- sweep(blocks, 1L, low, `+`)
            type <- as.integer(colSums(blocks[-1, ] * c(1L, 2L)))
            p <- apply(blocks, 2L, dmultinom, self$n, self$prob)
            pNorm <- tapply(p, type, function(x) x / sum(x))
            wt <- split(t(self$mass) %*% blocks, type)
            out <- setNames(vector("list", length(wt)), names(wt))
            for (nn in names(out)) {
                out[[nn]] <- sum(wt[[nn]] * pNorm[[nn]])
            }
            self$wt <- unlist(out, use.names = TRUE)
        }
    ),
    private = list(
        findVecLength = function() {
            super$findIndex(2L, 1e-6) + 2L * super$findIndex(3L, 1e-6)
        },
        atomGroupDist = function() {
            len <- private$findVecLength()
            out <- numeric(len)
            out[seq.int(self$type)] <- self$prob
            out <- fft(out)^self$n / len
            out <- Re(fft(out, inverse = TRUE))
            out <- private$find_index(out, self$peak)
            self$p <- private$setNms(out)
            self$resetLen()
        }
    )
)

# class that takes input chemical formula and initializes atom objects within
atomDataList <- R6Class(
    classname = "atomDataList",
    inherit = AtomData,
    public = list(
        elem = NULL,
        nElem = NULL,
        formula = NULL,
        atoms = NULL,
        initialize = function(elem, n, purity = NULL) {
            private$List <- IsotopeData$new()$List[elem]
            n <- setNames(as.integer(n), elem)
            self[["formula"]] <- n
            self[["elem"]] <- elem
            self[["nElem"]] <- length(elem)
            if (is.null(purity)) purity <- vector("list", self$nElem)
            names(purity) <- elem
            self[["atoms"]] <- setNames(vector("list", self$nElem), elem)
            for (el in elem) {
                tmp <- as.list(private$List[[el]])
                tmp[["type"]] <- tmp[["type"]][[1L]]
                if (tmp[["type"]] == 2L) {
                    self[["atoms"]][[el]] <- AtomData$new(tmp, n[[el]], purity[[el]])
                } else {
                    self[["atoms"]][[el]] <- AtomDataType3$new(tmp, n[[el]], purity[[el]])
                }
            }
        }
    ),
    private = list(
        List = list()
    )
)
