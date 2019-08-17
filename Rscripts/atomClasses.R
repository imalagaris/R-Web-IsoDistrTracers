# class containing isotope info table
IsotopeData <- R6Class(
    classname = "IsotopeData",
    public = list(
        List = NULL,
        filepath = "../data/iso4.RDS",
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
        reduceProbVec = function(cutoff = 1e-6) {
            id <- which(self$p > cutoff)
            self$p <- self$p[min(id):max(id)]
            self$low <- as.integer(names(self$p)[1])
            self$high <- as.integer(names(self$p)[length(self$p)])
            self$resetLen()
        },
        adjustVec = function() {
            tmp <- private$setNms(numeric(self$high + 1L))
            id <- names(self$p)
            tmp[id] <- self$p
            self$p <- tmp
            tmp[id] <- self$wt
            self$wt <- tmp
        },
        findIndex = function(id, allow = 1e-6) {
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
            nams <- as.character(List$neuEx)
            for (nn in names(List)) {
                self[[nn]] <- setNames(List[[nn]], nams)
            }
            if (!is.null(purity)) {
                self$prob <- c("0" = 1 - purity, "1" = purity)
            }
            self$n <- n
            self$peak<- ceiling(sum(self$prob * self$neuEx * n))
            self$type <- self$type[[1L]]
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
            if (length(vec) < 15L) {
                id <- max(which(vec > 1e-6))
                return(vec[seq(id)])
            } else {
            ProbableVec$new(vec, peak)$one2HighVec()
            }
        },
        atomGroupDist = function() {
            k <- self$findIndex(2L)
            out <- dbinom(0L:k, self$n, self$prob[2L])
            self$p <- private$setNms(out)
            self[["len"]] <- k + 1L
        },
        fft_unit = function() {
            # if (self$n == 1L) {
            #     self$p <- private$setNms(self$prob)
            #     self[["len"]] <- self[["type"]]
            # } else {
                private$atomGroupDist()
            # }
        }
    )
)
AtomDataType3 <- R6Class(
    classname = "AtomDataType3",
    inherit = AtomData,
    public = list(
        findIndex = function(id, allow = 1e-6) {
            low <- qbinom(allow, self$n, self$prob[id])
            high <- qbinom(1 - allow, self$n, self$prob[id])
            c(low, high)
        },
        getWt = function() {
            ranges <- matrix(0L, nrow = 2L, ncol = length(self$neuEx))
            colnames(ranges) <- names(self$neuEx)
            for (neutron in names(self$neuEx)) {
                ranges[, neutron] <- self$findIndex(neutron)
            }
            low <- ranges[1L, ]
            high <- ranges[2L, ]
            blocks <- blockparts(high - low, self$n - sum(low))
            blocks <- sweep(blocks, 1L, low, `+`)
            type <- as.integer(colSums(blocks[-1L, ] * self$neuEx[-1L]))
            p <- apply(blocks, 2L, dmultinom, self$n, self$prob[names(self$neuEx)])
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
            tmp <- lapply(names(self$neuEx), function(x) super$findIndex(x, 1e-5) * self$neuEx[x])
            Reduce("+", tmp)
        },
        atomGroupDist = function() {
            prob <- private$setNms(numeric(self$type))
            prob[names(self$prob)] <- self$prob
            self$prob <- prob
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
                if (tmp[["type"]][[1L]] == 2L) {
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
