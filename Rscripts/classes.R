require(R6)
IsotopeDat <- R6Class(
    classname = "IsotopeData",
    public = list(
        List = NULL,
        filepath = "../data/iso6.RDS",
        keepThem = c("type", "int_mass", "mass", "neuEx", "massEx", "prob"),
        initialize = function() {
            dat <- readRDS(self$filepath)
            self[["List"]] <- split(dat[self$keepThem], dat[["E"]])
        }
    )
)

atom <- R6Class(
    classname = "atom",
    public = list(
        n = NULL,
        prob = NULL,
        mass = NULL,
        massEx = NULL,
        type = NULL,
        meanWt = 0,
        baseWt = 0,
        print = function(...) {
            out <- list()
            for (nn in names(self)) {
                tmp <- class(self[[nn]])
                cond <- tmp != "function" && tmp != "environment"
                if (cond) out[[nn]] <- self[[nn]]
            }
            print(out)
        },
        findIndex = function(allow = 1e-15) {
            index <- 1L:self$type
            tmp <- sapply(
                index,
                function(id) {
                    qbinom(1 - allow, self$n, self$prob[id])
                }
            )
            max((index - 1L) * tmp) + 1L
        },
        getPwrOf2 = function(int) {
            int <- as.integer(int)
            expon <- as.integer(ceiling(log(int, 2L)))
            2L^expon
        },
        adjustVec = function(int) {
            tmp <- complex(int)
            tmp[1L:self$type] <- self$prob * exp(2i * pi * self$massEx)
            fft(tmp)^self$n
        },
        initialize = function(dat, n, purity = NULL) {
            self$n = n
            if (!is.null(purity)) {
                self$prob <- c("0" = 1 - purity, "1" = purity)
            } else {
                self$prob <- dat$prob
            }
            self$mass <- dat$mass
            self$massEx <- dat$massEx
            self$type <- dat$type[[1L]]
            self$baseWt <- self$mass[[1L]] * n
            self$meanWt <- sum(self$mass * self$prob) * n
        }
    )
)

atomList <- R6Class(
    classname = "atomList",
    inherit = atom,
    public = list(
        elem = NULL,
        atoms = NULL,
        baseWt = 0,
        meanWt = 0,
        initialize = function(elem, n, purity = NULL) {
            private$List <- IsotopeDat$new()$List[elem]
            n <- setNames(as.integer(n), elem)
            self[["elem"]] <- elem
            len <- length(elem)
            if (is.null(purity)) purity <- vector("list", len)
            names(purity) <- elem
            self[["atoms"]] <- setNames(vector("list", len), elem)
            for (el in elem) {
                tmp <- as.list(private$List[[el]])
                self[["atoms"]][[el]] <- atom$new(tmp, n[[el]], purity[[el]])
                self$baseWt <- self$baseWt + self$atoms[[el]]$baseWt
                self$meanWt <- self$meanWt + self$atoms[[el]]$meanWt
            }
        }
    ),
    private = list(
        List = list()
    )
)


MolData <- R6Class(
    classname = "MolData",
    inherit = atomList,
    public = list(
        tracee = NULL,
        tracer = NULL,
        len = NULL,
        p = NULL,
        pNorm = NULL,
        wt = NULL,
        MassDomain = NULL,
        result = NULL,
        initialize = function() {},
        setTracee = function(elem, n_atoms) {
            self[["tracee"]] = atomList$new(elem, n_atoms)
            invisible(self)
        },
        setTracer = function(elem, n_atoms, prob) {
            self[["tracer"]] = atomList$new(elem, n_atoms, prob)
            nn <- self$tracer$elem
            self$tracer$elem = names(self$tracer$atoms) <- paste0(nn, "1")
            invisible(self)
        },
        mergeType = function(tracer = FALSE) {
            type <- if (tracer) "tracer" else "tracee"
            for (el in self[[type]]$elem) {
                self$atoms[[el]] <- self[[type]]$atoms[[el]]
            }
            self$elem <- c(self$elem, self[[type]]$elem)
            self$baseWt <- self$baseWt + self[[type]]$baseWt
            self$meanWt <- self$meanWt + self[[type]]$meanWt
        },
        getMol = function () {
            if (!is.null(self$tracee)) {
                self$mergeType(tracer = FALSE)
            }
            if (!is.null(self$tracer)) {
                self$mergeType(tracer = TRUE)
            }
        },
        findVecLen = function() {
            tmp <- lapply(self$atoms, function(el) el$findIndex(1e-12))
            tmp <- Reduce("+", tmp)
            self$len <- self$getPwrOf2(tmp)
        },
        getMassDomain = function() {
            tmp <- lapply(self$atoms, function(el) el$adjustVec(self$len))
            tmp <- Reduce(`*`, tmp)
            self$MassDomain <- fft(tmp / self$len, inverse = TRUE)
        },
        getProb = function() {
            self$p <- Mod(self$MassDomain)
            self$pNorm <- self$p / max(self$p) * 100
        },
        getWt = function() {
            phi <- Arg(self$MassDomain) / (2 * pi)
            shiftedMass <- seq(0L, self$len - 1L)
            self$wt <- self$baseWt + shiftedMass + phi
        },
        getResultMat = function(digits = 6) {
            id <- which(self$p > 1e-5)
            index <- min(id):max(id)
            weight <- round(self$wt[index], 4)
            probability <- round(self$p[index], digits)
            normalized <- round(self$pNorm[index], digits - 2)
            self$result <- cbind(weight, probability, normalized)
            rownames(self$result) <- index - 1L
            colnames(self$result) <- c("Weight", "Distribution", "Normalized")
        },
        run = function(digits = 6) {
            self$getMol()
            self$findVecLen()
            self$getMassDomain()
            self$getProb()
            self$getWt()
            self$getResultMat(digits)
            invisible(self)
        }
    )
)

wrapperClass <- R6Class(
    classname = "wrapperClass",
    inherit = MolData,
    public = list(
        result = NULL,
        initialize = function(tracee, tracer = NULL) {
            at1 <- names(tracee)
            num1 <- unlist(tracee, use.names = FALSE)
            obj <- MolData$new()
            obj$setTracee(at1, num1)
            if (!is.null(tracer)) {
                at2 <- names(tracer)
                num2 <- sapply(tracer, `[`, 1)
                pr <- sapply(tracer, `[`, 2)
                obj$setTracer(at2, num2, pr)
            }
            obj$run()
            self$result <- obj$result
        }
    )
)

# at <- c("C", "H", "N", "O", "Si")
# num <-c(11, 30, 1, 2, 2)
#
#
# a <- MolData$new()$setTracee(at, num)
# a$setTracer("C", 6, 0.99)
# a$run()$result
#
# defaultPar <- par()
# do.call(par, defaultPar)
# par(mar = c(6, 6, 2, 6))
#
# plotAxis <- function(alist) {
#     axisList <- alist[c("side", "las", "cex.axis")]
#     mTextList <- alist[c("text", "side", "line", "cex", "font")]
#     do.call(axis, axisList)
#     do.call(mtext, mTextList)
# }
#
# argList <- list(
#     x <- NULL,
#     pch = 1,
#     yaxt = "n",
#     xaxt = "n",
#     xlab = "",
#     ylab = "",
#     lwd = 3L,
#     bty = "n",
#     type = "h"
# )
#
# axisList <- list(
#         side = 1L,
#         las = 1L,
#         cex.axis = 1.5,
#         text = "",
#         line = 4L,
#         cex = 2,
#         font = 2L
# )
#
# argList$x <- quote(a$result[, 1:2])
# do.call(plot, argList)
#
# axisList[c("text", "side")] <- list("Peaks", 1)
# do.call(plotAxis, list(axisList))
#
# axisList[c("text", "side")] <- list("Probability", 2)
# do.call(plotAxis, list(axisList))
#
# par(new = TRUE)
# argList[c("x", "y", "pch", "type")] <- list(quote(a$result[, 3L]),  NULL, "", NULL)
# axisList[c("text", "side")] <- list("Normalized", 4)
# do.call(plot, argList)
# do.call(plotAxis, list(axisList))
