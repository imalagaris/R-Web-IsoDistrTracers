MolData <- R6Class(
    classname = "MolData",
    inherit = atomDataList,
    public = list(
        tracee = NULL,
        tracer = NULL,
        mol = NULL,
        p = NULL,
        pNorm = NULL,
        wt = NULL,
        probsTest = NULL,
        result = NULL,
        initialize = function() {},
        importTracee = function(elem, n_atoms) {
            self[["tracee"]] = atomDataList$new(elem, n_atoms)
            invisible(self)
        },
        importTracer = function(elem, n_atoms, prob) {
            self[["tracer"]] = atomDataList$new(elem, n_atoms, prob)
            invisible(self)
        },
        calculateDist = function() {
            vec <- setNames(vector("list", self$nElem), self$elem)
            len <- Reduce(`+`, lapply(self$mol, `[[`, "len"))
            for (el in self$elem) {
                vec[[el]] <- vector("numeric", len)
                vec[[el]][seq.int(self$mol[[el]]$len)] <- self$mol[[el]]$p
                vec[[el]] <- fft(vec[[el]])
            }
            vec <- Reduce(`*`, vec) / len
            vec <- Re(fft(vec, inverse = TRUE))
            self$p <- private$setNms(vec, len)
            self$reduceProbVec(1e-5)
            self$pNorm <- self$p / max(self$p) * 100
        },
        mergeAtoms = function() {
            if (is.null(self$tracee)) {
                private$mergeInfo(self$tracer)
            } else {
                private$mergeInfo(self$tracee)
                if (!is.null(self$tracer)) {
                    nams <- paste0(names(self$tracer$atoms), "1")
                    names(self$tracer$atoms) = names(self$tracer$formula) = self$tracer$elem <- nams
                    private$mergeInfo(self$tracer)
                }
            }
            self$calculateDist()
            invisible(self)
        },
        reduceAll = function() {
            for (el in self$elem) {
                self$mol[[el]]$reduceProbVec()
                self$mol[[el]]$resetLen()
            }
            invisible(self)
        },
        getGroupWts = function() {
            for (el in self$elem) {
                self$mol[[el]]$getWt()
                id <- intersect(names(self$mol[[el]]$wt), names(self$mol[[el]]$p))
                self$mol[[el]]$wt <- self$mol[[el]]$wt[id]
                self$mol[[el]]$p <- self$mol[[el]]$p[id]
                self$mol[[el]]$resetLen()
            }
            invisible(self)
        },
        adjustVecs = function() {
            for (el in self$elem) {
                self$mol[[el]]$adjustVec()
            }
        },
        getMolWeight = function() {
            self$adjustVecs()
            low <- unlist(lapply(self$mol, `[[`, "low"), use.names = TRUE)
            high <- unlist(lapply(self$mol, `[[`, "high"), use.names = TRUE)

            nameIndex <- names(self$p)
            numIndex <- setNames(as.integer(nameIndex), nameIndex)
            wts = pTest = tmp <- setNames(vector("list", self$len), nameIndex)
            for (nn in nameIndex) {
                p = wt = tmp[[nn]] <- blockparts(high - low, numIndex[nn] - sum(low))
                tmp[[nn]] <- sweep(tmp[[nn]], 1L, low, `+`)
                tmp[[nn]][] <- as.character(tmp[[nn]])
                for (el in self$elem) {
                    index <- tmp[[nn]][el, ]
                    p[el, ] <- self$mol[[el]][["p"]][index]
                    wt[el, ] <- self$mol[[el]][["wt"]][index]
                }
                p <- apply(p, 2, prod)
                wt <- colSums(wt)
                pTest[[nn]] <- sum(p)
                wts[[nn]] <- sum(wt * p / pTest[[nn]])
            }
            self$probsTest <- unlist(pTest)
            self$wt <- unlist(wts)
            invisible(self)
        },
        getResult = function() {
            self$result <- cbind(self$wt, self$p, self$pNorm)
            colnames(self$result) <- c("Weight", "Probability", "Normalized")
            invisible(self)
        },
        computeFinal = function() {
            self$mergeAtoms()
            self$reduceAll()
            self$getGroupWts()
            self$getMolWeight()
            self$getResult()
            invisible(self)
        }
    ),
    private = list(
        mergeInfo = function(List) {
            for (el in List$elem) {
                self[["mol"]][[el]] <- List[["atoms"]][[el]]
            }
            self[["formula"]] <- c(self$formula, List$formula)
            self[["elem"]] <- c(self$elem, List$elem)
            self[["nElem"]] <- sum(self$nElem, List$nElem, na.rm = TRUE)
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
            obj$importTracee(at1, num1)
            if (!is.null(tracer)) {
                at2 <- names(tracer)
                num2 <- sapply(tracer, `[`, 1)
                pr <- sapply(tracer, `[`, 2)
                obj$importTracer(at2, num2, pr)
            }
            obj$computeFinal()
            self$result <- obj$result
        }
    )
)

