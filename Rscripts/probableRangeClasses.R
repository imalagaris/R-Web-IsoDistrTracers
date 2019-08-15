ProbableRange <- R6Class(
    classname = "ProbableRange",
    public = list(
        highLim = 1e-6,
        lowLim = 1e-8,
        minIndex = NULL,
        maxIndex = NULL,
        initialize = function(n, p) {
            private$p <- p
            private$n <- as.integer(n)
            private$mu <- n * p
        },
        getHigh = function() {
            private$low <- as.integer(ceiling(private$mu))
            private$high <- private$n
            self$maxIndex <- private$getMax()
            invisible(self)
        },
        getLow = function() {
            private$low <- 0L
            private$high <- as.integer(floor(private$mu))
            self$minIndex <- private$getMin()
            invisible(self)
        },
        getIndex = function() {
            self$getLow()
            self$getHigh()
            c(self$minIndex, self$maxIndex)
        }
    ),
    private = list(
        n = NULL,
        p = NULL,
        mu = NULL,
        prob = 0,
        low = NULL,
        middle = NULL,
        high = NULL,
        exitOn1 = function() {
            p <- private$prob
            p < self$highLim && p > self$lowLim
        },
        exitOn2 = function() {
            (private$high - private$low) < 2L
        },
        on2Low = function() {
            private$low
        },
        on2High = function() {
            private$high
        },
        getIntervalMin = function() {
            if (private$prob > self$highLim) {
                private$high <- private$middle
            } else if (private$prob < self$lowLim) {
                private$low <- private$middle
            }
        },
        getIntervalMax = function() {
            if (private$prob > self$highLim) {
                private$low <- private$middle
            } else if (private$prob < self$lowLim) {
                private$high <- private$middle
            }
        },
        getMiddle = function() {
            private$middle <- as.integer(round((private$low + private$high) * 0.5))
            private$prob <- dbinom(private$middle, private$n, private$p)
        },
        getMin = function() {
            private$getMiddle()
            if (private$exitOn1()) {
                return(private$middle)
            }
            private$getIntervalMin()
            if (private$exitOn2()) {
                return(private$on2Low())
            }
            private$getMin()
        },
        getMax = function() {
            private$getMiddle()
            if (private$exitOn1()) {
                return(private$middle)
            }
            private$getIntervalMax()
            if (private$exitOn2()) {
                return(private$on2High())
            }
            private$getMax()
        }
    )
)

ProbableVec <- R6Class(
    classname = "ProbableVec",
    inherit = ProbableRange,
    public = list(
        initialize = function(vec, peak = NULL) {
            self$highLim <- 1e-6
            self$lowLim <- 1e-8
            if (is.null(peak)) {
                private$mu <- which.max(vec)
            } else {
                private$mu <- as.integer(peak)
            }
            private$vec <- vec
            private$n <- length(vec)
        },
        getLow = function() {
            private$low <- 1L
            private$high <- private$mu
            self$minIndex <- private$getMin()
            invisible(self)
        },
        one2HighVec = function() {
            self$getHigh()
            private$vec[seq(self$maxIndex)]
        },
        low2HighVec = function() {
            self$getLow()
            self$getHigh()
            list(
                "low" = self$minIndex - 1L,
                "high" = self$maxIndex - 1L,
                "vec" = private$vec[seq(self$minIndex, self$maxIndex)]
            )
        }
    ),
    private = list(
        vec = NULL,
        n = NULL,
        mu = NULL,
        getMiddle = function() {
            private$middle <- as.integer(round((private$low + private$high) * 0.5))
            private$prob <- private$vec[private$middle]
        },
        on2Low = function() {
            tmp <- c(private$low, private$middle, private$high)
            p <- private$vec[tmp]
            tmp <- tmp[p > self$lowLim]
            min(tmp)
        },
        on2High = function() {
            tmp <- c(private$low, private$middle, private$high)
            p <- private$vec[tmp]
            tmp <- tmp[p > self$lowLim]
            min(tmp)
        }
    )
)
