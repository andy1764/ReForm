#' Constructing Reference Intervals
#'
#' `refint` is a generic function used to construct `pct`% reference intervals
#' from a model object. Currently supports both distributional regression models
#' (\link[stats]{lm} and \link[gamlss2]{gamlss2}) and quantile regression models
#' (\link[quantreg]{rq} and \link[qgam]{qgam}).
#'
#' @param object a model object used to construct the reference interval
#' @param upper for quantile regression, the upper quantile fit
#' @param pct for distributional regression, the desired percentage
#' @param ... additional arguments passed internally to \link[stats]{predict}
#'
#' @returns `refint` returns a list containing:
#' \item{fit}{original model fit for \link[stats]{lm} and \link[gamlss2]{gamlss2}}
#' \item{fits}{lower and upper fits for \link[quantreg]{rq} and \link[qgam]{qgam}}
#' \item{terms}{variables used in the original model fit}
#' \item{get.ri}{function used to get reference interval from model fit(s)}
#' \item{pct}{percentage reference interval (i.e. 95 refers to 95% reference interval}
#'
#' @import stats
#' @export
#'
#' @seealso
#' \link[ReForm]{predict.refint} for applying reference intervals to new data
#' \link[ReForm]{reform} for calibration to a new sample
#'
#' @rdname refint
#' @examples
#' # example for lm (not recommended)
#' ri <- refint(
#'   lm(Sepal.Length ~ Petal.Width + Species,
#'     data = iris[-(1:20),])
#' )
#' predict(ri, newdata = iris[1:20,])
#'
#' # example for quantreg::rq
#' if (require("quantreg")) {
#'   ri <- refint(
#'     rq(Sepal.Length ~ Petal.Length + Species,
#'       data = iris[-(1:20),], tau = 0.025),
#'     rq(Sepal.Length ~ Petal.Length + Species,
#'       data = iris[-(1:20),], tau = 0.975)
#'   )
#'   predict(ri, newdata = iris[1:20,])
#' }
#'
#' # example for qgam
#' if (require("quantreg")) {
#'   ri <- refint(
#'     qgam(Sepal.Length ~ s(Petal.Length) + Species,
#'       data = iris[-(1:20),], qu = 0.025),
#'     qgam(Sepal.Length ~ s(Petal.Length) + Species,
#'       data = iris[-(1:20),], qu = 0.975)
#'   )
#'   predict(ri, newdata = iris[1:20,])
#' }
#'
#' if (require("gamlss2")) {
#'   ri <- refint(
#'     gamlss2(Sepal.Length ~ s(Petal.Length) + Species |
#'       s(Petal.Length) + Species, data = iris[-(1:20),],
#'       family = BCT))
#'   predict(ri, newdata = iris[1:20,])
#' }
refint <- function(object, ...) {
  UseMethod("refint")
}

#' @rdname refint
#' @export
refint.lm <- function(object, pct = 95, ...) {
  get.ri <- function(newdata, ...) {
    alpha <- 1-pct/100
    params <- predict(object, newdata = newdata, ...)
    pred <- list(qnorm(alpha/2, mean = params, sd = summary(object)$sigma),
                 qnorm(1-alpha/2, mean = params, sd = summary(object)$sigma))
    names(pred) <- c("lower", "upper")
    pred
  }

  out <- list(fit = object, terms = names(object$model), get.ri = get.ri,
              pct = pct)
  class(out) <- "refint"
  out
}

#' @rdname refint
#' @export
refint.gamlss2 <- function(object, pct = 95, ...) {
  get.ri <- function(newdata, ...) {
    alpha <- 1-pct/100
    params = predict(object, newdata = newdata, ...)
    pred <- list(object$family$quantile(alpha/2, par = params),
                 object$family$quantile(1-alpha/2, par = params))
    names(pred) <- c("lower", "upper")
    pred
  }

  out <- list(fit = object, terms = names(object$model), get.ri = get.ri,
              pct = pct)
  class(out) <- "refint"
  out
}

#' @rdname refint
#' @export
refint.rq <- function(object, upper, ...) {
  if (!all.equal(object$tau, 1-upper$tau)) {
    stop("Lower and upper fits not matched")
  }
  if (!identical(names(object$model), names(upper$model))) {
    stop("Lower and upper fits have different terms")
  }

  get.ri <- function(newdata, ...) {
    list(lower = predict(object, newdata = newdata, ...),
         upper = predict(upper, newdata = newdata, ...))
  }
  pct <- 100*(1-2*object$tau)
  out <- list(fits = list(lower = object, upper = upper),
              terms = names(object$model), get.ri = get.ri, pct = pct)
  class(out) <- "refint"
  out
}

#' @rdname refint
#' @export
refint.qgam <- function(object, upper, ...) {
  if (!all.equal(object$family$getQu(), 1-upper$family$getQu())) {
    stop("Lower and upper fits not matched")
  }
  if (!identical(names(object$model), names(upper$model))) {
    stop("Lower and upper fits have different terms")
  }

  get.ri <- function(newdata, ...) {
    list(lower = predict(object, newdata = newdata, ...),
         upper = predict(upper, newdata = newdata, ...))
  }
  pct <- 100*(1-2*object$family$getQu())
  out <- list(fits = list(lower = object, upper = upper),
              terms = names(object$model), get.ri = get.ri, pct = pct)
  class(out) <- "refint"
  out
}

#' Predict Method for Reference Intervals
#'
#' Apply fitted reference interval from \link[ReForm]{refint} to new data
#'
#' @param object `refint` object obtained from \link[ReForm]{refint}
#' @param newdata input variables used to get reference intervals
#' @param ... additional arguments passed to \link[stats]{predict}
#'
#' @returns a data frame containing lower and upper bounds of the reference
#' interval, the test observation, and logicals for whether it is above or below
#'
#' @export
predict.refint <- function(object, newdata, ...) {
  if (any(!(object$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  out <- object$get.ri(newdata, ...)
  out[[object$terms[1]]] <-  newdata[,object$terms[1]]
  out$below <- newdata[,object$terms[1]] < out$lower
  out$above <- newdata[,object$terms[1]] > out$upper

  data.frame(out)
}

#' @export
print.refint <- function(x) {
  if (is.null(x$fits)) {
    cat(x$pct, "% reference interval using model of class ", class(x$fit)[1],
        sep = "")
  } else {
    cat(x$pct, "% reference interval using models of class ",
        class(x$fits[[1]])[1], sep = "")
  }
}
