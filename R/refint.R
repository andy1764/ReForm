#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # example for gamlss2
#' library(gamlss2)
#' fit <- gamlss2(Sepal.Length ~ Petal.Width + Species, data = iris[1:110,])
#' ri <- refint(fit, 95)
#' predict(ri, newdata = iris[111:150,])
refint <- function(x, ...) {
  UseMethod("refint")
}

#' @export
refint.gamlss <- function(fit, pct = 95, ...) {
  # helper function to get refint from model
  get.ri <- function(alpha, newdata, ...) {
    params = list(
      mu = predict(fit, newdata = newdata, what = "mu", type = "response", ...),
      sigma = predict(fit, newdata = newdata, what = "sigma", type = "response", ...))
    pred <- list(
      do.call(get(paste("q", fit$family[1], sep = "")), c(list(alpha/2), params)),
      do.call(get(paste("q", fit$family[1], sep = "")), c(list(1-alpha/2), params)))
    names(pred) <- c("lower", "upper")
    pred
  }

  out <- list(fit = fit, terms = as.character(attr(terms(fit), "variables"))[-1],
              get.ri = get.ri, pct = pct)
  class(out) <- "refint"
  out
}

#' @export
refint.gamlss2 <- function(x, pct = 95, ...) {
  # remove reference data
  x$model[] <- NA
  x$y[] <- NA
  x$x[] <- NA

  # helper function to get refint from model
  get.ri <- function(x, alpha, newdata, ...) {
    params = predict(x, newdata = newdata, ...)
    pred <- list(x$family$quantile(alpha/2, par = params),
                x$family$quantile(1-alpha/2, par = params))
    names(pred) <- c("lower", "upper")
    pred
  }

  out <- list(fit = x, terms = names(x$model), get.ri = get.ri, pct = pct)
  class(out) <- "refint"
  out
}

#' @export
refint.mqgam <- function(x, pct = 95, ...) {
  alpha = 1 - pct/100
  tryCatch(qdo(x, c(alpha/2, 1-alpha/2)),
           error = function(e) stop("Quantiles not found in model fit"))

  # helper function to get refint from model
  get.ri <- function(x, alpha, newdata, ...) {
    pred <- qgam::qdo(x, c(alpha/2, 1-alpha/2), predict, ...)
    names(pred) <- c("lower", "upper")
    pred
  }

  out <- list(fit = fit, terms = names(fit$model), get.ri = get.ri, pct = pct)
  class(out) <- "refint"
  out
}

#' @param x
#'
#' @param newdata
#' @param ...
#'
#' @export
predict.refint <- function(x, newdata, ...) {
  if (any(!(x$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  out <- x$get.ri(1 - x$pct/100, newdata)
  out[[x$terms[1]]] <-  newdata[,x$terms[1]]
  out$above <- newdata[,x$terms[1]] < out[[1]]
  out$below <- newdata[,x$terms[1]] > out[[2]]

  data.frame(out)
}

print.refint <- function(x) {
  cat(x$pct, "% reference interval using model of class ", class(x$fit)[1],
      sep = "")
}
