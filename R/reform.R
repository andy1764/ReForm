#' Reference interval calibration via conFormal prediction (ReForm)
#'
#' `reform` takes a fitted reference interval and calibrates to a new sample
#' using split conformalized quantile regression (Romano, et al., 2019).
#' `reform` requires a calibration set `caldata` from the new sample. For more
#' details on calibration set size and ReForm properties, refer to our preprint
#' Chen et al., 2025.
#'
#' @param object `refint` object obtained from \link[ReForm]{refint}
#' @param caldata calibration data used to fit ReForm
#' @param ... additional arguments passed internally to \link[stats]{predict}
#'
#' @returns `reform` returns a list containing all elements of `object` (see \link[ReForm]{refint}) and:
#' \item{cali}{lower and upper calibration constants used to adjust the interval}
#' \item{cali.df}{`caldata` stored in the output}
#' @export
#'
#' @seealso
#' \link[ReForm]{plot.reformint} for assessing fit of reference interval in
#' calibration set via diagnostic plots
#'
#' \link[ReForm]{predict.reformint} for applying a ReFormed reference
#' interval to new observations
#'
#' @examples
#' # example for lm (not recommended)
#' ri <- refint(
#'   lm(Sepal.Length ~ Sepal.Width,
#'     data = iris[51:100,]), pct = 90
#' )
#' predict(ri, newdata = iris[1:10,])
#' reformed_ri <- reform(ri, caldata = iris[11:50,])
#' predict(reformed_ri, newdata = iris[1:10,])
#'
#' # example for quantreg::rq
#' if (require("quantreg")) {
#'   ri <- refint(
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.05),
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.95)
#'   )
#'   reformed_ri <- reform(ri, caldata = iris[11:50,])
#'   predict(reformed_ri, newdata = iris[1:10,])
#' }
#'
#' @references
#' Romano, Y., Patterson, E., & Candes, E. (2019). Conformalized quantile regression. Advances in neural information processing systems, 32.
reform <- function(object, caldata, ...) {
  if (any(!(object$terms %in% names(caldata)))) {
    stop("Terms missing from caldata")
  }
  alpha = 1 - object$pct/100
  nc <- nrow(caldata)

  # calibrate using CQR
  pred <- object$get.ri(caldata)
  res <- data.frame(lower = pred[[1]] - caldata[,object$terms[1]],
                    upper = caldata[,object$terms[1]] - pred[[2]])
  if (ceiling((1 - alpha/2)*(nc+1)) > nc) {
    warning("Calibration set size is inadequate, results may be suboptimal")
  }
  q_cut <- min(ceiling((1 - alpha/2)*(nc+1)), nc)
  q <- apply(res, 2, function(x) sort(x)[q_cut])

  out <- object
  out$cali <- q
  out$cali.df <- caldata
  class(out) <- "reformint"
  out
}

#' Predict Method for ReFormed Reference Intervals
#'
#' Apply ReFormed reference interval from \link[ReForm]{reform} to new data
#'
#' @param object `reformint` object obtained from \link[ReForm]{reform}
#' @param newdata input variables used to get reference intervals
#' @param ... additional arguments passed to \link[stats]{predict}
#'
#' @returns a data frame containing lower and upper bounds of the reference
#' interval, the test observation, and logicals for whether it is above or below
#'
#' @export
#'
#' @examples
#' # example for lm (not recommended)
#' ri <- refint(
#'   lm(Sepal.Length ~ Sepal.Width,
#'     data = iris[51:100,]), pct = 90
#' )
#' reformed_ri <- reform(ri, caldata = iris[11:50,])
#' predict(reformed_ri, newdata = iris[1:10,])
#'
#' # example for quantreg::rq
#' if (require("quantreg")) {
#'   ri <- refint(
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.05),
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.95)
#'   )
#'   reformed_ri <- reform(ri, caldata = iris[11:50,])
#'   predict(reformed_ri, newdata = iris[1:10,])
#' }
predict.reformint <- function(object, newdata, ...) {
  if (any(!(object$terms %in% names(newdata)))) {
    stop("Terms missing from newdata")
  }

  # apply ReForm calibration
  out <- object$get.ri(newdata, ...)
  out$lower <- out$lower - object$cali[1]
  out$upper <- out$upper + object$cali[2]

  out[[object$terms[1]]] <- newdata[,object$terms[1]]
  out$below <- newdata[,object$terms[1]] < out$lower
  out$above <- newdata[,object$terms[1]] > out$upper

  data.frame(out)
}

#' Diagnostic Plots for ReForm
#'
#' Plot upper and lower residuals from the calibration set in
#' \link[ReForm]{reform} against a chosen variable `var`. Ideally, both plots
#' should display a roughly flat relationship.
#'
#' @param x `reformint` object obtained from \link[ReForm]{reform}
#' @param var string indicating to plot against variable to plot against (should
#'   be in `x$terms`). If none supplied, defaults to the first covariate (second
#'   value) in `x$terms`.
#' @param ... additional arguments passed to \link[base]{plot}
#'
#' @export
#'
#' @examples
#' # example for lm (not recommended)
#' ri <- refint(
#'   lm(Sepal.Length ~ Sepal.Width,
#'     data = iris[51:100,]), pct = 90
#' )
#' reformed_ri <- reform(ri, caldata = iris[11:50,])
#' plot(reformed_ri)
#'
#' # example for quantreg::rq
#' if (require("quantreg")) {
#'   ri <- refint(
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.05),
#'     rq(Sepal.Length ~ Sepal.Width,
#'       data = iris[51:100,], tau = 0.95)
#'   )
#'   reformed_ri <- reform(ri, caldata = iris[11:50,])
#'   plot(reformed_ri)
#' }
plot.reformint <- function(x, var = NULL, ...) {
  if (is.null(var)) {
    var <- x$terms[2]
  }
  pred <- x$get.ri(x$cali.df)
  res <- data.frame(lower = pred[[1]] - x$cali.df[,x$terms[1]],
                    upper = x$cali.df[,x$terms[1]] - pred[[2]])

  plot(x$cali.df[,var], res$lower, xlab = var, ylab = "lower residual")

  op <- par(ask=TRUE)
  plot(x$cali.df[,var], res$upper, xlab = var, ylab = "upper residual")
  par(op)
}

#' @export
print.reformint <- function(x) {
  if (is.null(x$fits)) {
    cat(x$pct, "% ReFormed reference interval using model of class ",
        class(x$fit), ", calibrated on ", nrow(x$cali.df), " observations\n",
        sep = "")
  } else {
    cat(x$pct, "% ReFormed reference interval using models of class ",
        class(x$fits[[1]]), ", calibrated on ", nrow(x$cali.df),
        " observations\n", sep = "")
  }
}
