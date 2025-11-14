#' Reference interval calibration via conFormal prediction (ReForm)
#'
#' `reform` takes a fitted reference interval and calibrates to a new sample
#' using split conformalized quantile regression (Romano, et al., 2019).
#' `reform` requires a calibration set from the new sample (`newdata`). For
#' more details on calibration set size and ReForm properties, refer to our
#' preprint Chen et al., 2025.
#'
#' @param x `refint` object obtained from \link[ReForm]{refint}
#' @param caldata calibration data used to fit ReForm
#' @param ... additional arguments passed internally to \link[stats]{predict}
#'
#' @return
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
#'   lm(Sepal.Length ~ Petal.Width + Species,
#'     data = iris[-(1:40),])
#' )
#' reformed_ri <- reform(ri, caldata = iris[1:30,])
#' predict(reformed_ri, newdata = iris[31:40,])
#'
#' # example for quantreg::rq
#' if (require("quantreg")) {
#'   ri <- refint(
#'     rq(Sepal.Length ~ Petal.Length,
#'       data = iris[1:50,], tau = 0.025),
#'     rq(Sepal.Length ~ Petal.Length,
#'       data = iris[1:50,], tau = 0.975)
#'   )
#'   reformed_ri <- reform(ri, caldata = iris[1:30,])
#'   predict(reformed_ri, newdata = iris[31:40,])
#' }
#'
#' @references
#' Romano, Y., Patterson, E., & Candes, E. (2019). Conformalized quantile regression. Advances in neural information processing systems, 32.
reform <- function(x, caldata, ...) {
  if (any(!(x$terms %in% names(caldata)))) {stop("Terms missing from caldata")}
  alpha = 1 - x$pct/100
  nc <- nrow(caldata)

  # calibrate using CQR
  pred <- x$get.ri(caldata)
  res <- data.frame(lower = pred[[1]] - caldata[,x$terms[1]],
                    upper = caldata[,x$terms[1]] - pred[[2]])
  q_cut <- min(ceiling((1 - alpha/2)*(nc+1)), nc)
  q <- apply(res, 2, function(x) sort(x)[q_cut])

  out <- x
  out$nc <- nc
  out$res <- res
  out$cali <- q
  out$cali.df <- caldata
  class(out) <- "reformint"
  out
}

#' @param x
#'
#' @param newdata
#' @param ...
#'
#' @export
predict.reformint <- function(x, newdata, ...) {
  if (any(!(x$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  # apply ReForm calibration
  out <- x$get.ri(newdata, ...)
  out$lower <- out$lower - x$cali[1]
  out$upper <- out$upper + x$cali[2]

  out[[x$terms[1]]] <- newdata[,x$terms[1]]
  out$above <- newdata[,x$terms[1]] < out[[1]]
  out$below <- newdata[,x$terms[1]] > out[[2]]

  data.frame(out)
}

#' @param x
#'
#' @param var
#' @param ...
#'
#' @export
plot.reformint <- function(x, var, ...) {
  plot(x$cali.df[,var], x$res$lower, xlab = var, ylab = "lower residual")

  op <- par(ask=TRUE)
  plot(x$cali.df[,var], x$res$upper, xlab = var, ylab = "upper residual")
  par(op)
}

#' @export
print.reformint <- function(x) {
  if (is.null(x$fits)) {
    cat(x$pct, "% ReFormed reference interval using model of class ", class(x$fit),
        ", calibrated on ", x$nc, " observations", sep = "")
  } else {
    cat(x$pct, "% ReFormed reference interval using models of class ", class(x$fits[[1]]),
        ", calibrated on ", x$nc, " observations", sep = "")
  }
}
