#' Reference interval calibration via conFormal prediction (ReForm)
#'
#'
#'
#' @param x `refint` object, typically result of \link[ReForm]{refint}
#' @param newdata calibration data used to fit ReForm
#' @param ... other parameters passed internally to `predict`
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
#' # example for gamlss2
#' library(gamlss2)
#' fit <- gamlss2(Sepal.Length ~ Petal.Width + Species, data = iris[1:110,])
#' ri <- refint(fit, 95)
#' reformint <- reform(ri, newdata = iris[111:130,])
#' predict(reformint, newdata = iris[131:150,])
reform <- function(x, newdata, ...) {
  if (any(!(x$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  alpha = 1 - x$pct/100
  nc <- nrow(newdata)

  # calibrate using ReForm
  pred <- x$get.ri(x$fit, alpha, newdata)
  res <- data.frame(lower = pred[[1]] - newdata[,x$terms[1]],
                    upper = newdata[,x$terms[1]] - pred[[2]])
  q_cut <- min(ceiling((1 - alpha/2)*(nc+1)), nc)
  q <- apply(res, 2, function(x) sort(x)[q_cut])

  out <- x
  out$nc <- nc
  out$res <- res
  out$cali <- q
  out$cali.df <- newdata
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
  out <- x$get.ri(x$fit, 1 - x$pct/100, newdata)
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

print.reformint <- function(x) {
  cat(x$pct, "% ReFormed reference interval using model of class ", class(x$fit),
      ", calibrated on ", x$nc, " observations", sep = "")
}
