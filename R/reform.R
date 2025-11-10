reform <- function(x, ...) {
  UseMethod("reform")
}

reform.refint <- function(x, newdata, ...) {
  if (any(!(x$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  alpha = 1 - x$pct/100
  nc <- nrow(newdata)

  # get intervals using predict for model class
  switch(
    class(x$fit),
    "gamlss2" = {
      params = predict(x$fit, newdata = newdata, ...)
      pred <- list(x$fit$family$quantile(alpha/2, par = params),
                  x$fit$family$quantile(1-alpha/2, par = params))
    },
    "mqgam" = {
      pred <- qgam::qdo(x$fit, c(alpha/2, 1-alpha/2), predict, ...)
    })

  # calibrate using ReForm
  res <- data.frame(lower = pred[[1]] - newdata[,x$terms[1]],
                    upper = newdata[,x$terms[1]] - pred[[2]])
  q_cut <- min(ceiling((1 - alpha/2)*(nc+1)), nc)
  q <- apply(res, 2, function(x) sort(x)[q_cut])

  out <- x
  out$nc <- nc
  out$res <- res
  out$cali <- q
  class(out) <- "reformint"
  out
}

predict.reformint <- function(x, newdata, ...) {
  if (any(!(x$terms %in% names(newdata)))) {stop("Terms missing from newdata")}

  alpha = 1 - x$pct/100

  # get intervals using predict for model class
  switch(
    class(x$fit),
    "gamlss2" = {
      params = predict(x$fit, newdata = newdata, ...)
      out <- list(x$fit$family$quantile(alpha/2, par = params),
                  x$fit$family$quantile(1-alpha/2, par = params))
    },
    "mqgam" = {
      out <- qgam::qdo(x$fit, c(alpha/2, 1-alpha/2), predict, ...)
    })
  names(out) <- c("lower", "upper")

  # apply ReForm calibration
  out$lower <- out$lower - x$cali[1]
  out$upper <- out$upper + x$cali[2]

  out[[x$terms[1]]] <- newdata[,x$terms[1]]
  out$above <- newdata[,x$terms[1]] < out[[1]]
  out$below <- newdata[,x$terms[1]] > out[[2]]

  data.frame(out)
}

print.reformint <- function(x) {
  cat(x$pct, "% ReFormed reference interval using model of class ", class(x$fit),
      ", calibrated on ", x$nc, " observations", sep = "")
}
