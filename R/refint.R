refint <- function(x, ...) {
  UseMethod("refint")
}

refint.gamlss2 <- function(x, pct = 95, ...) {
  fit <- x

  # remove reference data
  fit$model[] <- NA
  fit$y[] <- NA
  fit$x[] <- NA

  out <- list(fit = fit, terms = names(fit$model), pct = pct)
  class(out) <- "refint"
  out
}

refint.mqgam <- function(x, pct = 95, ...) {
  alpha = 1 - pct/100
  tryCatch(qdo(x, c(alpha/2, 1-alpha/2)),
           error = function(e) stop("Quantiles not found in model fit"))

  out <- list(fit = fit, terms = names(fit$model), pct = pct)
  class(out) <- "refint"
  out
}

predict.refint <- function(x, newdata, ...) {
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
  out[[x$terms[1]]] <-  newdata[,x$terms[1]]
  out$above <- newdata[,x$terms[1]] < out[[1]]
  out$below <- newdata[,x$terms[1]] > out[[2]]

  data.frame(out)
}

print.refint <- function(x) {
  cat(x$pct, "% reference interval using model of class ", class(x$fit),
      sep = "")
}
