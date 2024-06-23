spline.regression <- function(x, y, df = 6,
                              eval.points = seq(min(x), max(x), length = 100),
                              deriv = 0, display = TRUE, add = FALSE, ...) {
  m <- smooth.spline(x, y, df = df)
  estimate <- predict(m, eval.points, deriv = deriv)$y
  if (add == TRUE) {
    lines(eval.points, estimate, ...)
  } else if (display != "none") {
      if (deriv == 0) {
        plot(x, y, main = "spline.regression", xlab = substitute(x),
             ylab = substitute(y), cex = 0.5)
        lines(eval.points, estimate, col = "red", ...)
      }
      else {
        plot(eval.points, estimate, type = "l", main =
                 "spline.regression", xlab = substitute(x), ylab =
                 paste("D(", substitute(y), ",", deriv, ")", sep = ""))
    }
  }
  invisible(list(eval.points = eval.points, estimate = estimate, df = df, m))
}
