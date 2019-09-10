func.plot <- function(bb, coef) {
  
}

bb <- create.bspline.basis(nbasis = 10000, norder = 5)
coef <- rnorm(10000)
f. <- function(x) {
  return(eval.basis(x, bb) %*% coef)
}
plot(f.)
