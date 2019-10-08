source("source.R")
sim <- function(seed, vari = 1, p = 200, D = 3, lambda. = 10^(-1:2)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda", "j", "nl")
  ipath <- paste0("../4_Output/", p, "/", vari, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  
  n <- 250
  set.seed(seed)
  loc <- data.frame(PTID = rep(1:n, each = 20), loc = runif(20*n))
  Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, "polyno")
  # Data <- gData(seed, noise = 0, n, p, loc, "polyno")
  list2env(Data, envir = environment())
  rownames(G) <- 1:n
  Y$loc <- round(Y$loc, digits = 3)
  groups <- divide(Y$PTID, seed, "name")
  list2env(split(mget(c('Y', 'G')), groups), envir = environment())
  
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  method2 <- paste0("FN", D - 1)
  # if (!(method2 %in% (read.csv(ipath)$method))) {
  la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
  lz <- 20
  l <- round(ncol(G.train)^(1/3) * 20^(2/3))
  # if (!(method2 %in% (read.csv(ipath)$method))) {
  bb0 <- cbb(pos, nbasis = la)
  bb1 <- create.bspline.basis(norder = round(sqrt(l)), nbasis = l)
  
  bbz <- cbb(c(0, unlist(Y.train$loc), 1), nbasis = lz)
  Bases.fn <- c(list(bb0), rep(list(bb1), D - 1), list(bbz))
  error.fn <- Error.fnn(Y.train, X.train = NULL, G.train, Y.test, X.test = NULL,
                        G.test, pos, Bases.fn, A = A, lambda. = lambda.)
  prt(cbind(error.fn, method = method2, nl = paste(la, l, lz))[, vari.], 
      ipath)
  # }
  print(read.csv(ipath), row.names = FALSE, digits = 4)
}
integ <- function(D, B, int = T){
  if (is.null(D))
    return(NULL)
  # if (int) return(as.matrix(D) %*% B / nrow(B) * sqrt(ncol(B)))
  if (int) return(as.matrix(D) %*% B / nrow(B))
  else return(as.matrix(D) %*% B)
}
pen.fun <- function(bb, lambda = 1, order = 2) {
  if (is.basis(bb)) {
    if (order == 2)
      # return(bb$nbasis*bsplinepen(bb, order)/1e7)
      # return(bsplinepen(bb, order)/(sum(diag(bsplinepen(bb, 0)))^2)/10^2)
      # return(bb$nbasis*bsplinepen(bb, order)/sum(diag(bsplinepen(bb, 0)))/1e8)
      return(lambda*bb$nbasis*bsplinepen(bb, order) /
               sum(diag(bsplinepen(bb, 0))) /1e6)
    else return(bb$nbasis*bsplinepen(bb, order)/sum(diag(bsplinepen(bb, 0))))
  } else return(1)
}
f2m. <- function(bb, loc) {
  if (!is.data.frame(loc)) {
    # coef <- bb$nbasis^(1/2)
    # coef <- 1/sum(diag(bsplinepen(bb, 0)))
    coef <- (bb$nbasis/sum(diag(bsplinepen(bb, 0))))^(1/2)
    B <- eval.basis(loc, bb) * coef
  } else if (length(loc) == 1 && is.data.frame(loc))
    B <- f2m.(bb, unlist(loc))
  else B <- (f2m.(bb, loc$loc) - f2m.(bb, loc$loc0))
  return(B)
}