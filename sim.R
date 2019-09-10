source("source.R")
source("LM/flm.R")
source("LM/flm2.R")
sim <- function(seed, vari = 1, p = 200, D = 2, lambda. = 10^(-2:1)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda")
  n <- 250
  ipath <- paste0("../4_Output/", p, "/", vari, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  set.seed(seed)
  loc <- runif(20)
  Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, "polyno")
  list2env(Data, envir = environment())
  groups <- divide(Y, seed, names = c("train", "test"))
  list2env(split(mget(c("Y", "G")), groups), envir = environment())
  nz <- ceiling(ncol(G.train)*(1-exp(-nrow(G.train)/ncol(G.train))))
  # bb0 <- cbb(norder = 5, pos, nbasis = nz)
  bb0 <- cbb(norder = ceiling(sqrt(nz)), pos, nbasis = nz)
  bbz <- cbb(norder = 5, loc, nbasis = length(loc))
  if (!("FLM1" %in% (read.csv(ipath)$method)))
    # apply(data[ , cols ], 1 , paste0)
    prt(cbind(Error.flm1(Y.train, G.train, Y.test, G.test, bb0, pos), 
              method = "FLM1")[, vari.], ipath)
  if (!("FLM2" %in% (read.csv(ipath)$method)))
    prt(cbind(Error.flm2(Y.train, G.train, Y.test, G.test, pos, loc), 
              method = "FLM2")[, vari.], ipath)
  E.train. <- list(X = G.train, G = NULL)
  E.test. <- list(X = G.test, G = NULL)
  E.train <- list(X = NULL, G = G.train)
  E.test <- list(X = NULL, G = G.test)
  # bb0. <- cbb(norder = ceiling(sqrt(nz)), pos, nbasis = nz)
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  method1 <- paste0("NN", D-1)
  method2 <- paste0("FF", D-1)
  K <- round(ncol(G.train)^(1/3) * length(loc)^(2/3))
  if (!(method1 %in% (read.csv(ipath)$method))) {
    Bases.nn <- c(list(NULL), as.list(rep(K, D - 1)), list(length(loc)))
    error.nn <- Error.fnn(Y.train, E.train., NULL, Y.test, E.test., NULL,
                          pos, Bases.nn, A = A, lambda. = lambda.)
    prt(cbind(error.nn, method = method1)[, vari.], ipath)
  }
  if (!(method2 %in% (read.csv(ipath)$method))) {
    bb1 <- create.bspline.basis(norder = 5, nbasis = K)
    Bases.fn <- c(list(bb0), rep(list(bb1), D - 1), list(bbz))
    error.fn <- Error.fnn(Y.train, E.train, loc, Y.test, E.test, loc, 
                          pos, Bases.fn, A = A, lambda. = lambda./1e5)
    prt(cbind(error.fn, method = method2)[, vari.], ipath)
  }
  print(read.csv(ipath), row.names = FALSE, digits = 4)
}
# error.lm <- Error.lm(Y.train, G.train, Y.test, G.test, pos)
# prt(cbind(error.lm, K = 0, valid = 0, lambda = 0)[, vari.], ipath)
# error.flm <- Error.flm(Y.train, G.train, Y.test, G.test, pos, loc)
# prt(cbind(error.flm, K = 0, valid = 0, lambda = 0)[, vari.], ipath)
