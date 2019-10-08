source("source.R")
source('Data/fun_output2.R')
sim2 <- function(seed, vari = 1, D = 2, K. = c(1, 2), 
                 lambda. = c(0, 1e-5, 1e-4, 1e-3, 1e-2)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda")
  ipath <- paste0("../4_Output/matrix/", vari, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  set.seed(seed)
  loc1 <- runif(400)
  loc2 <- runif(400)
  Data <- gData.2(seed, noise = switch(vari, 0.3, 0.6, 1.2), 250, 200, 
                  loc1, loc2)
  list2env(Data, envir = environment())
  groups <- divide(Y, seed, names = c("train", "test"))
  list2env(split(mget(c("Y", "G")), groups), envir = environment())
  nz <- ceiling(ncol(G.train)*(1-exp(-rankMatrix(G.train)/ncol(G.train))))
  # bb0 <- cbb(norder = 5, pos, nbasis = nz)
  bb0 <- cbb(norder = ceiling(sqrt(nz)), pos, nbasis = nz)
  if (!("FLM1" %in% (read.csv(ipath)$method)))
    prt(cbind(Error.flm1(Y.train, G.train, Y.test, G.test, bb0, pos), 
              method = "FLM1")[, vari.], ipath)
  E.train. <- list(X = G.train, G = NULL)
  E.test. <- list(X = G.test, G = NULL)
  E.train <- list(X = NULL, G = G.train)
  E.test <- list(X = NULL, G = G.test)
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  method1 <- paste0("NN", D-1)
  method2 <- paste0("FN", D-1)
  K <- round(ncol(G.train)^(1/3) * length(loc1)^(2/3))
  if (!(method1 %in% (read.csv(ipath)$method))) {
    Bases.nn <- c(list(NULL), as.list(rep(K, D - 1)), list(length(loc1)))
    error.nn <- Error.fnn(Y.train, E.train., NULL, Y.test, E.test., NULL,
                          pos, Bases.nn, A = A, lambda. = lambda.)
    prt(cbind(error.nn, method = method1)[, vari.], ipath)
  }
  if (!(method2 %in% (read.csv(ipath)$method))) {
    bb1 <- create.bspline.basis(norder = 5, nbasis = K)
    bbz <- create.bspline.basis(norder = 5, nbasis = 20)
    Bases.z <- f2m(bbz, loc1)[, rep(1:20, 20)] *
      f2m(bbz, loc2)[, rep(1:20, each = 20)]
    Bases.fn <- c(list(bb0), rep(list(bb1), D-1), list(Bases.z))
    error.fn <- Error.fnn(Y.train, E.train, loc, Y.test, E.test, loc, 
                          pos, Bases.fn, A = A, lambda. = lambda.)
    prt(cbind(error.fn, method = method2)[, vari.], ipath)
  }
  # 
  # Bases.nn <- c(list(NULL), as.list(rep(nb*K, D - 1)), list(length(loc1)))
  # error.nn <- Error.fnn(Y.train, E.train., NULL, Y.test, E.test., NULL,
  #                       pos, Bases.nn, A = A, lambda. = lambda.)
  # prt(cbind(error.nn, K = K, method = paste0("NN", D-1))[, vari.], ipath)
  # ## ----------Basis----------
  # bb0 <- cbb(norder = 5, pos, ratio = 0.5)
  # bb1 <- create.bspline.basis(norder = 5, nbasis = nb*K)
  # bbz <- create.bspline.basis(norder = 5, nbasis = nb*K)
  # Bases.z <- f2m(bbz, loc1)[, rep(1:nb, nb*K)] *
  #   f2m(bbz, loc2)[, rep(1:nb, each = nb*K)]
  # Bases.ff <- c(list(bb0), rep(list(bb1), D-1), list(Bases.z))
  # error.fn <- Error.fnn(Y.train, E.train, NULL, Y.test, E.test, NULL,
  #                       pos, Bases.ff, A = A, lambda. = lambda.)
  # prt(cbind(error.fn, K = K, method = paste0("FN", D-1))[, vari.], ipath)
  print(read.csv(ipath), row.names = FALSE)
}