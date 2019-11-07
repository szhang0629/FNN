source("source.R")
source('Data/fun_output2.R')
## aligned condition without loc
sim2 <- function(seed, vari = 1, D = 3, lambda. = 0) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda", "j")
  ipath <- paste0("../4_Output/matrix/", vari, "/", seed, ".csv")
  set.seed(seed)
  n <- 250
  # loc <- data.frame(PTID = rep(1:n, each = 400), loc1 = rep(runif(1:400), n), 
  #                   loc2 = rep(runif(1:400), n))
  loc <- data.frame(loc1 = runif(1:400), loc2 = runif(1:400))
  Data <- gData.2(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, 200, loc)
  list2env(Data, envir = environment())
  # rownames(G) <- 1:n
  X <- NULL
  groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  Output.(rbind(vari.), ipath)
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  ## ----------Basis----------
  la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
  bb0 <- create.bspline.basis(norder = 4, nbasis = la)
  bb1 <- create.bspline.basis(norder = 4, nbasis = 15)
  bbz <- create.bspline.basis(norder = 4, nbasis = 15)
  Bases.z <- f2m(bbz, loc$loc1)[, rep(1:15, 15)] * 
    f2m(bbz, loc$loc2)[, rep(1:15, each = 15)]
  Bases <- c(list(bb0), rep(list(bb1), D - 1), list(Bases.z))
  # Bases <- c(list(bb0), rep(list(bb1), D - 1), list(list(bbz, bbz)))
  fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda. = lambda., pos)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos)
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos)
  
  error <- Error(c(Y.train), c(Y.test), c(Y.train.), c(Y.test.))
  error.fn <- cbind(error, data.frame(lambda = fnn.p$lambda, j = fnn.p$j))
  
  # train <- mean((Y.train - Y.train.)^2) /
  #   mean((Y.train - rep.row(colMeans(Y.train), nrow(Y.train)))^2)
  # test <- mean((Y.test - Y.test.)^2) /
  #   mean((Y.test - rep.row(colMeans(Y.test), nrow(Y.test)))^2)
  # cor1 <- cor(c(Y.train), c(Y.train.))
  # cor2 <- cor(c(Y.test), c(Y.test.))
  # error.fn <- data.frame(train = train, test = test, cor1 = cor1, cor2 = cor2, 
  #                        lambda = fnn.p$lambda, j = fnn.p$j)
  output <- format(cbind(error.fn, method = paste0("FN", D - 1))[, vari.]
                   , digits = 4)
  prt(output, ipath)
  print(read.csv(ipath), row.names = FALSE)
}
# nz <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
# bb0 <- cbb(norder = 5, pos, nbasis = nz)
# bb0 <- cbb(norder = ceiling(sqrt(nz)), pos, nbasis = nz)
# if (!("FLM1" %in% (read.csv(ipath)$method)))
#   prt(cbind(Error.flm1(Y.train, G.train, Y.test, G.test, bb0, pos), 
#             method = "FLM1")[, vari.], ipath)
# E.train. <- list(X = G.train, G = NULL)
# E.test. <- list(X = G.test, G = NULL)
# E.train <- list(X = NULL, G = G.train)
# E.test <- list(X = NULL, G = G.test)
# method1 <- paste0("NN", D - 1)
# method2 <- paste0("FN", D - 1)
# K <- round(ncol(G.train)^(1/3) * length(loc1)^(2/3))
# if (!(method1 %in% (read.csv(ipath)$method))) {
#   Bases.nn <- c(list(NULL), as.list(rep(K, D - 1)), list(length(loc1)))
#   error.nn <- Error.fnn(Y.train, E.train., NULL, Y.test, E.test., NULL,
#                         pos, Bases.nn, A = A, lambda. = lambda.)
#   prt(cbind(error.nn, method = method1)[, vari.], ipath)
# }
# if (!(method2 %in% (read.csv(ipath)$method))) {
#   bb1 <- create.bspline.basis(norder = 5, nbasis = K)
#   bbz <- create.bspline.basis(norder = 5, nbasis = 20)
#   Bases.z <- f2m(bbz, loc1)[, rep(1:20, 20)] *
#     f2m(bbz, loc2)[, rep(1:20, each = 20)]
#   Bases.fn <- c(list(bb0), rep(list(bb1), D - 1), list(Bases.z))
#   error.fn <- Error.fnn(Y.train, E.train, loc, Y.test, E.test, loc, 
#                         pos, Bases.fn, A = A, lambda. = lambda.)
#   prt(cbind(error.fn, method = method2)[, vari.], ipath)
# }
# Bases.nn <- c(list(NULL), as.list(rep(nb*K, D - 1)), list(length(loc1)))
# error.nn <- Error.fnn(Y.train, E.train., NULL, Y.test, E.test., NULL,
#                       pos, Bases.nn, A = A, lambda. = lambda.)
# prt(cbind(error.nn, K = K, method = paste0("NN", D-1))[, vari.], ipath)