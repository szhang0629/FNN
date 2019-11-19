source("source.R")
source("fnn_tl.R")
source("fnn_ml.R")
source("LM/flm.R")
## ------Bspline------
.bb <- create.bspline.basis(norder = 4, nbasis = 15)
# .P0 <- bsplinepen(.bb, 0)
# .P2 <- bsplinepen(.bb, 2)/1e5
real <- function(seed, snp.num = 4, lambda. = 10^(0:4)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda", "j")
  snp.name <- c("CHRNA5" ,"CHRNA3", "CHRNB4", "CHRNB3", "CHRNA6")[snp.num]
  ipath <- paste0("../4_Output/", snp.name, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  data <- readRDS(paste0("../2_Data/", snp.name, "r.rds"))
  list2env(data, envir = environment())
  # Y.n$PTID <- rownames(Y.n)
  # Y.o$PTID <- rownames(Y.o)
  Y.n$Y <- as.numeric(Y.n$Y)
  Y.o$Y <- as.numeric(Y.o$Y)
  Y.n <- as.matrix(log(Y.n$Y + 1))
  Y.o <- as.matrix(log(Y.o$Y + 1))
  G.n <- as.matrix(X.n[, -1])
  G.o <- as.matrix(X.o[, -1])
  # X.n <- as.matrix(X.n[, 1])
  # X.o <- as.matrix(X.o[, 1])
  # rownames(X.n) <- rownames(G.n)
  # rownames(X.o) <- rownames(G.o)
  X.n <- NULL
  X.o <- NULL
  library(tidyverse)
  groups <- divide(1:nrow(Y.n), seed)
  list2env(sep(mget(c('Y.n', 'X.n', 'G.n')), groups), envir = environment())
  A <- c(rep(list(sigmoid), 1), list(linear))
  
  flm1.p <- flm1(Y.train, G.train, Y.test, G.test, .bb, pos, X.train,
                 X.test, lambda. = c(10^(-2:2)))
  Y.train. <- pred.flm1(G.train, .bb, pos, flm1.p$para, X.train)
  Y.test. <- pred.flm1(G.test, .bb, pos, flm1.p$para, X.test)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  prt(format(cbind(error, method = "FLM1", j = 0,
                   lambda = flm1.p$lambda)[, vari.], digits = 4), ipath)
  print(read.csv(ipath), row.names = FALSE)
  D <- 2
  Bases <- c(list(.bb), rep(list(.bb), D - 1), list(1))
  fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A = A, lambda. = lambda., pos)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos)
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases,  pos)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  output <- format(cbind(error, data.frame(lambda = fnn.p$lambda, j = fnn.p$j,
                                           method = paste0("FN", D-1)))[, vari.], digits = 4)
  prt(output, ipath)

  lambda <- fnn.p$lambda
  
  fnn.p.o <- FNN.p.(Y.o, X.o, G.o, Bases, A = A, lambda, pos)
  fnn.tl <- FNN.tl(Y.train, X.train, G.train, Bases, A = A,
                   lambda = lambda, pos, fnn.p.o$Ac)
  Y.train. <- pred(fnn.tl$Ac, X.train, G.train, A, Bases, pos = pos)
  Y.test. <- pred(fnn.tl$Ac, X.test, G.test, A, Bases, pos = pos)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  output <- format(cbind(error, data.frame(lambda = lambda, j = fnn.tl$j,
                                           method = paste0("TL", D-1)))[, vari.], digits = 4)
  prt(output, ipath)
  
  fnn.ml <- FNN.ml(Y.o, X.o, G.o, Y.train, X.train, G.train, Bases, pos, A, lambda)
  Y.train. <- pred(fnn.ml$Ac2, X.train, G.train, A, Bases, pos = pos)
  Y.test. <- pred(fnn.ml$Ac2, X.test, G.test, A, Bases, pos = pos)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  output <- format(cbind(error, data.frame(lambda =  lambda, j = fnn.ml$j,
                                           method = paste0("ML", D-1)))[, vari.], digits = 4)
  prt(output, ipath)
  
  print(read.csv(ipath), row.names = FALSE)
}