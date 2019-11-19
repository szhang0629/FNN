source("source.R")
source("fnn_tl.R")
source("fnn_ml.R")
source("LM/flm.R")
## ------Bspline------
.bb <- create.bspline.basis(norder = 4, nbasis = 15)
.P0 <- bsplinepen(.bb, 0)
.P2 <- bsplinepen(.bb, 2)/1e5
real <- function(seed, lambda. = 10^(0:4)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda", "j")
  snp.name <- c("CHRNA5" ,"CHRNA3", "CHRNB4", "CHRNB3", "CHRNA6")[1]
  ipath <- paste0("../4_Output/", snp.name, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  data <- readRDS(paste0("../2_Data/", snp.name, ".rds"))
  list2env(data, envir = environment())
  Y.n$Y <- as.numeric(Y.n$Y)
  Y.o$Y <- as.numeric(Y.o$Y)
  Y.n <- as.matrix(log(Y.n$Y + 1))
  Y.o <- as.matrix(log(Y.o$Y + 1))
  G.n <- as.matrix(X.n[, -1])
  G.o <- as.matrix(X.o[, -1])
  X.n <- NULL
  X.o <- NULL
  snp.name <- c("CHRNA5" ,"CHRNA3", "CHRNB4", "CHRNB3", "CHRNA6")[4]
  data <- readRDS(paste0("../2_Data/", snp.name, ".rds"))
  G2.n <- as.matrix(data$X.n[, -1])
  G2.o <- as.matrix(data$X.o[, -1])
  pos2 <- data$pos
  pos <- list(pos1 = pos, pos2 = pos2)
  G.n <- cbind(G.n, G2.n)
  G.o <- cbind(G.o, G2.o)
  library(tidyverse)
  groups <- divide(1:nrow(Y.n), seed)
  list2env(sep(mget(c('Y.n', 'X.n', 'G.n')), groups), envir = environment())
  A <- c(rep(list(sigmoid), 1), list(linear))
  
  D <- 2
  Bases <- c(list(list(bb0 = .bb, bb0.1 = .bb)), rep(list(.bb), D - 1), list(1))
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