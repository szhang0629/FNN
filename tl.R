library(fda)
library(Deriv)
library(MASS)
source("fnn_tl.R")
source("fnn_ml.R")
source("LM/flm.R")
source('Common/common.R')
source('Common/divide.R')
## ------Bspline------
tl <- function(seed, snp.num = 1, Type = "r") {
  vari. <- c("method", "train", "test", "mse1", "mse2", "cor1", "cor2")
  snp.name <- c("CHRNA5" ,"CHRNA3", "CHRNB4", "CHRNB3", "CHRNA6")[snp.num]
  ipath <- paste0("../4_Output/", snp.name, Type, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  G <- read.csv(paste0("../2_Data/", snp.name, "/g.csv"), row.names = 1)
  X <- read.csv(paste0("../2_Data/", snp.name, "/x.csv"), row.names = 1)
  Y <- read.csv(paste0("../2_Data/", snp.name, "/y.csv"), row.names = 1)
  pos <- colnames(G)
  pos <- substring(pos, 2)
  pos <- as.numeric(pos)
  pos <- (pos - min(pos))/(max(pos) - min(pos))
  Y <- as.matrix(log(as.numeric(Y$Y) + 1))
  Y <- as.matrix((Y - mean(Y))/sd(Y))
  # X$sex <- as.numeric(X$sex)
  groups <- list(o = list((1:nrow(X))[(X$race == 1)]),
                 n = list((1:nrow(X))[(X$race == 0)]))
  # groups <- list(o = list((1:nrow(X))[(X$sex == 2)]),
  #                n = list((1:nrow(X))[(X$sex == 1)]))
  # # X <- as.matrix(X[, "sex", drop = FALSE])
  # # X <- NULL
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  # Y.o <- as.matrix(lm(Y~sex, cbind(Y = Y.o, X.o))$residual)
  # Y.n <- as.matrix(lm(Y~sex, cbind(Y = Y.n, X.n))$residual)
  # library(kernlab)
  # knn.model.o <- gausspr(X.o$age_int, Y.o)
  # Y.o <- Y.o - predict(knn.model.o, X.o$age_int)
  # knn.model.n <- gausspr(X.n$age_int, Y.n)
  # Y.n <- Y.n - predict(knn.model.n, X.n$age_int)
  library(tidyverse)
  groups <- divide(1:nrow(Y.n), seed)
  list2env(sep(mget(c('Y.n', 'X.n', 'G.n')), groups), envir = environment())
  
  X.o <- NULL
  X.train <- NULL
  X.test <- NULL
  loc <- NULL
  HL <- 1
  A <- c(rep(list(sigmoid), HL), list(linear))
  lambda. = 10^((6:10)/2)
  lambda. = 0
  .bb <- create.bspline.basis(norder = 4, nbasis = 20)
  bb0 <- cbb(norder = 4, pos, .bb$nbasis)
  # bb0 <- .bb
  Bases <- c(list(bb0), rep(list(.bb), HL), list(1))
  if (!("FLM1" %in% (read.csv(ipath)$method))) {
    flm1.p <- flm1(Y.train, G.train, bb0, pos, X.train, lambda.)
    Y.train. <- pred.flm1(G.train, bb0, pos, flm1.p$para, X.train)
    Y.test. <- pred.flm1(G.test, bb0, pos, flm1.p$para, X.test)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, method = "FLM1")[, vari.], digits = 4), ipath)
  }
  # if (!(paste0("FN", HL) %in% (read.csv(ipath)$method))) {
  source('FNN/f2m.R')
  source('FNN/cost.R')
  source('FNN/fnn_p.R')
  source('FNN/init.R')
  source('FNN/prop.R')
  source('FNN/pen.R')
    fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, loc)
    Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos)
    Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases,  pos)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, method = paste0("FN", HL))[, vari.], 
               digits = 4), ipath)
  # }
  # if (!(paste0("TL", HL) %in% (read.csv(ipath)$method))) {
    fnn.p.o <- FNN.p(Y.o, X.o, G.o, Bases, A = A, lambda., pos)
    res.train <- forw.prop(fnn.p.o$Ac, X.train, G.train, Bf2m(Bases, pos, NULL), 
                           A, lapply(Bases, is.basis))
    G..train <- res.train[[length(res.train) - 1]]$G
    res.test <- forw.prop(fnn.p.o$Ac, X.test, G.test, Bf2m(Bases, pos, NULL), 
                          A, lapply(Bases, is.basis))
    G..test <- res.test[[length(res.test) - 1]]$G
    pos. <- (1:ncol(G..train) - 0.5) / ncol(G..train)
    flm1..p <- flm1(Y.train, G..train, .bb, pos., X.train, lambda.)
    
    Y.train. <- pred.flm1(G..train, .bb, pos., flm1..p$para, X.train)
    Y.test. <- pred.flm1(G..test, .bb, pos., flm1..p$para, X.test)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, method = paste0("TL", HL))[, vari.], 
               digits = 4), ipath)
  # }
  
  print(read.csv(ipath), row.names = FALSE)
}

# fnn.tl <- FNN.tl(Y.train, X.train, G.train, Bases, A = A,
#                  lambda = lambda, pos, fnn.p.o$Ac)
# Y.train. <- pred(fnn.tl$Ac, X.train, G.train, A, Bases, pos = pos)
# Y.test. <- pred(fnn.tl$Ac, X.test, G.test, A, Bases, pos = pos)
# error <- Error(Y.train, Y.test, Y.train., Y.test.)
# output <- format(cbind(error, data.frame(method = paste0("TL", HL)))
#                  [, vari.], digits = 4)
# prt(output, ipath)

# fnn.ml <- FNN.ml(Y.o, X.o, G.o, Y.train, X.train, G.train, Bases, pos, A, lambda)
# Y.train. <- pred(fnn.ml$Ac2, X.train, G.train, A, Bases, pos = pos)
# Y.test. <- pred(fnn.ml$Ac2, X.test, G.test, A, Bases, pos = pos)
# error <- Error(Y.train, Y.test, Y.train., Y.test.)
# output <- format(cbind(error, data.frame(method = paste0("ML", HL)))
#                  [, vari.], digits = 4)
# prt(output, ipath)