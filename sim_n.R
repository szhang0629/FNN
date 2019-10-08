source("NN/cost.R")
source("NN/nn_p.R")
source("NN/prop.R")
source("NN/init.R")
source("NN/update.R")
source("Common/common.R")
source("Common/divide.R")
library(fda)
library(Deriv)
source('Data/gdata.R')
source('Data/data_input.R')
source('Data/fun_output.R')
source('Data/data_output.R')
sim <- function(seed, vari = 1, p = 200, D = 3, lambda. = 10^(-1:2)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda")
  n <- 250
  ipath <- paste0("../4_Output/", p, "/", vari, "/", seed, ".csv")
  Output.(rbind(vari.), ipath)
  set.seed(seed)
  # loc <- data.frame(PTID = rep(1:n, each = 20), loc = runif(20*n))
  loc <- data.frame(PTID = rep(1:n), loc = runif(n))
  # loc <- data.frame(PTID = rep(1:n), loc = rep(runif(1), n))
  Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, "polyno")
  list2env(Data, envir = environment())
  rownames(G) <- 1:n
  Y$loc <- round(Y$loc, digits = 3)
  X <- G[Y$PTID, ]
  
  # data <- cbind(Y, X)
  # Y <- subset(data, select = c(PTID, Y))
  # X <- subset(data, select = -Y)
  
  X <- as.matrix(cbind(Y$loc, X))
  # X <- as.matrix(X)
  Y <- subset(Y, select = c(PTID, Y))
  
  groups <- divide(Y$PTID, seed)
  list2env(split(mget(c('Y', 'X')), groups), envir = environment())
  # nb <- 20
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  method1 <- paste0("NN", D - 1)
  # K <- round(ncol(X.train)^(1/3) * nb^(2/3))
  K <- 5
  # if (!(method1 %in% (read.csv(ipath)$method))) {
    Bases.nn <- c(list(ncol(X.train)), as.list(rep(K, D - 1)), list(1))
    error.nn <- Error.nn(Y.train, X.train, Y.test, X.test, Bases.nn, 
                          A = A, lambda. = lambda.)
    prt(cbind(error.nn, method = method1)[, vari.], ipath)
  # }
  print(read.csv(ipath), row.names = FALSE, digits = 4)
}