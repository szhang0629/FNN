library(fda)
library(Deriv)
library(MASS)
source('Common/common.R')
source('Common/divide.R')
source('Data/gdata.R')
main <- function(seed, vari = 0, HL = 1, p = 200, fct = "polyno") {
  vari. <- c("index", "mse1", "mse2")
  Data <- gData(seed, vari, p = p, fct. = fct, method = paste0("FN", HL))
  list2env(Data, envir = environment())
  if (!is.null(rownames(Y))) groups <- divide(rownames(Y), seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  Output.(rbind(vari.), ipath)
  # if (!(seed %in% (read.csv(ipath)$index))) {
  sapply(paste0("FNN/", list.files("FNN")), source)
  A <- c(rep(list(sigmoid), HL), list(linear))
  fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, loc)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, 
                   if (is.null(loc)) Y.train[, -1, drop = FALSE] else loc)
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, 
                  if (is.null(loc)) Y.test[, -1, drop = FALSE] else loc)
  cols <- if (is.null(loc)) 1 else 1:ncol(Y.train)
  error <- Error(Y.train[, cols, drop = FALSE], 
                 Y.test[, cols, drop = FALSE], Y.train., Y.test.)
  prt(format(cbind(error, index = seed)[, vari.], digits = 4), ipath, 
      pscn = TRUE)
  # }
}