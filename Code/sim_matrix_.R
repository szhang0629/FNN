source("source_.R")
source('Data/fun_output2.R')
## aligned condition without loc
sim2 <- function(seed, vari = 1, HL = 1) {
  vari. <- c("index", "train", "test", "mse1", "mse2", "cor1", "cor2")
  ipath. <- paste0("../4_Output/matrix/", vari, "/")
  set.seed(seed)
  n <- 250
  # loc <- data.frame(PTID = rep(1:n, each = 400), loc1 = rep(runif(1:400), n), 
  #                   loc2 = rep(runif(1:400), n))
  loc <- list(loc1 = runif(1:400), loc2 = runif(1:400))
  Data <- gData.2(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, 200, loc)
  list2env(Data, envir = environment())
  # rownames(G) <- 1:n
  X <- NULL
  groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  ipath <- paste0(ipath., "FLM1.csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    bb0 <- create.bspline.basis(norder = 4, nbasis = 60)
    flm1.p <- flm1(Y.train, G.train, bb0, pos, X.train, lambda. = c(10^(-3:1)))
    Y.train. <- pred.flm1(G.train, bb0, pos, flm1.p$para, X.train)
    Y.test. <- pred.flm1(G.test, bb0, pos, flm1.p$para, X.test)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed, j = 0, lambda = flm1.p$lambda)
               [, vari.], digits = 4), ipath, pscn = TRUE)
  }
  ipath <- paste0(ipath., "NN", HL, ".csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    A <- c(rep(list(sigmoid), HL), list(linear))
    E.train <- cbind(X.train, G.train)
    E.test <- cbind(X.test, G.test)
    Bases <- c(list(ncol(E.train)), as.list(rep(60, HL)), list(ncol(Y)))
    nn.p <- NN.p(Y.train, E.train, Bases, A = A, lambda. = 10^(-2:2))
    Y.train. <- pred(nn.p$Ac, E.train, A)
    Y.test. <- pred(nn.p$Ac, E.test, A)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed, lambda = nn.p$lambda, 
                     j = nn.p$j)[, vari.], digits = 4), ipath, pscn = TRUE)
  }
}
