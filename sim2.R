source("source.R")
source('Data/fun_output2.R')
## aligned condition without loc
sim2 <- function(seed, vari = 1, D = 3, lambda. = 10^(-2:2)) {
  vari. <- c("index", "train", "test", "cor1", "cor2", "lambda", "j")
  ipath <- paste0("../4_Output/matrix/", vari, "/FN", D - 1, ".csv")
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
  Output.(rbind(vari.), ipath)
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  ## ----------Basis----------
  la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
  lb <- floor(sqrt(la * 15 * 15))
  bb0 <- create.bspline.basis(norder = 4, nbasis = la)
  bb1 <- create.bspline.basis(norder = 4, nbasis = lb)
  bbz <- create.bspline.basis(norder = 4, nbasis = 15)
  Bases <- c(list(bb0), rep(list(bb1), D - 1), list(list(bbz, bbz)))
  fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, loc)
  Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, loc)
  Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, loc)
  
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  prt(format(cbind(error, index = seed, lambda = fnn.p$lambda, 
                   j = fnn.p$j)[, vari.], digits = 4), ipath)
  print(read.csv(ipath), row.names = FALSE)
}