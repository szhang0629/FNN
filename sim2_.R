source("source_.R")
source('Data/fun_output2.R')
## aligned condition without loc
sim2 <- function(seed, vari = 1, D = 3) {
  vari. <- c("index", "train", "test", "cor1", "cor2", "lambda", "j")
  ipath <- paste0("../4_Output/matrix/", vari, "/", seed, ".csv")
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
    la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
    bb0 <- create.bspline.basis(norder = 4, nbasis = la)
    error <- Error.flm1(Y.train, G.train, Y.test, G.test, bb0, pos, X.train,
                        X.test, lambda. = c(10^(-3:1)))
    output <- format(cbind(error, method = "FLM1", j = 0)[, vari.], digits = 4)
    prt(output, ipath)
  }
  ipath <- paste0(ipath., "NN", D - 1, ".csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    A <- c(rep(list(sigmoid), D - 1), list(linear))
    la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
    lb <- floor(sqrt(la * 15 * 15))
    Bases <- c(list(ncol(cbind(X.train, G.train))), as.list(rep(lb, D - 1)), 
               list(ncol(Y)))
    error.nn <- Error.nn(Y.train, cbind(X.train, G.train), Y.test, 
                         cbind(X.test, G.test), Bases, A = A, lambda. = 10^(1))
    output <- format(cbind(error.nn, method = method1)[, vari.], digits = 4)
    prt(output, ipath)
  }
  print(read.csv(ipath), row.names = FALSE)
}
