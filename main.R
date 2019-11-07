source("source.R")
main <- function(seed, vari = 0, p = 200, D = 3, lambda. = 10^(0:1)) {
  vari. <- c("method", "train", "test", "cor1", "cor2", "lambda", "j")
  if (vari == 0) {
    ipath <- paste0("../4_Output/Real2/", seed, ".csv")
    Y <- read.csv("../2_Data/data_y.csv", as.is = TRUE)
    colnames(Y)[colnames(Y) == "AGE"] <- "loc"
    
    bbz <- create.bspline.basis(norder = 4, nbasis = 15)
    Y$loc <- loc.scale(Y$loc, bbz, 4)
    
    library(tidyverse)
    Y <- Y %>% group_by(PTID) %>% filter(n() > 1) %>% 
      arrange(loc, .by_group = TRUE) %>% mutate(loc0 = min(loc)) %>% 
      mutate(diff = Y - first(Y)) %>% slice(2:n()) %>% ungroup()
    Y <- as.data.frame(Y) 
    Y <- subset(Y, select = -Y)
    colnames(Y)[which(names(Y) == "diff")] <- "Y"
    Y$Y <- Y$Y/sd(Y$Y)
    # X <- read.csv("../2_Data/data_x_APOE.csv", row.names = 1)
    # X <- subset(X, select = -APOE4)
    # X <- as.matrix(cbind(subset(Y, select = -c(PTID, Y)), X[Y$PTID, ]))
    X <- as.matrix(read.csv("../2_Data/data_x_APOE.csv", row.names = 1))
    # X <- subset(X, select = -APOE4)
    G <- as.matrix(read.csv("../2_Data/data_g.csv", row.names = 1))
    pos <- read.csv("../2_Data/data_pos.csv")$POS
  } else {
    ipath <- paste0("../4_Output/", p, "/", vari, "/", seed, ".csv")
    n <- 250
    set.seed(seed)
    # loc <- data.frame(PTID = rep(1:n, each = 20), loc = runif(20*n))
    # loc <- data.frame(PTID = rep(1:n, each = 20), loc = rep(runif(20), n))
    loc <- runif(20)
    Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, "polyno")
    list2env(Data, envir = environment())
    rownames(G) <- 1:n
    
    bbz <- create.bspline.basis(norder = 4, nbasis = 20)
    if (is.data.frame(loc)) Y$loc <- loc.scale(Y$loc, bbz, 4)
    else loc <- loc.scale(loc, bbz, 4)
    
    X <- NULL
  }
  if (is.data.frame(Y)) groups <- divide(Y$PTID, seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  # groups <- list(train = unique(Y$PTID)[-seed], test = unique(Y$PTID)[seed])
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  Output.(rbind(vari.), ipath)
  method <- paste0("FN", D - 1)
  # if (!(method %in% (read.csv(ipath)$method))) {
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
  # l <- round(ncol(G.train)^(1/3) * 20^(2/3))
  # la <- 15
  l <- 15
  bb0 <- create.bspline.basis(norder = 4, nbasis = la)
  bb1 <- create.bspline.basis(norder = 4, nbasis = l)
  Bases <- c(list(bb0), rep(list(bb1), D - 1), list(bbz))
  
  fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, 
                 (if (is.list(Y.train)) NULL else loc))
  if (is.data.frame(Y.train)) {
    Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, 
                     loc = subset(Y.train, select = -Y))
    Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, 
                    loc = subset(Y.test, select = -Y))
    error <- Error(Y.train$Y, Y.test$Y, Y.train., Y.test.)
  } else {
    Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, loc)
    Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, loc)
    error <- Error(c(Y.train), c(Y.test), c(Y.train.), c(Y.test.))
  }
  error.fn <- cbind(error, data.frame(lambda = fnn.p$lambda, j = fnn.p$j))
  
  output <- format(cbind(error.fn, method = method)[, vari.], digits = 4)
  prt(output, ipath)
  # }
  print(read.csv(ipath), row.names = FALSE)
}